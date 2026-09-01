#!/bin/bash
# submit_pipeline.sh
#
# Drives submission of a (possibly huge) tasks.tsv against the work
# partition without ever exceeding MaxSubmit=1024 (per-user, per-partition
# live job count on pawsey0933/work), and without needing you to babysit
# it manually.
#
# Because it polls in a loop for potentially hours/days, this is exactly
# the kind of long-running "master task" Pawsey's Setonix WORKFLOW NODES
# exist for (setonix-07/setonix-08, dedicated to launching long-running
# workflow-manager processes like Nextflow/Snakemake controllers), rather
# than an ordinary login node -- see:
#   https://pawsey.atlassian.net/wiki/spaces/US/pages/286097469/How+to+Run+Workflows+on+the+Workflow+Nodes
# Use screen/tmux there to survive disconnects (both are modules on the
# workflow nodes, not preinstalled -- module load first):
#   module load tmux
#   tmux new -s deepvariant
#   submit_pipeline.sh tasks.tsv
#   # Ctrl-b d to detach; `tmux attach -t deepvariant` to check back in
# or:
#   nohup bash submit_pipeline.sh tasks.tsv > pipeline.log 2>&1 &
#   disown
#
# Usage:
#   mkdir -p logs
#   bash submit_pipeline.sh tasks.tsv
#
# Tune the knobs below if your situation changes (e.g. Pawsey adjusts your
# account limits, or the work partition gets much busier/quieter).

set -uo pipefail

MANIFEST_ARG="${1:?Usage: bash submit_pipeline.sh tasks.tsv}"
MANIFEST="$(cd "$(dirname "${MANIFEST_ARG}")" && pwd)/$(basename "${MANIFEST_ARG}")"
[[ -f "${MANIFEST}" ]] || { echo "ERROR: manifest not found: ${MANIFEST}" >&2; exit 1; }

SLICE_SIZE=200     # tasks per sbatch array call. Deliberately small
                   # relative to SAFE_CAP so several slices can be in
                   # flight at once (overlap) instead of one slice having
                   # to nearly fully drain before the next can start.
SAFE_CAP=1000      # don't let your total LIVE (pending+running) job count
                   # on the work partition exceed this before submitting
                   # the next slice -- leaves ~24 headroom under the 1024
                   # hard MaxSubmit limit for any other work-partition jobs
                   # you run. No per-slice %N throttle is used: MaxJobs=256
                   # (account-wide) and the fact each task needs ~1 full
                   # node's memory (220gb of ~239gb/node) already bound how
                   # many can actually run at once -- a %N throttle inside
                   # a slice does NOT reduce how much of MaxSubmit it uses
                   # (all N array elements count as "submitted" the moment
                   # you sbatch, throttled or not), so it was dropped.
POLL_SECONDS=120   # how often to re-check headroom while waiting
RETRY_SECONDS=120  # backoff before retrying a failed sbatch call

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ARRAY_SCRIPT="${SCRIPT_DIR}/deepvariant_array.slurm"

mkdir -p logs

TOTAL=$(wc -l < "${MANIFEST}")
echo "[$(date '+%F %T')] ${TOTAL} total tasks in ${MANIFEST}"
echo "[$(date '+%F %T')] slice size ${SLICE_SIZE}, safe cap ${SAFE_CAP} live jobs"

# live (pending+running) job count for this user on the work partition.
# IMPORTANT: use -r/--array so Slurm unrolls pending array jobs into one
# line per array element (e.g. "47268092_[1-900]" -> 900 lines). Without
# -r, squeue collapses an entire pending array into a single summary line
# and this count is wildly wrong (this is what caused the first run to
# blow straight through MaxSubmit).
live_count() {
  squeue -u "${USER}" -h -p work -r | wc -l
}

# submit an array over manifest line range [start, end], retrying with
# backoff on transient sbatch failures instead of killing the pipeline.
submit_slice() {
  local start="$1" end="$2" out ec
  while true; do
    out=$(sbatch --parsable --array="${start}-${end}" "${ARRAY_SCRIPT}" "${MANIFEST}" 2>&1)
    ec=$?
    if [[ ${ec} -eq 0 ]]; then
      echo "${out}"
      return 0
    fi
    echo "[$(date '+%F %T')] sbatch failed (exit ${ec}) for lines ${start}-${end}: ${out}" >&2
    echo "[$(date '+%F %T')] retrying in ${RETRY_SECONDS}s" >&2
    sleep "${RETRY_SECONDS}"
  done
}

SUBMITTED_JOBIDS=()

start=1
while (( start <= TOTAL )); do
  end=$(( start + SLICE_SIZE - 1 ))
  (( end > TOTAL )) && end=${TOTAL}
  n=$(( end - start + 1 ))

  # wait until submitting this slice won't push us over SAFE_CAP
  while true; do
    LIVE=$(live_count)
    if (( LIVE + n <= SAFE_CAP )); then
      break
    fi
    echo "[$(date '+%F %T')] ${LIVE} live jobs on work, next slice needs ${n} (cap ${SAFE_CAP}) -- waiting ${POLL_SECONDS}s"
    sleep "${POLL_SECONDS}"
  done

  echo "[$(date '+%F %T')] submitting lines ${start}-${end} (${n} tasks)"
  JOBID=$(submit_slice "${start}" "${end}")
  echo "[$(date '+%F %T')]   -> array job ${JOBID}"
  SUBMITTED_JOBIDS+=("${JOBID}")

  start=$(( end + 1 ))
done

echo "[$(date '+%F %T')] all slices submitted: ${SUBMITTED_JOBIDS[*]}"
echo "[$(date '+%F %T')] waiting for everything to finish..."

for JOBID in "${SUBMITTED_JOBIDS[@]}"; do
  while squeue -h -j "${JOBID}" 2>/dev/null | grep -q .; do
    sleep "${POLL_SECONDS}"
  done
done

echo "[$(date '+%F %T')] done. Check for stragglers with:"
echo "  bash make_missing.sh ${MANIFEST} | wc -l"