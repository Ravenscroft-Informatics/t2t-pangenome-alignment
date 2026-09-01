#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH -c 32
#SBATCH --job-name=deepvariant
#SBATCH --partition=work
#SBATCH --account=pawsey0933
#SBATCH --mem=200gb
#SBATCH --time=0-2:30:00
#SBATCH --export=NONE
#SBATCH --mail-type=END
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#
# One Slurm ARRAY TASK = one (sample, region) DeepVariant run.
# Each task reads its own row out of the manifest file using
# SLURM_ARRAY_TASK_ID as a 1-based line number.
#
# --- how to submit -------------------------------------------------------
#   mkdir -p logs
#   N=$(wc -l < tasks.tsv)
#   sbatch --array=1-${N}%25 deepvariant_array.slurm tasks.tsv
#
# %25 caps concurrency at 25 simultaneously-running tasks (25 x 32 cores =
# 800 cores, 25 x 220gb = 5.5TB mem in flight at once). Check your account's
# actual limits with `sacctmgr show assoc account=pawsey0933 format=account,
# maxsubmitjobs,maxsubmitjobsperuser` (or ask Pawsey support) and raise/lower
# the %N accordingly. If unsure, start conservative (e.g. %10-%20) and raise
# it once you see jobs completing cleanly.
#
# --- resuming after partial failures --------------------------------------
# This script SKIPS any (sample, region) whose output vcf+gvcf already
# exist, so it is always safe to just re-run the same sbatch command again
# on the same tasks.tsv. To avoid re-queuing thousands of instant no-op
# "skip" jobs, instead build a manifest of only what's missing:
#   make_missing.sh tasks.tsv > tasks_missing.tsv
#   N=$(wc -l < tasks_missing.tsv)
#   sbatch --array=1-${N}%25 deepvariant_array.slurm tasks_missing.tsv
# ---------------------------------------------------------------------------

set -euo pipefail

module load singularity/4.1.0-slurm

MANIFEST="${1:?Usage: sbatch --array=1-N deepvariant_array.slurm tasks.tsv}"

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${MANIFEST}")
if [[ -z "${LINE}" ]]; then
  echo "No manifest line for array task ${SLURM_ARRAY_TASK_ID} in ${MANIFEST}" >&2
  exit 1
fi
SAMPLE=$(cut -f1 <<< "${LINE}")
CRAM=$(cut -f2 <<< "${LINE}")
REGION=$(cut -f3 <<< "${LINE}")

SIF=/software/projects/pawsey0933/deepvariant/deepvariant_pangenome_aware_deepvariant-1.10.0.sif
OUTPUT_DIR=/scratch/pawsey0933/cfolland/deepvariant/hprc/pg_dv
REF=/software/projects/pawsey0933/pangenome/refs/hprc-v2.0-mc-grch38.fa
PANGENOME=/scratch/pawsey0933/cfolland/vg_snakemake/refs/hprcv2.0/hprc-v2.0-mc-grch38.gbz
REGIONS_DIR=/software/projects/pawsey0933/deepvariant/regions

OUT_VCF="${OUTPUT_DIR}/${SAMPLE}.deepvariant.${REGION}.vcf.gz"
OUT_GVCF="${OUTPUT_DIR}/${SAMPLE}.deepvariant.${REGION}.g.vcf.gz"

mkdir -p "${OUTPUT_DIR}"

if [[ -s "${OUT_VCF}" && -s "${OUT_GVCF}" ]]; then
  echo "${SAMPLE} / ${REGION}: outputs already exist, skipping"
  exit 0
fi

if [[ ! -f "${CRAM}" ]]; then
  echo "ERROR: CRAM not found for sample ${SAMPLE}: ${CRAM}" >&2
  exit 1
fi

TMP_DIR="${OUTPUT_DIR}/tmp_${SLURM_JOB_ID}"
mkdir -p "${TMP_DIR}/shm"

cleanup() {
  local ec=$?
  if [[ $ec -eq 0 ]]; then
    rm -rf "${TMP_DIR}"
  else
    echo "Job exited with status ${ec}; leaving ${TMP_DIR} in place for debugging" >&2
  fi
}
trap cleanup EXIT

export XDG_CONFIG_HOME=/scratch/pawsey0933/cfolland/tmp/config
export XDG_CACHE_HOME=/scratch/pawsey0933/cfolland/tmp/cache
mkdir -p "${XDG_CONFIG_HOME}" "${XDG_CACHE_HOME}"

echo "Running DeepVariant: sample=${SAMPLE} region=${REGION}"
echo "  cram=${CRAM}"
echo "  output=${OUT_VCF}"

singularity exec -B /scratch/pawsey0933 -B /software/projects/pawsey0933 -B "${TMP_DIR}/shm:/dev/shm" \
    "${SIF}" \
    /opt/deepvariant/bin/run_pangenome_aware_deepvariant \
    --model_type WGS \
    --ref "${REF}" \
    --reads "${CRAM}" \
    --output_vcf "${OUT_VCF}" \
    --output_gvcf "${OUT_GVCF}" \
    --regions "${REGIONS_DIR}/${REGION}_"*.bed \
    --pangenome "${PANGENOME}" \
    --intermediate_results_dir "${TMP_DIR}" \
    --gbz_shared_memory_size_gb 40 \
    --num_shards 32 \
    --ref_name_pangenome GRCh38