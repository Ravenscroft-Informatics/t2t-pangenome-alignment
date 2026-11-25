#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --time=0-04:00:00
#SBATCH --job-name=ref_vs_alt
#SBATCH --partition=work
#SBATCH --account=pawsey0933
#SBATCH --export=NONE
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/scratch/pawsey0933/cfolland/gaf/output/logs/%x_%A_%a.out
#SBATCH --error=/scratch/pawsey0933/cfolland/gaf/output/logs/%x_%A_%a.err

# ---------------------------------------------------------------------------
# Classify giraffe GAF reads as reference-only vs non-reference (one GAF per task).
#
# USAGE (run with bash, NOT sbatch -- it builds the to-do list and submits itself):
#     bash ref_vs_alt_array.sh            # submit all GAFs without finished output
#     MAX_PARALLEL=20 bash ref_vs_alt_array.sh
#
# A GAF counts as done when <id>.ref_vs_alt.tsv exists, is non-empty, and its
# header has all columns in required_cols below (current script version).
# Outputs from older versions are rerun and overwritten.
# Outputs are written to .tmp files and renamed only on success, so a task
# that times out or crashes leaves nothing behind and is picked up next time.
# ---------------------------------------------------------------------------
set -euo pipefail

py=/software/projects/pawsey0933/benchmarking/hprc/gaf_ref_vs_alt.py
ref_gaf=/scratch/pawsey0933/cfolland/vg_snakemake/refs/hprcv2.0/grch38.paths.gaf
gafdir=/scratch/pawsey0933/cfolland/gaf
outdir=${gafdir}/output
suffix=.hprc-v2.0-mc-grch38.gaf.gz
min_mapq=0

# Done = non-empty TSV whose header contains every required column and that has
# none of the retired GRCh38-orientation lines (from the short-lived previous version)
required_cols=(perfect_reads frac_perfect_in_category mapq0_reads proper_pair_reads mismatch_rate error_rate)
is_done() {
    local tsv=$1 col header
    [[ -s "$tsv" ]] || return 1
    header=$(head -n 1 "$tsv")
    for col in "${required_cols[@]}"; do
        [[ $'\t'"$header"$'\t' == *$'\t'"$col"$'\t'* ]] || return 1
    done
    grep -q '^## pair_fr_inward' "$tsv" && return 1
    return 0
}

# ============================ SUBMIT MODE ==================================
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    if [[ -n "${SLURM_JOB_ID:-}" ]]; then
        echo "ERROR: run this script with 'bash', not 'sbatch'." >&2
        exit 1
    fi

    mkdir -p "${outdir}/logs"
    list=${outdir}/logs/todo_gafs.$(date +%Y%m%d_%H%M%S).txt
    : > "$list"

    n_total=0; n_old=0
    for f in "${gafdir}"/*"${suffix}"; do
        [[ -e "$f" ]] || continue
        n_total=$((n_total + 1))
        id=$(basename "$f" "$suffix")
        tsv=${outdir}/${id}.ref_vs_alt.tsv
        is_done "$tsv" && continue
        if [[ -s "$tsv" ]]; then
            n_old=$((n_old + 1))
            echo "  outdated (missing new columns): ${id}"
        fi
        echo "$f" >> "$list"
    done

    n=$(wc -l < "$list")
    echo "GAFs found: ${n_total}   done: $((n_total - n))   to submit: ${n}" \
         "(${n_old} outdated, $((n - n_old)) new)"
    if (( n == 0 )); then
        rm -f "$list"
        echo "Nothing to do."
        exit 0
    fi

    sbatch --array=1-${n}%${MAX_PARALLEL:-50} "$(realpath "$0")" "$list"
    echo "To-do list: $list"
    exit 0
fi

# ============================= TASK MODE ===================================
list=$1
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$list")
id=$(basename "$sample" "$suffix")

ref_alt=${outdir}/${id}.ref_vs_alt.tsv
non_ref=${outdir}/${id}.nonref_reads.txt.gz

echo "[$(date)] task ${SLURM_ARRAY_TASK_ID}: ${id}"
echo "  input : ${sample}"
echo "  output: ${ref_alt}"

# Re-check in case another submission finished this sample in the meantime
if is_done "$ref_alt"; then
    echo "Up-to-date output already exists, skipping."
    exit 0
elif [[ -s "$ref_alt" ]]; then
    echo "Existing output is from an older script version; rerunning and overwriting."
fi

# Decompress in a separate process (2nd CPU) -- faster than Python's gzip
if command -v pigz > /dev/null; then decomp=(pigz -dc -p 1); else decomp=(gzip -dc); fi

run() {   # run [extra args...]
    "${decomp[@]}" "$sample" | python3 "$py" \
        --ref-gaf "$ref_gaf" \
        --reads - \
        --min-mapq "$min_mapq" \
        --nonref-names "${non_ref%.gz}.tmp.gz" \
        "$@" \
        > "${ref_alt}.tmp"
}

# Reuse the existing nonref_reads file (from an earlier run on this GAF, same
# min_mapq) to skip classifying ref_only reads. If it doesn't match the GAF
# (e.g. the GAF was regenerated), fall back to a full run.
if [[ -s "$non_ref" ]]; then
    echo "  mode  : fast (--known-nonref ${non_ref})"
    if ! run --known-nonref "$non_ref"; then
        echo "Fast mode failed (see .err); retrying full classification."
        run
    fi
else
    echo "  mode  : full"
    run
fi

mv "${non_ref%.gz}.tmp.gz" "$non_ref"
mv "${ref_alt}.tmp" "$ref_alt"

if ! is_done "$ref_alt"; then
    echo "ERROR: new output still lacks ${required_cols[*]} -- is ${py} the updated script?" >&2
    exit 1
fi

echo "[$(date)] done: ${id}"