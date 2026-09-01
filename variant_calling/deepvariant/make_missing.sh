#!/bin/bash
# make_missing.sh
#
# Filters a tasks.tsv manifest down to just the (sample, region) rows whose
# output vcf/gvcf don't exist yet (e.g. because they failed, timed out, or
# were never submitted). Use this to cheaply resubmit only what's left
# instead of re-running the full array.
#
#   bash make_missing.sh tasks.tsv > tasks_missing.tsv
#   N=$(wc -l < tasks_missing.tsv)
#   echo "${N} jobs left"
#   sbatch --array=1-${N}%25 deepvariant_array.slurm tasks_missing.tsv

set -euo pipefail

MANIFEST="${1:-tasks.tsv}"
OUTPUT_DIR=/scratch/pawsey0933/cfolland/deepvariant/hprc/pg_dv

while IFS=$'\t' read -r SAMPLE CRAM REGION; do
  VCF="${OUTPUT_DIR}/${SAMPLE}.deepvariant.${REGION}.vcf.gz"
  GVCF="${OUTPUT_DIR}/${SAMPLE}.deepvariant.${REGION}.g.vcf.gz"
  if [[ ! -s "${VCF}" || ! -s "${GVCF}" ]]; then
    printf '%s\t%s\t%s\n' "$SAMPLE" "$CRAM" "$REGION"
  fi
done < "${MANIFEST}"