#!/bin/bash
# make_manifest.sh
#
# Builds tasks.tsv: one line per (sample, region) DeepVariant job, as
# SAMPLE<TAB>CRAM_PATH<TAB>REGION
# make_manifest.sh > tasks.tsv

set -euo pipefail

CRAM_DIR="/scratch/pawsey0933/cfolland/crams_hprc"
CRAM_SUFFIX=".hprc-v2.0-mc-grch38.surj_realn.cram"
REGION_DIR="/software/projects/pawsey0933/deepvariant/regions"

region_from_filename() { basename "$1" | cut -d'_' -f1; }

mapfile -t SAMPLES < <(
  for f in "${CRAM_DIR}"/*"${CRAM_SUFFIX}"; do
    [[ -e "$f" ]] || continue
    b=$(basename "$f")
    echo "${b%"${CRAM_SUFFIX}"}"
  done | sort -u
)

mapfile -t REGIONS < <(
  for f in "${REGION_DIR}"/*.bed; do
    [[ -e "$f" ]] || continue
    region_from_filename "$f"
  done | sort -u
)

if [[ ${#SAMPLES[@]} -eq 0 ]]; then
  echo "ERROR: no CRAMs found under ${CRAM_DIR} matching *${CRAM_SUFFIX}" >&2
  exit 1
fi
if [[ ${#REGIONS[@]} -eq 0 ]]; then
  echo "ERROR: no region beds found under ${REGION_DIR}" >&2
  exit 1
fi

echo "Found ${#SAMPLES[@]} samples x ${#REGIONS[@]} regions = $(( ${#SAMPLES[@]} * ${#REGIONS[@]} )) total jobs" >&2

for s in "${SAMPLES[@]}"; do
  cram="${CRAM_DIR}/${s}${CRAM_SUFFIX}"
  for r in "${REGIONS[@]}"; do
    printf '%s\t%s\t%s\n' "$s" "$cram" "$r"
  done
done