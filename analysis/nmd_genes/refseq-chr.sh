#!/usr/bin/env bash
set -euo pipefail
# Usage: ./replace_refseq_in_bed.sh refseq_map.tsv input.bed output.bed

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 refseq_map.tsv input.bed output.bed" >&2
  exit 1
fi

mapfile="$1"
bedfile="$2"
outfile="$3"

awk -v OFS='\t' -v mapf="$mapfile" '
BEGIN {
  # Build lookup: ref2chr["NC_060925.1"] = "1"
  while ((getline line < mapf) > 0) {
    if (line ~ /^[[:space:]]*$/ || line ~ /^#/) continue
    n = split(line, a, /[[:space:]]+/)
    if (n >= 2) ref2chr[a[2]] = a[1]
  }
  close(mapf)
}
{
  # Cases to support:
  # 1) chrNC_060925.1  -> chr1
  # 2) NC_060925.1     -> chr1    (no leading "chr" in input)
  if ($1 ~ /^chrNC_[0-9]+\.[0-9]+$/) {
    ref = substr($1, 4)              # drop "chr" -> NC_...
    if (ref in ref2chr) $1 = "chr" ref2chr[ref]; else fprintf(stderr, "[WARN] No map for %s (line %d)\n", $1, NR)
  } else if ($1 ~ /^NC_[0-9]+\.[0-9]+$/) {
    ref = $1                          # already NC_...
    if (ref in ref2chr) $1 = "chr" ref2chr[ref]; else fprintf(stderr, "[WARN] No map for %s (line %d)\n", $1, NR)
  }
  print
}' "$bedfile" > "$outfile"

echo "✅ Wrote: $outfile"
