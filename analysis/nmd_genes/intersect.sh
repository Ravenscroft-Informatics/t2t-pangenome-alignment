#!/bin/bash
set -euo pipefail

# =============================
# Configuration
# =============================
base_dir=/software/projects/pawsey0933/benchmarking/nmd_genes
bed_chm13=${base_dir}/GCF_009914755.1_T2T-CHM13v2.0.transcripts.bed
bed_hg38=${base_dir}/GCF_000001405.40_GRCh38.p14.transcripts.bed

out_ids=${base_dir}/shared.transcript_ids.txt
out_chm13=${base_dir}/GCF_009914755.1_T2T-CHM13v2.0.transcripts.shared.bed
out_hg38=${base_dir}/GCF_000001405.40_GRCh38.p14.transcripts.shared.bed
out_join=${base_dir}/transcripts.chm13_vs_hg38.join.tsv

# Temporary sorted ID lists (kept in benchmark dir, not /tmp)
ids_hg38=${base_dir}/hg38.ids.txt
ids_chm13=${base_dir}/chm13.ids.txt

# =============================
# Optional flag
# =============================
IGNORE_VERSION=0  # default off

usage() {
  cat <<EOF
Usage: $(basename "$0") [--ignore-version]

Find transcripts shared between CHM13 and GRCh38 by transcript_id (col 5).

Outputs:
  - Shared IDs:           $out_ids
  - CHM13 shared BED:     $out_chm13
  - GRCh38 shared BED:    $out_hg38
  - Joined crosswalk TSV: $out_join

Options:
  --ignore-version   Ignore transcript version suffix (e.g. NM_001347931.2 -> NM_001347931)
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --ignore-version) IGNORE_VERSION=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1"; usage; exit 1 ;;
  esac
done

for f in "$bed_chm13" "$bed_hg38"; do
  [[ -s "$f" ]] || { echo "ERROR: Missing or empty file: $f" >&2; exit 1; }
done

# =============================
# Step 1: Collect unique transcript IDs
# =============================
echo "Extracting transcript IDs..."
if [[ $IGNORE_VERSION -eq 1 ]]; then
  awk -v OFS='\t' '{id=$5; sub(/\..*$/, "", id); print id}' "$bed_hg38" | LC_ALL=C sort -u > "$ids_hg38"
  awk -v OFS='\t' '{id=$5; sub(/\..*$/, "", id); print id}' "$bed_chm13" | LC_ALL=C sort -u > "$ids_chm13"
else
  cut -f5 "$bed_hg38"  | LC_ALL=C sort -u > "$ids_hg38"
  cut -f5 "$bed_chm13" | LC_ALL=C sort -u > "$ids_chm13"
fi

# =============================
# Step 2: Intersect the IDs
# =============================
echo "Finding shared transcript IDs..."
awk 'NR==FNR {a[$0]; next} ($0 in a)' "$ids_hg38" "$ids_chm13" | LC_ALL=C sort -u > "$out_ids"
echo "Shared transcript count: $(wc -l < "$out_ids")"

# =============================
# Step 3: Filter each BED by shared IDs
# =============================
echo "Filtering shared transcripts..."
if [[ $IGNORE_VERSION -eq 1 ]]; then
  awk 'NR==FNR {keep[$0]; next}
       {tid=$5; sub(/\..*$/, "", tid); if (tid in keep) print $0}' "$out_ids" "$bed_chm13" > "$out_chm13"

  awk 'NR==FNR {keep[$0]; next}
       {tid=$5; sub(/\..*$/, "", tid); if (tid in keep) print $0}' "$out_ids" "$bed_hg38" > "$out_hg38"
else
  awk 'NR==FNR {keep[$0]; next} ($5 in keep)' "$out_ids" "$bed_chm13" > "$out_chm13"
  awk 'NR==FNR {keep[$0]; next} ($5 in keep)' "$out_ids" "$bed_hg38"  > "$out_hg38"
fi

echo "CHM13 shared BED rows: $(wc -l < "$out_chm13")"
echo "GRCh38 shared BED rows: $(wc -l < "$out_hg38")"

# =============================
# Step 4: Create joined crosswalk table
# =============================
echo "Building crosswalk..."
printf "transcript_id\tgene_name_chm13\tchrom_chm13\tstart_chm13\tend_chm13\tgene_name_hg38\tchrom_hg38\tstart_hg38\tend_hg38\n" > "$out_join"

if [[ $IGNORE_VERSION -eq 1 ]]; then
  awk '
    NR==FNR {
      tid=$5; sub(/\..*$/, "", tid);
      a[tid]=$0;
      next
    }
    {
      tid=$5; sub(/\..*$/, "", tid);
      if (tid in a) {
        split(a[tid], L, "\t");
        split($0, R, "\t");
        tid_out=R[5]; if (tid_out=="") tid_out=L[5];
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
               tid_out, L[4], L[1], L[2], L[3], R[4], R[1], R[2], R[3];
      }
    }' "$out_chm13" "$out_hg38" >> "$out_join"
else
  awk '
    NR==FNR {a[$5]=$0; next}
    ($5 in a) {
      split(a[$5], L, "\t");
      split($0, R, "\t");
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
             $5, L[4], L[1], L[2], L[3], R[4], R[1], R[2], R[3];
    }' "$out_chm13" "$out_hg38" >> "$out_join"
fi

# =============================
# Cleanup & Summary
# =============================
rm -f "$ids_hg38" "$ids_chm13"
echo "Joined crosswalk written to: $out_join"
echo "All outputs saved under: $base_dir"
