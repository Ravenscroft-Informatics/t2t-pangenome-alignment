#!/bin/bash
set -euo pipefail

# Inputs
gtf_chm13=/software/projects/pawsey0933/benchmarking/nmd_genes/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz
gtf_hg38=/software/projects/pawsey0933/benchmarking/nmd_genes/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz

# Outputs (transcripts, not gene bodies)
out_chm13=/software/projects/pawsey0933/benchmarking/nmd_genes/GCF_009914755.1_T2T-CHM13v2.0.transcripts.bed
out_hg38=/software/projects/pawsey0933/benchmarking/nmd_genes/GCF_000001405.40_GRCh38.p14.transcripts.bed

make_transcripts_bed () {
  local in_gtf="$1"
  local out_bed="$2"

  # Output columns: chrom (as in file, with optional chr-ification for simple names), 0-based start, end, gene_name, transcript_id
  ( [[ "$in_gtf" = *.gz ]] && zcat "$in_gtf" || cat "$in_gtf" ) \
  | awk -v FS='\t' -v OFS='\t' '
      BEGIN {
        # regex for simple chroms that may need "chr" prefix
        simple="^(?:[0-9]+|[XY]|M|MT)$"
      }
      $0 ~ /^#/ { next }
      $3 == "transcript" {
        info = $9

        # --- transcript_id (keep version if present) ---
        tid = ""
        if (match(info, /transcript_id[= ]"([^"]+)"/, m)) { tid = m[1] }
        else if (match(info, /transcript_id=([^;]+)/, m))  { tid = m[1] }

        # --- gene_name ---
        # Prefer gene_name, then gene, then gene_id (NCBI files often put symbol in gene/gene_id)
        gnm = ""
        if (match(info, /gene_name[= ]"([^"]+)"/, n))      { gnm = n[1] }
        else if (match(info, /gene_name=([^;]+)/, n))      { gnm = n[1] }
        else if (match(info, /gene[= ]"([^"]+)"/, n))      { gnm = n[1] }
        else if (match(info, /gene=([^;]+)/, n))           { gnm = n[1] }
        else if (match(info, /gene_id[= ]"([^"]+)"/, n))   { gnm = n[1] }
        else if (match(info, /gene_id=([^;]+)/, n))        { gnm = n[1] }

        if (tid == "" || gnm == "") next

        chrom = $1
        # Only add "chr" for simple names; leave NC_..., GL..., KI..., etc. untouched
        if (chrom ~ simple) {
          if (chrom == "MT") chrom = "chrM"
          else if (chrom == "M") chrom = "chrM"
          else chrom = "chr" chrom
        }

        start = $4 - 1   # BED 0-based start
        end   = $5

        print chrom, start, end, gnm, tid
      }
    ' > "$out_bed"
}

# Build both transcript beds
make_transcripts_bed "$gtf_chm13" "$out_chm13"
make_transcripts_bed "$gtf_hg38"  "$out_hg38"

echo "Wrote:"
echo "  $out_chm13"
echo "  $out_hg38"
