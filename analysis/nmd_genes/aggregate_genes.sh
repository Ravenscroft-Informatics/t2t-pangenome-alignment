#!/usr/bin/env bash
set -euo pipefail

# Usage: ./aggregate_gene_coords.sh input.bed output.bed
# Input columns expected: chrom  start  end  gene  transcript_id
# Output columns: chrom  gene_start  gene_end  gene  n_transcripts

in="${1:?BED input required}"
out="${2:?BED output required}"

# Sort by chrom, gene, then start for stable aggregation
# - ensure bytewise sort for speed and predictability on HPC
LC_ALL=C sort -k1,1 -k4,4 -k2,2n "$in" \
| awk -v OFS='\t' '
  BEGIN {
    print "chrom","gene_start","gene_end","gene","n_transcripts"
  }
  {
    chrom=$1; start=$2; end=$3; gene=$4;
    key = chrom FS gene
    if (!(key in seen)) {
      # first time we see this chrom+gene
      min[key]=start
      max[key]=end
      cnt[key]=1
      seen[key]=1
      order[++i]=key
    } else {
      if (start < min[key]) min[key]=start
      if (end   > max[key]) max[key]=end
      cnt[key]++
    }
  }
  END {
    for (j=1; j<=i; j++) {
      k = order[j]
      split(k, a, FS)   # a[1]=chrom, a[2]=gene
      print a[1], min[k], max[k], a[2], cnt[k]
    }
  }' \
| LC_ALL=C sort -k1,1 -k2,2n -k3,3n -k4,4 > "$out"
