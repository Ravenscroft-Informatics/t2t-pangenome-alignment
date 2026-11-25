#!/bin/bash
# Build the 2-sample (giraffe, bwa) VCF for plot_allele_bias.py from existing
# DeepVariant VCFs, instead of re-calling with bcftools mpileup.
#
# Usage:
#   prep_allele_bias_deepvariant.sh <sample> <giraffe_dv.vcf.gz> <bwa_dv.vcf.gz> <ref.fasta> <out_dir>
#
# Output sample order is giraffe then bwa (what the python script assumes).
set -euo pipefail

sample="${1:?sample}"
giraffe_vcf="${2:?giraffe DeepVariant VCF (not gVCF)}"
bwa_vcf="${3:?bwa DeepVariant VCF (not gVCF)}"
ref="${4:?GRCh38 fasta used for left-alignment}"
out="${5:?output dir}/${sample}"
threads="${THREADS:-8}"

# Filter thresholds. DeepVariant leaves INFO empty (no INFO/DP or INFO/MQ),
# so per-sample FORMAT/DP and FORMAT/GQ are used instead.
MIN_DP="${MIN_DP:-25}"
MIN_GQ="${MIN_GQ:-20}"

mkdir -p "$out"/{giraffe,bwa,merged}
[[ -f "${ref}.fai" ]] || samtools faidx "$ref"

prep_one () {
  local mapper="$1" invcf="$2" wd="$out/$1"
  echo "== Preparing ${mapper}: ${invcf} =="

  # 1. Strip pangenome-style contig prefixes (e.g. GRCh38#0#chr1 -> chr1) so
  #    both VCFs and the reference share names. No-op if names are already plain.
  bcftools view -h "$invcf" \
    | sed -n 's/^##contig=<ID=\([^,>]*\).*/\1/p' \
    | awk '{n=$1; sub(/^GRCh38#[0-9]+#/, "", n); sub(/^GRCh38\./, "", n); print $1"\t"n}' \
    > "$wd/chr_rename.txt"

  # 2. Rename the sample to exactly "<mapper>"
  printf "%s\n" "$mapper" > "$wd/sample_name.txt"

  # 3. PASS only (drops RefCall / ./. rows), biallelic only (the plot uses the
  #    first ALT for indel length and sums ALT depths), left-align vs ref so the
  #    same indel gets the same representation in both VCFs.
  bcftools reheader -s "$wd/sample_name.txt" -o "$wd/${mapper}.renamed.vcf.gz" "$invcf"
  bcftools annotate --threads "$threads" --rename-chrs "$wd/chr_rename.txt" -Ou "$wd/${mapper}.renamed.vcf.gz" \
    | bcftools view -f PASS -m2 -M2 -Ou \
    | bcftools norm -f "$ref" -c w -Ou \
    | bcftools sort -Oz -T "$wd/sort_tmp" -o "$wd/${mapper}.prepped.vcf.gz"
  bcftools index -f "$wd/${mapper}.prepped.vcf.gz"
  rm -f "$wd/${mapper}.renamed.vcf.gz"
}

prep_one giraffe "$giraffe_vcf"
prep_one bwa     "$bwa_vcf"

# Merge; -m none keeps records with different ALTs separate, so a site only
# counts as "het in both" when both callers agree on the exact allele.
bcftools merge -m none --threads "$threads" -Oz -o "$out/merged/merged.giraffe_bwa.vcf.gz" \
  "$out/giraffe/giraffe.prepped.vcf.gz" \
  "$out/bwa/bwa.prepped.vcf.gz"
bcftools index -f "$out/merged/merged.giraffe_bwa.vcf.gz"

# High-confidence het in BOTH (per-sample DP/GQ since there's no INFO/DP or INFO/MQ)
expr="GT[0]=\"het\" && GT[1]=\"het\" && FMT/DP[0]>=${MIN_DP} && FMT/DP[1]>=${MIN_DP} && FMT/GQ[0]>=${MIN_GQ} && FMT/GQ[1]>=${MIN_GQ}"

bcftools view -i "$expr" -Oz -o "$out/merged/common.het.filtered.vcf.gz" \
  "$out/merged/merged.giraffe_bwa.vcf.gz"
bcftools index -f "$out/merged/common.het.filtered.vcf.gz"

echo "Sample order: $(bcftools query -l "$out/merged/common.het.filtered.vcf.gz" | paste -sd, -)"
echo "Common het sites: $(bcftools index -n "$out/merged/common.het.filtered.vcf.gz")"
echo "Merged VCF:     $out/merged/merged.giraffe_bwa.vcf.gz"
echo "Common het VCF: $out/merged/common.het.filtered.vcf.gz"
echo "Plot with: python plot_allele_bias.py $out/merged/common.het.filtered.vcf.gz ${sample}.allele_balance.png"
