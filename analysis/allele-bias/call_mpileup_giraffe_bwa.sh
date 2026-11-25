#!/bin/bash
set -euo pipefail

sample="D09-468"

ref="/Volumes/PERKINS-LL-001/Sequencing/wgs/secondary/T2T_realignment/references/Homo_sapiens_assembly38_masked.fasta"

giraffe_bam="/Volumes/PERKINS-LL-001/Sequencing/wgs/secondary/Pangenome_realignment/batch_1/vg_giraffe/bams/${sample}.hprc-v1.1-mc-grch38.surj_realn.bam"
bwa_cram="/Volumes/PERKINS-LL-001/Sequencing/wgs/secondary/T2T_realignment/batch_1/grch38/preprocessing/mapped/${sample}/${sample}.sorted.cram"

out="/Users/00104561/Library/CloudStorage/OneDrive-UWA/Research/Projects/AIM1_DATA-REANALYSIS/benchmark/batch1/pangenome_linear/mpileup_calls/${sample}"
threads=16

mkdir -p "$out"/{giraffe,bwa,merged}

# Reference index required for CRAM & mpileup
[[ -f "${ref}.fai" ]] || samtools faidx "$ref"

# Ensure alignment indexes exist
[[ -f "${giraffe_bam}.bai" ]] || samtools index -@ "$threads" "$giraffe_bam"
[[ -f "${bwa_cram}.crai"  ]] || samtools index -@ "$threads" "$bwa_cram"

# mpileup options (paper-like)
# -E enables BAQ adjustment
# -a adds tags we may filter/plot with
# -d avoids depth truncation surprises
MPILEUP_OPTS=(-f "$ref" -E -a DP -a SP -a ADF -a ADR -a AD -O u -d 100000)

call_one_mapper () {
  local mapper="$1"     # giraffe or bwa
  local inbam="$2"      # bam/cram path
  local workdir="$3"

  echo "== Calling whole-genome for ${mapper} =="

  bcftools mpileup "${MPILEUP_OPTS[@]}" "$inbam" \
    | bcftools call -m -v -O z -o "${workdir}/${mapper}.unsorted.vcf.gz"

  # Sort + index to be safe for merge/queries
  bcftools sort "${workdir}/${mapper}.unsorted.vcf.gz" -O z -o "${workdir}/${mapper}.vcf.gz"
  bcftools index -f "${workdir}/${mapper}.vcf.gz"
  rm -f "${workdir}/${mapper}.unsorted.vcf.gz"
}

call_one_mapper "giraffe" "$giraffe_bam" "$out/giraffe"
call_one_mapper "bwa"     "$bwa_cram"   "$out/bwa"

# Rename sample names inside each VCF to exactly "giraffe" and "bwa"
giraffe_sm=$(bcftools query -l "$out/giraffe/giraffe.vcf.gz")
bwa_sm=$(bcftools query -l "$out/bwa/bwa.vcf.gz")

printf "%s\tgiraffe\n" "$giraffe_sm" > "$out/merged/rename.giraffe.txt"
printf "%s\tbwa\n"     "$bwa_sm"     > "$out/merged/rename.bwa.txt"
bcftools reheader -s "$out/merged/rename.giraffe.txt" -o "$out/merged/giraffe.renamed.vcf.gz" "$out/giraffe/giraffe.vcf.gz"
bcftools reheader -s "$out/merged/rename.bwa.txt"     -o "$out/merged/bwa.renamed.vcf.gz"     "$out/bwa/bwa.vcf.gz"
bcftools index -f "$out/merged/giraffe.renamed.vcf.gz"
bcftools index -f "$out/merged/bwa.renamed.vcf.gz"

# Merge into a 2-sample VCF
bcftools merge -O z -o "$out/merged/merged.giraffe_bwa.vcf.gz" \
  "$out/merged/giraffe.renamed.vcf.gz" \
  "$out/merged/bwa.renamed.vcf.gz"
bcftools index -f "$out/merged/merged.giraffe_bwa.vcf.gz"

# Filter to "high-confidence het in BOTH"
# INFO/DP>=25
# INFO/MQ>=40
expr='GT[0]="het" && GT[1]="het" && INFO/DP>=25 && INFO/MQ>=40'

bcftools view -i "$expr" -O z -o "$out/merged/common.het.filtered.vcf.gz" \
  "$out/merged/merged.giraffe_bwa.vcf.gz"
bcftools index -f "$out/merged/common.het.filtered.vcf.gz"

echo "DONE"
echo "Merged VCF:     $out/merged/merged.giraffe_bwa.vcf.gz"
echo "Common het VCF: $out/merged/common.het.filtered.vcf.gz"
