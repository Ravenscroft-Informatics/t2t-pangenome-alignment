#!/bin/bash -l

#SBATCH --job-name=rename_vcf
#SBATCH --account=pawsey0933
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out
#SBATCH --export=ALL

module load bcftools/1.15--haf5b3da_0

cd /software/projects/pawsey0933/t2t/scripts

# 0) Inputs
VCF=/scratch/pawsey0933/cfolland/t2t/batch2/grch38/dv/linear-grch38_glnexus_dv_VEP.ann.vcf.gz
MAP=/software/projects/pawsey0933/t2t/scripts/vcf_mapping.tsv          # 2 columns: original<TAB>target
OUTPUT=/scratch/pawsey0933/cfolland/t2t/batch2/grch38/dv/linear-grch38_glnexus_dv_VEP.ann.renamed.vcf.gz

# 1) Extract sample names in order
bcftools query -l "$VCF" > samples.txt

# 2) Build the new ordered sample list:
#    If a sample is in the mapping TSV, replace it; otherwise keep original.
awk 'BEGIN{FS=OFS="\t"}
     NR==FNR {map[$1]=$2; next}
     { if ($1 in map) print map[$1]; else print $1 }' "$MAP" samples.txt > new_samples.txt

# 3) Safety: check for duplicates in the *new* names (should be zero)
if [ "$(sort new_samples.txt | uniq -d | wc -l)" -ne 0 ]; then
  echo "ERROR: duplicate sample names detected in new_samples.txt"; 
  echo "Duplicates:"
  sort new_samples.txt | uniq -d
  exit 1
fi

# 4) Reheader and reindex
bcftools reheader -s new_samples.txt -o $OUTPUT "$VCF"
bcftools index -t $OUTPUT
``