#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH -c 6                          
#SBATCH --job-name=concat_t2t
#SBATCH --partition=work
#SBATCH --account=pawsey0933
#SBATCH --mem=16G
#SBATCH --time=3:00:00
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out
#SBATCH --export=ALL

# Load modules
module load bcftools/1.15--haf5b3da_0

# Define variables
ref=grch38
dv=/scratch/pawsey0933/cfolland/t2t/batch2/grch38/dv/linear-${ref}_glnexus_dv_VEP.ann.sorted.vcf.gz
dysgu=/scratch/pawsey0933/cfolland/t2t/batch2/grch38/dysgu/batch2_linear_hg38_dysgu_merge_sorted_VEP.ann.vcf.gz
outdir=/scratch/pawsey0933/cfolland/t2t/batch2/${ref}/

# Concatenate, sort and index
bcftools concat -a -Oz --threads 12 -o $outdir/batch2_linear_${ref}_dv_dysgu_VEP.vcf.gz $dv $dysgu
bcftools index -t $outdir/batch2_linear_${ref}_dv_dysgu_VEP_concat.vcf.gz
bcftools sort -T $outdir -o $outdir/batch2_linear_${ref}_dv_dysgu_VEP.sorted.vcf.gz $outdir/batch2_linear_${ref}_dv_dysgu_VEP.vcf.gz
bcftools index -t $outdir/batch2_linear_${ref}_dv_dysgu_VEP_sorted.vcf.gz
bcftools view -e 'CSQ[*]="."' $outdir/batch2_linear_${ref}_dv_dysgu_VEP.sorted.vcf.gz -Oz -o $outdir/batch2_linear_${ref}_dv_dysgu_VEP.sorted.noCSQ.vcf.gz
bcftools index -t $outdir/batch2_linear_${ref}_dv_dysgu_VEP.sorted.noCSQ.vcf.gz