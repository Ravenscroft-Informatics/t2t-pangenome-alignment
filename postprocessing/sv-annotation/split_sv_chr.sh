#!/bin/bash -l

#SBATCH --job-name=subset_vcf
#SBATCH --account=pawsey0933
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --time=0:15:00
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out
#SBATCH --export=ALL

module load bcftools/1.15--haf5b3da_0

#VCF=/scratch/pawsey0933/cfolland/t2t/batch2/grch38/dysgu/batch2_linear_hg38_dysgu_merge_sorted.vcf.gz
#SUBSET=/scratch/pawsey0933/cfolland/t2t/batch2/grch38/dysgu/batch2_linear_hg38_dysgu_merge_sorted.chr4.vcf.gz

CHR=$1
VCF=/scratch/pawsey0933/cfolland/t2t/batch2/chm13/genmod/batch2_linear_chm13_dysgu_pass_VEP_genmod.vcf.gz
SUBSET=/scratch/pawsey0933/cfolland/t2t/batch2/chm13/genmod/sv/batch2_linear_chm13_dysgu_pass_VEP_genmod.${CHR}.vcf.gz

bcftools view --regions ${CHR} -Oz -o $SUBSET $VCF 
tabix -p vcf $SUBSET