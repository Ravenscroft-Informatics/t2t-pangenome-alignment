#!/bin/bash -l

#SBATCH --job-name=svfotate
#SBATCH --account=pawsey0933
# SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --time=0:30:00
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out
#SBATCH --export=ALL

module load bcftools/1.15--haf5b3da_0

conda activate /software/projects/pawsey0933/cfolland/miniforge3/envs/svafotate-env/

cd /scratch/pawsey0933/cfolland/t2t/batch2/chm13/genmod/sv/

#VCF=/scratch/pawsey0933/cfolland/t2t/batch2/grch38/dysgu/batch2_linear_hg38_dysgu_merge_sorted.vcf.gz
#SUBSET=/scratch/pawsey0933/cfolland/t2t/batch2/grch38/dysgu/batch2_linear_hg38_dysgu_merge_sorted.chr4.vcf.gz
#OUT=/scratch/pawsey0933/cfolland/t2t/batch2/grch38/dysgu/batch2_linear_hg38_dysgu_merge_sorted_chr4_svafotate.vcf.gz
#BED=/scratch/pawsey0933/cfolland/t2t/batch2/annotation/SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz

CHR=$1
SUBSET=/scratch/pawsey0933/cfolland/t2t/batch2/chm13/genmod/sv/batch2_linear_chm13_dysgu_pass_VEP_genmod.${CHR}.vcf.gz
OUT=/scratch/pawsey0933/cfolland/t2t/batch2/chm13/genmod/sv/batch2_linear_chm13_dysgu_pass_VEP_genmod.${CHR}.svafotate.vcf.gz
BED=/scratch/pawsey0933/cfolland/t2t/batch2/annotation/SVAFotate_core_SV_popAFs.CHM13.v4.2.bed.gz

svafotate annotate \
	-f 0.5 \
	-l 1000000 \
	-a mis best \
	-u 0.0001 \
	-v $SUBSET \
	-o $OUT \
	-b $BED