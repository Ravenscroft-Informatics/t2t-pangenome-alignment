#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH -c 12                          
#SBATCH --job-name=concat_genmod
#SBATCH --partition=work
#SBATCH --account=pawsey0933
#SBATCH --mem=16G
#SBATCH --time=10:00:00
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out
#SBATCH --export=ALL

# Load modules
module load bcftools/1.15--haf5b3da_0

# Define variables
ref=grch38
outdir=/scratch/pawsey0933/cfolland/t2t/batch2/$ref/genmod

# Concatenate, sort and index
bcftools concat -a -Oz --threads 24 -o $outdir/batch2_linear_${ref}_dv_dysgu_VEP.genmod.vcf.gz $outdir/batch2_linear_${ref}_dv_dysgu_VEP.sorted.noCSQ.genmod_*.genmod.vcf.gz
bcftools sort -T $outdir -o $outdir/batch2_linear_${ref}_dv_dysgu_VEP.genmod.sorted.vcf.gz $outdir/batch2_linear_${ref}_dv_dysgu_VEP.genmod.vcf.gz
bcftools index -t $outdir/batch2_linear_${ref}_dv_dysgu_VEP.genmod.sorted.vcf.gz