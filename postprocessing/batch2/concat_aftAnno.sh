#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH -c 6                          
#SBATCH --job-name=concat_t2t
#SBATCH --partition=work
#SBATCH --account=pawsey0933
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out
#SBATCH --export=ALL

# Load modules
module load bcftools/1.15--haf5b3da_0

# Define variables
ref=chm13
dv=/scratch/pawsey0933/cfolland/t2t/batch2/${ref}/genmod/sv/batch2_linear_${ref}_dysgu_pass_VEP_genmod.svafotate.sorted.vcf.gz
dysgu=/scratch/pawsey0933/cfolland/t2t/batch2/${ref}/genmod/batch2_linear_${ref}_dv_VEP.genmod.af.vcf.gz
outdir=/scratch/pawsey0933/cfolland/t2t/batch2/${ref}/

# Concatenate, sort and index
bcftools concat -a -Oz --threads 12 -o $outdir/batch2_linear_${ref}_dv_dysgu.annot.vcf.gz $dv $dysgu
bcftools sort -T $outdir -o $outdir/batch2_linear_${ref}_dv_dysgu.annot.sorted.vcf.gz $outdir/batch2_linear_${ref}_dv_dysgu.annot.vcf.gz
bcftools index -t $outdir/batch2_linear_${ref}_dv_dysgu.annot.sorted.vcf.gz