#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH -c 4                          
#SBATCH --job-name=bgzip_genmod
#SBATCH --partition=work
#SBATCH --account=pawsey0933
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out
#SBATCH --export=ALL

module load bcftools/1.15--haf5b3da_0

chr=$1
genome=chm13
dir=/scratch/pawsey0933/cfolland/t2t/batch2/$genome/genmod

bgzip -ki -@ 8 $dir/batch2_linear_${genome}_dv_dysgu_VEP.sorted.noCSQ.genmod_${chr}.genmod.vcf