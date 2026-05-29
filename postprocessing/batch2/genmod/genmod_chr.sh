#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH -c 1
#SBATCH --job-name=genmod_chr
#SBATCH --partition=long
#SBATCH --account=pawsey0933
#SBATCH --mem=200G
#SBATCH --time=6:00:00
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out
#SBATCH --export=ALL

source /scratch/pawsey0933/cfolland/anaconda3/etc/profile.d/mamba.sh
mamba activate /software/projects/pawsey0933/cfolland/miniforge3/envs/genmod

chr=$1
genome=chm13
dir=/scratch/pawsey0933/cfolland/t2t/batch2/$genome/genmod
ped=/software/projects/pawsey0933/t2t/scripts/batch2.new.ped

genmod models --vep -p 2 -f $ped -o $dir/batch2_linear_${genome}_dv_dysgu_VEP.sorted.noCSQ.genmod_${chr}.genmod.vcf $dir/batch2_linear_${genome}_dv_dysgu_VEP.sorted.noCSQ_${chr}.vcf.gz