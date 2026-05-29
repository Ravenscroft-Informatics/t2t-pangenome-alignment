#!/bin/bash -l

#SBATCH --job-name=dv_dysgu
#SBATCH --account=pawsey0933
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --time=8:00:00
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out
#SBATCH --export=ALL

#load modules
module load vcftools/0.1.16--pl5321hd03093a_7
module load bcftools/1.15--haf5b3da_0

input=/scratch/pawsey0933/cfolland/t2t/batch2/chm13/genmod/batch2_linear_chm13_dv_dysgu_VEP.genmod.sorted.vcf.gz
output=/scratch/pawsey0933/cfolland/t2t/batch2/chm13/genmod/
GENOME=chm13
PREFIX=batch2_linear_$GENOME

#Filter based on qc 
bcftools filter -e 'INFO/SVMETHOD == "DYSGUv1.8.7"' -Oz -o $output/${PREFIX}_dv_VEP_genmod.vcf.gz ${input}

#tabix
tabix -p vcf $output/${PREFIX}_dv_VEP_genmod.vcf.gz