#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH -c 10
#SBATCH --job-name=genmod_split_vcf
#SBATCH --partition=work
#SBATCH --account=pawsey0933
#SBATCH --mem=30gb
#SBATCH --time=0-03:00:00
#SBATCH --export=NONE
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out

module load bcftools/1.15--haf5b3da_0

set -ueo pipefail

#chm13
#input=/scratch/pawsey0933/cfolland/t2t/batch2/chm13/batch2_linear_chm13_dv_dysgu_VEP.sorted.noCSQ.vcf.gz
#output=/scratch/pawsey0933/cfolland/t2t/batch2/chm13/genmod

#grch38
input=/scratch/pawsey0933/cfolland/t2t/batch2/grch38/batch2_linear_grch38_dv_dysgu_VEP.sorted.noCSQ.vcf.gz
output=/scratch/pawsey0933/cfolland/t2t/batch2/grch38/genmod

base=`basename $input .vcf.gz`

for chrom in {1..22} X Y M ; do
    bcftools view --threads 20 -r chr${chrom} -Oz -o $output/${base}_chr${chrom}.vcf.gz $input
    bcftools index -f "$output/${base}_chr${chrom}.vcf.gz"
done