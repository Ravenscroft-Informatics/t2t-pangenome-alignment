#!/bin/bash -l

#SBATCH --job-name=reorder
#SBATCH --account=pawsey0933
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out
#SBATCH --export=ALL

module load bcftools/1.15--haf5b3da_0

cd /software/projects/pawsey0933/t2t/scripts

bcftools view -S samples_order.txt -O z -o /scratch/pawsey0933/cfolland/t2t/batch2/grch38/dv/linear-grch38_glnexus_dv_VEP.ann.sorted.vcf.gz /scratch/pawsey0933/cfolland/t2t/batch2/grch38/dv/linear-grch38_glnexus_dv_VEP.ann.renamed.vcf.gz
bcftools index /scratch/pawsey0933/cfolland/t2t/batch2/grch38/dv/linear-grch38_glnexus_dv_VEP.ann.sorted.vcf.gz