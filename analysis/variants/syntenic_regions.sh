#!/bin/bash -l

#SBATCH --job-name=synteny
#SBATCH --account=pawsey0933
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out
#SBATCH --export=ALL

#load modules
module load bedtools/2.30.0--h468198e_3
module load bcftools/1.15--haf5b3da_0

cd /software/projects/pawsey0933/benchmarking/mask

GENOME=chm13
prefix=${GENOME}_linear
cohort=/scratch/pawsey0933/cfolland/benchmark/vcfs/deepvariant/linear/cohort/T2T_dv_glnexus_VEP.ann.vcf.gz
#cohort=/scratch/pawsey0933/cfolland/vep/annotation/glnexus/merged-dv/dv_glnexus_VEP.ann.vcf.gz
vcf_output=/scratch/pawsey0933/cfolland/benchmark/vcfs/deepvariant/linear/cohort/
mask=chm13v2-unique_to_hg38.bed

# Output files
syntenic_vcf=${vcf_output}/${prefix}_dv_glnexus_syntenic.vcf.gz

# Preserve VCF header and intersect
bcftools view -T ^$mask ${cohort} -Oz -o ${syntenic_vcf}

# Compress with bgzip and index
tabix -p vcf ${syntenic_vcf}
