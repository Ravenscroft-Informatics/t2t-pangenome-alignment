#!/bin/bash -l

#SBATCH --job-name=mask
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

GENOME=grch38
prefix=${GENOME}_linear
cohort=/scratch/pawsey0933/cfolland/benchmark/vcfs/deepvariant/linear/cohort/T2T_dv_glnexus_VEP.ann.vcf.gz
#cohort=/scratch/pawsey0933/cfolland/vep/annotation/glnexus/merged-dv/dv_glnexus_VEP.ann.vcf.gz
vcf_output=/scratch/pawsey0933/cfolland/benchmark/vcfs/deepvariant/linear/cohort/
mask=/software/projects/pawsey0933/benchmarking/mask/hg38.combined_mask.bed
#mask=/software/projects/pawsey0933/benchmarking/mask/combined_mask.bed

# Output files
masked_vcf=${vcf_output}/${prefix}_dv_glnexus_mask.vcf
masked_vcf_gz=${masked_vcf}.gz

# Preserve VCF header and intersect
bcftools view -h ${cohort} > ${masked_vcf}
bedtools intersect -a ${cohort} -b ${mask} >> ${masked_vcf}

# Compress with bgzip and index
bgzip -f ${masked_vcf}  # Overwrites with ${masked_vcf}.gz
tabix -p vcf ${masked_vcf_gz}

echo "Masked VCF written to ${masked_vcf_gz} and indexed"
