#!/bin/bash -l

#SBATCH --job-name=ann_small
#SBATCH --account=pawsey0933
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out
#SBATCH --export=ALL

module load bcftools/1.15--haf5b3da_0

vcf_dir=/scratch/pawsey0933/cfolland/t2t/batch2/chm13/genmod
annotation_dir=/scratch/pawsey0933/cfolland/t2t/batch2/annotation

input=batch2_linear_chm13_dv_VEP_genmod.vcf.gz
temp=batch2_linear_chm13_dv_VEP.genmod.tmp.vcf.gz
output=batch2_linear_chm13_dv_VEP.genmod.af.vcf.gz

hprc=hprc-v2.0-mc-chm13.wave.biallelic.sorted.vcf
hgsvc=hgsvc3-hprc-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.biallelic.af.sorted.vcf

# make sure annotation files are bgzipped
#gunzip ${annotation_dir}/${hprc}.gz
#bgzip ${annotation_dir}/${hprc}
#tabix -p ${annotation_dir}/${hprc}.gz

#gunzip ${annotation_dir}/${hgsvc}.gz
#bgzip ${annotation_dir}/${hgsvc}
#tabix -p ${annotation_dir}/${hgsvc}.gz

bcftools annotate \
  -a ${annotation_dir}/${hprc}.gz \
  -c INFO/AF_HPRCwave:=INFO/AF \
  -Oz \
  -o ${vcf_dir}/${temp} \
  ${vcf_dir}/${input}

tabix -p vcf ${vcf_dir}/${temp}

bcftools annotate \
  -a ${annotation_dir}/${hgsvc}.gz \
  -c INFO/AF_HGSVC:=INFO/AF \
  -O z \
  -o ${vcf_dir}/${output} \
  ${vcf_dir}/${temp}

tabix -p vcf ${vcf_dir}/${output}