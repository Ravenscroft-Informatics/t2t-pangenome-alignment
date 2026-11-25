#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH -c 1
#SBATCH --job-name=mendel
#SBATCH --partition=work
#SBATCH --account=pawsey0933
#SBATCH --mem=16G
#SBATCH --time=1:10:00
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out
#SBATCH --export=ALL

GENOME=chm13
prefix=${GENOME}_linear
sif=/software/projects/pawsey0933/benchmarking/mendel/gatk_latest.sif
#ref=/software/projects/pawsey0933/sv/references/hg38_masked/Homo_sapiens_assembly38_masked.fasta
ref=/software/projects/pawsey0933/pangenome/refs/chm13v2.0.maskedY.rCRS.EBV.fasta
#input=/scratch/pawsey0933/cfolland/vep/annotation/glnexus/merged-dv/dv_glnexus_VEP.ann.vcf.gz
input=/scratch/pawsey0933/cfolland/benchmark/vcfs/deepvariant/linear/cohort/T2T_dv_glnexus_VEP.ann.vcf.gz
#input=/scratch/pawsey0933/cfolland/benchmark/vcfs/deepvariant/${approach}/cohort/${prefix}_dv_glnexus_VEP.gff.exons.vcf.gz 
output=/scratch/pawsey0933/cfolland/benchmark/mendel_viol/${approach}/$prefix.MVs.byFamily.table
ped=/software/projects/pawsey0933/pangenome/genmod/trios_trios.ped

module load singularity/4.1.0-slurm

singularity exec $sif \
    gatk \
    --java-options -Xmx16G \
    VariantEval \
    -R $ref \
    -O $output \
    --eval $input \
    -no-ev -no-st --lenient \
    -ST Family \
    -EV MendelianViolationEvaluator \
    -ped $ped -pedValidationType SILENT
