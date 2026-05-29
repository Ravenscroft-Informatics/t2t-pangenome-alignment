#!/bin/bash -l

#SBATCH --job-name=svfotate
#SBATCH --account=pawsey0933
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=80G
#SBATCH --nodes=1
#SBATCH --time=0:30:00
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out
#SBATCH --export=ALL

# --- Safety checks ---
if [ -z "$1" ]; then
    echo "Error: No chromosome provided"
    exit 1
fi

module load bcftools/1.15--haf5b3da_0

conda activate /software/projects/pawsey0933/cfolland/miniforge3/envs/svafotate-env/

cd /scratch/pawsey0933/cfolland/t2t/batch2/chm13/genmod/sv/

CHR=$1
SUBSET=/scratch/pawsey0933/cfolland/t2t/batch2/chm13/genmod/sv/batch2_linear_chm13_dysgu_pass_VEP_genmod.chr${CHR}.vcf.gz
OUT=/scratch/pawsey0933/cfolland/t2t/batch2/chm13/genmod/sv/batch2_linear_chm13_dysgu_pass_VEP_genmod.chr${CHR}.svafotate.vcf.gz
BED=/scratch/pawsey0933/cfolland/t2t/batch2/annotation/SVAFotate_core_SV_popAFs.CHM13.v4.2.bed.gz

svafotate annotate \
	-f 0.5 \
	-l 1000000 \
	-a mis best \
	-u 0.0001 \
	-v $SUBSET \
	-o $OUT \
	-b $BED

mkdir -p unique
mv uniques.bed ./unique/unique_$CHR.bed


#submit another job for the following chromosome
next_chr=$((CHR + 1))

if [ "$next_chr" -le 22 ]; then
    sbatch --dependency=afterok:$SLURM_JOB_ID \
        /software/projects/pawsey0933/t2t/sv_anno/svfotate.sh $next_chr
fi