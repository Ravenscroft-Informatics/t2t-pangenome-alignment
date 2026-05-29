#!/bin/bash -l

#SBATCH --job-name=liftover
#SBATCH --account=pawsey0933
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out
#SBATCH --export=ALL

# prepare environment
conda activate /software/projects/pawsey0933/cfolland/miniforge3/envs/annotsv/
module load bcftools/1.15--haf5b3da_0
module load bedtools/2.30.0--h468198e_3

# paths
ChainDir=/software/projects/pawsey0933/t2t/sv_anno
CHAIN=$ChainDir/hg38ToGCA_009914755.4.over.chain
LIFTOVER_SCRIPT=/software/projects/pawsey0933/t2t/sv_anno/liftover_script.sh
IN_BED=/scratch/pawsey0933/cfolland/t2t/batch2/annotation/SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz
OUT_BED=/scratch/pawsey0933/cfolland/t2t/batch2/annotation/SVAFotate_core_SV_popAFs.CHM13.v4.1.bed

# run liftover from TMPDIR (important for its own temp files)
echo "Running liftover:"
echo "  Input:  $IN_BED"
echo "  Output: $OUT_BED"
echo "  Chain:  $CHAIN"

"$LIFTOVER_SCRIPT" "$IN_BED" "$OUT_BED" "$CHAIN"

echo "Done: $OUT"