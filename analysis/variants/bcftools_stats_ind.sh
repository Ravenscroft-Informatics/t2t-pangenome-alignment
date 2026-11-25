#!/bin/bash -l

#SBATCH --job-name=bcftools_stats
#SBATCH --account=pawsey0933
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --time=3:00:00
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out
#SBATCH --export=NONE

# use this script to calculate bcftools stats on individual sample vcfs 

#load modules
module load bcftools/1.15--haf5b3da_0

#Define variables
GENOME=hg38
REF=/software/projects/pawsey0933/sv/references/hg38_masked/Homo_sapiens_assembly38_masked.fasta
#REF=/scratch/pawsey0933/cfolland/benchmark/refs/hprc-v1.1-mc-$GENOME.ref.fa 
#REF=/scratch/pawsey0933/pangenome/refs/chm13v2.0.maskedY.rCRS.EBV.fasta
INPUT_DIR=/scratch/pawsey0933/cfolland/benchmark/vcfs/deepvariant/linear/individual/hg38 #path to the individual vcf files (not g.vcf)
OUTPUT_DIR=/scratch/pawsey0933/cfolland/benchmark/output/linear/bcftools_stats/individual

ls "$INPUT_DIR" | grep '*.sorted.cram.vcf.gz' | grep -v 'tbi' | while read -r file; do

    basename=${file%%.*}
    OUTPUT="$OUTPUT_DIR/${basename}.${GENOME}.deepvariant.bcftools_stats.txt"

    if [ -f "$OUTPUT" ]; then
        echo "Output already exists for $file, skipping: $OUTPUT"
    else
        echo "Calculating bcftools stats for $basename"
        VCF="${basename}.sorted.cram.vcf.gz"
        bcftools stats "$INPUT_DIR/$VCF" > "$OUTPUT"
    fi

done
