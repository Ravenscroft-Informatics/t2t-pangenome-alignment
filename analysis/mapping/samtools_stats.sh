#!/bin/bash -l

#SBATCH --job-name=samtools_stats
#SBATCH --account=pawsey0933
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --time=1:50:00
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out
#SBATCH --export=ALL

#use this script to compute mapping statistics on bam files

set -euo pipefail
module load samtools/1.15--h3843a85_0

GENOME=grch38

REF=/software/projects/pawsey0933/sv/references/hg38_masked/Homo_sapiens_assembly38_masked.fasta #grch38-linear
REF=/software/projects/pawsey0933/pangenome/refs/chm13v2.0.maskedY.rCRS.EBV.fasta #chm13-linear
#REF=/scratch/pawsey0933/cfolland/benchmark/refs/hprc-v1.1-mc-$GENOME.ref.fa #pangenome
INPUT_DIR="$2" #input directory containing bam files
INPUT="$1"
OUTPUT_DIR=/scratch/pawsey0933/cfolland/benchmark/output/samtools_stats/linear

mkdir -p "$OUTPUT_DIR"

basename=$(basename "$INPUT")
basename=${basename%%.*}
OUTPUT_FILE="$OUTPUT_DIR/${basename}.${GENOME}.samtools.stats.txt"

if [ -f "$OUTPUT_FILE" ]; then
    echo "✅ Output already exists for $INPUT, skipping: $OUTPUT_FILE"
else
    echo "📊 Running samtools stats for $INPUT..."
    samtools view -b -T "$REF" "$INPUT_DIR/$basename/$INPUT" | \
    samtools stats -@ 16 - > "$OUTPUT_FILE"
    #samtools stats -r "$REF" -@ 16 "$INPUT_DIR/$basename/$INPUT" > "$OUTPUT_FILE"
    echo "✅ Done: $OUTPUT_FILE"
fi
