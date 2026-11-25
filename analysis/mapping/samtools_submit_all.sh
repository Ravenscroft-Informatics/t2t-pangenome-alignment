#!/bin/bash

GENOME=grch38 #chm13
INPUT_DIR=/scratch/pawsey0933/T2T/hg38_realignment/sarek_bwamem2_run2/preprocessing/mapped/
#/scratch/pawsey0933/cfolland/benchmark/bams/$GENOME/

ls $INPUT_DIR | grep 'bam' | grep -v 'md5\|crai\|bai' | while read -r file ; do
    echo "Submitting job for $file"
    sbatch samtools_stats.sh "$file" "$INPUT_DIR"
done
