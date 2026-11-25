#!/bin/bash

input_dir=/scratch/pawsey0933/cfolland/t2t/crams/chm13

# Find BAM/CRAM files, excluding index/md5, and submit one sbatch per file
find "$input_dir" -type f \( -name "*.cram" -o -name "*.bam" \) \
    ! -name "*.bai" ! -name "*.crai" ! -name "*.md5" \
| sort | while read -r file; do
    echo "Submitting job for $file"
    sbatch --export=ALL mosdepth.sh "$file"
    sleep 0.2  # slight delay to avoid overloading Slurm controller
done
