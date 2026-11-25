#!/bin/bash

#SBATCH --job-name=merge_regions_nmd
#SBATCH --account=pawsey0933
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out
#SBATCH --export=ALL

ref=chm13
mode=hprc
cd /scratch/pawsey0933/cfolland/mosdepth/results/$mode

module load samtools/1.15--h3843a85_0
source /software/projects/pawsey0933/cfolland/miniconda3/etc/profile.d/mamba.sh

python3 collect_coverage_mosdepth_regions.py *.$ref.$mode.regions.bed.gz -o mosdepth_nmd.$ref.$mode.merged.tsv
