#!/bin/bash

#SBATCH --job-name=merge
#SBATCH --account=pawsey0933
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=90G
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out
#SBATCH --export=ALL

#shared
ref=chm13
method=linear

#single gene/interval
gene='DMXL1'
interval="chr19:1040106-1065572"
single_py=/software/projects/pawsey0933/benchmarking/nmd_genes/scripts/collect_coverage_perbase.py

#mutliple genes/intervals
bed_py=/software/projects/pawsey0933/benchmarking/nmd_genes/scripts/collect_coverage_perbase_bed.py
bed=/software/projects/pawsey0933/benchmarking/nmd_genes/gene_list/nmd_gene_list_chm13.bed
merged=batch2.perbase_mosdepth_${method}_${ref}.tsv

cd /scratch/pawsey0933/cfolland/mosdepth/results/$method/$ref

module load samtools/1.15--h3843a85_0
source /software/projects/pawsey0933/cfolland/miniconda3/etc/profile.d/mamba.sh

#single gene/interval
python3 $single_py $gene $interval *.$ref.$method.per-base.bed.gz --index --threads 8 --end-inclusive -o $gene.perbase_mosdepth_${method}_${ref}.tsv

#bed file witrh multiple genes/intervals
python3 $bed_py $bed *.$ref.$method.per-base.bed.gz --index --threads 8 --bed-end-inclusive -o $merged
