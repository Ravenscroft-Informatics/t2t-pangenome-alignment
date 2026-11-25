#!/bin/bash

#SBATCH --job-name=merge
#SBATCH --account=pawsey0933
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --nodes=1
#SBATCH --time=6:00:00
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out
#SBATCH --export=ALL

#load modules
module load samtools/1.15--h3843a85_0
source /software/projects/pawsey0933/cfolland/miniconda3/etc/profile.d/mamba.sh

#define variables
ref=hg38
method=linear
agg_py=/software/projects/pawsey0933/benchmarking/nmd_genes/perbase_agg_summary.py
merged=all.perbase_mosdepth_${method}_${ref}.tsv
sorted=all.perbase_mosdepth_${method}_${ref}.sorted.tsv
final=all.perbase_mosdepth.summary.${method}_${ref}.tsv
tmp=/scratch/pawsey0933/cfolland/sort_tmp_${method}_${ref}

mkdir -p "$tmp"
cd /scratch/pawsey0933/cfolland/mosdepth/results/$method/$ref

# keep header, sort body only
{ head -n1 "$merged"; \
  tail -n +2 "$merged" | \
  LC_ALL=C sort \
    -T "$tmp" \
    --parallel=4 \
    -S 3G \
    --batch-size=32 \
    -k1,1 -k5,5 -k6,6n; } > "$sorted"

#calculate the mean, median, range and sd for each gene/sample using the file generated above
python3 $agg_py -i $sorted -o $final --stdev sample

