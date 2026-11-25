#!/bin/bash -l
#SBATCH --job-name=bcftools_stats
#SBATCH --account=pawsey0933
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out
#SBATCH --export=NONE

#use this script to generate a bcftools stats summary for a cohort vcf

set -euo pipefail

# Load modules
module load bcftools/1.15--haf5b3da_0

# config
approach=pangenome
input="/scratch/pawsey0933/cfolland/benchmark/vcfs/deepvariant/${approach}/cohort/fixed"
output_dir="/scratch/pawsey0933/cfolland/benchmark/output/${approach}/bcftools_stats/cohort/fixed"

mkdir -p "$output_dir"

process_one() {
  local vcf="$1"
  local filename
  filename=$(basename "$vcf" .vcf.gz)
  local output="${output_dir}/${filename}.bcftools_stats.txt"

  if [[ -s "$output" ]]; then
    echo "Exists, skipping: $output"
  else
    echo "Calculating bcftools stats for $(basename "$vcf")"
    bcftools stats "$vcf" > "$output"
  fi
}

if [[ -d "$input" ]]; then
  shopt -s nullglob
  files=( "$input"/*.vcf.gz )
  if (( ${#files[@]} == 0 )); then
    echo "No .vcf.gz files found in: $input"
    exit 1
  fi
  for vcf in "${files[@]}"; do
    process_one "$vcf"
  done
elif [[ -f "$input" ]]; then
  # If 'input' is a single VCF file
  process_one "$input"
else
  echo "Input not found: $input"
  exit 1
fi
