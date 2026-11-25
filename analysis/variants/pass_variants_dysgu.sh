#!/bin/bash -l

#SBATCH --job-name=pass_filter_dysgu
#SBATCH --account=pawsey0933
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out
#SBATCH --export=NONE

# Load modules
module load bcftools/1.15--haf5b3da_0

# Define variables
INPUT_DIR=/scratch/pawsey0933/cfolland/benchmark/vcfs/dysgu/linear/
OUTPUT_DIR=/software/projects/pawsey0933/benchmarking/bcftools
OUTFILE=${OUTPUT_DIR}/pass_dysgu_linear_cohort.tsv

# Initialize output file
echo -e "Sample\tPASS_Variants" > "$OUTFILE"

# Loop through VCF files
ls "$INPUT_DIR" | grep -E 'vcf$|vcf.gz$' | grep -v 'tbi' | while read -r file; do
    basename=${file%.vcf*}
    VCF="${INPUT_DIR}/${file}"

    # Count PASS variants
    count=$(bcftools view -f PASS "$VCF" | bcftools view -H | wc -l)

    # Write to output file
    echo -e "${basename}\t${count}" >> "$OUTFILE"
done
