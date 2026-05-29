#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH -c 2
#SBATCH --job-name=sarek_vep
#SBATCH --partition=work
#SBATCH --account=pawsey0933
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out
#SBATCH --export=ALL

export NXF_WORK=/scratch/pawsey0933/cfolland/vep
export SINGULARITY_CACHEDIR=/scratch/pawsey0933/cfolland/vep/.nextflow_singularity
export NXF_SINGULARITY_CACHEDIR=/scratch/pawsey0933/cfolland/vep/.nextflow_singularity
export NXF_HOME=/scratch/pawsey0933/cfolland/vep/.nextflow
export SINGULARITY_TMPDIR=/scratch/pawsey0933/cfolland/vep/tempdir
export NXF_OPTS='-Xms1g -Xmx19g'

module load singularity/4.1.0-slurm
module load nextflow/23.10.0 #nextflow/24.04.3

nextflow run nf-core/sarek \
	-r 3.4.0 \
	-profile singularity,pawsey_setonix \
	-config /software/projects/pawsey0933/gmonahan/T2T/bam2fastq/pawsey_setonix.config \
        --input /software/projects/pawsey0933/t2t/hg38/grch38-dysgu.csv \
	--step annotate \
        --genome null \
        --igenomes_ignore \
	--tools vep \
	--vep_custom_args "--everything --filter_common --per_gene --total_length --offline --format vcf --max_sv_size 100000000" \
	--vep_cache /scratch/pawsey0933/cfolland/vep/cache \
	--vep_species homo_sapiens \
	--vep_genome GRCh38 \
	--vep_cache_version 114 \
	--vep_include_fasta \
	--fasta /software/projects/pawsey0933/sv/references/hg38_masked/Homo_sapiens_assembly38_masked.fasta \
	--fasta_fai /software/projects/pawsey0933/sv/references/hg38_masked/Homo_sapiens_assembly38_masked.fasta.fai \
	--outdir /scratch/pawsey0933/cfolland/vep
