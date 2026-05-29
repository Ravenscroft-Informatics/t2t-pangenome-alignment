#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH -c 2
#SBATCH --job-name=sarek_vep_chm13
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
module load nextflow/24.10.0 #module load nextflow/23.10.0 #nextflow/24.04.3

nextflow run nf-core/sarek \
	-r 3.4.0 \
	-profile singularity,pawsey_setonix \
	-config /software/projects/pawsey0933/gmonahan/T2T/bam2fastq/pawsey_setonix.config \
        -resume \
        --input /software/projects/pawsey0933/t2t/vep/chm13-dv.csv \
	--step annotate \
        --genome null \
        --igenomes_ignore \
	--tools vep \
	--vep_custom_args "--everything --filter_common --per_gene --total_length --offline --format vcf --max_sv_size 100000000 --custom file=/software/projects/pawsey0933/gmonahan/T2T/vep/chm13v2.0_ClinVar20220313.vcf.gz,short_name=clinvar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN" \
	--vep_cache /software/projects/pawsey0933/gmonahan/T2T/vep/ \
	--vep_species homo_sapiens_gca009914755v4 \
	--vep_genome T2T-CHM13v2.0 \
	--vep_cache_version 107 \
	--vep_include_fasta \
	--fasta /software/projects/pawsey0933/pangenome/refs/chm13v2.0.maskedY.rCRS.EBV.fasta \
	--fasta_fai /software/projects/pawsey0933/pangenome/refs/chm13v2.0.maskedY.rCRS.EBV.fasta.fai \
	--outdir /scratch/pawsey0933/cfolland/vep
