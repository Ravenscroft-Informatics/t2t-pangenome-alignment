#!/bin/bash -l

#SBATCH --job-name=mosdepth
#SBATCH --account=pawsey0933
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out
#SBATCH --export=ALL
 
file=$1
filename=$(basename "$file")
basename=${filename%%.*}
input_dir=/scratch/pawsey0933/cfolland/t2t/crams/chm13
hg38_ref=/software/projects/pawsey0933/sv/references/hg38_masked/Homo_sapiens_assembly38_masked.fasta
chm13_ref=/software/projects/pawsey0933/pangenome/refs/chm13v2.0.maskedY.rCRS.EBV.fasta
hg38_bed=/software/projects/pawsey0933/benchmarking/nmd_genes/nmd_gene_list_hg38.bed
chm13_bed=/software/projects/pawsey0933/benchmarking/nmd_genes/nmd_gene_list_chm13.bed
sif=/software/projects/pawsey0933/benchmarking/nmd_genes/mosdepth_0.3.3--h37c5b7d_2.sif

#singularity pull quay.io/biocontainers/mosdepth:0.3.3--h37c5b7d_2

module load singularity/4.1.0-slurm

cd /scratch/pawsey0933/cfolland/mosdepth/results

singularity exec $sif \
    mosdepth \
    --fast-mode \
    --by $chm13_bed \
    --fasta $chm13_ref \
    $basename.chm13 \
    $file
