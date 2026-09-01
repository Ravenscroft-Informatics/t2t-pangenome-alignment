#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH -c 48
#SBATCH --job-name=dysgu
#SBATCH --partition=work
#SBATCH --account=pawsey0933
#SBATCH --mem=230gb
#SBATCH --time=0-06:00:00
#SBATCH --export=NONE
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out

set -ueo pipefail

conda activate /software/projects/pawsey0933/gmonahan/envs/dysgu

rundir=/scratch/pawsey0933/cfolland/benchmark/vcfs/dysgu
ref=/software/projects/pawsey0933/sv/references/hg38_masked/Homo_sapiens_assembly38_masked.fasta
samples=/software/projects/pawsey0933/dysgu/bams.txt

for bam in `cat ${samples}`
do
	base=`basename $bam`
	dysgu run --low-mem -p96 -x -c -v2 ${ref} ${rundir}/tmp ${bam} > ${rundir}/${base}_dysgu.vcf
done
