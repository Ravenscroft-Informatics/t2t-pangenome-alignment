#!/bin/bash -l

#SBATCH --job-name=merge_bed
#SBATCH --account=pawsey0933
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out
#SBATCH --export=ALL

dir=/scratch/pawsey0933/cfolland/t2t/batch2/annotation
svaf=${dir}/SVAFotate_core_SV_popAFs.CHM13.v4.1.bed.gz
hprc=${dir}/hprc-v2.0-mc-chm13.wave.biallelic.sorted.len50.bed.gz
hgsvc=${dir}/hgsvc3-hprc-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.biallelic.af.sorted.len50.bed.gz
out=${dir}/SVAFotate_core_SV_popAFs.CHM13.v4.2.bed.gz


{
zcat "${svaf}" | head -n 1

(
  zcat "${svaf}" | tail -n +2
  zcat "${hprc}" | grep -v '^#'
  zcat "${hgsvc}" | grep -v '^#'
) | sort -k1,1 -k2,2n

} | gzip > "${out}"
