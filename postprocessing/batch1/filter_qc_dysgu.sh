#!/bin/bash -l

#SBATCH --job-name=filter_qc
#SBATCH --account=pawsey0933
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out
#SBATCH --export=ALL

#load modules
module load vcftools/0.1.16--pl5321hd03093a_7
module load bcftools/1.15--haf5b3da_0
conda activate vcf_tools

cd /scratch/pawsey0933/cfolland/t2t/
ref=/software/projects/pawsey0933/t2t/chm13v2.0.maskedY.rCRS.EBV.fasta

#convert from multiallelic to bi-allelic
bcftools norm -m -any --check-ref skip -f $ref -o T2T_snps_dysgu_VEP_genmod_annotated_biallelic.vcf.gz ../annotated/T2T_snps_dysgu_VEP_genmod_annotated.vcf.gz
tabix -p vcf T2T_snps_dysgu_VEP_genmod_annotated_biallelic.vcf.gz

#Filter based on qc 
bcftools view -v snps,indels T2T_snps_dysgu_VEP_genmod_annotated_biallelic.vcf.gz | bcftools filter -e 'INFO/SVMETHOD == "DYSGUv1.7.0"' | bcftools view -i 'FILTER="PASS" && QUAL > 10' -Oz -o T2T_snps_dysgu_VEP_genmod_annotated_dv_qc.vcf.gz
bcftools view -v other,indels T2T_snps_dysgu_VEP_genmod_annotated_biallelic.vcf.gz | bcftools filter -i 'INFO/SVMETHOD == "DYSGUv1.7.0"' | bcftools view -i 'FILTER="PASS"' -Oz -o T2T_snps_dysgu_VEP_genmod_annotated_dysgu_qc.vcf.gz

#tabix
tabix -p vcf T2T_snps_dysgu_VEP_genmod_annotated_dysgu_qc.vcf.gz
tabix -p vcf T2T_snps_dysgu_VEP_genmod_annotated_dv_qc.vcf.gz

#concatenate files back
bcftools concat -a T2T_snps_dysgu_VEP_genmod_annotated_dv_qc.vcf.gz T2T_snps_dysgu_VEP_genmod_annotated_dysgu_qc.vcf.gz -Oz -o T2T_snps_dysgu_VEP_genmod_annotated_qc.vcf.gz

#Index the filtered VCF
tabix -p vcf T2T_snps_dysgu_VEP_genmod_annotated_qc.vcf.gz
