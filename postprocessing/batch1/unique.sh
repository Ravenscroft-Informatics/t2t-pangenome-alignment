#!/bin/bash -l

#SBATCH --job-name=annotate_unique
#SBATCH --account=pawsey0933
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --time=6:00:00
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out
#SBATCH --export=ALL

# Load modules
module load vcftools/0.1.16--pl5321hd03093a_7
module load bcftools/1.15--haf5b3da_0
conda activate vcf_tools

cd /scratch/pawsey0933/cfolland/t2t/annotation

# Sort the BED file
sort -k1,1 -k2,2n /software/projects/pawsey0933/t2t/refs/chm13v2-unique_to_hg38.bed > chm13v2-unique_to_hg38_sorted.bed

#add annotation column
awk '{print $0"\tUniqueT2T"}' chm13v2-unique_to_hg38_sorted.bed > chm13v2-unique_to_hg38_with_label.bed

#gzip and index
bgzip chm13v2-unique_to_hg38_with_label.bed
tabix -p bed chm13v2-unique_to_hg38_with_label.bed.gz

# Create a configuration file for vcfanno
cat <<EOF > unique.conf
[[annotation]]
file="chm13v2-unique_to_hg38_with_label.bed.gz"
columns=[4]
ops=["flag"]
names=["UniqueT2T"]
EOF

# Annotate VCF with unique regions using vcfanno
vcfanno -p 4 unique.conf T2T_snps_dysgu_VEP_genmod_annotated.vcf.gz > T2T_snps_dysgu_VEP_genmod_annotated_unique.vcf

# Compress and index the annotated VCF
bgzip -c T2T_snps_dysgu_VEP_genmod_annotated_unique.vcf > T2T_snps_dysgu_VEP_genmod_annotated_unique.vcf.gz
tabix -p vcf T2T_snps_dysgu_VEP_genmod_annotated_unique.vcf.gz
