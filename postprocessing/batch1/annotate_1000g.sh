#!/bin/bash -l

#SBATCH --job-name=annotate_vcfs_1000g
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

#load modules
module load vcftools/0.1.16--pl5321hd03093a_7
module load bcftools/1.15--haf5b3da_0
conda activate vcf_tools

cd /scratch/pawsey0933/cfolland/t2t/

#Create a configuration file for vcfanno
cat <<EOF > 1000G_anno.conf
[[annotation]]
file="1000g_af_chr.vcf.gz"
fields=["AF"]
ops=["self"]
names=["AF_1000G"]
type="Float"  # Explicitly define it as a Float
EOF

#Ensure the allele frequency (AF) field is present in the INFO column of the 1000 Genomes VCFs.
for chr in {1..22} X Y; do
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' ./vcfs/Homo_sapiens-GCA_009914755.4-2022_10-1000Genomes_chr${chr}.vcf.gz >> 1000g_af.txt
done

#Convert 1000g_af.txt into a VCF-like format:
awk 'BEGIN {
    print "##fileformat=VCFv4.2"
    print "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">"
    print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
} {
    print $1"\t"$2"\t.\t"$3"\t"$4"\t.\t.\tAF="$5
}' 1000g_af.txt > 1000g_af.vcf

Compress and index:
bgzip 1000g_af.vcf
tabix -p vcf 1000g_af.vcf.gz

#change chr names
bcftools annotate --rename-chrs chr_map.txt 1000g_af.vcf.gz -o 1000g_af_chr.vcf.gz -O z
tabix -p vcf 1000g_af_chr.vcf.gz

#Annotate VCF with allele frequencies using vcfanno
vcfanno -p 4 1000G_anno.conf sorted_toref_genmod_final_hprc_gnomad.vcf.gz > sorted_toref_genmod_final_hprc_gnomad_1000g.vcf

#Compress and index the annotated VCF
bgzip -c sorted_toref_genmod_final_hprc_gnomad_1000g.vcf > sorted_toref_genmod_final_hprc_gnomad_1000g.vcf.gz
tabix -p vcf sorted_toref_genmod_final_hprc_gnomad_1000g.vcf.gz
