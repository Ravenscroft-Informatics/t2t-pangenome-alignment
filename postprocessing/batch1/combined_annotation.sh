#!/bin/bash -l

#SBATCH --job-name=annotate_vcfs
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

cd /scratch/pawsey0933/cfolland/t2t/annotation

#Ensure the allele frequency (AF) field is present in the INFO column of the VCFs.
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' /software/projects/pawsey0933/pangenie/refs/chm13_cactus_filtered_ids_biallelic.vcf.gz > hprc_af.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' /scratch/pawsey0933/cfolland/t2t/vcfs/hgsvc3-hprc-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz > hgsvc3_af.txt
for chr in {1..22} X Y; do
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' ./vcfs/Homo_sapiens-GCA_009914755.4-2022_10-1000Genomes_chr${chr}.vcf.gz >> 1000g_af.txt
done
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' ./vcfs/Homo_sapiens-GCA_009914755.4-2022_10-gnomad.vcf.gz > gnomad_af.txt

#Convert txt files into VCF-like format:
awk 'BEGIN {
    print "##fileformat=VCFv4.2"
    print "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">"
    print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
} {
    print $1"\t"$2"\t.\t"$3"\t"$4"\t.\t.\tAF="$5
}' hprc_af.txt > hprc_af.vcf

awk 'BEGIN {
    print "##fileformat=VCFv4.2"
    print "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">"
    print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
} {
    print $1"\t"$2"\t.\t"$3"\t"$4"\t.\t.\tAF="$5
}' hgsvc3_af.txt > hgsvc3_af.vcf

awk 'BEGIN {
    print "##fileformat=VCFv4.2"
    print "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">"
    print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
} {
    print $1"\t"$2"\t.\t"$3"\t"$4"\t.\t.\tAF="$5
}' 1000g_af.txt > 1000g_af.vcf

awk 'BEGIN {
    print "##fileformat=VCFv4.2"
    print "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">"
    print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
} {
    print $1"\t"$2"\t.\t"$3"\t"$4"\t.\t.\tAF="$5
}' gnomad_af.txt > gnomad_af.vcf

#Ensure VCFs are bgzipped and indexed
bgzip -c hprc_af.vcf > hprc_af.vcf.gz
tabix -p vcf hprc_af.vcf.gz
bgzip -c hgsvc3_af.vcf > hgsvc3_af.vcf.gz
tabix -p vcf hgsvc3_af.vcf.gz
bgzip -c 1000g_af.vcf
tabix -p vcf 1000g_af.vcf.gz
bgzip -c gnomad_af.vcf > gnomad_af.vcf.gz
tabix -p vcf gnomad_af.vcf.gz

#change chr names
bcftools annotate --rename-chrs chr_map.txt 1000g_af.vcf.gz -o 1000g_af_chr.vcf.gz -O z
tabix -p vcf 1000g_af_chr.vcf.gz

#Create a configuration file for vcfanno
cat <<EOF > combined_anno.conf
[[annotation]]
file="hprc_af.vcf.gz"
fields=["AF"]
ops=["self"]
names=["AF_HPRC"]
type="Float"  # Explicitly define it as a Float

[[annotation]]
file="gnomad_af.vcf.gz"
fields=["AF"]
ops=["self"]
names=["AF_gnomad"]
type="Float"  # Explicitly define it as a Float

[[annotation]]
file="1000g_af_chr.vcf.gz"
fields=["AF"]
ops=["self"]
names=["AF_1000G"]
type="Float"  # Explicitly define it as a Float

[[annotation]]
file="hgsvc3_af.vcf.gz"
fields=["AF"]
ops=["self"]
names=["AF_hgsvc3"]
type="Float"  # Explicitly define it as a Float
EOF

#Annotate VCF with allele frequencies using vcfanno
vcfanno -p 4 combined_anno.conf ../cohort/T2T_snps_dysgu_VEP_ac10_CSQonly_rehead_genmod.vcf.gz > T2T_snps_dysgu_VEP_genmod_annotated.vcf

#Compress and index the annotated VCF
bgzip -c T2T_snps_dysgu_VEP_genmod_annotated.vcf > T2T_snps_dysgu_VEP_genmod_annotated.vcf.gz
tabix -p vcf T2T_snps_dysgu_VEP_genmod_annotated.vcf.gz
