#!/bin/bash -l

#SBATCH --job-name=annotate_vcfs_gnomad
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

#Ensure the allele frequency (AF) field is present in the INFO column of the 1000 Genomes VCFs.
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' ./vcfs/Homo_sapiens-GCA_009914755.4-2022_10-gnomad.vcf.gz > gnomad_af.txt

#Convert 1000g_af.txt into a VCF-like format:
awk 'BEGIN {
    print "##fileformat=VCFv4.2"
    print "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">"
    print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
} {
    print $1"\t"$2"\t.\t"$3"\t"$4"\t.\t.\tAF="$5
}' gnomad_af.txt > gnomad_af.vcf

#Ensure gnomad VCF is bgzipped and indexed
bgzip -c gnomad_af.vcf > gnomad_af.vcf.gz
tabix -p vcf gnomad_af.vcf.gz

#Create a configuration file for vcfanno
cat <<EOF > gnomad_anno.conf
[[annotation]]
file="gnomad_af.vcf.gz"
fields=["AF"]
ops=["self"]
names=["AF_gnomad"]
type="Float"  # Explicitly define it as a Float
EOF

#Annotate VCF with allele frequencies using vcfanno
vcfanno -p 4 gnomad_anno.conf sorted_toref_genmod_final_hprc.vcf.gz > sorted_toref_genmod_final_hprc_gnomad.vcf

#Compress and index the annotated VCF
bgzip -c sorted_toref_genmod_final_hprc_gnomad.vcf > sorted_toref_genmod_final_hprc_gnomad.vcf.gz
tabix -p vcf sorted_toref_genmod_final_hprc_gnomad.vcf.gz
