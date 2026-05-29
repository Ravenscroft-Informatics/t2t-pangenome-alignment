#!/bin/bash

cd /path/to/pangenie.vcfs

#index pangenie vcfs
for vcf in *genotyping_biallelic.vcf; do
    bgzip $vcf
    tabix -p vcf ${vcf}.gz
done

#merge VCF
bcftools merge -O v -o pangenie-t2t-merged.vcf.gz --write-index *genotyping_biallelic.vcf.gz

#index merged VCF
tabix -p vcf pangenie-t2t-merged.vcf.gz

#add in allele count
bcftools +fill-tags pangenie-t2t-merged.vcf.gz -o pangenie-t2t-merged-ac.vcf.gz -- -t AC,AN

#index merged VCF with  AC
tabix -p vcf pangenie-t2t-merged-AC.vcf.gz

#keep variants with allele count greater than 5
bcftools view -i 'INFO/AC >= 5' pangenie-t2t-merged-ac.vcf.gz -o pangenie-t2t-merged-ac-5.vcf.gz

#index filtered VCF
tabix -p vcf pangenie-t2t-merged-ac-5.vcf.gz

#use isec to filter for variants that are not in the pangenome
bcftools isec -C sorted_toref_genmod.vcf.gz pangenie-t2t-merged-ac-5.vcf.gz -p filtered
