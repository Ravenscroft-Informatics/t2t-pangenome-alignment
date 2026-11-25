# Allele mapping bias
Assess the mapping bias of linear (GRCh38) alignments (BWA-MEM) and GRCh38-surjected graph alignments (Giraffe and vg surject). 
Recreate the analysis of allele bias from paper ["Pangenomics enables genotyping of known structural variants in 5202 diverse genomes" (Siren et al. 2021)](https://www.science.org/doi/10.1126/science.abg8871). See [Figure 4A](https://www.science.org/doi/10.1126/science.abg8871#F4).

Use BAM/CRAM files generated previously by `nfcore/sarek` (linear, bwa-mem) and `vg_snakemake` (pangenome, giraffe). 

## Variant calling
Use `call_mpileup_giraffe_bwa.sh` to:
1. Call variants using `bcftools mpileup` and `call` in linear and pangenome alignments.
2. Merge VCFs to identify heterozygous variants common between the callsets.
3. Filter for high quality variants (`INFO/DP>=25 && INFO/MQ>=40`).

### DeepVariant
Prepare exisitng DeepVariant callsets using `prep_dfv_calls.sh`
```
bash prep_dv_calls.sh hg002 $giraffe $bwa $ref ./deepvariant_pg/
```

## Plotting
Use `plot_ab_mapers.py` to generate a plot showing the fraction of alternate alleles across reads for high-quality heterozygous variants, divided by allele length.
```
python plot_ab_mappers.py ./mpileup_calls/common.het.filtered.vcf.gz hg002_ab.png
```

### Note
DeepVariant derives the alelle depths from model-informed representations of aligned reads rather than raw pileup counts. Therefore, DeepVariant-derived allele depths reflect both mapping and model-specific filtering and interpretation of reads rather than an unbiased assessment of read counts. This is why we used `bcftools mpileup` for an unbiased assessment of allele depth.
