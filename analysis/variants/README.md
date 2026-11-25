## VARIANT QUALITY
### **`bcftools stats`**

Run `bcftools stats` on the `deepvariant` and `dysgu` VCFs for each sample (not the g.vcf files) using `bcftools_stats_ind.sh`. This is to generate per-sample summary statistics and also a plot of the number of variants of a certain variant quality score. See attached plot for example.

Run `bcftools stats` on the deepvariant and dysgu cohort VCFs using `bcftools_stats_cohort.sh`.

**Inputs**

- Use the individual deepvariant VCFs, i.e. `D09-468.deepvariant.vcf.gz`, and dysgu VCFs, i.e. `D09-468.sorted.cram_dysgu.vcf`.

***Linear-Alignment***
  - CHM13 available here: `/Volumes/PERKINS-LL-001/Sequencing/wgs/secondary/T2T_realignment/batch_1/chm13/variant_calling/`
  - GRCh38 available here: `/Volumes/PERKINS-LL-001/Sequencing/wgs/secondary/T2T_realignment/batch_1/grch38/variant_calling/`

***Pangenome-Alignment***
  - CHM13 + GRCh38 available here: `/Volumes/PERKINS-LL-001/Sequencing/wgs/secondary/Pangenome_realignment/batch_1/vg_giraffe`
     
- Use the cohort deepvariant callsets, i.e. `T2T_dv_glnexus_VEP.ann.vcf.gz`, and cohort dysgu variants, i.e. `dysgu_merge_T2T_VEP.ann_fixed.vcf.gz`
  
***Linear-Alignment***
  - CHM13 available here: `/Volumes/PERKINS-LL-001/Sequencing/wgs/secondary/T2T_realignment/batch_1/chm13/merge`
  - GRCh38 available here: `/Volumes/PERKINS-LL-001/Sequencing/wgs/secondary/T2T_realignment/batch_1/grch38/postprocess/vep` 

***Pangenome-Alignment***
  - CHM13 + GRCh38 available here: `/Volumes/PERKINS-LL-001/Sequencing/wgs/secondary/Pangenome_realignment/batch_1/postprocess/vep`
    - Note, the VCFs for different variant subsets are not available on the IRDS.
    
**Outputs**

- bcftools files are available on the IRDS.  
  - CHM13: `/Volumes/PERKINS-LL-001/Sequencing/wgs/secondary/T2T_realignment/batch_1/chm13/reports/bcftools_v1.15/`
  - GRCh38 `/Volumes/PERKINS-LL-001/Sequencing/wgs/secondary/T2T_realignment/batch_1/grch38/reports/bcftools/`

**Results**

*1. Frequency plot of variant quality*
- Use the R script, `plot_variant_quality_cohort.R` to generate a multifaceted plot of variant quality, subsetted according to particular variant callsets:
  - Genome-wide: genome variants (no filters)
  - Exome: exomic variants
  - Syntenic: variants in regions of CHM13 that are syntenic to GRCh38
  - Genome mask: variants in technically reliable regions accessible using short read mapping
  - ClinVar: variants described in ClinVar 
- Instructions to generate filtered VCFs are available [below](https://github.com/RAVING-Informatics/T2T-alignment/blob/main/benchmarking/README.md#filter-vcfs).
- Once the filtered VCFs are available, run `bcftools_stats_cohort.sh` to generate `bcftools stats` on each. Use `multiqc` to generate a report combining all files, and download the `bcftools_stats_vqc.tsv` data. For example, see the file [`linear_bcftools_stats_vqc_cohort.tsv`](https://github.com/RAVING-Informatics/T2T-alignment/blob/main/benchmarking/variants/linear_bcftools_stats_vqc_cohort.tsv)
- Use the functions in R script: `parse_bcftools_stats_vqc.R` to convert the data into a format that is easier to plot.
- Plot data with quality scores on x-axis and number of variants on y-axis using `plot_variant_quality_cohort.R`

*2. Plot of variant quality metrics*
- Use the bcftools stats individual multiqc reports to export out the general variant stats table.
  - E.g. linear deepvariant: [`linear_grch38-chm13_dv_general_stats_table.tsv`](https://github.com/RAVING-Informatics/T2T-alignment/blob/main/benchmarking/variants/linear_grch38-chm13_dv_general_stats_table.tsv)
  - E.g. linear dysgu: [`linear_grch38-chm13_dysgu_general_stats_table.tsv`](https://github.com/RAVING-Informatics/T2T-alignment/blob/main/benchmarking/variants/linear_grch38-chm13_dysgu_general_stats_table.tsv)
- Calculate the number of per sample pass variants from dysgu VCFs using `pass_variants_dysgu.sh` to generate [`pass_dysgu_linear.tsv`](https://github.com/RAVING-Informatics/T2T-alignment/blob/main/benchmarking/variants/pass_dysgu_linear.tsv)
- Combine data into single csv: [`hprc-linear_variant_stats.csv`](https://github.com/RAVING-Informatics/T2T-alignment/blob/main/benchmarking/variants/hprc-linear_variant_stats.csv)
- Use the R script, [`plot_variant_stats_linear_hprc.R`](https://github.com/RAVING-Informatics/T2T-alignment/blob/main/benchmarking/variants/plot_variant_stats_linear_hprc.R) to generate a multifaceted plot of variant quality metrics, including the following sample-wise metrics derived above:
  - The number of variants called by DeepVariant
  - The number of SVs called by Dysgu
  - The proportion of SVs called by Dysgu that pass the cut-off quality filter
- In addition, include the following family-wise metrics calculated using `VariantEval MendelianViolationEvaluator` as described [below](https://github.com/RAVING-Informatics/T2T-alignment/blob/main/benchmarking/README.md#mendelian-violations).
- The required input file is [`hprc-linear_mendel.csv`](https://github.com/RAVING-Informatics/T2T-alignment/blob/main/benchmarking/variants/hprc-linear_mendel.csv)
  - The Mendelian Violation rate
  - The number of *de novo* variants
  - The number of low-quality variants
  
### Filter VCFs

**Coding regions**

To view variant quality scores for the exome, subset the cohort VCF file to include only exons using a genome annotation file:

- T2T-CHM13:
```
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13.draft_v2.0.gene_annotation.gff3
```
- GRCh38:
```
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.annotation.gff3.gz
```

```
#create bed file to include only exons
zgrep "protein_coding" $gtf | awk '$3 == "exon" {print $1"\t"$4-1"\t"$5}' > protein_coding_{genome}.bed
# use bcftools view to subset vcf to only exons and index output
bcftools view -R protein_coding_chm13.bed -Oz -o hprc-v1.1-mc-chm13_dv_glnexus_VEP.gff.exons.vcf.gz hprc-v1.1-mc-chm13_dv_glnexus_VEP.ann.vcf.gz
tabix -p vcf hprc-v1.1-mc-chm13_dv_glnexus_VEP.gff.exons.vcf.gz
```

**Clinvar Variants**

To view variant quality scores for clinvar variants, intersect the cohort VCF file with clinvar variants:

- T2T-CHM13:
```
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/liftover/chm13v2.0_ClinVar20220313.vcf.gz
```
- GRCh38:
```
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_20220313.vcf.gz
```
Modify the GRCh38 Clinvar VCF to include only variants that could be lifted over to T2T-CHM13:
```
bcftools query -f '%ID\n' chm13v2.0_ClinVar20220313.vcf.gz | sort | uniq > chm13v2.0_ClinVar.ids
bcftools query -f '%ID\n' clinvar_20220313.vcf.gz | sort | uniq > grch38_ClinVar.ids
comm -12 grch38_ClinVar.ids chm13v2.0_ClinVar > common.ids
bcftools view -i 'ID=@common.ids' clinvar_20220313.vcf.gz -Oz -o common_clinvar_20220313.vcf.gz
bcftools view -H common_clinvar_20220313.vcf.gz | wc -l
>1113862
bcftools view -H chm13v2.0_ClinVar20220313.vcf.gz | wc -l
1113862
```
Intersect the clinvar variants with the variant callsets using `intersect_clinvar.sh`

**Masked Regions**

To view the variant quality scores for variants in technically reliable regions accessible using short read mapping, download the short-read accessibility mask:

```
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/accessibility/combined_mask.bed
https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/accessibility/hg38.combined_mask.bed.gz
gzip -d hg38.combined_mask.bed.gz
```
Intesect the cohort VCF file with the mask using the `intersect_mask.sh` script. 

**Syntenic Regions**

To view the variant quality scores for variants in regions of CHM13 that are syntenic to GRCh38, download the regions that are non-syntenic between GRCh38 and CHM13:

```
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-unique_to_hg38.bed
```

Use the script `syntenic_regions.sh` to generate a cohort VCF containing only syntenic variants. 

### **Mendelian Violations**

Calculate mendelian-violation rate using GATK `VariantEval MendelianViolationEvaluator`

**Inputs**

Use the cohort deepvariant callsets annotated with VEP: 

***Linear-Alignment***
- CHM13: `/Volumes/PERKINS-LL-001/Sequencing/wgs/secondary/T2T_realignment/batch_1/chm13/merge/glnexus+dysgu/T2T_dv_glnexus_VEP.ann.vcf.gz`
- GRCh38: `/Volumes/PERKINS-LL-001/Sequencing/wgs/secondary/T2T_realignment/batch_1/grch38/postprocess/vep/dv_glnexus_VEP.ann.vcf.gz`

***Pangenome-Alignment***
- CHM13:`/Volumes/PERKINS-LL-001/Sequencing/wgs/secondary/Pangenome_realignment/batch_1/postprocess/vep/hprc-v1.1-mc-chm13_dv_glnexus_VEP.ann.vcf.gz`
- GRCh38: `/Volumes/PERKINS-LL-001/Sequencing/wgs/secondary/Pangenome_realignment/batch_1/postprocess/vep/hprc-v1.1-mc-grch38_dv_glnexus_VEP.ann.vcf.gz`

**Setup Instructions**

Pull the singularity image from docker hub
```
singularity pull docker://broadinstitute/gatk@sha256:71b17ee42d149e8ec112603f5305c873ab60d93949ef8bb62a4fff85427f56fb
```

**Run Script**

[`gatk_mendel.sh`](https://github.com/RAVING-Informatics/T2T-alignment/blob/main/benchmarking/variants/gatk_mendel.sh)

**Outputs**
- Outputs are available for genome-wide and exome analysis:
  - chm13: [`chm13_linear.MVs.byFamily.fixed.table`](https://github.com/RAVING-Informatics/T2T-alignment/blob/main/benchmarking/variants/chm13_linear.MVs.byFamily.fixed.table) and [`chm13_linear.MVs.byFamily.exons.fixed.table`](https://github.com/RAVING-Informatics/T2T-alignment/blob/main/benchmarking/variants/chm13_linear.MVs.byFamily.exons.fixed.table)
  - grch38: [`grch38_linear.MVs.byFamily.fixed.table`](https://github.com/RAVING-Informatics/T2T-alignment/blob/main/benchmarking/variants/grch38_linear.MVs.byFamily.fixed.table) and [`grch38_linear.MVs.byFamily.exons.fixed.table`](https://github.com/RAVING-Informatics/T2T-alignment/blob/main/benchmarking/variants/grch38_linear.MVs.byFamily.exons.fixed.table)

**Results**
- Calculate the key results using the following spreadsheet [`mendel_viol_linear_calculations.xlsx`](https://github.com/RAVING-Informatics/T2T-alignment/blob/main/benchmarking/variants/mendel_viol_linear_calculations.xlsx).
