  
# Run VEP annotation on joint-called VCFs
`vep/`

To run `vep` on hg38 data, first download and unzip the cache:
```
wget https://ftp.ensembl.org/pub/release-114/variation/indexed_vep_cache/homo_sapiens_vep_114_GRCh38.tar.gz
tar -xvf homo_sapiens_vep_114_GRCh38.tar.gz
```
See run scripts for specifications.

# Combine deepvariant and dysgu VCFs
Run `concat_b4Anno.sh` on the VEP-annotated VCFs from Dysgu and DeepVariant. 
This produces a merged cohort VCF file that is sorted and does not contain any variants without CSQ annotation.

# Annotate with Genmod
`genmod/`

`Genmod` is memory-intensive to run, so we first split the VCF by chromosome using `genmod_split.sh`, before running genmod on each VCF using `genmod_chr.sh`. The `.ped` file describes the family relationships required for genmod to annotation inheritance.

After running Genmod, use `bgzip` on each vcf to compress and index the VCFs (`bgzip.sh`) before concatenating them with `bcftools concat` and sorting / indexing the cohort VCF (`concat-genmod.sh`).

# Add control population allele frequencies
- SVs and small variants need to be annotated separately. Use script `extract_dv.sh` and `extract_sv.sh` to separate VCF into two VCFs, one containing SNPs/Indels (`dv`) and the other containing SVs (`dysgu`).
- For annotation of SVs, see `sv-annotation/`.
- To annotate the small variants, use the script `annotate_smallvar.sh`

- 3202 srGS samples from 1000 Genomes Phase 3 recalled on T2T-CHM13
- The Human Pangenome Reference Consortium (HPRC) v1.0 T2T-CHM13 callset, including 44 high quality diploid human assemblies, as described in [Liao et al. 2023](https://www.nature.com/articles/s41586-023-05896-x). Sourced from https://zenodo.org/records/7839719 
- The combined HPRC and HGSVC3 T2T-CHM13 callset, including 42 HPRC assemblies + 65 HGVCS3, as described in [Logsdon et al. 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11451754/). VCF available from [here](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Graph_Genomes/1.0/2024_02_23_minigraph_cactus_hgsvc3_hprc/) (with [README](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Graph_Genomes/1.0/2024_02_23_minigraph_cactus_hgsvc3_hprc/20240223_Graph_Genomes_hgsvchprc_README.txt))
- 76,156 srGS aligned to GRCh38 from gnomAD v3.1.2 and lifted over to T2T-CHM13

See scripts: `combined_annotation.sh`

# Filter out low quality variants
- Filter for PASS variants
- Filter for variants in QUAL scores >10 (for SNPs/indels)
See script: `filter_qc_dysgu.sh`

Note, there is also a script for filtering callsets that include Manta and Tiddit:
- Filter for variants in QUAL scores >10 (for SVs called by Tiddit)
- Filter for variants with QUAL scores > 250 (for SVs called using Manta)
See script: `filter_qc_manta_tiddit.sh`

OPTIONAL - Filter out variants with allele frequencies greater than 0.01 in gnomAD and 1000G, and variants with AC > 10 in the callset.
See script: `filter_af.sh`

# Create family specific VCFs
- Subset the VCF to include variants only present in a given family
- Further filter this VCF to subset into VCFs containing only variants annotated as AD or AR by genmod
See script: `family.sh`

# Parse the VCF files to make them readable
- This step converts the VCF files into either TSVs or CSVs which can then be analysed in Excel or another platform
- The script also extracts out useful information from the INFO field in the VCFs and assigns them as their own columns for ease of filtering
See scripts: `process_vcfs.sh`, `parse.py`




