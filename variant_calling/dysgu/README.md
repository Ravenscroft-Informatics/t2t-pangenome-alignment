# Dysgu

## Installation
- Install Dysgu using conda:
```
conda create -n dysgu -c bioconda dysgu=1.7.0
```

## Run
- Run dysgu on cohort using `run_dysgu.sh`
- Split individual VCFs by contig and perform cohort merge for each contig - `dysgu_merge.sh`
- Concatenate cohort-wide VCFs across all contigs - `concat_dysgu.sh`

## Reference
Kez Cleal, Duncan M Baird, Dysgu: efficient structural variant calling using short or long reads, Nucleic Acids Research, Volume 50, Issue 9, 20 May 2022, Page e53, https://doi.org/10.1093/nar/gkac039