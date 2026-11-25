# Generating data

## Generate a .BED file containing coordinates for NMD genes

1) Source a list of genes associated with neuromuscular diseases. I used the list of genes run on the [PathWest neuromuscular gene panel](https://pathwest.health.wa.gov.au/Our-Services/Clinical-Services/Diagnostic-Genomics/Neurogenetics)
- This encompasses 912 genes, and FRG1 and FRG2 were added manually : `nmd_gene_list.tsv`

2) Source the appropriate annotation files
- I used the RefSeq annotation from 2025_08.
  - GRCh38 available [here](https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/)
  - CHM13 available [here](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/)
```
# chm13
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/GCF_009914755.1-RS_2025_08/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz
# this corresponds to ensembl v110 release for hg38
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/GCF_000001405.40-RS_2025_08/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz
```

3) Run script `gtf-bed_transcripts.sh` to extract out transcripts from gtf file
- Extracting out the coordinates from the “gene” entry in the gtf file was fallible as each “gene” is defined differently across the reference builds.
- It was best therefore to extract out the transcripts, find the intersection of the transcripts described across the two reference builds, then aggregate the coordinates across the transcripts to essentially define new gene regions for each reference build.

4) Find the intersection of transcripts described across the two reference builds with `intersect.sh`

5)  Map the chromosome RefSeq IDs (NC_00#) to chromosome names (chr#)
- The new GTF files contain Refseq IDs `NC_00*` instead of chromosome names i.e. `chr1`
    - Create a mapping file for the RefSeq Ids to chromosome names: `refseq-chr_chm13.tsv`  `refseq-chr_grch38.tsv`
    - Create a new script to change the refseq ids to chromsome names in the bed file: `refseq-chr.sh`
```
$ ./refseq-chr.sh refseq-chr_grch38.tsv GCF_000001405.40_GRCh38.p14.transcripts.shared.bed GCF_000001405.40_GRCh38.p14.transcripts.shared.fixed.bed
```

6) Use `aggregate_genes.sh` to aggregate across transcripts and define a new gene coordinate set for each reference genome bed file

```
$ ./aggregate_genes.sh GCF_000001405.40_GRCh38.p14.transcripts.shared.fixed.bed GCF_000001405.40_GRCh38.p14.genes_aggregate.bed
```
7) Filter bed file to include only NMD genes: `filter.py`
- Output file: `nmd_gene_list_chm13.bed` or `nmd_gene_list_grch38.bed`
```
$ python3 filter.py
```
8) Check that the genes match across all three files:
- `nmd_gene_list_chm13.bed`
- `nmd_gene_list_grch38.bed`
- `nmd_gene_list.tsv`

**1 gene missing NDUFA4 = COXFA4 - manually copy over gene coordinates into the .bed files**

## Run Mosdepth using .BED file as input with `--by` flag
1) Download mosdepth singularity image
```
$ singularity pull quay.io/biocontainers/mosdepth:0.3.3--h37c5b7d_2
```
2) Run Mosdepth on all samples
- Batch mosdepth using `run_mosdepth.sh` script which generates individual jobs for each `.cram` file using `mosdepth.sh`

## Collect region depth results across all samples for NMD regions.
   
The python script `collect_coverage.py` merges the files `*.regions.bed.gz` to create a single merged results file for export. Run the script using `merge_regions_nmd.sh`.
Output files are provided:
- CHM13: [mosdepth_nmd.chm13.linear.merged.tsv](https://github.com/RAVING-Informatics/T2T-alignment/blob/main/benchmarking/nmd_genes/mosdepth_nmd.chm13.linear.merged.tsv)
- GRCh38: [mosdepth_nmd.grch38.linear.merged.tsv](https://github.com/RAVING-Informatics/T2T-alignment/blob/main/benchmarking/nmd_genes/mosdepth_nmd.grch38.linear.merged.tsv)

## Collect per-base coverage for a specific gene/interval of interest.
   
The script `collect_coverage_perbase.py` extracts the per-gene coverage for a single genomic interval from many per-base mosdepth outputs and writes a single long-format table:
`gene	chr	start	end	sample	depth`
The script can be run easy using `merge_files_gene.sh` under `#single gene/interval`. The first time this script is run, lots of resources are needed as the python script will index the per-base.bed.gz files before extracting the depth metrics. After the first run, the time and memory allocations can be reduced significantly.
The output format is `$gene.perbase_mosdepth_${method}_${ref}.tsv`:
```
gene	chr	start	end	sample	depth
ACTN2	chr1	236087383	236087388	D09-468	38
ACTN2	chr1	236087388	236087394	D09-468	39
ACTN2	chr1	236087394	236087398	D09-468	43
```

## Collect per-base coverage for multiple gene/intervals of interest.
The script `collect_coverage_perbase_bed.py` achieves the same result as `collect_coverage_perbase.py`, however instead of taking a single gene/interval as input, it accepts a `.bed` file with multiple gene/intervals as input. This script can be run using `merge_files_gene.sh` under `#bed file with multiple genes/intervals`. Ensure to increase resources as needed. The output format is `all.perbase_mosdepth_${method}_${ref}.tsv`:

## Generate TSV files with exon and intron specifications for a gene
1) The python script `gtf_region_to_tables.py` can be used to extract out the co-ordinates of exons and introns for a gene of interest from a gtf file input.
The python script requires the following inputs:
- gtf: use the original gtf file (that matches the ref)
- region: coordinates of the gene `chr:start-end`
2) Use the bash script `gtf-exons.sh` to run the python script
- ref: either `chm13` or `grch38`
- gtf: same as above
- region_file: this is .bed file generated above (`nmd_gene_list_chm13.bed` or `nmd_gene_list_grch38.bed`) - this is used to collect the coordinates of the gene
- mapfile: use the same file used to map the chromosome RefSeq IDs (NC_00#) to chromosome names (chr#) - i.e. `refseq-chr_hg38.tsv`
- gene: gene name

## Calculate coverage uniformity
1) The output file `all.perbase_mosdepth_${method}_${ref}.tsv` generated [above](https://github.com/RAVING-Informatics/T2T-alignment/blob/main/benchmarking/nmd_genes/README.md#collect-per-base-coverage-for-a-specific-geneinterval-of-interest), can be used to calculate summary statstics for coverage across each gene per sample. If many genes were used as input, this will be a very large file, so it need to be sorted according to gene, sample and depth. 
2) The python script `perbase_agg_summary.sh` calculates the min/max/range/mean/median/SD per gene/sample from a pre-sorted file. Run the sorting and the python script using [`summarise_perbase_mosdepth.sh`](https://github.com/RAVING-Informatics/T2T-alignment/blob/main/benchmarking/nmd_genes/summarise_perbase_mosdepth.sh). The output is in the format: `all.perbase_mosdepth.summary.linear_chm13.tsv`
```
gene	sample	min_depth	max_depth	range_depth	mean_depth	median_depth	sd_depth
AARS1	D09-468.chm13	17	91	74	54.3797	54	10.7496
AARS1	D12-880.chm13	2	46	44	23.4331	24	7.4143
AARS1	D12-881.chm13	6	65	59	33.0484	33	8.47207
```

# Plotting

## Compare per-base gene coverage across specific gene in R
1) Use the script: [`single-gene_perbase_coverage_chm13-hg38.R`](https://github.com/RAVING-Informatics/T2T-alignment/blob/main/benchmarking/nmd_genes/single-gene_perbase_coverage_chm13-hg38.r) to generate a plot of the per-base coverage across a gene of interest.
- Ensure to specify the gene of interest: @line166
- As input, you need the exons and merged mosdepth results, for example:
  - `ACTN2.perbase_mosdepth_hg38.tsv`
  - `ACTN2.perbase_mosdepth_chm13.tsv`
  - `ACTN2.grch38.exons.tsv`
  - `ACTN2.chm13.exons.tsv`
- Example outputs are provided, see `ACTN2.grch38-chm13_perbase_cov.pdf`

## Create plots comparing the coverage of NMD genes between GRCh38 and CHM13
3) Use the files `all.perbase_mosdepth.summary.${method}_${ref}.tsv` as input into [`perbase_coverage_nmd_comparison.r`](https://github.com/RAVING-Informatics/T2T-alignment/blob/main/benchmarking/nmd_genes/perbase_coverage_nmd_comparison.r) to compare the average coverage and coverage uniformity for each NMD gene between CHM13 and GRCh38.
