# Pangenome-aware DeepVariant
Run pangenome-aware DeepVariant on cohort (234 samples) using a job array to submit jobs for each sample/contig combination.

## Requirements
- Pangenome-aware DeepVariant v1.10.0
```
singularity pull deepvariant_1.10.0.sif https://hub.docker.com/layers/google/deepvariant/1.10.0/images/sha256-962e5a83b1d76aae6990625d47102785f791603f2138aa1fa9aa4fb6a2eecbe6
singularity pull deepvariant_pangenome_aware_deepvariant-1.10.0.sif https://hub.docker.com/layers/google/deepvariant/pangenome_aware_deepvariant-1.10.0/images/sha256-996e456cd839e29abecd0fd3a9f552e8e09b59f945acfeb13022a2b88f9a2ede
```
- Bed files containing interval regions: see `./regions`. See "Notes" below. 
- HPRCv2.0 pangenome graph file (gbz) - `hprc-v2.0-mc-grch38.gbz`
```
wget https://human-pangenomics.s3.amazonaws.com/pangenomes/freeze/release2/minigraph-cactus/hprc-v2.0-mc-grch38.gbz
```
- Compatible reference fasta file: `hprc-v2.0-mc-grch38.fa`. See "Notes" below.

## Run using the following scripts:
- `make_manifest.sh` - generate a manifest specifying each sample / contig combination
- `deepvariant_array.slurm` - array script that submits one job for each sample / contig combination
- `submit_pipeline.sh`- used to execute `deepvariant_array.slurm` in a way that monitors the total number of jobs as to not overwhelm the slurm scheduler
- `make_missing.sh` - generate a manifest of sample / contig combinations that remain for completion
- `generate_sample_list.sh` - generate a list of samples ready for merging (VCFs for all contigs available)
- `merge_dv_array.sh` - array script to merge each samples contigs into one VCF

**How to use**
```
# Run deepvariant
$ bash make_manifest.sh > 01_tasks.tsv
$ bash submit_pipeline.sh 01_tasks.tsv
$ bash make_missing.sh 01_tasks.tsv > 02_tasks.tsv
$ bash submit_pipeline.sh > 02_tasks.tsv

# Merge deepvariant contigus / sample
$ bash generate_sample_list.sh
$ sbatch merge_dv_array.sh
```

## Notes
- Regions bed files derived from `sarek/nfcore` workflow [`create_intervals_bed`](https://github.com/nf-core/sarek/blob/3.4.2/modules/local/create_intervals_bed/main.nf)
- Generate FASTA file using `vg_snakemake` rule `extract_ref_fasta` which leverages `vg paths` to extract GRCh38 reference paths from the graph: https://github.com/vgteam/vg_snakemake/blob/master/workflow/rules/index.smk.
    - Need to add notes on correct ordering of dict

# GLNexus
Generate a joint-called cohort VCF using GLNexus.

## Scripts
- `prep_glnexus.sh` - Generates per-chromosome BED files for scattering GLnexus joint calling across a SLURM array.
- `glnexus_array.sh` - Run script for launching SLURM array task for Glnexus for each contig.