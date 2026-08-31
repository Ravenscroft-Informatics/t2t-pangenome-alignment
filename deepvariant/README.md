# Pangenome-aware DeepVariant
Run pangenome-aware DeepVariant on cohort (234 samples) using a job array to submit jobs for each sample/contig combination.

## Requirements
- Pangenome-aware DeepVariant v1.10.0
```
singularity pull deepvariant_1.10.0.sif https://hub.docker.com/layers/google/deepvariant/1.10.0/images/sha256-962e5a83b1d76aae6990625d47102785f791603f2138aa1fa9aa4fb6a2eecbe6
singularity pull deepvariant_pangenome_aware_deepvariant-1.10.0.sif https://hub.docker.com/layers/google/deepvariant/pangenome_aware_deepvariant-1.10.0/images/sha256-996e456cd839e29abecd0fd3a9f552e8e09b59f945acfeb13022a2b88f9a2ede
```
- Bed files containing interval regions: see `./regions`. See "Notes" below. 
- Pangenome graph file (gbz) - `hprc-v2.0-mc-grch38.gbz`
```
wget https://human-pangenomics.s3.amazonaws.com/pangenomes/freeze/release2/minigraph-cactus/hprc-v2.0-mc-grch38.gbz
```

## Run using the following scripts:
- `make_manifest.sh` - generate a manifest specifying each sample / contig combination
- `deepvariant_array.slurm` - array script that submits one job for each sample / contig combination
- `submit_pipeline.sh`- used to execute `deepvariant_array.slurm` in a way that monitors the total number of jobs as to not overwhelm the slurm scheduler
- `make_missing.sh` - generate a manifest of sample / contig combinations that remain for completion

**How to use**
```
$ bash make_manifest.sh > 01_tasks.tsv
$ bash submit_pipeline.sh 01_tasks.tsv
$ bash make_missing.sh 01_tasks.tsv > 02_tasks.tsv
```

## Notes
- Regions bed files derived from `sarek/nfcore` workflow [`create_intervals_bed`](https://github.com/nf-core/sarek/blob/3.4.2/modules/local/create_intervals_bed/main.nf)
- 
