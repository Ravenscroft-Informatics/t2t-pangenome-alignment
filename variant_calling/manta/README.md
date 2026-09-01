# Structural Variant Calling

## Manta
For pangenome-aligned / surjected BAMs, use the following version of Manta which has been made compatible with pangenome-based alignments (derived from [`vg_snakemake`](https://github.com/vgteam/vg_snakemake/blob/master/workflow/rules/common.smk)):
```
$ singularity pull manta-pg.1.6.sif docker://quay.io/jmonlong/manta:main
```
