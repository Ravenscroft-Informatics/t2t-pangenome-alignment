## MAPPING QUALITY 
### **`samtools stats`**

Run `samtools stats` on the mapped `BAM` files output from `sarek` or `vg surject`.
Use the `samtools_submit_all.sh` script to specify the location of the `BAM` files. This script will iterate through and submit individual jobs for each in combination with `samtools_stats.sh`.
- Alternatively, `samtools stats` is run as part of the `sarek/nfcore` workflow, with individual results stored in `./reports/samtools/`.

**Inputs**

***Linear-Alignment***
- CHM13 Mapped BAMs: `/Volumes/PERKINS-LL-001/Sequencing/wgs/secondary/T2T_realignment/batch_1/chm13/preprocessing/mapped`
- GRCh38 Mapped BAMs: `/Volumes/PERKINS-LL-001/Sequencing/wgs/secondary/T2T_realignment/batch_1/grch38/preprocessing/mapped`
  
***Pangenome-Alignment***
- Mapped BAMs: `/Volumes/PERKINS-LL-001/Sequencing/wgs/secondary/Pangenome_realignment/batch_1/vg_giraffe/bams`

**Outputs**
- Results for this project are available on the IRDS:
  - CHM13: `/Volumes/PERKINS-LL-001/Sequencing/wgs/secondary/T2T_realignment/batch_1/chm13/reports/samtools`
  - GRCh38: `/Volumes/PERKINS-LL-001/Sequencing/wgs/secondary/T2T_realignment/batch_1/grch38/reports/samtools`
- `multiqc` reports for this project available at: 
  - chm13: `/Volumes/PERKINS-LL-001/Sequencing/wgs/secondary/T2T_realignment/batch_1/chm13/multiqc/multiqc_report.html`
  - grch38: `/Volumes/PERKINS-LL-001/Sequencing/wgs/secondary/T2T_realignment/batch_1/grch38/multiqc/multiqc_report.html`

**Results**
- Use the `multiqc` report to extract out useful metric, including the number of mapped / unmapped reads.
  - Export the `samtools-stats-dp` data. See examples [`linear_sarek_chm13_samtools-stats-dp.tsv`](https://github.com/RAVING-Informatics/T2T-alignment/blob/main/benchmarking/mapping/linear_sarek_chm13_samtools-stats-dp.tsv) and [`linear_sarek_grch38_samtools-stats-dp.tsv`](https://github.com/RAVING-Informatics/T2T-alignment/blob/main/benchmarking/mapping/linear_sarek_grch38_samtools-stats-dp.tsv).
  - Export the `samtools_alignment_plot-1`. See examples [`linear_grch38_samtools_alignment_plot-1.tsv`](https://github.com/RAVING-Informatics/T2T-alignment/blob/main/benchmarking/mapping/linear_grch38_samtools_alignment_plot-1.tsv) and [`linear_sarek_chm13_samtools-stats-dp.tsv`](https://github.com/RAVING-Informatics/T2T-alignment/blob/main/benchmarking/mapping/linear_sarek_chm13_samtools-stats-dp.tsv)
- Use the script: `extract-mismatch-rate.sh` to create TSV file with the mismatch rates derived from the `samtools stats` files. 
- Manipulate data from all exports in excel to calculate mapping statistics for each sample. See excel template [`linear_t2t-hg38_mapping_calculations.xlsx`](https://github.com/RAVING-Informatics/T2T-alignment/blob/main/benchmarking/mapping/linear_t2t-hg38_mapping_calculations.xlsx) and [`linear_t2t-hg38_mapping.tsv`](https://github.com/RAVING-Informatics/T2T-alignment/blob/main/benchmarking/mapping/linear_t2t-hg38_mapping.csv) for result. 
  - Proportion reads unmapped = reads unmapped / total reads
  - Proportion reads mapped = reads mapped / total reads
  - Proportion reads MQ0 = MQ0 reads / total reads
  - Proportion properly paired = reads properly paired bit set / total reads
  - Proportion misoriented reads = read pairs with other orientation / total reads
- Use the `plot-mapping-stats-linear-hprc.R` script to generate a multifaceted plot of the key mapping statistics. Note, this script also contains data from pangenome realignment. Input data containing pangenome data is available [here](https://github.com/RAVING-Informatics/T2T-alignment/blob/main/benchmarking/mapping/linear-hprc_t2t-hg38_mapping.csv).
