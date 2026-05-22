# SV Annotation

## AnnotSV 
- As of Jul 2025, AnnotSV supports T2T-CHM13-based VCFs.
- To use AnnotSV with CHM13-based VCFs, ensure the version is >3.5.
- Once AnnotSV is installed as per the [documentation](https://github.com/lgmgeo/AnnotSV/blob/master/share/doc/AnnotSV/quickstart.md), the CHM13-based BED annotations need to be prepared via liftover. 
    - Source the chain file:
```
wget https://hgdownload.soe.ucsc.edu/gbdb/hg38/liftOver/hg38ToGCA_009914755.4.over.chain.gz
gunzip hg38ToGCA_009914755.4.over.chain.gz
```
    - `liftover_script.sh` - adapted from script from AnnotSV [`lift_Over_a_BED_file.tcl`](https://github.com/lgmgeo/AnnotSV/blob/master/share/tcl/AnnotSV/Scripts/lift_Over_a_BED_file.tcl)
    - Use `liftover_annotsv_chm13.sh` to loop over relevent files.

**Issues with AnnotSV**
- Does not appear to match SVs based on SVTYPE (i.e. INS, DEL, INV). Only matches based on co-ordinate overlaps.
- Cannot get AnnotSV to output to VCF (issue with python incompatibility).

## SVAFotate
- Tool designed to annotate SVs with population level AFs and other related metrics.
- Allows for matching with SVTYPE.
- Outputs naturally to VCF.
- Does not currently support CHM13
- Highly customisable

**Liftover SVAFotate BED**
A BED file of SV annotations is available [here](https://zenodo.org/records/11642574). To liftover the BED file to T2T-CHM13 co-ordinates, use the following scripts:
- `liftover_svfatotate_chm13.sh` - uses `liftover_script.sh` to liftover BED co-ordinates.

These are the results of the liftover:
```
----------------------------------
Total regions:    3455239
Lifted regions:   3307338
Unmapped regions: 147901
----------------------------------
```

## Other SV AF Annotation Sources
**HPRC-HGSV**
These files were constructed with minigraph-cactus v2.7.2 from the HGSVC samples, HPRCv1, CHM13 and GRCh38 (216 haplotypes total).
[Link to original paper](https://www.nature.com/articles/s41586-025-09140-6)
Source the VCF files: 
```
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Graph_Genomes/1.0/2024_02_23_minigraph_cactus_hgsvc3_hprc/hgsvc3-hprc-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Graph_Genomes/1.0/2024_02_23_minigraph_cactus_hgsvc3_hprc/hgsvc3-hprc-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz.tbi
```
**HPRCv2.0**
[Link to the original paper](https://www.nature.com/articles/s41586-023-05896-x) 
More information on the VCFs available [here](https://42basepairs.com/browse/s3/human-pangenomics/pangenomes/scratch/2025_02_28_minigraph_cactus).
Source of the VCF files:
```
wget https://42basepairs.com/download/s3/human-pangenomics/pangenomes/scratch/2025_02_28_minigraph_cactus/hprc-v2.0-mc-chm13/hprc-v2.0-mc-chm13.wave.vcf.gz
wget https://42basepairs.com/download/s3/human-pangenomics/pangenomes/scratch/2025_02_28_minigraph_cactus/hprc-v2.0-mc-chm13/hprc-v2.0-mc-chm13.wave.vcf.gz.tbi
```

**Generate BED file for SVAFotate**
To create a BED file compatible with running SVAFotate using the VCFs above, use the following script `vcf_bed.sh`. This will also normalise the chromosome prefix, and recompute the AF for the HGSVC data. 
To filter the BED file to include only unique SVs (>50bp), run `filter_bed.sh`.
To merge the BED files with the original SVAFotate BED file, run `merge_bed.sh`.

## Run SVAFotate
- First, extract out VCF file containing only the DYSGU variant calls: `extract_sv_chr.sh`. This further splits the VCF by chromosome to avoid OOM errors when running SVAFotate.
- Run SVAFotate on each chromosome using `svafotate.sh`:
```
for chr in chr{1..22}; do sbatch svafotate.sh $chr ; done
for chr in chrY chrX chrM; do sbatch svafotate.sh $chr ; done
```
