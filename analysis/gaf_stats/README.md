# Calculate GAF alignment statistics

The python script `gaf_alignment_stats.py` parses reads in a GAF file and counts the number of:
- Unmapped reads (path == '*')
- Low MAPQ reads (mapped, MAPQ < --min-mapq)
- MAPQ 0 reads (mapped)
- Properly paired reads (giraffe pd:b:1 tag)
- Perfect reads (whole read aligned (no soft-clip), no mismatches/indels vs its graph path: qstart==0, qend==qlen, matches==block==qlen))
- Reference-only reads (every node is on the reference path AND every node-to-node step follows the reference path (read is fully explained by GRCh38))
- Non-reference edge reads (every node is on the reference path but at least one step is not a reference edge (deletion / rearrangement vs GRCh38))
- Non-reference node reads (at least one node is absent from the reference path (insertion / SNV allele / divergent sequence from other haplotypes))

The script also calculates the 
- mismatch_rate (mismatches / (matches + mismatches), from the cs tag)
- error_rate (mismatches + inserted + deleted bases) / (matches + mismatches + inserted bases) -- same definition as samtools stats "error rate" (NM / bases mapped (cigar))