#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH -c 4                          
#SBATCH --job-name=filter_bed
#SBATCH --partition=work
#SBATCH --account=pawsey0933
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out
#SBATCH --export=ALL

#IN_BED=/scratch/pawsey0933/cfolland/t2t/batch2/annotation/hgsvc3-hprc-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.biallelic.af.sorted.bed.gz
#OUT_BED=/scratch/pawsey0933/cfolland/t2t/batch2/annotation/hgsvc3-hprc-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.biallelic.af.sorted.len50.bed.gz

IN_BED=/scratch/pawsey0933/cfolland/t2t/batch2/annotation/hprc-v2.0-mc-chm13.wave.biallelic.sorted.bed.gz
OUT_BED=/scratch/pawsey0933/cfolland/t2t/batch2/annotation/hprc-v2.0-mc-chm13.wave.biallelic.sorted.len50.bed.gz

zcat "$IN_BED" | awk '
BEGIN { 
    before=0; 
    after=0; 
    unique=0;
    OFS="\t" 
}

{
    before++

    # Apply your filters
    if ($4 > 50 && $8+0 > 0) {

        # Remove duplicate rows (based on full line)
        if (!seen[$0]++) {
            print
            unique++
        }

        after++
    }
}

END {
    print "Total variants before filtering: " before > "/dev/stderr"
    print "Total variants passing filters: " after > "/dev/stderr"
    print "Total unique variants after deduplication: " unique > "/dev/stderr"
}
' | gzip > "$OUT_BED"