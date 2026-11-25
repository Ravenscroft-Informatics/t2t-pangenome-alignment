#!/bin/bash

#SBATCH --job-name=exons
#SBATCH --account=pawsey0933
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --time=0:15:00
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out
#SBATCH --export=ALL

##chm13
ref=chm13
gtf=/software/projects/pawsey0933/benchmarking/nmd_genes/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz
region_file=/software/projects/pawsey0933/benchmarking/nmd_genes/nmd_gene_list_chm13.bed
mapfile=/software/projects/pawsey0933/benchmarking/nmd_genes/refseq-chr_chm13.tsv
##hg38
#ref=grch38
#gtf=/software/projects/pawsey0933/benchmarking/nmd_genes/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz
#region_file=/software/projects/pawsey0933/benchmarking/nmd_genes/nmd_gene_list_grch38.bed
#mapfile=/software/projects/pawsey0933/benchmarking/nmd_genes/refseq-chr_hg38.tsv

source /software/projects/pawsey0933/cfolland/miniconda3/etc/profile.d/mamba.sh

# --- Config (you can override via env or edit) ---
ps=${ps:-/software/projects/pawsey0933/benchmarking/nmd_genes/gtf_region_to_tables.py}
gene=${1:-FRG2}

die(){ echo "ERROR: $*" >&2; exit 1; }

[[ -f "$ps" ]] || die "Python script not found: $ps"
[[ -f "$gtf" ]] || die "GTF not found: $gtf"
[[ -f "$region_file" ]] || die "Region BED not found: $region_file"

# Detect naming scheme of the GTF (first non-comment seqname)
gtf_scheme=$(
  zgrep -v '^#' "$gtf" | head -n1 | awk -F'\t' '{print $1}' | awk '
    /^NC_/ {print "NC"; exit}
    /^chr/ {print "CHR"; exit}
    {print "PLAIN"; exit}
  '
)
echo "GTF naming scheme detected: $gtf_scheme"

# Pull BED interval for gene (accept gene name in col4 OR col5)
read -r bed_chr bed_start bed_end <<<"$(
  awk -v FS='\t' -v g="$gene" '
    $4==g || $5==g { print $1, $2, $3; found=1; exit }
    END{ if(!found) exit 1 }
  ' "$region_file" || true
)"
[[ -n "${bed_chr:-}" ]] || die "Gene '$gene' not found in $region_file"

# Convert BED (0-based, half-open) -> region (1-based, closed)
start_1b=$(( bed_start + 1 ))
end_1b=$(( bed_end ))
(( start_1b < end_1b )) || die "Start ($start_1b) must be < End ($end_1b)"

echo "BED interval:   ${bed_chr}:${bed_start}-${bed_end}   (0-based, half-open)"

# Normalize a chromosome token without 'chr'
nochr="${bed_chr#chr}"

# Determine target contig name to match the GTF
contig=""
case "$gtf_scheme" in
  NC)
    # Need NC_* contig; map from {1..22,X,Y} -> NC_* using mapfile
    [[ -f "$mapfile" ]] || die "Need $mapfile to map '$bed_chr' -> NC_* for NC-named GTF."
    contig=$(awk -v FS='\t' -v k="$nochr" '$1==k{print $2; exit}' "$mapfile")
    [[ -n "$contig" ]] || die "Could not map '$bed_chr' (key '$nochr') via $mapfile"
    ;;
  CHR)
    # Ensure chr prefix
    [[ "$bed_chr" =~ ^chr ]] && contig="$bed_chr" || contig="chr$nochr"
    ;;
  PLAIN)
    contig="$nochr"
    ;;
  *)
    die "Unrecognized GTF scheme: $gtf_scheme"
    ;;
esac

echo "Resolved contig: $contig"

# Construct region string (chrom:start-end)
region="${contig}:${start_1b}-${end_1b}"

# Validate region format strictly
if ! [[ "$region" =~ ^(NC_[0-9]+\.[0-9]+|chr([0-9]{1,2}|1[0-9]|2[0-2]|X|Y)|[0-9]{1,2}|X|Y):[0-9]+-[0-9]+$ ]]; then
  die "Region looks malformed: '$region'"
fi

echo "Region (1-based): $region"
echo "GTF:              $gtf"
echo "Running python…"

python3 "$ps" \
  --gtf "$gtf" \
  --region "$region" \
  --exons_out "${gene}.${ref}.exons.tsv" \
  --introns_out "${gene}.${ref}.introns.tsv"

echo "Done."
