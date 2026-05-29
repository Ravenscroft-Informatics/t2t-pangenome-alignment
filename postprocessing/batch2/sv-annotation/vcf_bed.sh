#!/bin/bash -l

#SBATCH --job-name=bed_vcf
#SBATCH --account=pawsey0933
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out
#SBATCH --export=ALL

module load bcftools/1.15--haf5b3da_0

annotate=/scratch/pawsey0933/cfolland/t2t/batch2/annotation/
ref=/software/projects/pawsey0933/pangenome/refs/chm13v2.0.maskedY.rCRS.EBV.fasta
prefix=hgsvc3-hprc-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm
#prefix=hprc-v2.0-mc-chm13.wave

cd $annotate

# PREPROCESSING 
# convert multiallelic to biallelic, recompute AF and sort vcf
bcftools norm -m -any --check-ref skip -f $ref -o $annotate/$prefix.biallelic.vcf.gz $annotate/$prefix.vcf.gz
bcftools annotate -x INFO/AF $annotate/$prefix.biallelic.vcf.gz -Oz -o $annotate/tmp.vcf.gz #only for hgsvc vcf
bcftools +fill-tags $annotate/tmp.vcf.gz -Oz -o $dir/$prefix.biallelic.af.vcf.gz -- -t AF #only for hgsvc vcf
bcftools sort -T $annotate -o $annotate/$prefix.biallelic.af.sorted.vcf.gz $annotate/$prefix.biallelic.af.vcf.gz
tabix -p vcf $annotate/$prefix.biallelic.af.sorted.vcf.gz

## CONVERT VCF TO BED
# NOTE, this script uses the following logic:
# CHROM → column 1
# POS → START (VCF POS is 1-based; BED is 0-based, so subtract 1 if you want strict BED)
# END → compute from POS + SVLEN (or from REF/ALT lengths if SVLEN not present)
# SVLEN → from INFO (LEN= or derived) (absolute value)
# SVTYPE → from INFO (TYPE= if present, otherwise infer)
# SOURCE → fixed label (e.g., your caller name) -> change "HPRCv2_wave" to fit the source of the VCF
# SV_ID → VCF ID column
# AF → from INFO (AF=)
# Other columns - added to be compatible with SVAFotate

#INPUT=hgsvc3-hprc-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.biallelic.af.sorted.vcf.gz
#OUTPUT=hgsvc3-hprc-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.biallelic.af.sorted.bed.gz
#SOURCE="HGVSC"

INPUT=hprc-v2.0-mc-chm13.wave.biallelic.sorted.vcf.gz
OUTPUT=hprc-v2.0-mc-chm13.wave.biallelic.sorted.bed.gz
SOURCE=HPRCv2_wave


LC_ALL=C zcat "$INPUT" | awk -v source="$SOURCE" '
BEGIN {
    OFS="\t";

    header = "#CHROM\tSTART\tEND\tSVLEN\tSVTYPE\tSOURCE\tSV_ID\tAF\tHomRef\tHet\tHomAlt\tMale_AF\tMale_HomRef\tMale_Het\tMale_HomAlt\tMale_HemiAlt\tMale_HemiAF\tFemale_AF\tFemale_HomRef\tFemale_Het\tFemale_HomAlt\tAFR_AF\tAFR_HomRef\tAFR_Het\tAFR_HomAlt\tAFR_Male_AF\tAFR_Male_HomRef\tAFR_Male_Het\tAFR_Male_HomAlt\tAFR_Male_HemiAlt\tAFR_Male_HemiAF\tAFR_Female_AF\tAFR_Female_HomRef\tAFR_Female_Het\tAFR_Female_HomAlt\tAMI_AF\tAMI_HomRef\tAMI_Het\tAMI_HomAlt\tAMI_Male_AF\tAMI_Male_HomRef\tAMI_Male_Het\tAMI_Male_HomAlt\tAMI_Male_HemiAlt\tAMI_Male_HemiAF\tAMI_Female_AF\tAMI_Female_HomRef\tAMI_Female_Het\tAMI_Female_HomAlt\tAMR_AF\tAMR_HomRef\tAMR_Het\tAMR_HomAlt\tAMR_Male_AF\tAMR_Male_HomRef\tAMR_Male_Het\tAMR_Male_HomAlt\tAMR_Male_HemiAlt\tAMR_Male_HemiAF\tAMR_Female_AF\tAMR_Female_HomRef\tAMR_Female_Het\tAMR_Female_HomAlt\tASJ_AF\tASJ_HomRef\tASJ_Het\tASJ_HomAlt\tASJ_Male_AF\tASJ_Male_HomRef\tASJ_Male_Het\tASJ_Male_HomAlt\tASJ_Male_HemiAlt\tASJ_Male_HemiAF\tASJ_Female_AF\tASJ_Female_HomRef\tASJ_Female_Het\tASJ_Female_HomAlt\tEAS_AF\tEAS_HomRef\tEAS_Het\tEAS_HomAlt\tEAS_Male_AF\tEAS_Male_HomRef\tEAS_Male_Het\tEAS_Male_HomAlt\tEAS_Male_HemiAlt\tEAS_Male_HemiAF\tEAS_Female_AF\tEAS_Female_HomRef\tEAS_Female_Het\tEAS_Female_HomAlt\tEUR_AF\tEUR_HomRef\tEUR_Het\tEUR_HomAlt\tEUR_Male_AF\tEUR_Male_HomRef\tEUR_Male_Het\tEUR_Male_HomAlt\tEUR_Male_HemiAlt\tEUR_Male_HemiAF\tEUR_Female_AF\tEUR_Female_HomRef\tEUR_Female_Het\tEUR_Female_HomAlt\tFIN_AF\tFIN_HomRef\tFIN_Het\tFIN_HomAlt\tFIN_Male_AF\tFIN_Male_HomRef\tFIN_Male_Het\tFIN_Male_HomAlt\tFIN_Male_HemiAlt\tFIN_Male_HemiAF\tFIN_Female_AF\tFIN_Female_HomRef\tFIN_Female_Het\tFIN_Female_HomAlt\tMID_AF\tMID_HomRef\tMID_Het\tMID_HomAlt\tMID_Male_AF\tMID_Male_HomRef\tMID_Male_Het\tMID_Male_HomAlt\tMID_Male_HemiAlt\tMID_Male_HemiAF\tMID_Female_AF\tMID_Female_HomRef\tMID_Female_Het\tMID_Female_HomAlt\tNFE_AF\tNFE_HomRef\tNFE_Het\tNFE_HomAlt\tNFE_Male_AF\tNFE_Male_HomRef\tNFE_Male_Het\tNFE_Male_HomAlt\tNFE_Male_HemiAlt\tNFE_Male_HemiAF\tNFE_Female_AF\tNFE_Female_HomRef\tNFE_Female_Het\tNFE_Female_HomAlt\tOTH_AF\tOTH_HomRef\tOTH_Het\tOTH_HomAlt\tOTH_Male_AF\tOTH_Male_HomRef\tOTH_Male_Het\tOTH_Male_HomAlt\tOTH_Male_HemiAlt\tOTH_Male_HemiAF\tOTH_Female_AF\tOTH_Female_HomRef\tOTH_Female_Het\tOTH_Female_HomAlt\tSAS_AF\tSAS_HomRef\tSAS_Het\tSAS_HomAlt\tSAS_Male_AF\tSAS_Male_HomRef\tSAS_Male_Het\tSAS_Male_HomAlt\tSAS_Male_HemiAlt\tSAS_Male_HemiAF\tSAS_Female_AF\tSAS_Female_HomRef\tSAS_Female_Het\tSAS_Female_HomAlt\tPopMax_AF\tInPop";

    print header;

    n_header = split(header, h, "\t");
    fixed_cols = 8;
    extra_cols = n_header - fixed_cols;
}

!/^#/ {
    chrom=$1; pos=$2; id=$3; ref=$4; alt=$5; info=$8;

    # Normalize chromosome
    sub(/^chr/, "", chrom);
    if (chrom == "M") chrom = "MT";

    start = pos - 1;

    # Fix ID: take first if multiple, replace > with -
    split(id, ids, ";");
    id = ids[1];
    gsub(/>/, "-", id);

    svlen=""; svtype=""; af="";
    n = split(info, fields, ";");

    for (i=1; i<=n; i++) {
        if (fields[i] ~ /^LEN=/)  { split(fields[i], a, "="); svlen=a[2]; }
        if (fields[i] ~ /^TYPE=/) { split(fields[i], a, "="); svtype=a[2]; }
        if (fields[i] ~ /^AF=/)   { split(fields[i], a, "="); af=a[2]; }
    }

    if (svlen == "") svlen = length(alt) - length(ref);
    if (svlen < 0) svlen = -svlen;

    end = start + svlen;

    if (svtype == "") {
        if (length(alt) > length(ref)) svtype="INS";
        else if (length(alt) < length(ref)) svtype="DEL";
        else svtype="MNP";
    }

    # Build output safely using array (robust fix)
    n_out = 0;

    out[++n_out] = chrom;
    out[++n_out] = start;
    out[++n_out] = end;
    out[++n_out] = svlen;
    out[++n_out] = svtype;
    out[++n_out] = source;
    out[++n_out] = id;
    out[++n_out] = (af == "" ? "NA" : af);

    for (i=1; i<=extra_cols; i++) {
        out[++n_out] = "NA";
    }

    for (i=1; i<=n_out; i++) {
        printf "%s%s", out[i], (i==n_out ? "\n" : OFS);
    }
}
' | gzip > "$OUTPUT"