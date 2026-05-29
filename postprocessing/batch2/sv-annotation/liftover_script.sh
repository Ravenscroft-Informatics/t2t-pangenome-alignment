#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 3 ]; then
    echo "USAGE: $0 input.bed(.gz) output.bed chainfile"
    exit 1
fi

IN="$1"
OUT="$2"
CHAIN="$3"

TMPDIR=${TMPDIR:-$MYSCRATCH/liftover_tmp}
mkdir -p "$TMPDIR"
WORK=$(mktemp -d "$TMPDIR/run.XXXXXX")

echo "Working in: $WORK"

# -----------------------------
# Step 0: extract header
# -----------------------------
echo "...extracting header"

if [[ "$IN" == *.gz ]]; then
    zcat "$IN" | grep '^#' > "$WORK/header.txt" || true
else
    grep '^#' "$IN" > "$WORK/header.txt" || true
fi

# -----------------------------
# Step 1: prepare input + IDs
# -----------------------------
echo "...preparing BED + IDs"

if [[ "$IN" == *.gz ]]; then
    zcat "$IN"
else
    cat "$IN"
fi \
| grep -v '^#' \
| awk -F'\t' 'BEGIN{OFS="\t"} NF>2 {
    id=NR
    coords = $1"\t"$2"\t"$3

    ann=""
    for(i=4;i<=NF;i++) ann = ann (i==4?"":"\t") $i

    print coords, id, ann
}' > "$WORK/with_id.bed"

# -----------------------------
# coords with chr prefix
# -----------------------------
cut -f1-4 "$WORK/with_id.bed" \
| awk 'BEGIN{OFS="\t"} {
    if ($1 !~ /^chr/) $1="chr"$1;
    if ($1=="chrMT") $1="chrM";
    print
}' > "$WORK/coords.bed"

# -----------------------------
# store annotations (FIXED)
# -----------------------------
SOURCE=$(basename "$CHAIN" | sed -E 's/^([^_-]+).*/\1/')

awk -F'\t' -v OFS='\t' -v src="$SOURCE" '
{
    id = $4

    ann=""
    for(i=5;i<=NF;i++) ann = ann (i==5?"":"\t") $i

    # prefix embedded coordinates
    gsub(/[0-9XYchrMT]+:[0-9]+-[0-9]+/, src"_" "&", ann)

    print id, ann
}
' "$WORK/with_id.bed" > "$WORK/annotations.tsv"

# -----------------------------
# Step 2: liftover
# -----------------------------
echo "...running liftOver"

liftOver \
    "$WORK/coords.bed" \
    "$CHAIN" \
    "$WORK/lifted.bed" \
    "$WORK/unmapped.bed"

# -----------------------------
# Step 3: merge annotations
# -----------------------------
echo "...re-attaching annotations"

sort -k4,4 "$WORK/lifted.bed" > "$WORK/lifted.sorted"
sort -k1,1 "$WORK/annotations.tsv" > "$WORK/annotations.sorted"

join -1 4 -2 1 -t $'\t' \
    "$WORK/lifted.sorted" \
    "$WORK/annotations.sorted" \
| cut -f2- \
> "$WORK/merged.bed"

# -----------------------------
# Step 4: clean + enforce columns
# -----------------------------
echo "...cleaning"

awk -v OFS="\t" '{$1=$1; print}' "$WORK/merged.bed" > "$WORK/clean.bed"

EXPECTED=$(head -1 "$WORK/clean.bed" | awk '{print NF}')

awk -v n="$EXPECTED" 'NF==n' "$WORK/clean.bed" > "$WORK/filtered.bed"

# -----------------------------
# Step 5: contig → chromosome
# -----------------------------
echo "...converting chromosomes"

awk 'BEGIN{OFS="\t"}
{
    if ($1=="CP068277.2") $1="1";
    else if ($1=="CP068276.2") $1="2";
    else if ($1=="CP068275.2") $1="3";
    else if ($1=="CP068274.2") $1="4";
    else if ($1=="CP068273.2") $1="5";
    else if ($1=="CP068272.2") $1="6";
    else if ($1=="CP068271.2") $1="7";
    else if ($1=="CP068270.2") $1="8";
    else if ($1=="CP068269.2") $1="9";
    else if ($1=="CP068268.2") $1="10";
    else if ($1=="CP068267.2") $1="11";
    else if ($1=="CP068266.2") $1="12";
    else if ($1=="CP068265.2") $1="13";
    else if ($1=="CP068264.2") $1="14";
    else if ($1=="CP068263.2") $1="15";
    else if ($1=="CP068262.2") $1="16";
    else if ($1=="CP068261.2") $1="17";
    else if ($1=="CP068260.2") $1="18";
    else if ($1=="CP068259.2") $1="19";
    else if ($1=="CP068258.2") $1="20";
    else if ($1=="CP068257.2") $1="21";
    else if ($1=="CP068256.2") $1="22";
    else if ($1=="CP068255.2") $1="X";
    else if ($1=="CP086569.2") $1="Y";
    else if ($1=="CP068254.1") $1="MT";
    print
}' "$WORK/filtered.bed" > "$WORK/final.bed"

# -----------------------------
# Step 6: final sort
# -----------------------------
sort -k1,1V -k2,2n "$WORK/final.bed" > "$WORK/final.sorted.bed"

# -----------------------------
# Step 7: header
# -----------------------------
if [ -s "$WORK/header.txt" ]; then
    cat "$WORK/header.txt" "$WORK/final.sorted.bed" > "$OUT"
else
    cp "$WORK/final.sorted.bed" "$OUT"
fi

# -----------------------------
# Step 8: stats
# -----------------------------
TOTAL=$(wc -l < "$WORK/coords.bed")
LIFTED=$(wc -l < "$WORK/lifted.bed")
UNMAPPED=$(grep -vc '^#' "$WORK/unmapped.bed" || true)

echo "----------------------------------"
echo "Total regions:    $TOTAL"
echo "Lifted regions:   $LIFTED"
echo "Unmapped regions: $UNMAPPED"
echo "----------------------------------"

cp "$WORK/unmapped.bed" "${OUT}.unmapped"

rm -rf "$WORK"

echo "✅ Liftover complete: $OUT"