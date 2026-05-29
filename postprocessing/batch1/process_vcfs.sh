#!/bin/bash

# Directory containing the VCF subdirectories
INPUT_DIR="/data/t2t/families"

# Directory where the output files will be stored
OUTPUT_DIR="/data/t2t/families_parsed"

# Path to the Python script
PARSE_SCRIPT="/data/t2t/scripts/parse.py"

# Output format (e.g., tsv or csv)
OUTPUT_FORMAT="csv"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through each subdirectory
find "$INPUT_DIR" -type f \( -name "*.vcf" -o -name "*.vcf.gz" \) | while read -r VCF_FILE; do
    # Get the relative path of the VCF file (relative to INPUT_DIR)
    RELATIVE_PATH=$(dirname "${VCF_FILE#$INPUT_DIR/}")
    
    # Create the corresponding subdirectory in the output directory
    OUTPUT_SUBDIR="$OUTPUT_DIR/$RELATIVE_PATH"
    mkdir -p "$OUTPUT_SUBDIR"
    
    # Get the prefix of the VCF file (without extension)
    PREFIX=$(basename "$VCF_FILE" | sed -E 's/\.(vcf|vcf\.gz)$//')
    
    # Construct the output file path
    OUTPUT_FILE="$OUTPUT_SUBDIR/${PREFIX}_parsed.$OUTPUT_FORMAT"
    
    # Check if the file is gzipped and decompress inline
    if [[ "$VCF_FILE" == *.vcf.gz ]]; then
        echo "Processing compressed file: $VCF_FILE -> $OUTPUT_FILE"
        zcat "$VCF_FILE" | python "$PARSE_SCRIPT" /dev/stdin "$OUTPUT_FILE" --output_format "$OUTPUT_FORMAT"
    else
        echo "Processing uncompressed file: $VCF_FILE -> $OUTPUT_FILE"
        python "$PARSE_SCRIPT" "$VCF_FILE" "$OUTPUT_FILE" --output_format "$OUTPUT_FORMAT"
    fi
done
