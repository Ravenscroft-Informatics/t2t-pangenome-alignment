import pandas as pd
import argparse
import re

def parse_vcf(input_vcf, output_file, output_format='tsv'):
    # VEP headers as provided
    vep_headers = [
        "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE",
        "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position", "Protein_position",
        "Amino_acids", "Codons", "Existing_variation", "DISTANCE", "STRAND", "FLAGS", "VARIANT_CLASS",
        "SYMBOL_SOURCE", "HGNC_ID", "CANONICAL", "MANE_SELECT", "MANE_PLUS_CLINICAL", "TSL", "APPRIS",
        "CCDS", "ENSP", "SWISSPROT", "TREMBL", "UNIPARC", "UNIPROT_ISOFORM", "SOURCE", "GENE_PHENO",
        "DOMAINS", "miRNA", "HGVS_OFFSET", "AF", "AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF",
        "gnomADe_AF", "gnomADe_AFR_AF", "gnomADe_AMR_AF", "gnomADe_ASJ_AF", "gnomADe_EAS_AF",
        "gnomADe_FIN_AF", "gnomADe_NFE_AF", "gnomADe_OTH_AF", "gnomADe_SAS_AF", "gnomADg_AF",
        "gnomADg_AFR_AF", "gnomADg_AMI_AF", "gnomADg_AMR_AF", "gnomADg_ASJ_AF", "gnomADg_EAS_AF",
        "gnomADg_FIN_AF", "gnomADg_MID_AF", "gnomADg_NFE_AF", "gnomADg_OTH_AF", "gnomADg_SAS_AF",
        "MAX_AF", "MAX_AF_POPS", "FREQS", "CLIN_SIG", "SOMATIC", "PHENO", "PUBMED", "MOTIF_NAME",
        "MOTIF_POS", "HIGH_INF_POS", "MOTIF_SCORE_CHANGE", "TRANSCRIPTION_FACTORS", "clinvar",
        "clinvar_CLNSIG", "clinvar_CLNREVSTAT", "clinvar_CLNDN"
    ]

    # Fields to extract from the INFO column
    info_fields = ["AF", "AQ", "AC", "AN", "AF_HPRC", "AF_1000G", "AF_gnomad", "AF_hgsvc3", "UniqueT2T", "SVTYPE", "END", "CHR2", "NGRP", "CT", "CIEND95", "CIPOS95", "SVLEN", "GC", "NEXP", "STRIDE", "RPOLY", "OL", "SU", "WR", "PE", "SR", "SC", "BND", "LPREC", "RT", "MeanPROB", "MaxPROB"]

    # Open the VCF file and process
    with open(input_vcf, 'r') as f:
        lines = f.readlines()

    # Filter out metadata lines
    metadata_lines = [line for line in lines if line.startswith('#')]
    vcf_data = [line.strip() for line in lines if not line.startswith('#')]

    # Extract sample column names from the VCF header
    headers = metadata_lines[-1].strip().split('\t')
    sample_columns = headers[9:]  # Sample names start at the 10th column

    # Initialize a list for parsed rows
    parsed_rows = []

    for row in vcf_data:
        fields = row.split('\t')

        # Ensure there are at least 9 fields
        if len(fields) < 9:
            raise ValueError(f"Invalid VCF row format. Expected at least 9 columns, got {len(fields)}: {row}")

        # Extract fixed fields
        chrom, pos, v_id, ref, alt, qual, filt, info, fmt = fields[:9]
        sample_data = fields[9:]  # Sample genotype data

        # Parse INFO field
        info_dict = {key.split('=')[0]: key.split('=', 1)[1] for key in info.split(';') if '=' in key}
        info_present = {key.split('=')[0]: True for key in info.split(';') if '=' not in key}  # Handle presence-only fields

        # Extract specified INFO fields, setting UniqueT2T to 'Yes' if present
        extracted_info = {field: info_dict.get(field, None) for field in info_fields}
        extracted_info['UniqueT2T'] = 'Yes' if 'UniqueT2T' in info_present or 'UniqueT2T' in info_dict else None

        # Parse CSQ field
        csq_data = info_dict.get('CSQ', '').split(',')
        genetic_models = info_dict.get('GeneticModels', '')  # Preserve full GeneticModels field

        # Expand CSQ entries
        for csq_entry in csq_data:
            csq_fields = csq_entry.split('|')
            csq_dict = dict(zip(vep_headers, csq_fields))

            # Add fixed fields to the row
            csq_dict['CHROM'] = chrom
            csq_dict['POS'] = pos
            csq_dict['ID'] = v_id
            csq_dict['REF'] = ref
            csq_dict['ALT'] = alt
            csq_dict['QUAL'] = qual
            csq_dict['FILTER'] = filt
            csq_dict['FORMAT'] = fmt  # Add FORMAT column

            # Add GeneticModels, INFO fields, and sample data
            csq_dict['GeneticModels'] = genetic_models
            csq_dict.update(extracted_info)
            for sample, data in zip(sample_columns, sample_data):
                csq_dict[sample] = data

            parsed_rows.append(csq_dict)

    # Convert to DataFrame
    df = pd.DataFrame(parsed_rows)

    # Reorder columns: Fixed fields, sample columns, INFO fields, then the rest
    fixed_columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'FORMAT']
    all_columns = fixed_columns + sample_columns + info_fields + [col for col in df.columns if col not in fixed_columns + sample_columns + info_fields ]
    df = df[all_columns]

    # Write to file
    if output_format == 'tsv':
        df.to_csv(output_file, sep='\t', index=False)
    else:
        df.to_csv(output_file, index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Parse a VCF file and extract INFO fields")
    parser.add_argument("input_vcf", help="Input VCF file")
    parser.add_argument("output_file", help="Output file (TSV or CSV)")
    parser.add_argument("--output_format", choices=['tsv', 'csv'], default='tsv', help="Output file format (default: TSV)")
    args = parser.parse_args()

    parse_vcf(args.input_vcf, args.output_file, output_format=args.output_format)
