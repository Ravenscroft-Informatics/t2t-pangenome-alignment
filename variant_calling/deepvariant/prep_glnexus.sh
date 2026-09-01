#!/bin/bash -l
# One-time prep: derive per-chromosome BED files + an ordered contig list
# for scattering GLnexus joint genotyping across a SLURM array, instead of
# one long serial glnexus_cli run over the whole genome (the previous
# attempt timed out at chr17 after 72h on --partition=highmem --mem=440G).
#
# Run ONCE (from a login or workflow node), after build_karyotypic_dict.sh
# and before submitting glnexus_array.sh.

set -euo pipefail

merged_dir=/scratch/pawsey0933/cfolland/deepvariant/hprc/pg_dv/merged
dict=${merged_dir}/hprc-v2.0-mc-grch38.karyotypic.dict
out_dir=${merged_dir}/glnexus_regions
contigs_file=${out_dir}/contigs.txt
other_bed=${out_dir}/other.bed

if [[ ! -s "${dict}" ]]; then
    echo "ERROR: dictionary not found: ${dict}" >&2
    echo "Run ./build_karyotypic_dict.sh first." >&2
    exit 1
fi

mkdir -p "${out_dir}"
: > "${contigs_file}"
: > "${other_bed}"

# sn<TAB>ln for every contig in the dictionary.
sq_file=$(mktemp)
trap 'rm -f "${sq_file}"' EXIT
awk -F'\t' '
    $1 == "@SQ" {
        sn = ""; ln = ""
        for (i = 2; i <= NF; i++) {
            if ($i ~ /^SN:/) sn = substr($i, 4)
            if ($i ~ /^LN:/) ln = substr($i, 4)
        }
        if (sn != "" && ln != "") print sn "\t" ln
    }
' "${dict}" > "${sq_file}"

declare -A length_of
while IFS=$'\t' read -r sn ln; do
    length_of["${sn}"]=${ln}
done < "${sq_file}"

# Primary chromosomes: own whole-chromosome BED + own line in contigs.txt,
# in chr1..chr22,chrX,chrY order (array task N = line N of contigs.txt).
for n in $(seq 1 22) X Y; do
    contig="chr${n}"
    if [[ -z "${length_of[${contig}]+x}" ]]; then
        echo "WARNING: ${contig} not found in ${dict}, skipping" >&2
        continue
    fi
    printf '%s\t0\t%s\n' "${contig}" "${length_of[${contig}]}" > "${out_dir}/${contig}.bed"
    echo "${contig}" >> "${contigs_file}"
done

# Everything else in the dictionary (alt/decoy/unplaced contigs) bundled
# into one "other" region so nothing in the merged gVCFs is silently
# dropped from joint genotyping.
while IFS=$'\t' read -r sn ln; do
    [[ -s "${out_dir}/${sn}.bed" ]] && continue
    printf '%s\t0\t%s\n' "${sn}" "${ln}" >> "${other_bed}"
done < "${sq_file}"
if [[ -s "${other_bed}" ]]; then
    echo "other" >> "${contigs_file}"
fi

echo "Wrote $(wc -l < "${contigs_file}") region(s) to ${out_dir}" >&2