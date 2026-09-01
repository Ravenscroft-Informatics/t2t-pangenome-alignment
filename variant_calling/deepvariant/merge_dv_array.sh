#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH -c 6
#SBATCH --job-name=merge_dv_array
#SBATCH --partition=work
#SBATCH --account=pawsey0933
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=END
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --array=1-234%20

# ---------------------------------------------------------------------------
# One-time prep (run ONCE before the first sbatch submission, from a login
# or workflow node):
#
#   ./build_karyotypic_dict.sh
#   ./generate_sample_list.sh
#
# build_karyotypic_dict.sh writes a reordered copy of the reference
# sequence dictionary (@SQ lines in chr1..chr22,chrX,chrY order) to
# ${output_dir}/hprc-v2.0-mc-grch38.karyotypic.dict. The ORIGINAL dict's
# contig order is NOT karyotypic (chr10 sorts before chr1 in it), and that
# same non-karyotypic order is baked into every per-region VCF's own header
# (the CRAMs/gVCFs were called against the original dict). Because of that:
#   - MergeVcfs is run against the ORIGINAL dict (src_dict below), since
#     that's what's actually embedded in the per-region VCF headers.
#   - The merged-but-still-wrong-order output is then reheadered with
#     UpdateVcfSequenceDictionary onto the karyotypic dict.
#   - Only THEN can SortVcf run against the karyotypic dict -- SortVcf hard
#     -asserts that every input's embedded dictionary exactly matches
#     --SEQUENCE_DICTIONARY (order included), so handing it the karyotypic
#     dict against a file still declaring the original order throws
#     "SAM dictionaries are not the same".
#
# generate_sample_list.sh writes sample_list.txt containing only samples
# that have BOTH a gVCF and a called-VCF for every region listed in
# /software/projects/pawsey0933/deepvariant/regions (any incomplete sample
# is logged to incomplete_samples.txt instead, not silently merged short).
# Re-run it any time more region VCFs land, to pick up newly-completed
# samples.
#
# Then submit with the array size matched to the sample count (safer than
# trusting the hardcoded 1-234 above if the sample count ever changes):
#
#   sbatch --array=1-$(wc -l < /scratch/pawsey0933/cfolland/deepvariant/hprc/pg_dv/merged/sample_list.txt)%20 \
#     merge_dv_array.sh
#
# Each array task merges ONE sample's per-region outputs into two
# per-sample files: a merged gVCF and a merged called VCF.
# ---------------------------------------------------------------------------

set -euo pipefail

module load singularity/4.1.0-slurm

input_dir=/scratch/pawsey0933/cfolland/deepvariant/hprc/pg_dv
output_dir=/scratch/pawsey0933/cfolland/deepvariant/hprc/pg_dv/merged
regions_dir=/software/projects/pawsey0933/deepvariant/regions
sample_list=${output_dir}/sample_list.txt
src_dict=/software/projects/pawsey0933/pangenome/refs/hprc-v2.0-mc-grch38.dict
dict=/software/projects/pawsey0933/pangenome/refs/hprc-v2.0-mc-grch38.karyotypic.dict
sif=/software/projects/pawsey0933/deepvariant/gatk4:4.4.0.0.sif
tag=hprc-v2.0-mc-grch38.pg-deepvariant.1-10-0

mkdir -p "${output_dir}"

if [[ ! -s "${sample_list}" ]]; then
    echo "ERROR: sample list not found or empty: ${sample_list}" >&2
    echo "Run the one-time prep steps at the top of this script first." >&2
    exit 1
fi

if [[ ! -s "${dict}" ]]; then
    echo "ERROR: karyotypic dictionary not found or empty: ${dict}" >&2
    echo "Run ./build_karyotypic_dict.sh first." >&2
    exit 1
fi

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${sample_list}")
if [[ -z "${sample}" ]]; then
    echo "ERROR: no sample found at line ${SLURM_ARRAY_TASK_ID} of ${sample_list}" >&2
    exit 1
fi

tmp_dir=${output_dir}/tmp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "${tmp_dir}"
trap 'rm -rf "${tmp_dir}"' EXIT

gvcf_list=${tmp_dir}/${sample}.gvcf.list
vcf_list=${tmp_dir}/${sample}.vcf.list
: > "${gvcf_list}"
: > "${vcf_list}"

# Build the per-sample --INPUT list from the regions dir (existence-checked
# against each expected contig). Order here doesn't matter for correctness
# any more (SortVcf does the real ordering downstream) -- `sort -V` just
# keeps the list/log output readable.
mapfile -t regions < <(
    find "${regions_dir}" -maxdepth 1 -type f -name '*.bed' -printf '%f\n' \
        | sed -E 's/^([^_]+)_.*$/\1/' \
        | sort -V -u
)
for contig in "${regions[@]}"; do
    g="${input_dir}/${sample}.deepvariant.${contig}.g.vcf.gz"
    v="${input_dir}/${sample}.deepvariant.${contig}.vcf.gz"
    [[ -f "${g}" ]] && echo "${g}" >> "${gvcf_list}"
    [[ -f "${v}" ]] && echo "${v}" >> "${vcf_list}"
done

n_gvcf=$(wc -l < "${gvcf_list}")
n_vcf=$(wc -l < "${vcf_list}")

if [[ "${n_gvcf}" -eq 0 && "${n_vcf}" -eq 0 ]]; then
    echo "ERROR: no per-region VCFs found for sample ${sample} in ${input_dir}" >&2
    exit 1
fi

# merge_sort <list_file> <n_files> <final_output> <label>
#
# Three-step pipeline, required because the per-region VCFs' embedded
# dictionary is the original (non-karyotypic) one:
#   1. MergeVcfs against src_dict (matches what's embedded in the inputs).
#      Output is complete but still in the original contig order.
#   2. UpdateVcfSequenceDictionary reheaders that output onto the
#      karyotypic dict, WITHOUT touching record order yet.
#   3. SortVcf against the karyotypic dict now succeeds (the file's
#      embedded dict matches) and physically reorders every record.
merge_sort() {
    local list_file=$1 n_files=$2 final_output=$3 label=$4
    local unsorted="${tmp_dir}/${sample}.${label}.unsorted.vcf.gz"
    local reheadered="${tmp_dir}/${sample}.${label}.reheadered.vcf.gz"

    echo "Merging ${n_files} ${label} region file(s) for ${sample}"
    singularity exec -B /scratch/pawsey0933 -B /software/projects/pawsey0933 \
        "${sif}" \
            gatk --java-options "-Xmx3276M -XX:-UsePerfData" \
                MergeVcfs \
                --INPUT "${list_file}" \
                --OUTPUT "${unsorted}" \
                --SEQUENCE_DICTIONARY "${src_dict}" \
                --TMP_DIR "${tmp_dir}"

    echo "Reheadering merged ${label} onto karyotypic dictionary for ${sample}"
    singularity exec -B /scratch/pawsey0933 -B /software/projects/pawsey0933 \
        "${sif}" \
            gatk --java-options "-Xmx3276M -XX:-UsePerfData" \
                UpdateVcfSequenceDictionary \
                --INPUT "${unsorted}" \
                --OUTPUT "${reheadered}" \
                --SEQUENCE_DICTIONARY "${dict}"

    echo "Sorting merged ${label} for ${sample}"
    singularity exec -B /scratch/pawsey0933 -B /software/projects/pawsey0933 \
        "${sif}" \
            gatk --java-options "-Xmx3276M -XX:-UsePerfData" \
                SortVcf \
                --INPUT "${reheadered}" \
                --OUTPUT "${final_output}" \
                --SEQUENCE_DICTIONARY "${dict}" \
                --TMP_DIR "${tmp_dir}"
}

# Skip work whose final output already exists, so a resubmitted array only
# (re)does samples that actually need it -- e.g. after a fix like this one,
# or when retrying failed tasks. Set FORCE_REMERGE=1 (sbatch --export=... or
# `export FORCE_REMERGE=1` before sbatch, since --export=NONE above strips
# the environment otherwise) to force regeneration regardless.
final_gvcf=${output_dir}/${sample}.${tag}.g.vcf.gz
final_vcf=${output_dir}/${sample}.${tag}.vcf.gz
force=${FORCE_REMERGE:-0}

if [[ "${n_gvcf}" -gt 0 ]]; then
    if [[ "${force}" != "1" && -s "${final_gvcf}" ]]; then
        echo "SKIP: ${final_gvcf} already exists (set FORCE_REMERGE=1 to redo)"
    else
        merge_sort "${gvcf_list}" "${n_gvcf}" "${final_gvcf}" gvcf
    fi
else
    echo "WARNING: no gVCF region files found for ${sample}, skipping gVCF merge" >&2
fi

if [[ "${n_vcf}" -gt 0 ]]; then
    if [[ "${force}" != "1" && -s "${final_vcf}" ]]; then
        echo "SKIP: ${final_vcf} already exists (set FORCE_REMERGE=1 to redo)"
    else
        merge_sort "${vcf_list}" "${n_vcf}" "${final_vcf}" vcf
    fi
else
    echo "WARNING: no called-VCF region files found for ${sample}, skipping VCF merge" >&2
fi

echo "Done: ${sample}"