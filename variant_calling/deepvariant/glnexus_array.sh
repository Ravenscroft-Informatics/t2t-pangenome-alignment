#!/bin/bash -l

#SBATCH --job-name=glnexus_array
#SBATCH --account=pawsey0933
#SBATCH --partition=highmem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=0
#SBATCH --nodes=1
#SBATCH --time=16:00:00
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=FAIL
#SBATCH --error=%x.%A_%a.err
#SBATCH --output=%x.%A_%a.out
#SBATCH --export=NONE
#SBATCH --array=1-25%4

# ---------------------------------------------------------------------------
# Joint-genotypes ONE contig per array task instead of the whole genome in a
# single serial glnexus_cli run. Restricting each task to one chromosome via 
# --bed turns a >100h serial job into ~25 independent, much shorter tasks run 
# in parallel (subject to the %4 concurrency cap # below).
#
# One-time prep (run ONCE before the first sbatch submission):
#   ./prep_glnexus_regions.sh
# writes glnexus_regions/contigs.txt (one contig per line, in the order
# array tasks should run) and glnexus_regions/<contig>.bed for each.
#
# GLnexus's own performance notes say its internal parallelism has only
# been tuned up to ~16 cores -- the original single-job run's
# --cpus-per-task=24 was already past the point of returns, and the failed
# array run's low CPU efficiency (16%, 6.7%) confirms memory, not CPU, is
# the binding constraint here. --cpus-per-task is left at 16 for that
# reason even though the node may grant far more cores than that anyway.
#
# Prefer to size the array off contigs.txt directly rather than trust the
# hardcoded 1-25 above (in case a contig was missing from the dict):
#
#   sbatch --array=1-$(wc -l < glnexus_regions/contigs.txt)%4 glnexus_array.sh
#
# ---------------------------------------------------------------------------

set -ueo pipefail

merged_dir=/scratch/pawsey0933/cfolland/deepvariant/hprc/pg_dv/merged
cd "${merged_dir}"

PREFIX=hprc-v2.0-mc-grch38.pg-deepvariant.1-10-0
regions_dir=${merged_dir}/glnexus_regions
contigs_file=${regions_dir}/contigs.txt
out_dir=${merged_dir}/glnexus_out
mkdir -p "${out_dir}"

if [[ ! -s "${contigs_file}" ]]; then
    echo "ERROR: ${contigs_file} not found. Run ./prep_glnexus_regions.sh first." >&2
    exit 1
fi

contig=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${contigs_file}")
if [[ -z "${contig}" ]]; then
    echo "ERROR: no contig found at line ${SLURM_ARRAY_TASK_ID} of ${contigs_file}" >&2
    exit 1
fi

bed=${regions_dir}/${contig}.bed
if [[ ! -s "${bed}" ]]; then
    echo "ERROR: BED file not found for ${contig}: ${bed}" >&2
    exit 1
fi

final_bcf=${out_dir}/${contig}.${PREFIX}.bcf
force=${FORCE_REDO:-0}
if [[ "${force}" != "1" && -s "${final_bcf}" ]]; then
    echo "SKIP: ${final_bcf} already exists (set FORCE_REDO=1 to redo)"
    exit 0
fi

mem_gb=$(( SLURM_MEM_PER_NODE / 1024 - 40 ))
if [[ "${mem_gb}" -lt 1 ]]; then
    echo "ERROR: SLURM_MEM_PER_NODE (${SLURM_MEM_PER_NODE:-unset} MB) too small to budget from -- is --mem=0 in effect?" >&2
    exit 1
fi

# Scoped to this array job + task index so concurrent tasks never collide
# on the same scratch DB path (glnexus_cli refuses to run if --dir already
# exists). Cleaned up on exit either way.
db_dir=${merged_dir}/GLnexus.DB.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
trap 'rm -rf "${db_dir}"' EXIT

n_gvcf=$(ls -1 ./*."${PREFIX}".g.vcf.gz 2>/dev/null | wc -l)
if [[ "${n_gvcf}" -eq 0 ]]; then
    echo "ERROR: no *.${PREFIX}.g.vcf.gz files found in ${merged_dir}" >&2
    exit 1
fi
echo "Joint-calling ${contig} across ${n_gvcf} gVCF(s)" >&2

# Write to a temp name and mv into place on success only, so a task killed
# mid-run (e.g. another timeout) never leaves a truncated file at
# final_bcf that a later resubmit's skip-if-exists check would mistake for
# a completed contig.
tmp_bcf=${out_dir}/${contig}.${PREFIX}.tmp.bcf
/software/projects/pawsey0933/gmonahan/T2T/glnexus/glnexus_cli \
    -m "${mem_gb}" \
    -c DeepVariantWGS \
    --dir "${db_dir}" \
    --bed "${bed}" \
    ./*."${PREFIX}".g.vcf.gz \
    > "${tmp_bcf}"
mv "${tmp_bcf}" "${final_bcf}"

echo "Done: ${contig} -> ${final_bcf}"