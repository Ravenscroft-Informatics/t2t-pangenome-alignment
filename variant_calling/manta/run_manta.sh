#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH -c 4
#SBATCH --job-name=manta
#SBATCH --partition=work
#SBATCH --account=pawsey0933
#SBATCH --mem=32GB
#SBATCH --time=4:00:00
#SBATCH --array=0-4%10
#SBATCH --export=NONE
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END,FAIL
#SBATCH --error=logs_2/%A_%a.err
#SBATCH --output=logs_2/%A_%a.out

set -euo pipefail

module load singularity/4.1.0-slurm
module load samtools/1.15--h3843a85_0
mamba activate /software/projects/pawsey0933/cfolland/miniforge3/envs/manta/

#REF=/software/projects/pawsey0933/sv/references/hg38_masked/Homo_sapiens_assembly38_masked.fasta
REF=/software/projects/pawsey0933/manta/hprc-v2.0-mc-grch38.fa
SIF=/software/projects/pawsey0933/manta/manta_main.sif

#CRAM_DIR=/scratch/pawsey0933/cfolland/hg38_crams
#CRAM_DIR=/scratch/pawsey0933/cfolland/bams_hprc
CRAM_DIR=/scratch/pawsey0933/cfolland/vg_snakemake/results/bams
#OUT=/scratch/pawsey0933/cfolland/manta/results_hg38_bwamem
OUT=/scratch/pawsey0933/cfolland/manta/results_hg38_giraffe

mkdir -p "${OUT}"

# Restrict to a fixed set of samples
SAMPLES=(D16-1946 D19-0335 D23-0015 D23-0033 D23-0055)

# Create array of CRAMs
#CRAMS=(${CRAM_DIR}/*.cram)
#CRAMS=(${CRAM_DIR}/*.bam)
CRAMS=()
for s in "${SAMPLES[@]}"; do
    CRAMS+=("${CRAM_DIR}/${s}.hprc-v2.0-mc-grch38.surj_realn.bam")
done

BAM=${CRAMS[$SLURM_ARRAY_TASK_ID]}

if [[ ! -e "${BAM}" ]]; then
    echo "ERROR: expected BAM not found: ${BAM}" >&2
    exit 1
fi

SAMPLE=$(basename "${BAM}" .hprc-v2.0-mc-grch38.surj_realn.bam)
#SAMPLE=$(basename "${BAM}" .sorted.cram)

WORKDIR=/scratch/pawsey0933/cfolland/manta/manta_runs/${SAMPLE}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
#PREFIX=${SAMPLE}.bwamem_grch38.manta
PREFIX=${SAMPLE}.giraffe_grch38.manta

echo "Processing ${SAMPLE}"
echo "CRAM: ${BAM}"

# Skip if outputs already exist
EXPECTED_OUTPUTS=(
    "${OUT}/${PREFIX}.candidate_small_indels.vcf.gz"
    "${OUT}/${PREFIX}.candidate_small_indels.vcf.gz.tbi"
    "${OUT}/${PREFIX}.candidate_sv.vcf.gz"
    "${OUT}/${PREFIX}.candidate_sv.vcf.gz.tbi"
    "${OUT}/${PREFIX}.diploid_sv.vcf.gz"
    "${OUT}/${PREFIX}.diploid_sv.vcf.gz.tbi"
)

ALL_EXIST=true
for f in "${EXPECTED_OUTPUTS[@]}"; do
    if [[ ! -s "${f}" ]]; then
        ALL_EXIST=false
        break
    fi
done

if [[ "${ALL_EXIST}" == true ]]; then
    echo "All outputs already exist for ${SAMPLE}, skipping."
    exit 0
fi

mkdir -p "${WORKDIR}"
cd "${WORKDIR}"

ulimit -s 65520

# Configure Manta
singularity exec -B /scratch/pawsey0933 -B /software/projects/pawsey0933 \
    "${SIF}" \
    configManta.py \
        --bam "${BAM}" \
        --referenceFasta "${REF}" \
        --runDir manta

# Run workflow
singularity exec -B /scratch/pawsey0933 -B /software/projects/pawsey0933 \
    "${SIF}" \
    "${WORKDIR}/manta/runWorkflow.py" \
        -m local \
        -j 8

# Rename outputs
mv manta/results/variants/candidateSmallIndels.vcf.gz \
   ${OUT}/${PREFIX}.candidate_small_indels.vcf.gz

mv manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
   ${OUT}/${PREFIX}.candidate_small_indels.vcf.gz.tbi

mv manta/results/variants/candidateSV.vcf.gz \
   ${OUT}/${PREFIX}.candidate_sv.vcf.gz

mv manta/results/variants/candidateSV.vcf.gz.tbi \
   ${OUT}/${PREFIX}.candidate_sv.vcf.gz.tbi

mv manta/results/variants/diploidSV.vcf.gz \
   ${OUT}/${PREFIX}.diploid_sv.vcf.gz

mv manta/results/variants/diploidSV.vcf.gz.tbi \
   ${OUT}/${PREFIX}.diploid_sv.vcf.gz.tbi