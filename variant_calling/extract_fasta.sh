#!/bin/bash -l
#SBATCH --job-name=build_reference
#SBATCH --account=pawsey0933
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --time=07:00:00
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out
#SBATCH --export=ALL

set -euo pipefail

module load singularity/4.1.0-slurm
module load samtools/1.15--h3843a85_0

assembly=grch38
base=hprc-v2.0-mc-${assembly}
sif=/software/projects/pawsey0933/manta/vg_v1.59.0.sif
gbz=/scratch/pawsey0933/cfolland/vg_snakemake/refs/hprcv2.0/${base}.gbz
paths_list=/software/projects/pawsey0933/manta/${base}_path_lists.txt
sorted_paths_list=/software/projects/pawsey0933/manta/temp.sorted_path_list.${SLURM_JOB_ID}.${base}.txt
tmp=/software/projects/pawsey0933/manta/temp.path_name.${SLURM_JOB_ID}.${base}.txt
fasta=/software/projects/pawsey0933/manta/${base}.fa
dict="${fasta%.*}.dict"
seqn_prefix="GRCh38#0#"

# --- 1. Generate the path list ----------------------------------------------
singularity exec -B /scratch/pawsey0933/cfolland -B /software/projects/pawsey0933 "${sif}" \
    vg paths \
    -L \
    -x "${gbz}" \
    | grep GRCh > "${paths_list}"

# --- 2. Sort the path list into karyotypic order ----------------------------
# chr1..chr22 (numeric order), then chrX, chrY, chrM/chrMT, then anything else
# (alt/unplaced contigs) alphabetically. Extracting in this order means the
# FASTA -- and therefore the .fai and .dict built from it -- come out in
# correct karyotypic order without needing a separate reorder step later.
awk -v prefix="${seqn_prefix}" '
{
    chrom = $0
    sub("^" prefix, "", chrom)
    name = chrom
    sub(/^chr/, "", name)
    if (name == "X")                      { rank = 1000; num = 0 }
    else if (name == "Y")                 { rank = 1001; num = 0 }
    else if (name == "M" || name == "MT") { rank = 1002; num = 0 }
    else if (name ~ /^[0-9]+$/)           { rank = 0; num = name + 0 }
    else                                   { rank = 2000; num = 0 }
    printf "%d\t%d\t%s\t%s\n", rank, num, name, $0
}' "${paths_list}" | sort -k1,1n -k2,2n -k3,3 | cut -f4- > "${sorted_paths_list}"

# --- 3. Extract the FASTA in karyotypic order -------------------------------
: > "${fasta}"

while read -r CHR; do
    echo "${CHR}" > "${tmp}"
    singularity exec -B /scratch/pawsey0933 -B /software/projects/pawsey0933 "${sif}" \
        vg paths --extract-fasta -p "${tmp}" --xg "${gbz}" \
        | sed -e "s/>${seqn_prefix}/>/g" >> "${fasta}"
done < "${sorted_paths_list}"

rm -f "${tmp}" "${sorted_paths_list}"

# --- 4. Index and build the sequence dictionary -----------------------------
# samtools dict preserves the input FASTA's sequence order, so it inherits
# the karyotypic order established above -- .fa and .dict stay consistent.
samtools faidx -o "${fasta}.fai" "${fasta}"
samtools dict -o "${dict}" "${fasta}"