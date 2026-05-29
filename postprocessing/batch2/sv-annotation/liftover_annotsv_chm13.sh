#!/bin/bash -l

#SBATCH --job-name=liftover
#SBATCH --account=pawsey0933
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out
#SBATCH --export=ALL

#source the chain file
#wget https://hgdownload.soe.ucsc.edu/gbdb/hg38/liftOver/hg38ToGCA_009914755.4.over.chain.gz
#gunzip hg38ToGCA_009914755.4.over.chain.gz

# prepare environment
conda activate /software/projects/pawsey0933/cfolland/miniforge3/envs/annotsv/
module load bcftools/1.15--haf5b3da_0
module load bedtools/2.30.0--h468198e_3

# directories
ANNOTSV=/scratch/pawsey0933/cfolland/AnnotSV
ChainDir=/software/projects/pawsey0933/t2t/sv_anno
chain=$ChainDir/hg38ToGCA_009914755.4.over.chain
liftover=/software/projects/pawsey0933/t2t/sv_anno/liftover_script.sh


# function
run_liftover () {
    local infile="$1"
    local outfile

    outfile=$(basename "$infile")
    outfile=${outfile/GRCh38/CHM13}

    echo "----------------------------------"
    echo "INPUT : $infile"
    echo "OUTPUT: $outfile"
    echo "----------------------------------"

    # ✅ call your robust liftover script directly
    "$liftover" "$infile" "$outfile" "$chain"
}


# BENIGN
cd "$ANNOTSV/share/AnnotSV/Annotations_Human/SVincludedInFt/BenignSV/CHM13/"

for f in ../GRCh38/benign*_GRCh38.sorted.bed; do
    run_liftover "$f"
done

# PATHOGENIC SV
cd "$ANNOTSV/share/AnnotSV/Annotations_Human/FtIncludedInSV/PathogenicSV/CHM13/"

for f in ../GRCh38/pathogenic_*_SV_GRCh38.sorted.bed; do
    run_liftover "$f"
done

# CLINVAR
cd "$ANNOTSV/share/AnnotSV/Annotations_Human/FtIncludedInSV/PathogenicSNVindel/CHM13/"

run_liftover "../GRCh38/pathogenic_SNVindel_GRCh38.sorted.bed"
