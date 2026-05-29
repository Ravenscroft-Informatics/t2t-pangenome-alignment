#!/bin/bash -l

#SBATCH --job-name=filter_family
#SBATCH --account=pawsey0933
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --time=6:00:00
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out
#SBATCH --export=ALL

# Load modules
module load vcftools/0.1.16--pl5321hd03093a_7
module load bcftools/1.15--haf5b3da_0
conda activate vcf_tools

cd /scratch/pawsey0933/cfolland/t2t/

prefix="T2T_VEP_genmod"

# Define suffixes for each group
AR_suffixes="AR_comp|AR_comp_dn|AR_hom|AR_hom_dn"
AD_suffixes="AD|AD_dn"
X_suffixes="XD|XR|XD_dn|XR_dn"

for family in $(awk '{print $1}' trios.ped | sort -nu); do
    mkdir -p $family
    cd $family
    rm *.vcf*
    samples=$(awk -v var="$family" '$1==var {print $2}' ../trios.ped | tr '\n' ',')
    echo "Processing family: $family with samples: $samples"
    
    # Extract family-specific samples
    bcftools view --threads 6 -s ${samples%,*} -c1 -Ov -o ${family}_${prefix}.vcf ../sorted_toref_genmod_filtered.vcf.gz
    bgzip -f ${family}_${prefix}.vcf
    tabix -p vcf ${family}_${prefix}.vcf.gz
    
    # Filter variants for each group
      for filter in AR AD X; do
        filt="^${family}:.*${filter}.*"
        echo "Filtering for: $filt"
        
        # Apply filter with bcftools
        bcftools view --threads 6 -i "INFO/GeneticModels~'${filt}'" -Ov -o ${family}_${filter}_${prefix}.vcf ${family}_${prefix}.vcf.gz
        
        # Check for output
        if [[ -f ${family}_${filter}_${prefix}.vcf ]]; then
            echo "Filtered VCF created: ${family}_${filter}_${prefix}.vcf"
        else
            echo "Warning: No variants found for ${filter} in family $family"
        fi
    done
    cd ..
done
