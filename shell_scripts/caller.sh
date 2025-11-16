#!/bin/bash

# Written by: tyak9569
# Date: 11_15_2025
# Purpose: caller script for tyak9569

#SBATCH --partition=amilan    # Summit partition
#SBATCH --qos=normal                 # Summit qos
#SBATCH --time=001:00:00           # Max wall time in HHH:MM:SS
#SBATCH --ntasks=4           # Number of tasks per job  
#SBATCH --nodes=1             # Number of nodes per job
#SBATCH --job-name=caller      # Job submission name
#SBATCH --output=caller%j.out   # Output file name with Job ID


# purge all existing modules
module purge

# load the module needed to run the software container, and set up temporary directories
module load gatk
outdirectory=/gpfs/alpine1/scratch/$USER/GTM/bqsr/out
indirectory=/gpfs/alpine1/scratch/$USER/GTM/bqsr
ref=/gpfs/alpine1/scratch/$USER/GTM/genome/GCA_000001405.15_GRCh38_full_analysis_set.fna

# Running gatk
gatk HaplotypeCaller -R ${ref} -I ${indirectory}/${filename}_recal.bam -O ${outdirectory}/${filename}_raw_variants.vcf
gatk SelectVariants -R ${ref} -V ${outdirectory}/${filename}_raw_variants.vcf --select-type-to-include SNP -O ${outdirectory}/${filename}_raw_snps.vcf
gatk SelectVariants -R ${ref} -V ${outdirectory}/${filename}_raw_variants.vcf --select-type-to-include INDEL -O ${outdirectory}/${filename}_raw_indels.vcf

gatk VariantFiltration \
-R ${ref} \
-V ${outdirectory}/${filename}_raw_snps.vcf \
--filter-name "QD_filter" \
--filter-expression "QD<2.0" \
--filter-name "FS_filter" \
--filter-expression "FS>60.0" \
--filter-name "MQ_filter" \
--filter-expression "MQ<40.0" \
--filter-name "SOR_filter" \
--filter-expression "SOR>10.0" \
-O ${outdirectory}/${filename}_filtered_snps.vcf

gatk VariantFiltration \
-R ${ref} \
-V ${outdirectory}/${filename}_raw_indels.vcf \
--filter-name "QD_filter" \
--filter-expression "QD<2.0" \
--filter-name "FS_filter" \
--filter-expression "FS>200.0" \
--filter-name "SOR_filter" \
--filter-expression "SOR>10.0" \
-O ${outdirectory}/${filename}_filtered_indels.vcf