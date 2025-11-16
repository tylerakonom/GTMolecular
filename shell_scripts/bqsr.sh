#!/bin/bash

# Written by: tyak9569
# Date: 11_15_2025
# Purpose: bqsr script for tyak9569

#SBATCH --partition=amilan    # Summit partition
#SBATCH --qos=normal                 # Summit qos
#SBATCH --time=001:00:00           # Max wall time in HHH:MM:SS
#SBATCH --ntasks=4           # Number of tasks per job  
#SBATCH --nodes=1             # Number of nodes per job
#SBATCH --job-name=bqsr      # Job submission name
#SBATCH --output=bqsr%j.out   # Output file name with Job ID


# purge all existing modules
module purge

# load the module needed to run the software container, and set up temporary directories
module load gatk
outdirectory=/gpfs/alpine1/scratch/$USER/GTM/bqsr
indirectory=/gpfs/alpine1/scratch/$USER/GTM/aligned
ks=/gpfs/alpine1/scratch/$USER/GTM/genome/resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf
ref=/gpfs/alpine1/scratch/$USER/GTM/genome/GCA_000001405.15_GRCh38_full_analysis_set.fna

# Running gatk
gatk BaseRecalibrator -R ${ref} -I ${indirectory}/${filename}_dedup.bam -known-sites ${ks} -O ${outdirectory}/${filename}_recal_data.table