#!/bin/bash

# Written by: tyak9569
# Date: 03_13_2023
# Purpose: BWA script for tyak9569

#SBATCH --partition=amilan    # Summit partition
#SBATCH --qos=normal                 # Summit qos
#SBATCH --time=004:00:00           # Max wall time in HHH:MM:SS
#SBATCH --ntasks=24           # Number of tasks per job  
#SBATCH --nodes=1             # Number of nodes per job
#SBATCH --job-name=bwa      # Job submission name
#SBATCH --output=bwa%j.out   # Output file name with Job ID


# purge all existing modules
module purge

# load the module needed to run the software container, and set up temporary directories
module load bwa
outdirectory=/gpfs/alpine1/scratch/$USER/GTM/aligned
indirectory=/gpfs/alpine1/scratch/$USER/GTM
genome=/gpfs/alpine1/scratch/$USER/GTM/genome/GCA_000001405.15_GRCh38_full_analysis_set.fna
rg=@RG\\tID:${filename}\\tLB:${filename}\\tPL:ILLUMINA\\tSM:${filename}


# Running bwa
bwa mem -t 8 -M -R ${rg} ${genome} ${indirectory}/${filename}.fastq > ${outdirectory}/${filename}.sam