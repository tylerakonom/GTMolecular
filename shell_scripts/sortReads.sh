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
module load picard
outdirectory=/gpfs/alpine1/scratch/$USER/GTM/aligned/sorted
indirectory=/gpfs/alpine1/scratch/$USER/GTM/aligned

# Running picard
java -jar $PICARD_JAR SortSam INPUT=CDTC33_1.sam OUTPUT=CDTC33_1_sorted.bam SORT_ORDER=coordinate