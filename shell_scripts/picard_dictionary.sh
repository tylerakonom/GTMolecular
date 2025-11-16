#!/bin/bash

# Written by: tyak9569
# Date: 11_15_2025
# Purpose: dictionary script for tyak9569

#SBATCH --partition=amilan    # Summit partition
#SBATCH --qos=normal                 # Summit qos
#SBATCH --time=001:00:00           # Max wall time in HHH:MM:SS
#SBATCH --ntasks=24           # Number of tasks per job  
#SBATCH --nodes=1             # Number of nodes per job
#SBATCH --job-name=dictionary      # Job submission name
#SBATCH --output=dictionary%j.out   # Output file name with Job ID


# purge all existing modules
module purge

# load the module needed to run the software
module load picard
# Running picard
java -jar $PICARD CreateSequenceDictionary R=GCA_000001405.15_GRCh38_full_analysis_set.fna O=GCA_000001405.15_GRCh38_full_analysis_set.fna.dict