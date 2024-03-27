#!/bin/bash
#SBATCH --partition=Orion
#SBATCH --time=12:00:00
#SBATCH --job-name=job_for_amova_test
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=4
#SBATCH --mem=128gb

# load R
module load R

# run the r file
Rscript /projects/cooper_research2/eric/data_landscape_genetics/amova_test/amova_test_cluster.R
 