#!/bin/bash
#SBATCH --partition=Orion
#SBATCH --time=24:00:00
#SBATCH --job-name=partial_RDA_significance 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=4
#SBATCH --mem=256gb

# load R
module load R/4.1.3

# run the r file
Rscript /projects/cooper_research2/eric/data_landscape_genetics/RDA_significant_test/partial_rda/Partial_RDA_geo_env.R
 
