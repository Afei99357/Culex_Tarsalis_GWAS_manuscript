#!/bin/bash
#SBATCH --partition=Orion
#SBATCH --time=24:00:00
#SBATCH --job-name=job_for_admixture
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1
#SBATCH --mem=256gb

FILE=culex_plink_new

cd /projects/cooper_research2/eric/data_landscape_genetics/plink_file_for_admixture

# load R
for i in {6..8}
do
	echo start no.$i
	/projects/cooper_research2/eric/data_landscape_genetics/admixture_linux-1.3.0/admixture --cv $FILE.bed $i > log${i}.out
	echo end no.$I
done

# identify the best value of k clusters 
grep "CV" *out | awk '{print $3,$4}' | cut -c 4,7-20 > $FILE.cv.error

# make table for plot
awk '{split($1,name,"."); print $1,name[2]}' $FILE.nosex > $FILE.list

