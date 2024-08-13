#!/bin/bash

#SBATCH --job-name="ct_filter"
#SBATCH --partition=Orion
#SBATCH --time=72:00:00
#SBATCH --nodes=4
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=12
#SBATCH --mem=48gb

module load vcftools

# No missing

vcftools --vcf ct_filtInd.recode.vcf --max-missing 1 --maf 0.05 --recode --recode-INFO-all --out ct_nomissing_filtSNP
vcftools --vcf ct_nomissing_filtSNP.recode.vcf --max-alleles 2 --minQ 30 --recode --recode-INFO-all --out bi_nomissing_filtSNP

# 10% missing

vcftools --vcf ct_filtInd.recode.vcf --max-missing 0.9 --maf 0.05 --recode --recode-INFO-all --out ct_10missing_filtSNP
vcftools --vcf ct_10missing_filtSNP.recode.vcf --max-alleles 2 --minQ 30 --recode --recode-INFO-all --out bi_10missing_filtSNP

# 20% missing

vcftools --vcf ct_filtInd.recode.vcf --max-missing 0.8 --maf 0.05 --recode --recode-INFO-all --out ct_20missing_filtSNP
vcftools --vcf ct_20missing_filtSNP.recode.vcf --max-alleles 2 --minQ 30 --recode --recode-INFO-all --out bi_20missing_filtSNP