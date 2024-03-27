library(vcfR)
library(hierfstat)
library(adegenet)
library(ggplot2)
library(radiator)

rm(list = ls())

setwd("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/R_code")

##Load genetic data and convert it to genlight objects, and prepare it for analysis
vcf <- read.vcfR( "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/original_vcf/bi_20missing_filtSNP_maf_005.recode.vcf", verbose = FALSE ) # 20% missing values

## sample information
sample_info <- read.csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Ctarsalis_sample_w_GPS_climate_average_new_filtered_id_region.csv", row.names = 1, stringsAsFactors = TRUE)

genind <- vcfR2genind(vcf)

genind@pop <- sample_info$region

hierfstat <- genind2hierfstat(genind)

write.bayescan(hierfstat)

