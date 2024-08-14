library("vcfR")     # for vcf file
library("adegenet") # for genlight obeject and tab() function
library("dplyr") ## for filtering table
library("hierfstat") ## for pairwise Fst

rm(list = ls())

# ## SNP data
vcf <- read.vcfR( "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/original_vcf/bi_20missing_filtSNP_maf_005.recode.vcf", verbose = FALSE )

## Get the names of individual samples
vcfIDs <- as.vector(colnames(vcf@gt))[-1]

## enviromental data
## read sample information file
env <- read.csv("/Users/ericliao/Desktop/WNV_project_files/Ctarsalis_sample_w_GPS_climate_average_new_filtered_id_region.csv")

sample_info_popid <- subset(env, select = c(vcfID, popID, region)) ## get popID and vcfID

matched_pop <- sample_info_popid %>% filter(vcfID %in% vcfIDs) ## find Population ID

population <- matched_pop$popID ## get the population info
regions <- matched_pop$region ## get the region infor
indivduals <- as.character(matched_pop$vcfID) ## get the individual sample name

SNP_genind <- vcfR2genind(vcf, return.alleles = TRUE, ind.names = matched_pop$vcfID) ## convert to genind object

# SNP_genind$pop <- as.factor(matched_pop$region) ## add region information as pop information to genind object
SNP_genind$pop <- as.factor(matched_pop$popID) ## between populations

SNP_genind@call$ind.names <- indivduals ## add sample individual names to genind objects

strata_df <- data.frame(regions, population) 
strata(SNP_genind) <- strata_df ## add strata

## Pairwise Fst method option: y “Nei87” and “WC84” return pairwise FSTs estimated following Nei (1987) pairwise.neifst and Weir & Cockerham (1984) pp.fst respectively
pairwise_fsts <- genet.dist(SNP_genind, method = "WC84")

# Convert distance object to a matrix
pairwise_fsts_matrix <- as.matrix(pairwise_fsts)

write.csv(pairwise_fsts_matrix, "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/IBD_IBE_test/pairwise_fsts_matrix_between_population.csv")
