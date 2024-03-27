rm(list = ls())

##Loading library
library(vcfR)     # for vcf file
library(adegenet) # for genlight obeject and tab() function
library(vegan)    # Used to run PCA & RDA

##Load genetic data and convert it to genlight objects, and prepare it for analysis
vcf <- read.vcfR( "/projects/cooper_research2/eric/data_landscape_genetics/RDA_significant_test/bi_20missing_filtSNP_maf_005.recode.vcf", verbose = FALSE ) # no missing values

SNP_genlight <- vcfR2genlight(vcf) # convert to genlight objects

allele_matrix_genlight <- tab(SNP_genlight, freq = FALSE, NA.method = "asis") # get the matrix of allele counts

dim(allele_matrix_genlight) # check the dimension of genetic data

sum(is.na(allele_matrix_genlight)) # check the data dimmension and number of NAs

gen.imp <- apply(allele_matrix_genlight, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x)))))) # impute the NA values into most common allele

sum(is.na(gen.imp)) # check if No NAs

##Load the enviromental data and get it ready for analysis

env <- read.csv("/projects/cooper_research2/eric/data_landscape_genetics/RDA_significant_test/Ctarsalis_sample_w_GPS_climate_average_new_filtered_id_region.csv", row.names = 1) # load enviromental data

env <- env[, !(names(env) %in% c("locID", "State", "City", "date"))] # remove unused columns

str(env) # check the structure of the enviromental data

env$vcfID <- as.character(env$vcfID) # Make individual names characters (not factors)

identical(rownames(gen.imp), env[,1]) # Confirm that genotypes and environmental data are in the same order

pred <- env[, 18:30] # subset enviromental predictors and shorten their names

## remove as little variables as we can based on any pair if the correlation is higher than 0.7
pred <- subset(pred, select=-c(avg_ssr, avg_evabs, avg_tp, avg_swvl1, avg_sf))

##Run Redundancy Analysis (RDA): a multivariate GEA

mos.rda <- rda(gen.imp ~ ., data=pred, scale=T) # ## Run RDA
mos.rda

RsquareAdj(mos.rda) # # R2 will be biased and should be adjusted based on the number of predictors

summary(mos.rda)$concont ## The eigenvalues for the constrained axes reflect the variance explained by each canonical axis

screeplot(mos.rda) # screeplot of the canonical eigenvalues

# full model
signif.full <- anova.cca(mos.rda, parallel=getOption("mc.cores")) # default is permutation=999
signif.full

# run a formal test of statistical significance of each constrained axis
# # The null hypothesis is that no linear relationship exists between the SNP data and the environmental predictors
signif.axis <- anova.cca(mos.rda, by="axis", parallel=getOption("mc.cores"))
signif.axis
