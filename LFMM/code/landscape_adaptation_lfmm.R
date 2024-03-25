rm(list = ls())

##Loading library
library(vcfR)     # for vcf file
library(adegenet) # for genlight obeject and tab() function
library(vegan)    # Used to run PCA & RDA
library(lfmm)     # Used to run LFMM
library(qvalue)   # Used to post-process LFMM output
library(psych)
library(VennDiagram)

##Load genetic data and convert it to genlight objects, and prepare it for analysis

vcf <- read.vcfR( "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/original_vcf/bi_20missing_filtSNP_maf_005.recode.vcf", verbose = FALSE ) # no missing values

SNP_genlight <- vcfR2genlight(vcf) # convert to genlight objects

allele_matrix_genlight <- tab(SNP_genlight, freq = FALSE, NA.method = "asis") # get the matrix of allele counts

dim(allele_matrix_genlight) # check the dimension of genetic data

sum(is.na(allele_matrix_genlight)) # check the data dimmension and number of NAs

gen.imp <- apply(allele_matrix_genlight, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x)))))) # impute the NA values into most common allele

sum(is.na(gen.imp)) # check if No NAs

##Load the enviromental data and get it ready for analysis
env <- read.csv("/Users/ericliao/Desktop/WNV_project_files/Ctarsalis_sample_w_GPS_climate_average_new_filtered_id_region.csv", row.names = 1) # load enviromental data

env <- env[, !(names(env) %in% c("locID", "State", "City", "date"))] # remove unused columns

str(env) # check the structure of the enviromental data

env$vcfID <- as.character(env$vcfID) # Make individual names characters (not factors)

identical(rownames(gen.imp), env[,1]) # Confirm that genotypes and environmental data are in the same order

pred <- env[, 18:30] # subset enviromental predictors and shorten their names

pred <- subset(pred, select=-c(avg_ssr, avg_evabs, avg_tp, avg_swvl1, avg_sf))
##Run Redundancy Analysis (RDA): a multivariate GEA

##Run PCA analysis

pred.pca <- rda(pred, scale=T)

summary(pred.pca)$cont # show the summary of PCA

screeplot(pred.pca, main = "Screeplot: Eigenvalues of Mosquitoes Predictor Variables", npcs = 15) # scree plot

round(scores(pred.pca, choices=1:8, display="species", scaling=0), digits=3) # ## correlations between the PC axis and predictors:

pred.PC1 <- scores(pred.pca, choices=1, display="sites", scaling=0) # store our synthetic PC1 axis predictor as pred.PC1 for use in LFMM.
pred.PC2 <- scores(pred.pca, choices = 2, display = "site", scaling = 0) # # store our synthetic PC2 axis predictor as pred.PC1 for use in LFMM.
pred.PC3 <- scores(pred.pca, choices = 3, display = "site", scaling = 0)

##Run LFMM analysis

# # estimate of the number of populations in the data (K).
screeplot(pred.pca, main = "Screeplot of Mosquitoes Predictor Variables with Broken Stick", bstick=TRUE, type="barplot", npcs = 15) # according to broken stick criterion to determine K

gen.pca <- rda(gen.imp, scale=T) # run a PCA with the genommic data and plot the eigenvalues with the broken stick criterion

screeplot(gen.pca, main = "Screeplot of Genetic Data with Broken Stick", bstick=TRUE, type="barplot", npcs = 15) # scree plot with broken stick

K <- 4 # aaccording to admixture analysis

mos.lfmm.PC1 <- lfmm_ridge(Y=gen.imp, X=pred.PC1, K=K) 
mos.lfmm.PC2 <- lfmm_ridge(Y=gen.imp, X=pred.PC2, K=K) 
mos.lfmm.PC3 <- lfmm_ridge(Y=gen.imp, X=pred.PC3, K=K) 

##Identify LFMM candidates using False Discovery Rate

#### for PC1
## pc1 step 1. Look at the genomic inflation factor (GIF), which gives us a sense for how well the model has accounted for confounding factors in the data.
mos.pv.PC1 <- lfmm_test(Y=gen.imp, X=pred.PC1, lfmm=mos.lfmm.PC1, calibrate="gif")
names(mos.pv.PC1) # this object includes raw z-scores and p-values, as well as GIF-calibrated scores and p-values
mos.pv.PC1$gif # Genomic inflation factor (GIF)  1.203221

## pc1 step 2. Plot the p-values to see how application of the GIF influences the p-value distribution.
hist(mos.pv.PC1$pvalue[,1], main="Unadjusted p-values")        
hist(mos.pv.PC1$calibrated.pvalue[,1], main="GIF-adjusted p-values")
# 
# # #####  manually adjust the GIF correction factor #################
# # # Let's change the GIF and readjust the p-values:
# zscore <- mos.pv.PC1$score[,1]   # zscores for first predictor, we only have one in our case...
# (gif <- mos.pv.PC1$gif[1])       ## d.fault GIF for this predictor
# # 
# new.gif1 <- 1.0               ## c.oose your new GIF
# # 
# # # Manual adjustment of the p-values:
# adj.pv1 <- pchisq(zscore^2/new.gif1, df=1, lower = FALSE)
# # 
# # # Plot the p-value histograms:
# hist(mos.pv.PC1$pvalue[,1], main="Unadjusted p-values")
# hist(mos.pv.PC1$calibrated.pvalue[,1], main="GIF-adjusted p-values (GIF=1.0)")
# hist(adj.pv1, main="REadjusted p-values (GIF=1.0)")
# ############## manually adjust the GIF correction factor ###########

## pc1 step 3. convert the adjusted p-values to q-values.
mos.qv.PC1 <- qvalue(mos.pv.PC1$calibrated.pvalue)$qvalues

length(which(mos.qv.PC1 < 0.1)) # how many SNPs have an FDR < 10%?

(mos.FDR.PC1 <- colnames(gen.imp)[which(mos.qv.PC1 < 0.1)]) # identify which SNPs these are

#### follow up with an LFMM model using the second axis as a predictor
## pc2 step 1. Look at the genomic inflation factor (GIF), which gives us a sense for how well the model has accounted for confounding factors in the data.
mos.pv.PC2 <- lfmm_test(Y=gen.imp, X=pred.PC2, lfmm=mos.lfmm.PC2, calibrate="gif")
names(mos.pv.PC2) # this object includes raw z-scores and p-values, as well as GIF-calibrated scores and p-values
mos.pv.PC2$gif # Genomic inflation factor (GIF)

## pc2 step 2. Plot the p-values to see how application of the GIF influences the p-value distribution.
hist(mos.pv.PC2$pvalue[,1], main="Unadjusted p-values")        
hist(mos.pv.PC2$calibrated.pvalue[,1], main="GIF-adjusted p-values")

## pc2 step 3. convert the adjusted p-values to q-values.
mos.qv.PC2 <- qvalue(mos.pv.PC2$calibrated.pvalue)$qvalues

length(which(mos.qv.PC2 < 0.1)) # how many SNPs have an FDR < 10%?

(mos.FDR.PC2 <- colnames(gen.imp)[which(mos.qv.PC2 < 0.1)]) # identify which SNPs these are

#### follow up with an LFMM model using the third axis as a predictor
## pc3 step 1. Look at the genomic inflation factor (GIF), which gives us a sense for how well the model has accounted for confounding factors in the data.
mos.pv.PC3 <- lfmm_test(Y=gen.imp, X=pred.PC3, lfmm=mos.lfmm.PC3, calibrate="gif")
names(mos.pv.PC3) # this object includes raw z-scores and p-values, as well as GIF-calibrated scores and p-values
mos.pv.PC3$gif # Genomic inflation factor (GIF)

## pc3 step 2. Plot the p-values to see how application of the GIF influences the p-value distribution.
hist(mos.pv.PC3$pvalue[,1], main="Unadjusted p-values")        
hist(mos.pv.PC3$calibrated.pvalue[,1], main="GIF-adjusted p-values")

## pc3 step 3. convert the adjusted p-values to q-values.
mos.qv.PC3 <- qvalue(mos.pv.PC3$calibrated.pvalue)$qvalues

length(which(mos.qv.PC3 < 0.1)) # how many SNPs have an FDR < 10%?

(mos.FDR.PC3 <- colnames(gen.imp)[which(mos.qv.PC3 < 0.1)]) # identify which SNPs these are


### output the result of LFMM candidates
LFMM_PC1_df <- as.data.frame(mos.FDR.PC1)
colnames(LFMM_PC1_df) <- c("adaptation_candidates")
LFMM_PC1_df$origin <-"LFMM_PC1"

LFMM_PC2_df <- as.data.frame(mos.FDR.PC2)
colnames(LFMM_PC2_df) <- c("adaptation_candidates")
LFMM_PC2_df$origin <-"LFMM_PC2"

LFMM_PC3_df <- as.data.frame(mos.FDR.PC3)
colnames(LFMM_PC3_df) <- c("adaptation_candidates")
LFMM_PC3_df$origin <-"LFMM_PC3"

combined_LFMM_candidates <- rbind(LFMM_PC1_df, LFMM_PC2_df, LFMM_PC3_df)

# Split the "adaptation_candidates" column into "seq_name" and "location"
split_columns <- strsplit(combined_LFMM_candidates$adaptation_candidates, "_")

# Create new columns "seq_name" and "location"
combined_LFMM_candidates$seq_name <- sapply(split_columns, function(x) x[1])
combined_LFMM_candidates$location <- sapply(split_columns, function(x) x[2])

# Remove the original "adaptation_candidates" column if no longer needed
combined_LFMM_candidates$adaptation_candidates <- NULL

write.csv(combined_LFMM_candidates[combined_LFMM_candidates$origin == "LFMM_PC1", ], "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/LFMM_LatentFactorMixedModels/lfmm_PC1_adaptation_candidates.csv", row.names = FALSE)
write.csv(combined_LFMM_candidates[combined_LFMM_candidates$origin == "LFMM_PC2", ], "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/LFMM_LatentFactorMixedModels/lfmm_PC2_adaptation_candidates.csv", row.names = FALSE)
write.csv(combined_LFMM_candidates[combined_LFMM_candidates$origin == "LFMM_PC3", ], "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/LFMM_LatentFactorMixedModels/lfmm_PC3_adaptation_candidates.csv", row.names = FALSE)

