library("adegenet")
library("hierfstat")
library("pegas")
library(dplyr) ## for filtering table
library(vcfR)

### Install required package devtools if not installed
if (!("devtools" %in% installed.packages())){install.packages("devtools")}
library(devtools)

rm(list = ls())

## #####Load genetic data and convert it to genlight objects, and prepare it for analysis
# vcf <- read.vcfR( "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/ct_filtInd.recode.vcf", verbose = FALSE, index_col= FALSE ) # read vcf file 
vcf <- read.vcfR( "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/original_vcf/bi_20missing_filtSNP_maf_005.recode.vcf", verbose = FALSE) # no missing values

## Get the names of individual samples
vcfIDs <- as.vector(colnames(vcf@gt))[-1]

## read sample information file
env <- read.csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Ctarsalis_sample_w_GPS_climate_average_new_filtered_id_region.csv", row.names = NULL)

## get popID and vcfID
sample_info_popid <- subset(env, select = c(vcfID, popID, region))

## find Population ID
matched_pop <- sample_info_popid %>% filter(vcfID %in% vcfIDs)

## get the population and region information
population <- matched_pop$popID
regions <- matched_pop$region
indivduals <- as.character(matched_pop$vcfID)

## convert to genind object
SNP_genind <- vcfR2genind(vcf, return.alleles = TRUE, ind.names = matched_pop$vcfID)

## add pop information to genind object
SNP_genind$pop <- as.factor(population)

# add sample individual names to genind objects
SNP_genind@call$ind.names <- indivduals

# add strata
strata_df <- data.frame(regions, population)
strata(SNP_genind) <- strata_df

SNP_genind ## Check that markers are polymorphic

nAll(SNP_genind) # Number of alleles per locus

# # Genetic diversity (observed and expected heterozygosity)
div <- summary(SNP_genind)
names(div)

plot(div$Hobs, xlab="Loci number", ylab="Observed Heterozygosity", main="Observed heterozygosity per locus")

# # displays expected heterozygosity as a function of observed heterozygosity per locus
plot(div$Hobs, div$Hexp, xlab="Observed Heterozygosity", ylab="Expect Heterozygosity under HWE", main="heterozygosity-observed heterozygosity plot")

bartlett.test(list(div$Hexp, div$Hobs)) # a test : H0: Hexp = Hobs. from the result, p-values < 2.2e -16 which means the variances of div$Hexp and div$Hobs are significantly different.

## Converts genind objects from adegenet into a hierfstat data frame
hfstat <- genind2hierfstat(SNP_genind)

basicstat <- basic.stats(hfstat, diploid = TRUE, digits = 2) 
names(basicstat) ## the observed heterozygosity (Ho), mean gene diversities within population (Hs), Fis, and Fst.

# # confidence interval for Fis
boot.ppfis(hfstat) 

# # PCA on a matrix of individuals genotypes frequencies
x <- indpca(hfstat) 
plot(x, cex = 0.7)

# # Testing for Hardy-Weinberg Equilibrium
# # hw.test(SNP_genind, B = 10)

# Weir and Cockerham's Fst estimate ()
wc(SNP_genind) ## FST : 0.05119577.  FIS:0.09179463

## Hierachical Fst tests (= AMOVA for SNP dataset)
# loci <- hfstat[, -1] # Remove the population column
# hierarchical_fst <- varcomp.glob(levels = data.frame(population, regions), loci, diploid = TRUE) # The function varcomp.glob() produces a Hierarchical Fst (=AMOVA for SNPs or bi-allelic markers) 
amova.result <- poppr::poppr.amova(SNP_genind, hier=~regions/population, within=FALSE)

## in the result, "Variations within samples" refer to " within populations within the same region",
## "Variations between samples" refer to "the genetic differentiation between populations within the same region"
## "Variations between regions" refers to "the genetic variation between different regions"
amova.test <- randtest(amova.result) # Test for significance
plot(amova.test)

## Pairwise Fst method option: y “Nei87” and “WC84” return pairwise FSTs estimated following Nei (1987) pairwise.neifst and Weir & Cockerham (1984) pp.fst respectively
genet.dist(SNP_genind, method = "WC84")
