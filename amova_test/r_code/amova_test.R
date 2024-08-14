library("vcfR")     # for vcf file
library("adegenet") # for genlight obeject and tab() function
library("dplyr") ## for filtering table
library("hierfstat")
library("poppr")
library("ade4")

rm(list = ls())

# ## SNP data
vcf <- read.vcfR( "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/original_vcf/bi_20missing_filtSNP_maf_005.recode.vcf", verbose = FALSE )

## Get the names of individual samples
vcfIDs <- as.vector(colnames(vcf@gt))[-1]

## enviromental data
## read sample information file
env <- read.csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Ctarsalis_sample_w_GPS_climate_average_new_filtered_id_region.csv")

sample_info_popid <- subset(env, select = c(vcfID, popID, region)) ## get popID and vcfID

matched_pop <- sample_info_popid %>% filter(vcfID %in% vcfIDs) ## find Population ID

population <- matched_pop$popID ## get the population info
regions <- matched_pop$region ## get the region infor
indivduals <- as.character(matched_pop$vcfID) ## get the individual sample name

SNP_genind <- vcfR2genind(vcf, return.alleles = TRUE, ind.names = matched_pop$vcfID) ## convert to genind object

# SNP_genind$pop <- as.factor(matched_pop$region) ## add region information as pop information to genind object
SNP_genind$pop <- as.factor(matched_pop$region) ## between populations

SNP_genind@call$ind.names <- indivduals ## add sample individual names to genind objects

strata_df <- data.frame(regions, population) 
strata(SNP_genind) <- strata_df ## add strata

# # Genetic diversity (observed and expected heterozygosity)
div <- summary(SNP_genind)

## the plot visualizes how observed heterozygosity varies as a function of the number of loci, providing insights into genetic diversity within data
png("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/basic_population_genetic_statistics/amova_test/observed_heterozygosity_plot.png", width = 6, height = 6, units = "in", res = 300)
plot(div$Hobs, xlab="Loci number", ylab="Observed Heterozygosity", main="Observed heterozygosity per locus") 
dev.off()

## displays expected heterozygosity as a function of observed heterozygosity per locus
png("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/basic_population_genetic_statistics/amova_test/heterozygosity_observed_heterozygosity_plot.png", width = 6, height = 6, units = "in", res = 300)
plot(div$Hobs, div$Hexp, xlab="Observed Heterozygosity", ylab="Expect Heterozygosity under HWE", main="heterozygosity-observed heterozygosity plot")
dev.off()


## Hierachical Fst tests (= AMOVA for SNP dataset)
# loci <- hfstat[, -1] # Remove the population column
# hierarchical_fst <- varcomp.glob(levels = data.frame(population, regions), loci, diploid = TRUE) # The function varcomp.glob() produces a Hierarchical Fst (=AMOVA for SNPs or bi-allelic markers) 
amova.result <- poppr::poppr.amova(SNP_genind, hier=~regions/population, within=FALSE)
amova.result

## in the result, "Variations within samples" refer to "the variance within populations,
## "Variations between samples" refer to "the variance between populations"
## "Variations between regions" refers to "the variance between different regions"
amova.test <- randtest(amova.result, nrepet = 999) # Test for significance
amova.test

## modify the default Test name
amova.test$names <- c("(A) Variations within Populations", "(B) Variations between Populations within the Same Region", "(C) Variations between Regions")

png("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/basic_population_genetic_statistics/amova_test/AMOVA_significant_test_plot.png", width = 6, height = 6, units = "in", res = 300)

# font size of plot title, no bold
par(family = "arial", font.lab = 2, font.axis = 2, font.main = 2, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.1)

# plot the result and change the font for title
plot(amova.test)

# Close the device to save the plot
dev.off()

