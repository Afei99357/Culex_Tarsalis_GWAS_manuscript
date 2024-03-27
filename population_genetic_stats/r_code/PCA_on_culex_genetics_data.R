library(vcfR)     # for vcf file
library(adegenet) # for genlight obeject and tab() function
library(vegan)    # Used to run PCA & RDA
library(lfmm)     # Used to run LFMM
library(qvalue)   # Used to post-process LFMM output
library(dplyr) ## for filtering table

rm(list = ls())

# ## read vcf file
vcf <- read.vcfR( "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/original_vcf/bi_20missing_filtSNP_maf_005.recode.vcf", verbose = FALSE )

## Get the names of individual samples
vcfIDs <- as.vector(colnames(vcf@gt))[-1]

## read sample information file
sample_info <- read.csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Ctarsalis_sample_w_GPS_climate_average_new_filtered_id_region.csv", row.names = NULL)

## get popID and vcfID
sample_info_popid <- subset(sample_info, select = c(vcfID, region))

## find Population ID
matched_pop <- sample_info_popid %>% filter(vcfID %in% vcfIDs)

# ## convert to genlight object, store large genomic data much more efficiently.
SNP_genlight <- vcfR2genlight(vcf)

## PCA analysis
pca <- glPca(SNP_genlight, nf=10) # choose 10 axes

# Create a vector to assign colors to populations
levels(sample_info$region) <- c("West Coast","Northwest", "Midwest", "Southwest")
color_vector <- c("goldenrod","skyblue","hotpink","forestgreen") # 4 nice colors for our regions
color_mapping <- color_vector[as.factor(sample_info$region)]

png("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/outlier_analysis/PCA_RESULTS/pca_region_pc1_pc2_new.png", width=20, height=10, units="in", res=300)

par(mar=c(5, 5, 4, 4))

# Calculate the percentage of variance explained by PC1
percentage_var_explained_PC1 <- (pca$eig[1] / sum(pca$eig)) * 100

# Calculate the percentage of variance explained by PC2
percentage_var_explained_PC2 <- (pca$eig[2] / sum(pca$eig)) * 100

# plot pca between pc1 and pc2 by region with filled dots
plot(pca$scores[, 1], pca$scores[, 2], bg = color_mapping, cex = 2, pch = 21, xlab = paste("PC1 (", sprintf("%.2f", percentage_var_explained_PC1), "%)", sep=""),
     ylab = paste("PC2 (", sprintf("%.2f", percentage_var_explained_PC2), "%)", sep=""),
     cex.lab = 1.5)

# add label text
# text(pca$scores[,1], pca$scores[,2] + 0.7, labels=rownames(pca$scores), cex= 0.7)

# add legend and make the legend bigger
legend("topright", legend=levels(sample_info$region), bty="n", col="gray32", pch=21, cex=3, pt.bg=color_vector)

dev.off()

