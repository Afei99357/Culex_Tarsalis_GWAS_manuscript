library(pcadapt)
library(vcfR)
library(ggplot2)
library(plotly)
library(dplyr)
library(viridis)
library(VennDiagram)

rm(list = ls())

data <- read.pcadapt("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/original_vcf/bi_20missing_filtSNP_maf_005.recode.vcf", type = "vcf")

##Load genetic data and convert it to genlight objects, and prepare it for analysis
vcf_pca <- read.vcfR( "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/original_vcf/bi_20missing_filtSNP_maf_005.recode.vcf", verbose = FALSE ) # no missing values

data3 <- vcfR2genlight(vcf_pca)
snp <- as.data.frame(data3@loc.names)
ind <- as.data.frame(data3@ind.names)
colnames(ind) <-"vcfID"

## sample information
pop_map <- read.csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Ctarsalis_sample_w_GPS_climate_average_new_filtered_id_region.csv", row.names = 1, stringsAsFactors = TRUE)

## Merge individuals and pop_map information.
ind_pop_map <- merge(ind, pop_map, by=c("vcfID"))

## First run pcadapt with a large number K of Principal Components, for instance K = 20
data_pcadapt_trial <- pcadapt(data, K = 20, min.maf = 0.01, ploidy=2) 

png("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/outlier_analysis/PCA_RESULTS/scree_plot_pca.png", width = 8, height = 6, units = "in", res = 300)
plot(data_pcadapt_trial, option = "screeplot", col="blue", snp.info = NULL, 
     plt.pkg = "ggplot")
dev.off()

## Check the plot to select the optimal number K of Principal Components
## First, check PC1 and PC2.
png("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/outlier_analysis/PCA_RESULTS/pca_score_plot_pc1_pc2.png", width = 8, height = 6, units = "in", res = 300)
plot(data_pcadapt_trial, option = "scores", pop=pop_map$popID, 
     snp.info = NULL, plt.pkg = "ggplot")
dev.off()

## Second, check PC3 and PC4.
png("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/outlier_analysis/PCA_RESULTS/pca_score_plot_pc3_pc4.png", width = 8, height = 6, units = "in", res = 300)
plot(data_pcadapt_trial, option = "scores", pop=pop_map$popID, 
     i = 3, j = 4, snp.info = NULL, plt.pkg = "ggplot")
dev.off()

## Third, check PC5 and PC6.
png("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/outlier_analysis/PCA_RESULTS/pca_score_plot_pc5_pc6.png", width = 8, height = 6, units = "in", res = 300)
plot(data_pcadapt_trial, option = "scores", pop=pop_map$popID, 
     i = 5, j = 6, snp.info = NULL, plt.pkg = "ggplot")
dev.off()

#### 1. according to the scree plot and the score plots where no significant 
####    population structure was found in the axes after 4. We choose K=4 here


### now Compute the statistic test,
## The test statistic for detecting outlier SNPs is the Mahalanobis distance, 
## which is a multi-dimensional approach that measures how distant a point is from the mean.

## By default, the parameter min.maf is set to 5%; here we are changing it to 1%. 
## p-values of SNPs with a minor allele frequency smaller than the threshold are not computed (NA is returned).
data_pcadapt <- pcadapt(data, K = 4, min.maf = 0.01) 

## Create a dataframe gathering pcadapt results and sampling sites.
pca_adapt_pop_map <-cbind(data_pcadapt$scores, ind_pop_map)
colnames(pca_adapt_pop_map) <- c("PC1","PC2","PC3","PC4","IND","SITES", "REGIONS")

pca_adapt_pop_map <- pca_adapt_pop_map[,1:7]

pca_adapt_pop_map$SITES <- ind_pop_map$popID

## get the list of p-values.
snps_pvalues <- cbind(snp, data_pcadapt$pvalues)
snps_pvalues_no_na <- na.omit(snps_pvalues)
write.table(snps_pvalues, "All_Pvalues.txt", sep="\t", quote=FALSE)

## Change the value of the P_VALUE column: 0 == 0.0001
attach(data_pcadapt)
class(data_pcadapt$pvalues)

data_pcadapt$pvalues <- as.numeric(data_pcadapt$pvalues) 
data_pcadapt$pvalues[data_pcadapt$pvalues == 0] <- 0.0001

## plot Q-Q plot
### Q-Q plot confirms that most of the p-values follow the expected uniform distribution.
png("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/outlier_analysis/PCA_RESULTS/qqplot_pca.png", width=15, height=10, units="in", res=300)
plot(data_pcadapt, option = "qqplot", col, snp.info = NULL, plt.pkg = "ggplot")
dev.off()


## Check the distribution of these p-values.
### A histogram of p-values confirms that most of the p-values follow a uniform distribution. 
### The excess of small P-values indicate the presence of outliers.
png("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/outlier_analysis/PCA_RESULTS/pvalue_histogram_pca.png", width=15, height=10, units="in", res=300)
hist(data_pcadapt$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
dev.off()


## Visualize the distribution of p-values.
quantile(snps_pvalues_no_na$`data_pcadapt$pvalues`, probs = c(0.01, 0.99))
###            1%          99% 
##   1.799655e-46 9.960129e-01 

## Get only the markers showing extreme p-values: the top 1%.
top_1percent <- subset(snps_pvalues_no_na, snps_pvalues_no_na$`data_pcadapt$pvalues` <= 1.799655e-46)
colnames(top_1percent) <- c("LOCUS","PVALUE")

write.table(top_1percent, "173outliers_pcadapt.txt", sep="\t", quote=FALSE, row.names = FALSE)

## Save pcadapt names in a vector.
pcadapt_outliers <- top_1percent$LOCUS
## Check the number of outliers found by pcadapt.
length(pcadapt_outliers)

## get the snps are not in the candidate list
non_cand <- colnames(allele_matrix_genlight)[!colnames(allele_matrix_genlight) %in% pcadapt_outliers]

## create a data frame to store the list of non-candidate SNPs
non_cand <- cbind.data.frame(non_cand)

## colnames is snp
colnames(non_cand) <- c("snp")

## output non_cand to a file 
write.csv(non_cand, file = "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/outlier_analysis/PCAdapt_RESULTS/pcadapt_non_candidates_snps.csv", row.names = FALSE)

