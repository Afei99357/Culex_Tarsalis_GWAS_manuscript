rm(list = ls())

##Loading library
library(vcfR)     # for vcf file
library(adegenet) # for genlight obeject and tab() function
library(vegan)    # Used to run PCA & RDA
library(qvalue)   # Used to post-process LFMM output
library(psych)
library(ggplot2)

##Load genetic data and convert it to genlight objects, and prepare it for analysis

vcf <- read.vcfR( "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/original_vcf/bi_20missing_filtSNP_maf_005.recode.vcf", verbose = FALSE ) # no missing values

SNP_genlight <- vcfR2genlight(vcf) # convert to genlight objects

allele_matrix_genlight <- tab(SNP_genlight, freq = FALSE, NA.method = "asis") # get the matrix of allele counts

dim(allele_matrix_genlight) # check the dimension of genetic data

sum(is.na(allele_matrix_genlight)) # check the data dimmension and number of NAs

gen.imp <- apply(allele_matrix_genlight, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x)))))) # impute the NA values into most common allele

sum(is.na(gen.imp)) # check if No NAs

##Load the enviromental data and get it ready for analysis
env <- read.csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Ctarsalis_sample_w_GPS_climate_average_new_filtered_id_region.csv", row.names = 1) # load enviromental data

env <- env[, !(names(env) %in% c("locID", "State", "City", "date"))] # remove unused columns

str(env) # check the structure of the enviromental data

env$vcfID <- as.character(env$vcfID) # Make individual names characters (not factors)

identical(rownames(gen.imp), env[,1]) # Confirm that genotypes and environmental data are in the same order

pred <- env[, 18:30] # subset enviromental predictors and shorten their names

# ## rename the predictors
new_names <- c("eastward_wind", "northward_wind", "temperature",
               "evaporation_from_bare_soil", "high_vegetation",
               "low_vegetation", "water_retention_capacity",
               "snowfall", "surface_net_solar_radiation", "surface_runoff",
               "evaporation", "total_precipitation", "volumetric_soil_water_layer1")

## renmae the predictors
colnames(pred) <- new_names

# check the Linear dependency
usdm::vif(pred)

png("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/variable_pair_correlation_panel.png", width = 8, height = 6, units = "in", res = 300)
pairs.panels(pred, scale = TRUE, 
             cex = 1.5)       # Controls overall text size

dev.off()

## remove as little variables as we can based on any pair if the correlation is higher than 0.7
pred <- subset(pred, select=-c(surface_net_solar_radiation, evaporation_from_bare_soil, total_precipitation, volumetric_soil_water_layer1, snowfall))

# check the Linear dependency
usdm::vif(pred)

##Run Redundancy Analysis (RDA): a multivariate GEA

png("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/variable_pair_correlation_panel_after_reduced.png", width = 8, height = 6, units = "in", res = 300)
pairs.panels(pred, scale = TRUE, 
             cex = 1.5)       # Controls overall text size

dev.off()

mos.rda <- rda(gen.imp ~ ., data=pred, scale=T) # ## Run RDA, Observe the shortcut (.) to tell the function to use all variables present in data gene, without having to enumerate them.
mos.rda

RsquareAdj(mos.rda) # # R2 will be biased and should be adjusted based on the number of predictors


# which environmental predictors are most strongly correlated with the first 4 RDA axes
intersetcor(mos.rda)[,1:4]

summary_rda <- summary(mos.rda)
summary_rda$concont ## The eigenvalues for the constrained axes reflect the variance explained by each canonical axis

## calculate each axis of the first 5 axes explains the variance of the predictors
rda1_var <- summary_rda$concon$importance[1,1]/ sum(summary_rda$concon$importance)
rda2_var <- summary_rda$concon$importance[1,2]/ sum(summary_rda$concon$importance)
rda3_var <- summary_rda$concon$importance[1,3]/ sum(summary_rda$concon$importance)
rda4_var <- summary_rda$concon$importance[1,4]/ sum(summary_rda$concon$importance)
rda5_var <- summary_rda$concon$importance[1,5]/ sum(summary_rda$concon$importance)


png("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/RDA_scree_plot.png", width = 8, height = 6, units = "in", res = 300)
screeplot(mos.rda) # screeplot of the canonical eigenvalues
dev.off()

# run a formal test of statistical significance of each constrained axis
# # The null hypothesis is that no linear relationship exists between the SNP data and the environmental predictors
# signif.full <- anova.cca(mos.rda, by="axis")
# signif.full

# Scaling 1
plot(mos.rda,
     scaling = 1,
     display = c("sp", "lc", "cn"),
     main = "Triplot RDA spe.hel ~ env3 - scaling 1 - lc scores"
)

# # checking Variance Inflation Factors for the predictor variables used in the model
vif.cca(mos.rda)

### plot the RDA
plot(mos.rda, scaling=3) ## plot of the RDA output, default is axes 1 and 2

# # see where these candidate SNPs are in the ordination space
# Set up the color scheme for plotting:
levels(env$region) <- levels(as.factor(env$region))
color_vector <- c("hotpink","skyblue","forestgreen","goldenrod") # 4 nice colors for our regions
# 9 nice colors for our regions
color_mapping <- color_vector[as.factor(env$region)]

cex_single_text <- 1
fon_single_text <- 2


# axes 1 & 2
png("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/RDA_region_ax_1_2.png", width = 8, height = 6, units = "in", res = 300)
plot(mos.rda, type="n", scaling=3,
     xlab=paste("RDA1 (", round(rda1_var * 100, 2), "%)"),
     ylab=paste("RDA2 (", round(rda2_var * 100, 2), "%)"))
points(mos.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(mos.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=color_mapping) # the mosquitoes
text(mos.rda, scaling=3, display="bp", col="purple", cex=cex_single_text, font=fon_single_text)                           # the predictors
legend("bottomright", legend=levels(env$region), bty="n", col="gray32", pch=21, cex=cex_single_text, pt.bg=color_vector)

dev.off()

# axes 1 & 3
png("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/RDA_region_ax_1_3.png", width = 8, height = 6, units = "in", res = 300)
plot(mos.rda, type="n", scaling=3, choices=c(1,3), 
     xlab=paste("RDA1 (", round(rda1_var * 100, 2), "%)"),
     ylab=paste("RDA3 (", round(rda3_var * 100, 2), "%)"))
points(mos.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3, choices=c(1,3))
points(mos.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=color_mapping, choices=c(1,3))
text(mos.rda, scaling=3, display="bp", col="purple", cex=cex_single_text, choices=c(1,3), font=fon_single_text)
legend("bottomright", legend=levels(env$region), bty="n", col="gray32", pch=21, cex=cex_single_text, pt.bg=color_vector)
dev.off()

# axes 2 & 3
png("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/RDA_region_ax_2_3.png", width = 8, height = 6, units = "in", res = 300)
plot(mos.rda, type="n", scaling=3, choices=c(2,3),
     xlab=paste("RDA2 (", round(rda2_var * 100, 2), "%)"),
     ylab=paste("RDA3 (", round(rda3_var * 100, 2), "%)"))
points(mos.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3, choices=c(2,3))
points(mos.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=color_mapping, choices=c(2,3))
text(mos.rda, scaling=3, display="bp", col="purple", cex=cex_single_text, font=fon_single_text, choices=c(2,3))
legend("bottomright", legend=levels(env$region), bty="n", col="gray32", pch=21, cex=cex_single_text, pt.bg=color_vector)
dev.off()

# axes 3 & 4
png("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/RDA_region_ax_3_4.png", width = 8, height = 6, units = "in", res = 300)
plot(mos.rda, type="n", scaling=3, choices=c(3,4),
     xlab=paste("RDA3 (", round(rda3_var * 100, 2), "%)"),
     ylab=paste("RDA4 (", round(rda4_var * 100, 2), "%)"))
points(mos.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3, choices=c(3,4))
points(mos.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=color_mapping, choices=c(3,4))
text(mos.rda, scaling=3, display="bp", col="purple", cex=cex_single_text, font=fon_single_text, choices=c(2,3))
legend("bottomright", legend=levels(env$region), bty="n", col="gray32", pch=21, cex=cex_single_text, pt.bg=color_vector)
dev.off()

# axes 4 & 5
png("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/RDA_region_ax_4_5.png", width = 8, height = 6, units = "in", res = 300)
plot(mos.rda, type="n", scaling=3, choices=c(4,5),
     xlab=paste("RDA4 (", round(rda4_var * 100, 2), "%)"),
     ylab=paste("RDA5 (", round(rda5_var * 100, 2), "%)"))
points(mos.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3, choices=c(4,5))
points(mos.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=color_mapping, choices=c(4,5))
text(mos.rda, scaling=3, display="bp", col="purple", cex=cex_single_text, font=fon_single_text, choices=c(2,3))
legend("bottomright", legend=levels(env$region), bty="n", col="gray32", pch=21, cex=cex_single_text, pt.bg=color_vector)
dev.off()

### Identify candidate SNPs involved in local adaptation
load.rda <- scores(mos.rda, choices=c(1:4), display="species") ## indentify RDA candidates, based on scree plot of the canonical eigenvalues, extract the SNP loadings from the first 4 constrained axes
png("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/RDA_ax_1_loading_histogram.png", width = 8, height = 6, units = "in", res = 300)
hist(load.rda[,1], main="Loadings on RDA1") # histograms of the loadings on RDA axis 1
dev.off()

png("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/RDA_ax_2_loading_histogram.png", width = 8, height = 6, units = "in", res = 300)
hist(load.rda[,2], main="Loadings on RDA2") # histograms of the loadings on RDA axis 2
dev.off()

png("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/RDA_ax_3_loading_histogram.png", width = 8, height = 6, units = "in", res = 300)
hist(load.rda[,3], main="Loadings on RDA3") # histograms of the loadings on RDA axis 3
dev.off()

png("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/RDA_ax_4_loading_histogram.png", width = 8, height = 6, units = "in", res = 300)
hist(load.rda[,4], main="Loadings on RDA4") # histograms of the loadings on RDA axis 4
dev.off()

## function to identify SNPs that load in the tails of these distributions
## define the function here as outliers, where x is the vector of loadings and z is the number of standard deviations to use
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x) ## f.nd loadings +/- z SD from mean loading     
  x[x < lims[1] | x > lims[2]]           # locus names in these tails
}

cand1 <- outliers(load.rda[,1], 3) 
cand2 <- outliers(load.rda[,2], 3) 
cand3 <- outliers(load.rda[,3], 3) 
cand4 <- outliers(load.rda[,4], 3) 


ncand <- length(cand1) + length(cand2) +  length(cand3) + length(cand4)

# organize our results by making one data frame with the axis, SNP name, loading, & correlation with each predictor:
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
cand4 <- cbind.data.frame(rep(4,times=length(cand4)), names(cand4), unname(cand4))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- colnames(cand4) <- c("axis","snp","loading")

cand <- rbind(cand1, cand2, cand3, cand4)
cand$snp <- as.character(cand$snp)

## get the snps are not in the candidate list
non_cand <- colnames(allele_matrix_genlight)[!colnames(allele_matrix_genlight) %in% cand$snp]

## create a data frame to store the list of non-candidate SNPs
non_cand <- cbind.data.frame(non_cand)

## colnames is snp
colnames(non_cand) <- c("snp")

## output non_cand to a file 
write.csv(non_cand, file = "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/rda_non_candidate_snps.csv", row.names = FALSE)

#  add in the correlations of each candidate SNP with the eight environmental predictors:
foo <- matrix(nrow=(ncand), ncol=8)  # 9 columns for 9 predictors
colnames(foo) <- colnames(pred)

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- gen.imp[,nam]
  foo[i,] <- apply(pred,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)
head(cand)

# Investigate the candidates
length(cand$snp[duplicated(cand$snp)])  # 10 duplicate detections

foo <- cbind(cand$axis, duplicated(cand$snp)) 
table(foo[foo[,1]==1,2]) # no duplicates on axis 1
table(foo[foo[,1]==2,2]) #  no duplicates on axis 2
table(foo[foo[,1]==3,2]) # 10 duplicates on axis 3
table(foo[foo[,1]==4,2]) # 8 duplicates on axis 4

cand <- cand[!duplicated(cand$snp),] # remove duplicate detections

# # see which of the predictors each candidate SNP is most strongly correlated with
for (i in 1:length(cand$snp)) {
  
  bar <- cand[i,]
  cand[i,12] <- names(which.max(abs(bar[4:11]))) # gives the variable
  cand[i,13] <- max(abs(bar[4:11]))              # gives the correlation
}

colnames(cand)[12] <- "predictor"
colnames(cand)[13] <- "correlation"

## output candidate SNPs
write.csv(cand, file = "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/rda_analysis_candidates_correlation_to_env.csv", row.names = FALSE)

## read groupsection csv file
group_section_df <- read.csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/rda_analysis_candidates_correlation_to_env.csv")

candidates_pred_sum <- table(cand$predictor)

## adding a column to cand using the column "axis" from group_section_df
cand$section_group <- group_section_df$axis[match(cand$snp, group_section_df$snp)]

# Convert the table to a data frame for plotting
plot_data <- as.data.frame(candidates_pred_sum)
# Name the columns appropriately
colnames(plot_data) <- c("predictor", "count")

png("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/RDA_candidate_counts_variable_cor.png", width = 8, height = 6, units = "in", res = 300)
ggplot(plot_data, aes(x = predictor, y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  xlab("Predictor Variable") +
  ylab("Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

## output result
write.csv(cand, file = "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/rda_analysis_candidates.csv", row.names = FALSE)

### Plot the SNPs
sel <- cand$snp
envi <- cand$predictor
envi[envi=="evaporation"] <- '#000000'
envi[envi=="eastward_wind"] <- '#E69F00'
envi[envi=="northward_wind"] <- '#F0E442'
envi[envi=="temperature"] <- '#56B4E9'
envi[envi=="high_vegetation"] <- '#009E73'
envi[envi=="low_vegetation"] <- '#D55E00'
envi[envi=="water_retention_capacity"] <- '#CC79A7'
envi[envi=="surface_runoff"] <- '#0072B2'

group_section <- cand$section_group
group_section[group_section==1] <- '#a6611a'
group_section[group_section==2] <- '#1E88E5'
group_section[group_section==3] <- '#FFC107'
group_section[group_section==4] <- '#80cdc1'

# pull the SNP names
col.pred <- rownames(mos.rda$CCA$v) 
col.group <- rownames(mos.rda$CCA$v)

# color code candidate SNPs by candidate snps
for (i in 1:length(sel)) {          
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- envi[i]
}

# color code candidate SNPs by group section
for (i in 1:length(sel)) {          
  foo <- match(sel[i],col.group)
  col.group[foo] <- group_section[i]
}



# non-candidate SNPs
col.pred[grep("tig",col.pred)] <- '#f1eef6'
empty <- col.pred

col.group[grep("tig",col.group)] <- '#f1eef6'
empty.group <- col.group

## make the non_candidate SNPs transparent
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.group[grep("#f1eef6",empty.group)] <- rgb(0,1,0, alpha=0) # transparent

## make the candidate SNPs outline color gray, non-candidate outline color transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
empty.group.outline <- ifelse(empty.group=="#00FF0000","#00FF0000","gray32")

bg <- c( '#000000','#E69F00','#F0E442','#56B4E9','#009E73', '#D55E00','#CC79A7','#0072B2')
bg.group <- c('#a6611a','#1E88E5','#FFC107','#80cdc1')


legends <- c("RDA1 SNP Group", "RDA2 SNP Group", "RDA3 SNP Group", "RDA4 SNP Group")


# plot the SNPs axes 1 & 2
png("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/RDA_snps_plot_ax_1_2.png", width = 8, height = 6, units = "in", res = 300)
plot(mos.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1),
     xlab=paste("RDA1 (", round(rda1_var * 100, 2), "%)"),
     ylab=paste("RDA2 (", round(rda2_var * 100, 2), "%)"))
points(mos.rda, display="species", pch=21, cex=1, col="white", bg=col.pred, scaling=3)
points(mos.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(mos.rda, scaling=3, display="bp", col="purple", cex=cex_single_text, font=fon_single_text)
legend("topright", legend=c("evaporation", "eastward_wind", "northward_wind", "temperature", "high_vegetation", "low_vegetation", "water_retention_capacity", "surface_runoff"), bty="n", col="gray32", pch=21, cex=cex_single_text, pt.bg=bg)
dev.off()

# plot the SNPs axes 1 and 2 with group section color
png("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/RDA_snps_plot_ax_1_2_group_section.png", width = 8, height = 6, units = "in", res = 300)
plot(mos.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1),
     xlab=paste("RDA1 (", round(rda1_var * 100, 2), "%)"),
     ylab=paste("RDA2 (", round(rda2_var * 100, 2), "%)"))
points(mos.rda, display="species", pch=21, cex=1, col="white", bg=col.group, scaling=3)
points(mos.rda, display="species", pch=21, cex=1, col=empty.group.outline, bg=empty.group, scaling=3)
text(mos.rda, scaling=3, display="bp", col="purple", cex=cex_single_text, font=fon_single_text)
legend("topright", legend=legends, bty="n", col="gray32", pch=21, cex=cex_single_text, pt.bg=bg.group)
dev.off()

# axes 1 & 3
png("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/RDA_snps_plot_ax_1_3.png", width = 8, height = 6, units = "in", res = 300)
plot(mos.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(1,3),
     xlab=paste("RDA1 (", round(rda1_var * 100, 2), "%)"),
     ylab=paste("RDA3 (", round(rda3_var * 100, 2), "%)"))
points(mos.rda, display="species", pch=21, cex=1, col="white", bg=col.pred, scaling=3, choices=c(1,3))
points(mos.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3, choices=c(1,3))
text(mos.rda, scaling=3, display="bp", col="purple", cex=cex_single_text, font=fon_single_text, choices=c(1,3))
legend("topright", legend=c("evaporation", "eastward_wind", "northward_wind", "temperature", "high_vegetation", "low_vegetation", "water_retention_capacity", "surface_runoff"), bty="n", col="gray32", pch=21, cex=cex_single_text, pt.bg=bg)
dev.off()

# plot the SNPs axes 1 and 3 with group section color
png("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/RDA_snps_plot_ax_1_3_group_section.png", width = 8, height = 6, units = "in", res = 300)
plot(mos.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(1,3),
     xlab=paste("RDA1 (", round(rda1_var * 100, 2), "%)"),
     ylab=paste("RDA3 (", round(rda3_var * 100, 2), "%)"))
points(mos.rda, display="species", pch=21, cex=1, col="white", bg=col.group, scaling=3, choices=c(1,3))
points(mos.rda, display="species", pch=21, cex=1, col=empty.group.outline, bg=empty.group, scaling=3, choices=c(1,3))
text(mos.rda, scaling=3, display="bp", col="purple", cex=cex_single_text, font=fon_single_text, choices=c(1,3))
legend("topright", legend=legends, bty="n", col="gray32", pch=21, cex=cex_single_text, pt.bg=bg.group)
dev.off()

# axes 2 & 3
png("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/RDA_snps_plot_ax_2_3.png", width = 8, height = 6, units = "in", res = 300)
plot(mos.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(2,3),
     xlab=paste("RDA2 (", round(rda2_var * 100, 2), "%)"),
     ylab=paste("RDA3 (", round(rda3_var * 100, 2), "%)"))
points(mos.rda, display="species", pch=21, cex=1, col="white", bg=col.pred, scaling=3, choices=c(2,3))
points(mos.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3, choices=c(2,3))
text(mos.rda, scaling=3, display="bp", col="purple", cex=cex_single_text, font=fon_single_text, choices=c(2,3))
legend("topright", legend=c("evaporation", "eastward_wind", "northward_wind", "temperature", "high_vegetation", "low_vegetation", "water_retention_capacity", "surface_runoff"), bty="n", col="gray32", pch=21, cex=cex_single_text, pt.bg=bg)
dev.off()

# plot the SNPs axes 2 and 3 with group section color
png("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/RDA_snps_plot_ax_2_3_group_section.png", width = 8, height = 6, units = "in", res = 300)
plot(mos.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(2,3),
     xlab=paste("RDA2 (", round(rda2_var * 100, 2), "%)"),
     ylab=paste("RDA3 (", round(rda3_var * 100, 2), "%)"))
points(mos.rda, display="species", pch=21, cex=1, col="white", bg=col.group, scaling=3, choices=c(2,3))
points(mos.rda, display="species", pch=21, cex=1, col=empty.group.outline, bg=empty.group, scaling=3, choices=c(2,3))
text(mos.rda, scaling=3, display="bp", col="purple", cex=cex_single_text, font=fon_single_text, choices=c(2,3))
legend("topright", legend=legends, bty="n", col="gray32", pch=21, cex=cex_single_text, pt.bg=bg.group)
dev.off()


# axes 3 & 4
png("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/RDA_snps_plot_ax_3_4.png", width = 8, height = 6, units = "in", res = 300)
plot(mos.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(3,4),
     xlab=paste("RDA3 (", round(rda3_var * 100, 2), "%)"),
     ylab=paste("RDA4 (", round(rda4_var * 100, 2), "%)"))
points(mos.rda, display="species", pch=21, cex=1, col="white", bg=col.pred, scaling=3, choices=c(3,4))
points(mos.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3, choices=c(3,4))
text(mos.rda, scaling=3, display="bp", col="purple", cex=cex_single_text, font=fon_single_text, choices=c(2,3))
legend("topright", legend=c("evaporation", "eastward_wind", "northward_wind", "temperature", "high_vegetation", "low_vegetation", "water_retention_capacity", "surface_runoff"), bty="n", col="gray32", pch=21, cex=cex_single_text, pt.bg=bg)
dev.off()

# plot the SNPs axes 3 and 4 with group section color
png("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/RDA_snps_plot_ax_3_4_group_section.png", width = 8, height = 6, units = "in", res = 300)
plot(mos.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(3,4),
     xlab=paste("RDA3 (", round(rda3_var * 100, 2), "%)"),
     ylab=paste("RDA4 (", round(rda4_var * 100, 2), "%)"))
points(mos.rda, display="species", pch=21, cex=1, col="white", bg=col.group, scaling=3, choices=c(3,4))
points(mos.rda, display="species", pch=21, cex=1, col=empty.group.outline, bg=empty.group, scaling=3, choices=c(3,4))
text(mos.rda, scaling=3, display="bp", col="purple", cex=cex_single_text, font=fon_single_text, choices=c(2,3))
legend("topright", legend=legends, bty="n", col="gray32", pch=21, cex=cex_single_text, pt.bg=bg.group)
dev.off()



###plot a 2*2 plots
cex_multi_text = 2.3
font_multi_text = 3
legend_cex_text = 2.5
vector_line_thick = 2

## create a panel 2 x 2 to store four plots
png("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/RDA_region_and_SNPs_section_group.png", width = 34, height = 22, units = "in", res = 300)
par(mfrow=c(2,2))
par(family='Arial', cex.lab=2, cex.axis=2, lwd=2)
par(mar=c(5,5,5,5))

plot(mos.rda, type="n", scaling=3, xlab=paste("RDA1 (", round(rda1_var * 100, 2), "%)"),
     ylab=paste("RDA2 (", round(rda2_var * 100, 2), "%)"))
points(mos.rda, display="species", pch=20, cex=2, col="gray32", scaling=3)           # the SNPs
points(mos.rda, display="sites", pch=21, cex=2, col="gray32", scaling=3, bg=color_mapping) # the mosquitoes
text(mos.rda, scaling=3, display="bp", col="purple", cex=cex_multi_text, font=font_multi_text, lwd=vector_line_thick)                           # the predictors
legend("topright", legend=levels(env$region), bty="n", col="gray32", pch=21, cex=legend_cex_text, pt.bg=color_vector)
mtext("(A)", side=3, line=0.5, adj=0, cex=2.5, font=1) # Add label to the first plot


plot(mos.rda, type="n", scaling=3, choices=c(3,4), xlab=paste("RDA3 (", round(rda3_var * 100, 2), "%)"),
     ylab=paste("RDA4 (", round(rda4_var * 100, 2), "%)"))
points(mos.rda, display="species", pch=20, cex=2, col="gray32", scaling=3, choices=c(3,4))
points(mos.rda, display="sites", pch=21, cex=2, col="gray32", scaling=3, bg=color_mapping, choices=c(3,4))
text(mos.rda, scaling=3, display="bp", col="purple",cex=cex_multi_text, font=font_multi_text, choices=c(3,4), lwd=vector_line_thick)
legend("topright", legend=levels(env$region), bty="n", col="gray32", pch=21, cex=legend_cex_text, pt.bg=color_vector)
mtext("(B)", side=3, line=0.5, cex=2.5, adj=0, font=1)

plot(mos.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), xlab=paste("RDA1 (", round(rda1_var * 100, 2), "%)"),
     ylab=paste("RDA2 (", round(rda2_var * 100, 2), "%)"))
points(mos.rda, display="species", pch=21, cex=2, col="white", bg=col.group, scaling=3)
points(mos.rda, display="species", pch=21, cex=2, col=empty.group.outline, bg=empty.group, scaling=3)
text(mos.rda, scaling=3, display="bp", col="purple", cex=cex_multi_text, font=font_multi_text, lwd=vector_line_thick)
legend("topright", legend=legends, bty="n", col="gray32", pch=21, cex=legend_cex_text, pt.bg=bg.group)
mtext("(C)", side=3, line=0.5, cex=2.5, adj=0, font=1)

plot(mos.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(3,4), xlab=paste("RDA3 (", round(rda3_var * 100, 2), "%)"),
     ylab=paste("RDA4 (", round(rda4_var * 100, 2), "%)"))
points(mos.rda, display="species", pch=21, cex=2, col="white", bg=col.group, scaling=3, choices=c(3,4))
points(mos.rda, display="species", pch=21, cex=2, col=empty.group.outline, bg=empty.group, scaling=3, choices=c(3,4))
text(mos.rda, scaling=3, display="bp", col="purple",cex=cex_multi_text, font=font_multi_text, choices=c(3,4), lwd=vector_line_thick)
legend("topright", legend=legends, bty="n", col="gray32", pch=21, cex=legend_cex_text, pt.bg=bg.group)
mtext("(D)", side=3, line=0.5, cex=2.5, adj=0, font=1)

dev.off()

