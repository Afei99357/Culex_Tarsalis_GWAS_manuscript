library(adegenet)
library(ggplot2)
library(radiator)
library(vcfR)

rm(list = ls())

# setwd("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/R_code")

bayescan=read.table("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/BayeScan/output_result/dat_fst.txt") 

##Load genetic data and convert it to genlight objects, and prepare it for analysis
vcf <- read.vcfR( "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/original_vcf/bi_20missing_filtSNP_maf_005.recode.vcf", verbose = FALSE ) # no missing values

## sample information
sample_info <- read.csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Ctarsalis_sample_w_GPS_climate_average_new_filtered_id_region.csv", row.names = 1, stringsAsFactors = TRUE)

## grab loci names info
genlight <- vcfR2genlight(vcf)
SNP_names <- genlight@loc.names

# Merge the names of the outliers with the results from the bayescan dataframe
bayescan=cbind(SNP_names, bayescan) 

## Rename columns
colnames(bayescan)=c("SNP","PROB","LOG_PO","Q_VALUE","ALPHA","FST") 

## Write the results
write.table(bayescan, "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/BayeScan/output_result/culex-bayescan-results.txt", quote=FALSE, sep="\t", row.names=FALSE) 

## Change the value of the Q_VALUE column: 0 == 0.0001
attach(bayescan)
class(bayescan$Q_VALUE)

bayescan$Q_VALUE <- as.numeric(bayescan$Q_VALUE) 
bayescan[bayescan$Q_VALUE<=0.0001,"Q_VALUE"]=0.0001 

## Round the values
bayescan$LOG_PO <- (round(bayescan$LOG_PO, 4)) 
bayescan$Q_VALUE <- (round(bayescan$Q_VALUE, 4)) 
bayescan$ALPHA <- (round(bayescan$ALPHA, 4)) 
bayescan$FST <- (round(bayescan$FST, 6))

## Add a column for the type of selection grouping based on a Q-VALUE < 0.05. You can also choose a Q-VALUE < 0.01 if you want to be more conservative.
bayescan$SELECTION <- ifelse(bayescan$ALPHA>=0&bayescan$Q_VALUE<=0.05,"diversifying",ifelse(bayescan$ALPHA>=0&bayescan$Q_VALUE>0.05,"neutral","balancing"))
# bayescan$SELECTION <- ifelse(bayescan$ALPHA>=0&bayescan$Q_VALUE<=0.01,"diversifying",ifelse(bayescan$ALPHA>=0&bayescan$Q_VALUE>0.01,"neutral","balancing")) 
bayescan$SELECTION<- factor(bayescan$SELECTION)
levels(bayescan$SELECTION) 

## Save the results of the SNPs potentially under positive (divergent) and balancing selection (qvalue < 0.05).
positive <- bayescan[bayescan$SELECTION=="diversifying",] 
neutral <- bayescan[bayescan$SELECTION=="neutral",] 
balancing <- bayescan[bayescan$SELECTION=="balancing",]

## comnbine neutral and balancing 
neutral_balancing <- rbind(neutral, balancing)

## output the neutral_balancing to csv as non candidates snps
write.csv(neutral_balancing, "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/BayeScan/output_result/non_candidate_list_bayescan.csv", row.names = FALSE)

## Check the number of SNPs belonging to each category.
xtabs(data=bayescan, ~SELECTION) 

## Write the results of the SNPs potentially under selection (qvalue < 0.05 or 0.01).
write.table(neutral, "neutral.txt", row.names=F, quote=F)
write.table(balancing, "balancing.txt", row.names=F, quote=F) 
write.table(positive, "positive.txt", row.names=F, quote=F) 

## Transformation Log of the Q value in order to create the ggplot graph.
range(bayescan$Q_VALUE) 
bayescan$LOG10_Q <- -log10(bayescan$Q_VALUE) 

## qqplot
x_title="Log(q-value)" 
y_title="Fst" 

ggplot(bayescan,aes(x=LOG10_Q,y=FST)) +
  geom_point(aes(fill=SELECTION), pch=21, size=2)+ 
  scale_fill_manual(name="Selection",values=c("white","red","orange"))+ 
  labs(x=x_title)+ 
  labs(y=y_title)+   
  theme_classic()

## Save the file in a pdf format
ggsave("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/outlier_analysis/bayescan/bayescan_seacumcumber.png", dpi=300, width=10, height=5) 
dev.off()

##  or plot can use the function plot_R.r already available in BayeScan
source("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/BayeScan/BayeScan2.1/R functions/plot_R.r")
plot_bayescan("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/BayeScan/output_result/dat_fst.txt")



