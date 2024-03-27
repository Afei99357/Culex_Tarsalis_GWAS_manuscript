library(VennDiagram)
library(dplyr)

rm(list = ls())

# load lfmm pc1 candidates
lfmm_pc1_results <- read.csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/LFMM_LatentFactorMixedModels/lfmm_pc1_gene_match_final.csv")

# load lfmm RDA candidates
rda_results <- read.csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/RDA_Redundancy_Analysis/rda_gene_match_final.csv")

# load bayescan candidates
bayescan_results <- read.csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/outlier_analysis/bayescan_gene_match_final.csv")

# load pcadapt candidates
pcadapt_results <- read.csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/outlier_analysis/pcadapt_gene_match_final.csv")

lfmm_pc1_unique <- unique(lfmm_pc1_results$Name)
rda_unique <- unique(rda_results$Name)
bayescan_unique <- unique(bayescan_results$Name)
pcadapt_unique <- unique(pcadapt_results$Name)

## lfmm pc1 rda common venn diagram
venn.diagram(
  x = list(lfmm_pc1_unique,rda_unique, bayescan_unique, pcadapt_unique),
  category.names = c("LFMM_PC1" , "RDA", "Bayescan_Outlier", "Pcadapt"),
  filename = '/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/common_landscape_and_outlier_analysis/Venn_diagramm_overlapping_landscape_outlier_analysis.png',
  output=TRUE,
  imagetype="png",
  resolution = 300,
  lty = 1,  # Set line type (solid line)
  lwd = 2,  # Set line width
  alpha = 0.5,  # Adjust transparency
  cex = 4,
  cat.cex = 3,  # Increase the size of category labels
  cat.pos = 0,
  cat.col = c("black", "blue", "purple", "green")
)


# check common
common_rows <- Reduce(intersect, list(lfmm_pc1_results, rda_results, bayescan_results, pcadapt_results))

## ouput the common results
write.csv(common_rows, "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/common_landscape_and_outlier_analysis/landscape_outlier_analysis_common_gene_match_final.csv", row.names = FALSE)

