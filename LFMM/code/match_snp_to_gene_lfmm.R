library(dplyr) # manipulate tables 

rm(list = ls())

# read table from GFF file
gff_table <- read.delim("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/original_vcf/Culex-tarsalis-v1.0.a1-merged-2019-08-30-4-45-01.gff3", header=F, comment.char="#")

# read candidate locus information
candidates_pc1 <- read.csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/LFMM_LatentFactorMixedModels/lfmm_PC1_adaptation_candidates.csv", sep = ",")
candidates_pc1 <- candidates_pc1[c("seq_name", "location")]

candidates_pc2 <- read.csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/LFMM_LatentFactorMixedModels/lfmm_PC2_adaptation_candidates.csv", sep = ",")
candidates_pc2 <- candidates_pc2[c("seq_name", "location")]

candidates_pc3 <- read.csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/LFMM_LatentFactorMixedModels/lfmm_PC3_adaptation_candidates.csv", sep = ",")
candidates_pc3 <- candidates_pc3[c("seq_name", "location")]


# based on the value in gff_table, get subset of gff_table where seq_name == v1, location between the begin and end position
filtered_table_pc1 <- gff_table %>%
  inner_join(candidates_pc1, by = c("V1" = "seq_name")) %>%
  filter(location >= V4, location <= V5) %>%
  select(seq_name = V1, location, source = V2, type = V3, start = V4, end = V5,  score = V6, strand = V7, phase = V8, attributes = V9) # add the candidates location reference

filtered_table_pc2 <- gff_table %>%
  inner_join(candidates_pc2, by = c("V1" = "seq_name")) %>%
  filter(location >= V4, location <= V5) %>%
  select(seq_name = V1, location, source = V2, type = V3, start = V4, end = V5,  score = V6, strand = V7, phase = V8, attributes = V9) # add the candidates location reference

filtered_table_pc3 <- gff_table %>%
  inner_join(candidates_pc3, by = c("V1" = "seq_name")) %>%
  filter(location >= V4, location <= V5) %>%
  select(seq_name = V1, location, source = V2, type = V3, start = V4, end = V5,  score = V6, strand = V7, phase = V8, attributes = V9) # add the candidates location reference


# output
write.csv(filtered_table_pc1, file="/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/LFMM_LatentFactorMixedModels/lfmm_pc1_matched_candidate_to_gene.csv")
write.csv(filtered_table_pc2, file="/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/LFMM_LatentFactorMixedModels/lfmm_pc2_matched_candidate_to_gene.csv")
write.csv(filtered_table_pc3, file="/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/LFMM_LatentFactorMixedModels/lfmm_pc3_matched_candidate_to_gene.csv")
