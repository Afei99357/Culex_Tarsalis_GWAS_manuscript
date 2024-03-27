library(dplyr) # manipulate tables 

rm(list = ls())

# read table from GFF file
gff_table <- read.delim("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/original_vcf/Culex-tarsalis-v1.0.a1-merged-2019-08-30-4-45-01.gff3", header=F, comment.char="#")

# read candidate locus information
candidates <- read.csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/rda_analysis_candidates.csv", sep = ",")
candidates <- candidates["snp"]

# Split the "adaptation_candidates" column into "seq_name" and "location"
split_columns <- strsplit(candidates$snp, "_")

# Create new columns "seq_name" and "location"
candidates$seq_name <- sapply(split_columns, function(x) x[1])
candidates$location <- sapply(split_columns, function(x) x[2])

# Remove the original "adaptation_candidates" column if no longer needed
candidates$snp <- NULL

# based on the value in gff_table, get subset of gff_table where seq_name == v1, location between the begin and end position
filtered_table <- gff_table %>%
  inner_join(candidates, by = c("V1" = "seq_name")) %>%
  filter(location >= V4, location <= V5) %>%
  dplyr::select(seq_name = V1, location, source = V2, type = V3, start = V4, end = V5,  score = V6, strand = V7, phase = V8, attributes = V9) # add the candidates location reference

# output
write.csv(filtered_table, file="/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/rda_matched_candidate_to_gene.csv")

