
rm(list = ls())

### Read in the pairwise Fst Calculations
fst_matrix = read.csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/IBD_IBE_test/pairwise_fsts_matrix_between_population.csv", header=TRUE, row.names = 1)

# Get the dimensions of the matrix
n <- nrow(fst_matrix)

# Create an empty data frame to store the results
result_df <- data.frame(
  Pop1 = character(0),
  Pop2 = character(0),
  PairwiseFst = numeric(0)
)

# Loop through the upper triangle of the matrix to extract pairs and distances
for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    result_df <- rbind(result_df, data.frame(
      Pop1 = rownames(fst_matrix)[i],
      Pop2 = colnames(fst_matrix)[j],
      PairwiseFst = fst_matrix[i, j]
    ))
  }
}

# Save the result_df data frame as a CSV file
write.csv(result_df, file = "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/IBD_IBE_test/pairwise_fst_table.csv", row.names = FALSE)
