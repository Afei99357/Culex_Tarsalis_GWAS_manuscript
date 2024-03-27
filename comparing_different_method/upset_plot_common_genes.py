import pandas as pd
from upsetplot import plot, UpSet
from matplotlib import pyplot as plt

lfmm_pc1_df = pd.read_csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/"
                          "Landscape_genetics_GEA/LFMM_LatentFactorMixedModels/lfmm_pc1_gene_match_final.csv",
                          header=0, sep=",", index_col=False)

lfmm_pc2_df = pd.read_csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/"
                          "Landscape_genetics_GEA/LFMM_LatentFactorMixedModels/lfmm_pc2_gene_match_final.csv",
                          header=0, sep=",", index_col=False)

lfmm_pc3_df = pd.read_csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/"
                          "Landscape_genetics_GEA/LFMM_LatentFactorMixedModels/lfmm_pc3_gene_match_final.csv",
                          header=0, sep=",", index_col=False)

rda_df = pd.read_csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/"
                     "Landscape_genetics_GEA/RDA_Redundancy_Analysis/rda_gene_match_final.csv",
                     header=0, sep=",", index_col=False)

bayescan_df = pd.read_csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/outlier_analysis/"
                          "bayescan/bayescan_gene_match_final.csv", header=0, sep=",", index_col=False)

pcadapt_df = pd.read_csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/outlier_analysis/"
                         "PCA_RESULTS/pcadapt_gene_match_final.csv", header=0, sep=",", index_col=False)

## get the unique name from each df in the Name column and store in a list
lfmm_pc1_gene_name = lfmm_pc1_df["Name"].unique().tolist()
# lfmm_pc2_gene_name = lfmm_pc2_df["Name"].unique().tolist()
# lfmm_pc3_gene_name = lfmm_pc3_df["Name"].unique().tolist()
rda_gene_name = rda_df["Name"].unique().tolist()
bayescan_gene_name = bayescan_df["Name"].unique().tolist()
pcadapt_gene_name = pcadapt_df["Name"].unique().tolist()

# Create a DataFrame from the lists for plotting upset plot
# all_genes = set(lfmm_pc1_gene_name + lfmm_pc2_gene_name + lfmm_pc3_gene_name + rda_gene_name + bayescan_gene_name + pcadapt_gene_name)
all_genes = set(lfmm_pc1_gene_name + rda_gene_name + bayescan_gene_name + pcadapt_gene_name)
# all_genes = set(lfmm_pc2_gene_name + rda_gene_name + bayescan_gene_name + pcadapt_gene_name)
# all_genes = set(lfmm_pc3_gene_name + rda_gene_name + bayescan_gene_name + pcadapt_gene_name)

# Convert set to list for DataFrame index
all_genes_list = list(all_genes)

data = pd.DataFrame(index=all_genes_list, columns=['LFMM_PC1', 'RDA', 'Bayescan', 'PCA'])
# data = pd.DataFrame(index=all_genes_list, columns=['LFMM_PC2', 'RDA', 'Bayescan', 'PCA'])
# data = pd.DataFrame(index=all_genes_list, columns=['LFMM_PC3', 'RDA', 'Bayescan', 'PCA'])
# data = pd.DataFrame(index=all_genes_list, columns=['LFMM_PC1', 'LFMM_PC2', 'LFMM_PC3', 'RDA', 'Bayescan', 'PCA'])

# Fill the DataFrame: 1 if gene is in the list, 0 otherwise
data['LFMM_PC1'] = data.index.isin(lfmm_pc1_gene_name).astype(int)
# data['LFMM_PC2'] = data.index.isin(lfmm_pc2_gene_name).astype(int)
# data['LFMM_PC3'] = data.index.isin(lfmm_pc3_gene_name).astype(int)
data['RDA'] = data.index.isin(rda_gene_name).astype(int)
data['Bayescan'] = data.index.isin(bayescan_gene_name).astype(int)
data['PCA'] = data.index.isin(pcadapt_gene_name).astype(int)

# Convert DataFrame format for UpSetPlot
upset_data = data.groupby(list(data.columns)).size()

# Create the UpSet plot
plot(upset_data, show_counts=True, sort_by='cardinality', show_percentages=True, element_size=50)

## save the plot to 300 dpi png file
plt.savefig("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/"
            "compare_landscape_and_outlier_analysis/upset_plot_lfmm_pc1_rda_bayescan_pca.png", dpi=300)

# plt.savefig("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/"
#             "compare_landscape_and_outlier_analysis/upset_plot_lfmm_pc2_rda_bayescan_pca.png", dpi=300)
#
# plt.savefig("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/"
#             "compare_landscape_and_outlier_analysis/upset_plot_lfmm_pc3_rda_bayescan_pca.png", dpi=300)

# plt.savefig("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/"
#             "compare_landscape_and_outlier_analysis/upset_plot_lfmm_pc1_pc2_pc3_rda_bayescan_pca.png", dpi=300)

# get the total common genes between all the methods
common_genes = set(lfmm_pc1_gene_name).intersection(rda_gene_name).intersection(
    bayescan_gene_name).intersection(pcadapt_gene_name)

# common_genes = set(lfmm_pc2_gene_name).intersection(rda_gene_name).intersection(
#     bayescan_gene_name).intersection(pcadapt_gene_name)
#
# common_genes = set(lfmm_pc3_gene_name).intersection(rda_gene_name).intersection(
#     bayescan_gene_name).intersection(pcadapt_gene_name)
#
# common_genes = set(lfmm_pc1_gene_name).intersection(lfmm_pc2_gene_name).intersection(lfmm_pc3_gene_name).intersection(
#     rda_gene_name).intersection(bayescan_gene_name).intersection(pcadapt_gene_name)

## get the gene information any one of the dfs
common_gene_info = bayescan_df[bayescan_df["Name"].isin(common_genes)]

# # output the common gene information to a csv file
common_gene_info.to_csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/"
                        "compare_landscape_and_outlier_analysis/common_genes_lfmm_pc1_rda_bayescan_pcadapt.csv", index=False)
#
# common_gene_info.to_csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/"
#                         "compare_landscape_and_outlier_analysis/common_genes_lfmm_pc2_rda_bayescan_pcadapt.csv", index=False)
#
# common_gene_info.to_csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/"
#                         "compare_landscape_and_outlier_analysis/common_genes_lfmm_pc3_rda_bayescan_pcadapt.csv", index=False)
#
# common_gene_info.to_csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/"
#                         "compare_landscape_and_outlier_analysis/common_genes_lfmm_pc1_pc2_pc3_rda_bayescan_pcadapt.csv", index=False)
#
