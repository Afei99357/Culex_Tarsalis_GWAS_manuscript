import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D  # For custom legend

# Define colors for each gene
gene_colors = {
    'Ct.00g025080': 'red',
    'Ct.00g026900': 'blue',
    'Ct.00g030230': 'green',
    'Ct.00g032480': 'c',
    'Ct.00g049290': 'magenta',
    'Ct.00g051300': 'goldenrod',
    'Ct.00g062900': 'black',
    'Ct.00g064410': 'orange',
    'Ct.00g095350': 'purple',
    'Ct.00g154760': 'brown',
    'Ct.00g176220': 'deeppink',
    'Ct.00g179740': 'gray',
    'Ct.00g237940': 'lime',
    'Ct.00g238000': 'olive',
    'Ct.00g280270': 'chocolate',
    'Ct.00g280280': 'teal',
    'Ct.00g290200': 'coral'
}

# Load data
correlation_data = pd.read_csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/rda_analysis_candidates_correlation_to_env.csv")
snp_gene_data = pd.read_csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/candidate_gene_map_to_RDA_SNPs.csv")

### for 17 genes snps, create a df to store only the genes are in the gene_colors' keys
gene_17_data = snp_gene_data[snp_gene_data['Gene'].isin(gene_colors.keys())]

## get axis and snp columns from the correlation_data
axis_snp = correlation_data[['axis', 'snp']]

## merge the axis_snp with snp_gene_data on snp
gene_17_data = pd.merge(gene_17_data, axis_snp, on='snp', how='left')

## order the snp_gene_data by axis
gene_17_data = gene_17_data.sort_values(by='axis')

## change the axis column values if 1 to "RDA1 SNP Group", 2 to "RDA2 SNP Group", 3 to "RDA3 SNP Group", 4 to "RDA4 SNP Group"
gene_17_data['axis'] = gene_17_data['axis'].replace({1: 'RDA1 SNP Group', 2: 'RDA2 SNP Group', 3: 'RDA3 SNP Group', 4: 'RDA4 SNP Group'})

## reset the index
gene_17_data.reset_index(drop=True, inplace=True)

## output the snp_gene_data to a csv file
gene_17_data.to_csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/17_gene_map_to_RDA_SNPs.csv", index=False)

# Drop the 'axis' column
correlation_data.drop(columns=['axis', 'loading', 'predictor', 'correlation'], inplace=True)

# Merge the datasets on SNP
full_data = pd.merge(correlation_data, snp_gene_data, on='snp', how='left')

# Set 'snp' as the index after ensuring it is not dropped in merging process
data = full_data.set_index('snp')

# Transpose data to get SNPs on the x-axis and variables on the y-axis
data_transposed = data.drop('Gene', axis=1).T  # Drop the 'Gene' column from the data used in the heatmap

# Create the heatmap
plt.figure(figsize=(200, 20))
ax = sns.heatmap(data_transposed, cmap='coolwarm', annot=False, cbar_kws={'shrink': 1, 'pad': 0.01, 'aspect': 30})

# Customize x-axis labels with gene colors
# This assumes that 'Gene' information is still in 'full_data' dataframe indexed by 'snp'
for ticklabel, (snp, row) in zip(ax.get_xticklabels(), full_data.iterrows()):
    ticklabel.set_color(gene_colors.get(row['Gene'], 'black'))  # Default to black if no gene is matched

# Rotate tick labels for better visibility
plt.xticks(rotation=90, fontsize=10)
plt.yticks(rotation=0, fontsize=20)

# Create a legend for the colors
legend_elements = [Line2D([0], [0], marker='o', color=gene_colors[gene], label=gene, markersize=10) for gene in gene_colors]
plt.legend(handles=legend_elements, title="Gene Colors", bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=20)

# Title and save
plt.title('Heatmap of the Correlation between SNPs and Environmental Variables with Gene Color Coded', fontsize=20)
plt.tight_layout()
plt.savefig('/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/rda_analysis_candidates_correlation_to_env_heatmap_annotated_gene.png', dpi=100)
