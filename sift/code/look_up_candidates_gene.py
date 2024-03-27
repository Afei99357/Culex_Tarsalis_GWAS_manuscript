import pandas as pd

# read tab delimited file
df_interproscan = pd.read_csv(
    "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/original_vcf/"
    "Culex-tarsalis-v1.0.a1.5d6405151b078-interproscan.tab",
    header=8,
    index_col=False,
    sep="\t",
)

# read candidate gene list
candidate_gene_df = pd.read_csv(
    "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/rda_matched_candidate_to_gene.csv",
    header=0,
    index_col=0,
    sep=",",
)

# get the subset of candidate genes df by type is "gene"
candidate_gene_df_gene = candidate_gene_df[candidate_gene_df["type"] == "gene"]

# ## get gene list and build gene_list value by append 20471- in front of each gene name and .m01 at the end
gene_list = []
for gene_entry in candidate_gene_df_gene["attributes"]:
    key_value_pairs = gene_entry.split(";")[0]
    gene = key_value_pairs.split("=")[1]
    gene_list.append("20471-" + gene + ".m01")

# from df_interproscan, the Name column contains the characters that we want to match in the gene_list. get the rows that contains gene_list values as part of the valu
df_interproscan_gene = df_interproscan[df_interproscan["Name"].isin(gene_list)]

# save the df_interproscan_gene as a csv file
df_interproscan_gene.to_csv(
    "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/rda_gene_match_final.csv",
    index=False,
)
