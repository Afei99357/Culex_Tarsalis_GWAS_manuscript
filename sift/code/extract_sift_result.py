import pandas as pd
import os

# read the candidate gene file
candidate_gene_df = pd.read_csv(
    "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/rda_gene_match_final.csv",
    header=0,
    index_col=None,
    sep=",",
)

# get the unique gene list from the Name column, remove 20471- and .m01 from the gene name
gene_list = []
for gene_entry in candidate_gene_df["Name"]:
    gene = gene_entry.split("-")[1].split(".m")[0]
    gene_list.append(gene)

# get the unique gene list
gene_list = list(set(gene_list))

## read the candidate snp match gene file
candidate_snp_df = pd.read_csv(
    "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/rda_matched_candidate_to_gene.csv",
    header=0,
    index_col=0,
    sep=",",
)

## for each gene in the gene_list, find all the rows that column attributes contains the gene name
## and save the rows to a new dataframe

# create a new dataframe
gene_df = pd.DataFrame()

# for each gene in the gene_list, find all the rows that column attributes contains the gene name,
# then add to the new dataframe
for gene in gene_list:
    gene_df = gene_df.append(
        candidate_snp_df[candidate_snp_df["attributes"].str.contains(gene)]
    )


## get the unique pair of seq_name and location from the gene_df
seq_name_location_list = []
for index, row in gene_df.iterrows():
    seq_name_location = row["seq_name"] + "_" + str(row["location"])
    seq_name_location_list.append(seq_name_location)

seq_name_location_list = list(set(seq_name_location_list))

# store the seq_name_location_list to a new dataframe and have two columns: seq_name and location
seq_name_location_df = pd.DataFrame()
seq_name_location_df["seq_name"] = [
    seq_name_location.split("_")[0] for seq_name_location in seq_name_location_list
]
seq_name_location_df["location"] = [
    seq_name_location.split("_")[1] for seq_name_location in seq_name_location_list
]

## create a new dataframe to store the sift result
sift_result_df = pd.DataFrame()

# read all the sift tsv files names in the directory
tsv_file_list = []
for file in os.listdir(
    "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/sift_prediciton_results/sift_prediction_results"
):
    if file.endswith(".tsv"):
        tsv_file_list.append(file)

# for each row in seq_name_location_df, find the tsv file which file name contains the seq_name,
# then find rows that second column value is the same as location, save the rows to a new dataframe
for index, row in seq_name_location_df.iterrows():
    for tsv_file in tsv_file_list:
        if row["seq_name"] in tsv_file:
            tsv_df = pd.read_csv(
                "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/sift_prediciton_results/sift_prediction_results/"
                + tsv_file,
                sep="\t",
                header=None,
                index_col=None,
            )
            sift_result_df = sift_result_df.append(
                tsv_df[tsv_df[1] == int(row["location"])]
            )

## column names:
sift_result_df.columns = [
    "CHROM",
    "POSITION",
    "REF_ALLELE",
    "ALT_ALLELE",
    "TRANSCRIPT_ID",
    "GENE_ID",
    "GENE_NAME",
    "REGION",
    "VARIANT_TYPE",
    "REF_AA",
    "ALT_AA",
    "AA_POS",
    "SIFT_SCORE",
    "SIFT_MEDIAN",
    "NUM_SEQs",
    "dbSNP",
    "PREDICTION",
]

## output the results to a csv file
sift_result_df.to_csv(
    "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/sift_prediciton_results/rda_candidates_sift_prediction.csv",
    index=False,
    header=True,
    sep=",",
)
