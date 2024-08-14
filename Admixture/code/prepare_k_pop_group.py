import pandas as pd

## READ TXT FILE WITH space and there is no collumn names
df_group = pd.read_csv(
    "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/admixture/"
    "plink_file_for_admixture/culex_plink_new.4.Q",
    sep=" ",
    header=None,
)


## read file with mosquito ID and population
df_pop = pd.read_csv(
    "/Users/ericliao/Desktop/WNV_project_files/"
    "Ctarsalis_sample_w_GPS_climate_average_new_filtered_id.csv",
    sep=",",
    header=0,
    index_col=0,
)

## get the vcfID list from the df_pop and add it to the df_group
vcfID_list = df_pop["vcfID"].tolist()
df_group["vcfID"] = vcfID_list

## add column names to the df_group which are group1, group2, group3, group4
df_group.columns = ["group1", "group2", "group3", "group4", "vcfID"]

## add a new empty column called region to the df_pop
df_pop["region"] = ""

## for each row in df_group, find the max value from group1, group2, group3, group4 columns
# and add the corresponding region to the df_pop region column based on the vcfID
for index, row in df_group.iterrows():
    max_value = max(row["group1"], row["group2"], row["group3"], row["group4"])
    if max_value == row["group1"]:
        df_pop.loc[df_pop["vcfID"] == row["vcfID"], "region"] = "Northwest"
    elif max_value == row["group2"]:
        df_pop.loc[df_pop["vcfID"] == row["vcfID"], "region"] = "Midwest"
    elif max_value == row["group3"]:
        df_pop.loc[df_pop["vcfID"] == row["vcfID"], "region"] = "West Coast"
    elif max_value == row["group4"]:
        df_pop.loc[df_pop["vcfID"] == row["vcfID"], "region"] = "Southwest"

## save the df_pop as a csv file
df_pop.to_csv(
    "/Users/ericliao/Desktop/WNV_project_files/Ctarsalis_sample_w_GPS_climate_average_new_filtered_id_region.csv"
)
