import pandas as pd
import matplotlib.pyplot as plt

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
    "Ctarsalis_sample_w_GPS_climate_average_new_filtered_id_region.csv",
    sep=",",
    header=0,
    index_col=0,
)

## get the State, City, region, vcfID, GPS.Lat, GPS.Lon from the df_pop and add it to the df_group
state_list = df_pop["State"].tolist()
city_list = df_pop["City"].tolist()
region_list = df_pop["region"].tolist()
vcfID_list = df_pop["vcfID"].tolist()
GPS_Lat_list = df_pop["GPS.Lat"].tolist()
GPS_Lon_list = df_pop["GPS.Lon"].tolist()
df_group["State"] = state_list
df_group["City"] = city_list
df_group["region"] = region_list
df_group["vcfID"] = vcfID_list
df_group["GPS.Lat"] = GPS_Lat_list
df_group["GPS.Lon"] = GPS_Lon_list

## add column names to the df_group which are group1, group2, group3, group4
df_group.columns = [
    "Northwest",
    "Midwest",
    "West Coast",
    "Southwest",
    "State",
    "City",
    "region",
    "vcfID",
    "GPS.Lat",
    "GPS.Lon",
]

## plot one stack bar plot with axis is each vcfID and the bar is the values of group1, group2, group3, group4, sorting the
## plot by GPS.lon first and then GPS.Lat, color the bar by the region, where northwest is skyblue,
## midwest is hotpink, west coast is goldenrod, southwest is forestgreen

## sort the df_group by region and GPS.Lon in descending way
df_group = df_group.sort_values(by=["region", "GPS.Lon"], ascending=False)

## plot the bar plot with x axis is vcfID descending order
df_group.plot.bar(
    x="vcfID",
    y=["Northwest", "Midwest", "West Coast", "Southwest"],
    stacked=True,
    figsize=(20, 5),
    color=["skyblue", "hotpink", "goldenrod", "forestgreen"],
)

## x axis label use the City name and state name, if the neibouring x axis label is the same, then only show once
## y axis label is the percentage
plt.xticks(
    range(len(df_group["vcfID"])),
    df_group["City"] + ", " + df_group["State"],
    rotation=90,
)

ax = plt.gca()
# if the current x axis label is the same as next, dont show the current label
preview_label = ""
for label in ax.xaxis.get_ticklabels():
    # if the previous label is not empty and the current label is the same as the previous label, then dont show the current label
    if preview_label != "" and preview_label == label.get_text():
        label.set_visible(False)
    preview_label = label.get_text()

## legend on the left upper corner
plt.legend(loc="upper left")

# plt.title("Admixture plot for Culex tarsalis")
plt.tight_layout()

## remove the x axis and y axis label
plt.xlabel("")
plt.ylabel("")

# plt.show()
plt.savefig(
    "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Admixture/stack_bar_plot_admixture_K_4.png",
    dpi=300,
)
