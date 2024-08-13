import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Load the data
file_path = "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/rda_analysis_candidates_correlation_to_env.csv"  # Update the path to your CSV file
data = pd.read_csv(file_path)

# Set the 'snp' column as the index
data.set_index('snp', inplace=True)

## change axis column values, if 1 to "Group 1", if 2 to "Group 2", if 3 to "Group 3", if 4 to "Group 4"
data['axis'] = data['axis'].replace({1: 'RDA1 SNP Group', 2: 'RDA2 SNP Group', 3: 'RDA3 SNP Group', 4: 'RDA4 SNP Group'})

## set the color for each group
group_colors = {
    'RDA1 SNP Group': '#a6611a',
    'RDA2 SNP Group': '#1E88E5',
    'RDA3 SNP Group': '#FFC107',
    'RDA4 SNP Group': '#80cdc1'
}

## map the group_colors to the axis column
data['axis'] = data['axis'].map(group_colors)

## get the axis column as group_info
group_info = data.pop('axis')

# Drop the 'axis' column
data.drop(columns=['loading', 'predictor', 'correlation'], inplace=True)

## transpose the data
data = data.T

## set the
fig, ax = plt.subplots(figsize=(20, 7))

sns.heatmap(data, cmap='coolwarm', annot=False)


## remove the x axis ticks labels
ax.set_xticks([])

## adding curly brace to the x axis for each group in the group_info as the x axis tick labels
## get the unique group_info
unique_group_info = group_info.unique()

## get the number of unique group_info
num_unique_group_info = len(unique_group_info)

## get the number of snps in each unique group
group_size_list = [len(group_info[group_info == group]) for group in unique_group_info]

group_labels = ["RDA1 SNP Group", "RDA2 SNP Group", "RDA3 SNP Group", "RDA4 SNP Group"]

## Position for each curly brace and labels
first_position = group_size_list[0]/2
second_position = group_size_list[0] + group_size_list[1]/2
third_position = group_size_list[0] + group_size_list[1] + group_size_list[2]/2
fourth_position = group_size_list[0] + group_size_list[1] + group_size_list[2] + group_size_list[3]/2

## store the x tick location
x_tick_locations = [0, group_size_list[0], sum(group_size_list[:2]), sum(group_size_list[:3]), sum(group_size_list[:4])]

## adding the x axis tick at the x_tick_locations
plt.xticks(x_tick_locations)

## change the x axis tick length
ax.tick_params(axis='x', length=10)

## y position for the curly brace using the min y - 0.5
y_brace = -0.01

## get the minimun value of y axis in the plot
y_min = ax.get_ylim()[0]

# adding text for each group
for i, group in enumerate(unique_group_info):
    if i == 0:
        start_index = 0
        end_index = group_size_list[i]
        ax.text(x=(start_index + end_index)/2, y=y_min + 0.3, s=group_labels[i], fontsize=15, ha='center', color=group)

    if i == 1:
        start_index = group_size_list[0]
        end_index = group_size_list[i] + group_size_list[0]
        ax.text(x=(start_index + end_index) / 2, y=y_min + 0.3, s=group_labels[i], fontsize=15, ha='center', color=group)

    elif i == 2 or i == 3:
        start_index = sum(group_size_list[:i])
        end_index = sum(group_size_list[:i+1])
        ax.text(x=(start_index + end_index) / 2, y=y_min + 0.3, s=group_labels[i], fontsize=15, ha='center', color=group)



# ## set the y axis tick labels size
plt.yticks(fontsize=12)

## remove x axis label
plt.xlabel('')


## adding a custom legend with the group colors
legend_elements = [Line2D([0], [0], marker='o', color=group_colors[group], label=group, markersize=10) for group in group_colors]
plt.legend(handles=legend_elements, title="SNP Group", bbox_to_anchor=(1.11, 1), loc='upper left', fontsize=12)


plt.title('SNPs and Variables Correlation Heatmap', fontsize=20)
plt.tight_layout()
plt.savefig('/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/rda_analysis_candidates_correlation_to_env_heatmap.png', dpi=300)
