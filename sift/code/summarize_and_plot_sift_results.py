import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

## read the sift results
sift_result_df = pd.read_csv(
    "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/sift_prediciton_results/overlapping_gene_snps_sift_info_combined.csv",
    header=0,
    index_col=False,
    sep=",",
)

## create a new empty table with column names are: Gene_ID, Total_SNPs, Nonsynonymous, nsyn_tolerated, nsyn_deleterious,
# Synonymous, syn_tolerated, syn_deleterious
sift_summary_df = pd.DataFrame(columns=["gene_ID", "total_SNPs", "nonsynonymous", "nsyn_tolerated", "nsyn_deleterious",
                                        "synonymous", "syn_tolerated", "syn_deleterious"])

## get the unique gene_ID from sift_result_df
gene_ID_list = sift_result_df["GENE_ID"].unique().tolist()

## for each gene_ID, get the unique rows from sift_result_df based on CHROM, POSITION, VARIANT_TYPE AND PREDICTION
for gene_ID in gene_ID_list:
    gene_sift_result_df = sift_result_df[sift_result_df["GENE_ID"] == gene_ID]
    ## DROP DUPLICATES BASED ON CHROM, POSITION, VARIANT_TYPE AND PREDICTION
    gene_sift_result_df = gene_sift_result_df.drop_duplicates(subset=["CHROM", "POSITION", "VARIANT_TYPE", "PREDICTION"])
    ## get the total number of SNPs
    total_SNPs = gene_sift_result_df.shape[0]
    ## get the total number of NONSYNONYMOUS and SYNONYMOUS SNPs
    nonsynonymous = gene_sift_result_df[gene_sift_result_df["VARIANT_TYPE"] == "NONSYNONYMOUS"].shape[0]
    synonymous = gene_sift_result_df[gene_sift_result_df["VARIANT_TYPE"] == "SYNONYMOUS"].shape[0]
    ## get the total number of TOLERATED and DELETERIOUS SNPs FOR NONSYNONYMOUS AND SYNONYMOUS
    nonsyn_tolerated = gene_sift_result_df[(gene_sift_result_df["VARIANT_TYPE"] == "NONSYNONYMOUS") &
                                           (gene_sift_result_df["PREDICTION"] == "TOLERATED")].shape[0]
    nonsyn_deleterious = gene_sift_result_df[(gene_sift_result_df["VARIANT_TYPE"] == "NONSYNONYMOUS") &
                                             (gene_sift_result_df["PREDICTION"] == "DELETERIOUS")].shape[0]
    syn_tolerated = gene_sift_result_df[(gene_sift_result_df["VARIANT_TYPE"] == "SYNONYMOUS") &
                                        (gene_sift_result_df["PREDICTION"] == "TOLERATED")].shape[0]
    syn_deleterious = gene_sift_result_df[(gene_sift_result_df["VARIANT_TYPE"] == "SYNONYMOUS") &
                                          (gene_sift_result_df["PREDICTION"] == "DELETERIOUS")].shape[0]

    ## append the summary to the sift_summary_df
    sift_summary_df = sift_summary_df.append({"gene_ID": gene_ID, "total_SNPs": total_SNPs,
                                              "nonsynonymous": nonsynonymous, "nsyn_tolerated": nonsyn_tolerated,
                                              "nsyn_deleterious": nonsyn_deleterious, "synonymous": synonymous,
                                              "syn_tolerated": syn_tolerated, "syn_deleterious": syn_deleterious},
                                             ignore_index=True)

## output the sift_summary_df to a csv file
sift_summary_df.to_csv(
    "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/sift_prediciton_results/overlapping_sift_result_summary.csv",
    header=True,
    index=False,
    sep=",",
)

# Create a figure and a 3x1 grid of subplots
fig, axs = plt.subplots(3, 1, figsize=(6, 12))

## for each gene_ID, get the unique rows from sift_result_df based on CHROM, POSITION, VARIANT_TYPE AND PREDICTION
axs[0].bar(sift_summary_df["gene_ID"], sift_summary_df["nonsynonymous"], color="blue", label="nonsynonymous")
axs[0].bar(sift_summary_df["gene_ID"], sift_summary_df["synonymous"], bottom=sift_summary_df["nonsynonymous"],
         color="red", label="synonymous")
## only show whole numbers on the y axis
axs[0].yaxis.set_major_locator(MaxNLocator(integer=True))
## no grid lines
axs[0].grid(False)
## remove the x axis ticks
axs[0].set_xticks([])
## add legend
legend = ["nonsynonymous", "synonymous"]
axs[0].legend(legend, loc="upper right", title="Variant type")

## for nonsynonymous plot stacked bar plot for the number of tolerated and deleterious SNPs for each gene_ID, in the plot, the bar for tolerated SNPs is on top of the bar for deleterious SNPs
axs[1].bar(sift_summary_df["gene_ID"], sift_summary_df["nsyn_tolerated"], color="#0000e6", label="tolerated")
axs[1].bar(sift_summary_df["gene_ID"], sift_summary_df["nsyn_deleterious"], bottom=sift_summary_df["nsyn_tolerated"],
       color="#9999ff", label="deleterious")

## only show whole numbers on the y axis
axs[1].yaxis.set_major_locator(MaxNLocator(integer=True))
## no grid lines
axs[1].grid(False)
## remove the x axis ticks
axs[1].set_xticks([])
## add legend
legend = ["tolerated", "deleterious"]
## add text under the legend
axs[1].legend(legend, loc="upper right", title="Nonsynonymous")

## for synonymous, plot stacked bar plot for the number of tolerated and deleterious SNPs for each gene_ID, in the plot, the bar for tolerated SNPs is on top of the bar for deleterious SNPs
axs[2].bar(sift_summary_df["gene_ID"], sift_summary_df["syn_tolerated"], color="red", label="tolerated")
axs[2].bar(sift_summary_df["gene_ID"], sift_summary_df["syn_deleterious"], bottom=sift_summary_df["syn_tolerated"],
       color="#800000", label="deleterious")
axs[2].set_xticklabels(axs[2].get_xticklabels(), rotation=90)
## only show whole numbers on the y axis
axs[2].yaxis.set_major_locator(MaxNLocator(integer=True))
## no grid lines
axs[2].grid(False)
## add legend
legend = ["tolerated", "deleterious"]
## add text before the legend
axs[2].legend(legend, loc="upper right", title="Synonymous")

plt.subplots_adjust(hspace=0.05, top=0.99, left=0.1)
# plt.show()

## save the plot to 300 dpi png file
plt.savefig("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/sift_prediciton_results/overlapping_sift_result_summary.png", dpi=300)