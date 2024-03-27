library(conStruct)
library(ggplot2)
library(scatterpie)
library(dplyr)

rm(list = ls())

admix_proportions <- read.csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Admixture/culex_plink_new.4.Q", sep=' ', header = FALSE)

mat <- as.matrix(admix_proportions)

env <- read.csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Ctarsalis_sample_w_GPS_climate_average_new_filtered_id_region.csv", row.names = 1)

coordinates <- env[c("region","GPS.Lon", "GPS.Lat")]

## c("Northwest", "Midwest", "West Coast", "Southwest") ==> c("skyblue", "hotpink", "goldenrod", "forestgreen")
layer_colors <- c("skyblue", "hotpink", "goldenrod", "forestgreen")

# Get USA map data
us_map <- map_data(map = "state", region = ".")

# Plot the US map
base_map <- ggplot(us_map, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill = "gray30", color = "black") + coord_fixed() 

admix_coords <- cbind(admix_proportions, coordinates)

# Group by longitude and latitude, then summarize to get sample size for each group
admix_coords_summary <- admix_coords %>% group_by(GPS.Lon, GPS.Lat) %>% summarize(SampleCount = n()) %>% ungroup()

result <- merge(admix_coords, admix_coords_summary, by = c("GPS.Lon", "GPS.Lat"), all = FALSE)

colnames(result) <- c("GPS.Lon", "GPS.Lat", "Northwest", "Midwest", "West Coast", "Southwest", "region", "Sample_Count")

## adjust the radius range based on sample count
result$Sample_Count <- sqrt(result$Sample_Count) / 3.5

admixture_pieplot <- base_map + geom_scatterpie(aes(x = GPS.Lon, y = GPS.Lat, r = Sample_Count), 
                                                data = result, 
                                                pie_scale=0.1,
                                                cols = c("Northwest", "Midwest", "West Coast", "Southwest"), 
                                                color=NA) + 
  coord_equal(expand = FALSE) + 
  scale_fill_manual(values=layer_colors) +
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = "Regions") + 
  theme(
        legend.position = "none",  # Adjust the legend position
        # legend.text = element_text(size = 30),  # Adjust the size of the legend labels
        # legend.title = element_text(size = 35),
        legend.title = element_blank(), # Remove legend title
        panel.background = element_blank(),  # Remove panel background
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.text.x = element_blank(),  # Remove x-axis tick labels
        axis.text.y = element_blank(),  # Remove y-axis tick labels
        axis.ticks = element_blank(),
        axis.title.x = element_blank(), #remove x axis label 
        axis.title.y = element_blank(), #remove y axis label
        ) 

admixture_pieplot
ggsave("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Admixture/admixtrue_pie_map_plot.png", 
       plot = admixture_pieplot, width = 20, height = 22, units = "in", dpi = 300)

