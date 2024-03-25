library(tidyverse)
library(ggplot2)

rm(list = ls())

cv_errors <- read.csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/admixture/plink_file_for_admixture/admixture_CV_error_list.csv")

# check for which K is the best based on Cross validation error
ggplot(data=cv_errors, aes(x=K, y=CV_error, group=1)) +
  geom_line(linetype = "dashed")+
  geom_point()+
  ylab("cross-validation error")+
  xlab("K")+
  scale_x_continuous(breaks = c(1:13))+
  theme_classic()

sample_list <- read.csv("/Users/ericliao/Desktop/WNV_project_files/sample_list.csv", header = FALSE)

pop_list <- as.data.frame(sample_list$V1)
sample_name_list <- as.data.frame(sample_list$V2)
latitude_list <- as.data.frame(sample_list$V3)
longitude_list <- as.data.frame(sample_list$V4)

colnames(pop_list) <- "popID"
colnames(sample_name_list) <- "sampleID"
colnames(latitude_list) <- "latitude"
colnames(longitude_list) <- "longitude"

read_delim("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/admixture/plink_file_for_admixture/culex_plink_new.4.Q",
           col_names = paste0("Q", seq(1:2)),
           delim = " ")

all_data <- tibble(sampleID = character(),
                   k = numeric(),
                   Q = character(),
                   value = numeric(),
                   latitude = numeric(),
                   longitude = numeric())

for (k in 1:13) {
  data <- read_delim(paste0("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/admixture/plink_file_for_admixture/culex_plink_new.", k, ".Q"),
                     col_names = paste0("Q", seq(1:k)),
                     delim = " ")
  
  data$popID <- pop_list$popID
  data$sampleID <- sample_name_list$sampleID
  data$latitude <- latitude_list$latitude
  data$longitude <- longitude_list$longitude
  data$k <- k
  
  
  # This step converts from wide to long.
  data %>% gather(Q, value,-popID, -sampleID,  -k, -latitude, -longitude) -> data
  all_data <- rbind(all_data, data)
}

all_data

# plottin!! abandon, plotted in python

# custom_colors <- c("skyblue", "hotpink", "goldenrod", "forestgreen")
# # only for K=4
# all_data %>% 
#   filter(k == 4) %>%
#   ggplot(.,aes(x=sampleID,y=value,fill=factor(Q))) + 
#   geom_bar(stat="identity",position="stack") +
#   xlab("Sample") + ylab("Ancestry") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
#   scale_fill_manual(values=custom_colors,name="K",
#                     labels=c("1", "2", "3", "4"))




