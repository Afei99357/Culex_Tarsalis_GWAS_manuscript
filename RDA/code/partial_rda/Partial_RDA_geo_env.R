rm(list = ls())

##Loading library
library(vcfR)     # for vcf file
library(adegenet) # for genlight obeject and tab() function
library(vegan)    # Used to run PCA & RDA
library(qvalue)   # Used to post-process LFMM output
library(psych)
library(ggplot2)

##Load genetic data and convert it to genlight objects, and prepare it for analysis

#vcf <- read.vcfR( "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/original_vcf/bi_20missing_filtSNP_maf_005.recode.vcf", verbose = FALSE ) # no missing values
vcf <- read.vcfR( "/projects/cooper_research2/eric/data_landscape_genetics/original_vcf/bi_20missing_filtSNP_maf_005.recode.vcf", verbose = FALSE ) # no missing values


SNP_genlight <- vcfR2genlight(vcf) # convert to genlight objects

allele_matrix_genlight <- tab(SNP_genlight, freq = FALSE, NA.method = "asis") # get the matrix of allele counts

dim(allele_matrix_genlight) # check the dimension of genetic data

sum(is.na(allele_matrix_genlight)) # check the data dimmension and number of NAs

gen.imp <- apply(allele_matrix_genlight, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x)))))) # impute the NA values into most common allele

sum(is.na(gen.imp)) # check if No NAs

## get the environmental data
##Load the enviromental data and get it ready for analysis
#env <- read.csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Ctarsalis_sample_w_GPS_climate_average_new_filtered_id_region.csv", row.names = 1) # load enviromental data
env <- read.csv("/projects/cooper_research2/eric/data_landscape_genetics/Ctarsalis_sample_w_GPS_climate_average_new_filtered_id_region.csv", row.names = 1) # load enviromental data

env <- env[, !(names(env) %in% c("locID", "State", "City", "date"))] # remove unused columns

geo_vars <- env[, c("GPS.Lat", "GPS.Lon")] # get the geographic variables

str(env) # check the structure of the enviromental data

env$vcfID <- as.character(env$vcfID) # Make individual names characters (not factors)

identical(rownames(gen.imp), env[,1]) # Confirm that genotypes and environmental data are in the same order

pred <- env[, 18:30] # subset enviromental predictors and shorten their names

# ## rename the predictors
new_names <- c("eastward_wind", "northward_wind", "temperature",
               "evaporation_from_bare_soil", "high_vegetation",
               "low_vegetation", "water_retention_capacity",
               "snowfall", "surface_net_solar_radiation", "surface_runoff",
               "evaporation", "total_precipitation", "volumetric_soil_water_layer1")

## renmae the predictors
colnames(pred) <- new_names

# check the Linear dependency
usdm::vif(pred)

#png("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/variable_pair_correlation_panel.png", width = 8, height = 6, units = "in", res = 300)
png("variable_pair_correlation_panel.png", width = 8, height = 6, units = "in", res = 300)

pairs.panels(pred, scale = TRUE, 
             cex = 1.5)       # Controls overall text size

dev.off()

## remove as little variables as we can based on any pair if the correlation is higher than 0.7
pred <- subset(pred, select=-c(surface_net_solar_radiation, evaporation_from_bare_soil, total_precipitation, volumetric_soil_water_layer1, snowfall))

# check the Linear dependency
usdm::vif(pred)

# Normalize the data and convert to data frames
pred <- as.data.frame(scale(pred))
geo_vars <- as.data.frame(scale(geo_vars))

### adding geographical variables to the data
pred_with_geo <- cbind(pred, geo_vars)

## Partial RDA for Environment Controlling for Geography
partial_rda_env_geo <- rda(gen.imp ~ . + Condition(GPS.Lat + GPS.Lon), data = pred_with_geo, scale = TRUE)

# Summarize results
summary(partial_rda_env_geo)

# Check adjusted R^2
RsquareAdj(partial_rda_env_geo)

# Test significance
paste("test significance -- Partial RDA for Environment Controlling for Geography")
anova.cca(partial_rda_env_geo, permutations = 999)


### Partial RDA for Geography Controlling for Environment
partial_rda_geo_env <- rda(gen.imp ~ GPS.Lat + GPS.Lon + Condition(eastward_wind + northward_wind + temperature + high_vegetation + low_vegetation + water_retention_capacity + surface_runoff + evaporation), 
                           data = pred_with_geo, scale = TRUE)

# Summarize results
summary(partial_rda_geo_env)

# Check adjusted R^2
RsquareAdj(partial_rda_geo_env)

# Test significance
paste("test significance -- Partial RDA for Geography Controlling for Environment")
anova.cca(partial_rda_geo_env, permutations = 999)

# plot the results
#png("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Landscape_genetics_GEA/RDA_Redundancy_Analysis/partial_rda_env_geo.png", width = 8, height = 6, units = "in", res = 300)
png("partial_rda_env_geo.png", width = 8, height = 6, units = "in", res = 300)
