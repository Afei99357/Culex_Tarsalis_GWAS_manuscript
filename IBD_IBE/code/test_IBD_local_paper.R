### Note tCtarsalis_sample_w_GPS_climate_new_filtered_idhat a lot of this code (especially the model testing)
### comes from bookdown.org/hhwagner1/LandGenCourse
### Set working directory and environment variables
#setwd("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/R_code")
#options(stringsAsFactors = FALSE)
library(geosphere)
library(vegan)
library(Fragman)
library(extrafont)
library(cowplot)
library(corMLPE)

rm(list = ls())

### Read in the pairwise Fst Calculations
fst = read.csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/IBD_IBE_test/pairwise_fst_table.csv", header=TRUE)

### Read in the Latitude and Longitudes
### Fix come of the locIDs
coords = read.csv("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Ctarsalis_sample_w_GPS_climate_average_new_filtered_id_region.csv")

### Simplify the coords table to only have 1 entry per population
pop.coords = coords[!duplicated(coords[,c('locID')]),]

### Get the Lats and Long for the first pop in each comparison
m = match(fst$Pop1, pop.coords$popID)
fst$Lat1 = pop.coords$GPS.Lat[m]
fst$Lon1 = pop.coords$GPS.Lon[m]

### Now get the coordinates for the second pop in each comparison
m2 = match(fst$Pop2, pop.coords$popID)
fst$Lat2 = pop.coords$GPS.Lat[m2]
fst$Lon2 = pop.coords$GPS.Lon[m2]

### Now calculate the Geographic Distances
Dgeo = rep(0, nrow(fst))
for (i in 1:nrow(fst)) {
  my.dist = distm(c(fst$Lon1[i], fst$Lat1[i]), c(fst$Lon2[i], fst$Lat2[i]), fun = distHaversine)
  Dgeo[i] = my.dist
}

### Pull out the desired Fst value and transform
Dgen = as.vector(fst$PairwiseFst)
Dgen = Dgen/(1-Dgen)

### Convert the table into 2 matrices (to use in Mantel test)
temp.df = data.frame(Pop1 = fst$Pop1, Pop2=fst$Pop2, Dgen=Dgen, Dgeo=Dgeo)
nams <- with(temp.df, unique(c(as.character(Pop1), as.character(Pop2))))
mat.gen = matrix(0, ncol=length(nams), nrow=length(nams))
mat.geo = matrix(0, ncol=length(nams), nrow=length(nams))

for (i in 1:length(nams)) {
  for (j in 1:length(nams)) {
    
    ## if the 2 population IDs are the same, then skip to the next ID
    if (nams[i]==nams[j]) {
      next
    } else {
      
      ## Figure out which ID should be pop1 and which should be pop2
      ## Then get the row that actually has data and save in the matrices
      x = which(temp.df$Pop1==nams[i] & temp.df$Pop2==nams[j])
      y = which(temp.df$Pop1==nams[j] & temp.df$Pop2==nams[i])
      
      if ((length(x)) > length(y)) {
        mat.gen[i,j] = temp.df$Dgen[x]
        mat.geo[i,j] = temp.df$Dgeo[x]
      } else {
        mat.gen[i,j] = temp.df$Dgen[y]
        mat.geo[i,j] = temp.df$Dgeo[y]
      }
    }
  }
}

### Mantel test of correlation with Physical Distance
IBD.geo <- vegan::mantel(mat.gen,mat.geo, method="spearman", permutations=9999)

IBD.geo

# get the population info
row_names <- pop.coords$popID

### adding environmental data, remove the highly correlated variables from RDA analysis r script (avg_ssr, avg_evabs, avg_e, avg_tp)
clim.table <- pop.coords[, c("avg_u10", "avg_v10", "avg_t2m", "avg_lai_hv", "avg_lai_lv", "avg_src", "avg_sro", "avg_e")]

# standardize the data before Euclildian distance
# standardizw_matrix <- round(scale(clim.table[, -1]), 4)

### Next, get the Canberra distances based on the mean environmental variables
dist.env <- dist(clim.table[, -1], method = 'canberra')
IBE.env <- vegan::mantel(mat.gen,dist.env, method="spearman", permutations=9999)

IBE.env

### Get the environmental distances into a column of the dataframe with the Genetic Distance
### and Geographic distance variables
env.mat = as.matrix(dist.env)
temp.df$Denv = rep(NA, nrow(temp.df))
for (i in 1:nrow(temp.df)) {
  my.row = which(nams==temp.df$Pop1[i])
  my.col = which(nams==temp.df$Pop2[i])
  temp.df$Denv[i] = env.mat[my.row,my.col]
}

## Denv
Denv <- temp.df$Denv

### Scale the data
temp.df$Dgenz <- scale(temp.df$Dgen, center=TRUE, scale=TRUE)
temp.df$Dgeoz <- scale(temp.df$Dgeo, center=TRUE, scale=TRUE)
temp.df$Denvz <- scale(temp.df$Denv, center=TRUE, scale=TRUE)

### Check for colinearity
CSF.df = temp.df[,6:8]
usdm::vif(CSF.df)
# Variables      VIF
# 1     Dgenz 1.314403
# 2     Dgeoz 1.310785
# 3     Denvz 1.006497

### Run a Mixed Model to look at Both Distance AND Environment together
mixMod <- nlme::gls(Dgenz ~ Dgeoz + Denvz, 
                    correlation=corMLPE::corMLPE(form=~Pop1+Pop2), 
                    data=temp.df, method="REML")

### Change predictors to just distance or just environment
modD = update(mixMod, ~Dgeoz)
modE = update(mixMod, ~Denvz)

### Compare evidence for the different models
mod1noREML = update(mixMod, method="ML")
mod2noREML = update(modD, method="ML")
mod3noREML = update(modE, method="ML")

Models <- list(Full=mod1noREML, Distance=mod2noREML, Environment=mod3noREML)
CSF.IC <- data.frame(AIC = sapply(Models, AIC),
                     BIC = sapply(Models, BIC)) 
CSF.IC <- data.frame(CSF.IC, k = sapply(Models, function(ls) attr(logLik(ls), "df")))

N = nrow(temp.df)  # Number of unique pairs
CSF.IC$AICc <- CSF.IC$AIC + 2*CSF.IC$k*(CSF.IC$k+1)/(N-CSF.IC$k-1)

AICcmin <- min(CSF.IC$AICc)
RL <- exp(-0.5*(CSF.IC$AICc - AICcmin))
sumRL <- sum(RL)
CSF.IC$AICcmin <- RL/sumRL

BICmin <- min(CSF.IC$BIC)
RL.B <- exp(-0.5*(CSF.IC$BIC - BICmin))
sumRL.B <- sum(RL.B)
CSF.IC$BICew <- RL.B/sumRL.B
round(CSF.IC,3)

library(ggplot2)
library(ggpubr)
library(MASS) 

dens <- MASS::kde2d(Dgeo, Dgen, n=300)
dens_2 <- MASS::kde2d(Denv, Dgen, n=300)
myPal <- colorRampPalette(c("white", "blue", "gold", "orange", "red"))

# Convert density data for plots 3 and 4
dens_df <- expand.grid(x = dens$x, y = dens$y)
dens_df$z <- as.vector(dens$z)

dens2_df <- expand.grid(x = dens_2$x, y = dens_2$y)
dens2_df$z <- as.vector(dens_2$z)

# Common theme adjustments for consistency
common_theme <- theme(text = element_text(size=24, family="Arial"),
                      plot.title = element_text(hjust = 0.5, size = 30, family="Arial"),
                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                      axis.title = element_text(size=16, family="Arial"))

# Plot 1: Geographic Distance vs. Genetic Distance
p1 <- ggscatter(temp.df, x="Dgeo", y="Dgen", color="#0072B2",
                xlab="Geographic Distance",
                ylab="Genetic Distance",
                add = "reg.line",
                add.params = list(color="#CC79A7", fill="lightpink"),
                conf.int = TRUE,
                cor.coef = FALSE) + 
  ggtitle("(A)") +
  common_theme

# Plot 2: Environmental Distance vs. Genetic Distance
p2 <- ggscatter(temp.df, x="Denv", y="Dgen", color="#0072B2",
                xlab="Environmental Distance",
                ylab="Genetic Distance",
                add = "reg.line",
                add.params = list(color="#CC79A7", fill="lightpink"),
                conf.int = TRUE,
                cor.coef = FALSE) + 
  ggtitle("(B)") +
  common_theme

# Plot 3: Geographic Distance with Density
p3 <- ggplot(temp.df, aes(x = Dgeo, y = Dgen)) +
  geom_point(color = "#0072B2", size = 0.5) +
  geom_raster(data = dens_df, aes(x = x, y = y, fill = z), alpha = 0.7) +
  scale_fill_gradientn(colours = myPal(300)) +
  geom_smooth(method = "lm", color = "#CC79A7", se = FALSE) +
  geom_smooth(color = "red", se = FALSE) +
  labs(x = "Geographic Distance", y = "Genetic Distance", title = "(C)") +
  theme_minimal() +
  common_theme +
  theme(legend.position = c(0.9, 0.2)) +
  theme(legend.title = element_blank())

# Plot 4: Environmental Distance with Density
p4 <- ggplot(temp.df, aes(x = Denv, y = Dgen)) +
  geom_point(color = "#0072B2", size = 0.5) +
  geom_raster(data = dens2_df, aes(x = x, y = y, fill = z), alpha = 0.7) +
  scale_fill_gradientn(colours = myPal(300)) +
  geom_smooth(method = "lm", color = "#CC79A7", se = FALSE) +
  geom_smooth(color = "red", se = FALSE) +
  labs(x = "Environmental Distance", y = "Genetic Distance", title = "(D)") +
  theme_minimal() +
  common_theme +
  theme(legend.position = c(0.9, 0.2)) +
  theme(legend.title = element_blank())

# Save plots to PNG
png("/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/IBD_IBE_test/r_code/IBD_IBE_Plots.png", width=20, height=20, units="in", res=300)

# Arrange and save the combined plot
plot_grid(p1, p2, p3, p4, ncol = 2)

dev.off()

