# ------------------------------------------------------------------------------
# Project: DETECTING SIGNALS OF SPECIES’ ECOLOGICAL NICHES IN RESULTS OF STUDIES
#          WITH DEFINED SAMPLING PROTOCOLS: EXAMPLE APPLICATION TO PATHOGEN NICHES
# Author: Marlon E. Cobos and A. Townsend Peterson
# Date: 10/26/2021
# ------------------------------------------------------------------------------

# Description ------------------------------------------------------------------
# This script contains code to reproduce analyses and figures for this project
# 
# All data required can be obtained using the code in this script. Make sure to 
# have internet connection to download environmental data.
#
# Note: some results will be written in your working directory.
# ------------------------------------------------------------------------------

# R packages needed ------------------------------------------------------------
# the package evniche is in a GitHub repository. To install it use the 
# following lines of code
if (!require(remotes)) {
  install.packages("remotes")
}

remotes::install_github("marlonecobos/evniche")

# all other packages are on CRAN and cab be installed as:
# install.packages("package_name") 

# loading packages
library(evniche)
library(raster)
library(ellipse)
library(scales)
library(maps)
library(biosurvey)
# ------------------------------------------------------------------------------


# Functions to run analyses of similarity --------------------------------------
source("https://raw.githubusercontent.com/marlonecobos/host-pathogen/main/Scripts/niche_signal.R")
source("https://raw.githubusercontent.com/marlonecobos/host-pathogen/main/Scripts/plot_niche_signal.R")
# ------------------------------------------------------------------------------


# Working directory ------------------------------------------------------------
setwd("YOUR/DIRECTORY")
# ------------------------------------------------------------------------------


# Getting and preparing data ---------------------------------------------------
# This part is to prepare environmental data to be used in further analyses

# new folder for data
dir.create("Data")

# downloading environmental variables
bios <- getData("worldclim", var = "bio", res = 10, path = "Data")

# keeping only temperature and precipitation and masking to SA (BIO 1 and 12)
## SA extent
SAext <- extent(c(-85, -30, -60, 15))

## cropping and renaming
bio_SA <- crop(bios[[c(1, 12)]], SAext)
names(bio_SA) <- c("Temperature", "Precipitation")

## dividing temperature by 10
bio_SA$Temperature <- bio_SA$Temperature / 10 

# preparing data for analyses (from raster to matrix)
data_T_P <- rasterToPoints(bio_SA)

# saving data (complete)
write.csv(data_T_P, "Data/environmental_data.csv", row.names = FALSE)

# sampling data to make it easier when simulating data
set.seed(1)
data_T_P <- data_T_P[sample(nrow(data_T_P), 5000), ]

# saving data (sample)
write.csv(data_T_P, "Data/environmental_data_sample.csv", row.names = FALSE)
# ------------------------------------------------------------------------------


# Characteristics of virtual niches --------------------------------------------
# Here we define the characteristics of virtual niches (limits)

# host niche (always the same niche)
## niche range (limits)
host_niche_range <- cbind(Temperature = c(12, 26), Precipitation = c(700, 2800))

# pathogen (7 scenarios)
## niche range (limits) (various examples of levels of similarity and 
## type of pathogen niche)
path_niche_ranges <- list(
  scenario_1 = cbind(Temperature = c(12, 26), # exactly the same niche
                     Precipitation = c(700, 2800)),
  scenario_2 = cbind(Temperature = c(14, 28), # the position changes 
                     Precipitation = c(800, 2900)),
  scenario_3 = cbind(Temperature = c(13, 22), # the limits and position are different 
                     Precipitation = c(900, 2600)),
  scenario_4  = cbind(Temperature = c(15, 23), # the limits change slightly 
                      Precipitation = c(1000, 2500)),
  scenario_5 = cbind(Temperature = c(18, 25), # position and limits changes
                     Precipitation = c(1000, 4000)),
  scenario_6 = cbind(Temperature = c(14, 29), # position, limits and direction changes
                     Precipitation = c(600, 4500)),
  scenario_7 = cbind(Temperature = c(8, 30), # a bigger, different niche for pathogen that 
                     Precipitation = c(200, 3200)) # could be detected as similar?
)
# ------------------------------------------------------------------------------

# Creating virtual niches ------------------------------------------------------
# This portion of the script creates virtual niches based on limits defined above

# host virtual niche
## variances
vars <- evniche:::var_from_range(range = host_niche_range)

## covariance limits
cov_lim <- evniche:::covariance_limits(range = host_niche_range)

## variance-covariance matrix
cov <- cov_lim$max_covariance * 0.2 # covariance selected
varcov <- evniche:::var_cov_matrix(variances = vars, covariances = cov) 

## centroid
cent <- evniche:::centroid(range = host_niche_range)

## ellipsoid characteristics (virtual niche)
host_niche <- ell_features(centroid = cent, covariance_matrix = varcov,
                           level = 0.99)


# pathogen virtual niches (for all scenarios)
path_niches <- lapply(1:length(path_niche_ranges), function(x) {
  ## variances
  vars <- evniche:::var_from_range(range = path_niche_ranges[[x]])
  
  ## covariance limits
  cov_lim <- evniche:::covariance_limits(range = path_niche_ranges[[x]])
  
  ## variance-covariance matrix
  cov <- cov_lim$max_covariance * ifelse(x == 6, 0.5, 0.2) # covariance selected
  varcov <- evniche:::var_cov_matrix(variances = vars, covariances = cov)
  
  ## centroid
  cent <- evniche:::centroid(range = path_niche_ranges[[x]])
  
  ## ellipsoid characteristics (virtual niche)
  ell_features(centroid = cent, covariance_matrix = varcov, level = 0.99)
})

names(path_niches) <- names(path_niche_ranges)
# ------------------------------------------------------------------------------


# Generating data for analyses derived from virtual niches ---------------------
# Now we generate host and pathogen records based on virtual niches and 
# available envrionmental conditions 

# host data
## predict suitability on bioclimatic data based on ellipsoids
pred_host <- evniche:::ell_predict(data = data_T_P, features = host_niche, 
                                   longitude = "x", latitude = "y")

## generate new data
## based on the ellipsoid and available conditions
set.seed(1)
vd_pre_host <- virtual_data(features = host_niche, from = "prediction",
                            data = data_T_P, prediction = pred_host, n = 200)

# pathogen data (for all scenarios)
vd_pre_path <- lapply(1:length(path_niches), function(x) {
  ## predict suitability
  predp <- evniche:::ell_predict(data = data_T_P, features = path_niches[[x]], 
                                 longitude = "x", latitude = "y")
  
  ## generate new data
  ## based on the ellipsoid and available conditions
  set.seed(1)
  vd <- virtual_data(features = path_niches[[x]], from = "prediction",
                     data = data_T_P, prediction = predp, n = 400)
  
  if (x == 1) { # sample of pathogen niche to avoid total coincidence with host
    vd <- vd[sample(nrow(vd), 200), ]
  }
  vd
})

names(vd_pre_path) <- names(path_niche_ranges)
# ------------------------------------------------------------------------------


# Visualization of host pathogen niche combinations ----------------------------
# This part generates plots to check how virtual niches look like (all scenarios)

# host ellipsoid niche
host_ell <- ellipse(x = host_niche$covariance_matrix,
                    centre = host_niche$centroid,
                    level = host_niche$level)

# visualizations considering host and pathogen niches (plots will show up in 
# new windows)
path_ells <- lapply(1:length(vd_pre_path), function(x) {
  ## ellipsoid
  path_ell <- ellipse(x = path_niches[[x]]$covariance_matrix,
                      centre = path_niches[[x]]$centroid,
                      level = path_niches[[x]]$level)
  
  # this figures are to check a each combination host - pathogen
  x11()
  plot(data_T_P[, 3:4], pch = 16, col = alpha("gray45", 0.4), # background
       main = "Host (black) and pathogen (red) niches") 
  lines(host_ell, col = "black", lwd = 2) # host ellipsoid
  points(vd_pre_host[, 3:4], pch = 16, cex = 1.3, col = "black") # host points
  
  lines(path_ell, col = "red") # pathogen ellipsoid
  points(vd_pre_path[[x]][, 3:4], pch = 16, cex = 0.6, col = "red") # pathogen points
  
  path_ell
}) 

names(path_ells) <- names(path_niche_ranges)
# ------------------------------------------------------------------------------


# Putting data together for analyses -------------------------------------------
# This generates the datasets needed for running the protocols for pathogen
# niche detection

data_host_path <- lapply(1:length(path_niches), function(x) {
  ## infected hosts (points in host that coincide with points in pathogen)
  infection <- paste(vd_pre_host[, 1], vd_pre_host[, 2]) %in% 
    paste(vd_pre_path[[x]][, 1], vd_pre_path[[x]][, 2])
  
  ## data for comparisons:
  ## columns are (longitude, latitude, infection-(0, 1), variable1, variable2)
  data_niches <- data.frame(vd_pre_host[, 1:2], 
                            infection = as.numeric(infection),
                            vd_pre_host[, 3:4])
  
  data_niches <- data_niches[order(data_niches$infection, decreasing = TRUE), ]
  
  ## writing data for comparisons
  fname <- paste0("Data/data_host_pathogen_", names(path_niches)[x], ".csv")
  
  write.csv(data_niches, file = fname, row.names = FALSE)
  
  data_niches
})

names(data_host_path) <- names(path_niche_ranges)

# saving all as RData just in case
save(bio_SA, data_T_P, path_niches, host_niche, vd_pre_path, vd_pre_host, 
     data_host_path, file = "Data/data_host_pathogen1.RData")
# ------------------------------------------------------------------------------


# Detecting pathogen niche signals ---------------------------------------------
# This part performs the analyses described in the protocols

# univariate approach
uni_comps <- lapply(data_host_path, function(x) {
  ## tenmperature
  u_t <- niche_signal(data = x, condition = "infection", 
                      variables = "Temperature", method = "univariate")
  
  ## precipitation
  u_p <- niche_signal(data = x, condition = "infection", 
                      variables = "Precipitation", method = "univariate")
  
  list(Temperature = u_t, Precipitation = u_p)
})

# multivariate approach (PERMANOVA)
multi_comps <- lapply(data_host_path, function(x) {
  m_comp <- niche_signal(data = x, condition = "infection", 
                         variables = c("Temperature", "Precipitation"), 
                         method = "permanova")
})
# ------------------------------------------------------------------------------


# Summarizing results ----------------------------------------------------------
# directory for results
dir.create("Results")

# univariate results in a table
## getting the  results
u_res <- lapply(uni_comps, function(x) {
  rbind(
    c(variable = "Temperature", x$Temperature$univariate_results$hypothesis_test),
    c(variable = "Precipitation", x$Precipitation$univariate_results$hypothesis_test)
  )
})

## preparing the table of results
u_res <- do.call(rbind, u_res)
u_res <- data.frame(comparison = paste(rep(names(uni_comps), each = 2), "vs host"),
                    u_res)

## writing results
write.csv(u_res, "Results/summary_univariate1.csv", row.names = FALSE)

# multivariate results
## getting the  results
m_res <- lapply(1:length(multi_comps), function(x) {
  mtab <- as.data.frame(multi_comps[[x]]$permanova_results)
  
  ## writing results
  fname <- paste0("Results/permanova_", names(multi_comps)[x], "_vs_host.csv")
  write.csv(mtab, fname, row.names = FALSE)
  
  mtab
})

names(m_res) <- names(multi_comps)

# saving results as RData
save(uni_comps, multi_comps, u_res, m_res,
     file = "Results/niniche_comparisons1.RData")
# ------------------------------------------------------------------------------



# Figure to explain univariate results -----------------------------------------
# directory for figures
dir.create("Figures")

# figure
jpeg("Figures/univariate_explain1.jpg", width = 80, height = 50, units = "mm", 
     res = 600)
par(cex = 0.7, mar = c(4.5, 4.5, 1, 0.5))
plot_niche_signal(uni_comps$scenario_5$Temperature, breaks = 50, lty = c(1, 2))
dev.off()
# ------------------------------------------------------------------------------



# Figure representing all scenarios --------------------------------------------
# exporting figure
jpeg("Figures/all_niches1.jpg", width = 120, height = 195, units = "mm", res = 600)

## layout
matlay <- matrix(c(1, 2, 9, 10, 3, 4, 11, 12, 5, 6, 13, 14, 7, 8, 15, 16, 17, 17), 
                 nrow = 9, byrow = TRUE)

layout(matlay, heights = c(rep(c(1.2, 10), 4), 3))

## titles
par(cex = 0.6)
par(mar = rep(0, 4))
plot.new(); text(0.53, 0.3, labels = "Host niche", cex = 0.8)
plot.new(); text(0.53, 0.3, labels = "Scenario 1", cex = 0.8)
plot.new(); text(0.53, 0.3, labels = "Scenario 2", cex = 0.8)
plot.new(); text(0.53, 0.3, labels = "Scenario 3", cex = 0.8)
plot.new(); text(0.53, 0.3, labels = "Scenario 4", cex = 0.8)
plot.new(); text(0.53, 0.3, labels = "Scenario 5", cex = 0.8)
plot.new(); text(0.53, 0.3, labels = "Scenario 6", cex = 0.8)
plot.new(); text(0.53, 0.3, labels = "Scenario 7", cex = 0.8)

## host niche
par(mar = c(2.1, 2.1, 0, 0.5))
plot(data_T_P[, 3:4], pch = 16, col = alpha("gray55", 0.3), axes = FALSE) 
lines(host_ell, col = "black", lwd = 1.5) 
points(vd_pre_host[, 3:4], pch = 16, col = "black", cex = 0.6) 

axis(1, line = 0.03, tcl = -0.2, cex.axis = 0.6, mgp = c(0, 0.2, 0)) 
axis(2, line = 0.03, tcl = -0.2, cex.axis = 0.6, mgp = c(0, 0.2, 0))
title(ylab = "Precipitation", line = 1.1, cex.lab = 0.8)
box()

## pathogen niches
for (i in 1:length(path_ells)) {
  ## virtual niche E space
  plot(data_T_P[, 3:4], pch = 16, col = alpha("gray55", 0.3), axes = FALSE) 
  lines(host_ell, col = "black", lwd = 1.5) 
  points(vd_pre_host[, 3:4], pch = 16, col = "black", cex = 0.6) 
  lines(path_ells[[i]], col = "red", lwd = 1) 
  points(vd_pre_path[[i]][, 3:4], pch = 3, col = "red", cex = 0.3) 
  
  axis(1, line = 0.03, tcl = -0.2, cex.axis = 0.6, mgp = c(0, 0.2, 0)) 
  axis(2, line = 0.03, tcl = -0.2, cex.axis = 0.6, mgp = c(0, 0.2, 0))
  
  if (i %in% c(2, 4, 6)) {
    title(ylab = "Precipitation", line = 1.1, cex.lab = 0.8)
  }
  
  if (i %in% 6:7) {
    title(xlab = "Temperature", line = 1.1, cex.lab = 0.8)
  }
  
  box()
}

## legend
par(cex = 0.5)
par(mar = rep(0, 4)); plot.new() 
legend(0.05, 1, legend = c("Available conditions", "Host records", 
                           "Pathogen records"), 
       pch = c(16, 16, 3), pt.cex = c(1, 0.8, 0.5), 
       col = c(alpha("gray5", 0.3), "black", "red"), bty = "n")
legend(0.75, 0.8, legend = c("Host niche", "Pathogen niche"), lty = 1, 
       lwd = c(1.5, 1), col = c("black", "red"), bty = "n")

dev.off()
# ------------------------------------------------------------------------------


# Visualization of prepared data -----------------------------------------------
# This is how the data of a host and a pathogen would look like in real life
# we know where the host is and which which individual hosts have the pathogens
jpeg("Figures/all_cases_in_reality1.jpg", width = 166, height = 95, 
     units = "mm", res = 600)

par(mfrow = c(2, 4), mar = c(4, 4, 2, 0.1))
for (i in 1:length(data_host_path)) {
  xlb <- ifelse(i %in% c(1:3), "", "Temperature") 
  ylb <- ifelse(i %in% c(2:4, 6:7), "", "Precipitation") 
  par(cex = 0.5)
  
  # all host records
  plot(data_host_path[[i]][, 4:5], pch = 16, cex = 1.5, col = "black", 
       xlab = xlb, ylab = ylb)
  title(main = paste("Scenario", i), line = 0.7, font.main = 1)
  
  ## records of infected hosts
  points(data_host_path[[i]][data_host_path[[i]][, 3] == 1, 4:5], 
         pch = 3, cex = 2, col = "red")
}

# legend
plot.new()
legend("center", legend = c("Host records", "Infected hosts"), 
       pch = c(16, 3), pt.cex = c(1.5, 2), col = c("black", "red"), 
       y.intersp = 2, bty = "n")

dev.off()
# ------------------------------------------------------------------------------


# Representation of multivariate results ---------------------------------------
# This is a representation of how the niches of host and pathogens will be 
# reconstructed, and the outcomes of multivariate tests
jpeg("Figures/all_multi_comps1.jpg", width = 166, height = 95, 
     units = "mm", res = 600)

par(mfrow = c(2, 4), mar = c(4, 4, 2, 0.1))
for (i in 1:length(multi_comps)) {
  xlb <- ifelse(i %in% c(1:3), "", "Temperature") 
  ylb <- ifelse(i %in% c(2:4, 6:7), "", "Precipitation") 
  par(cex = 0.5)
  
  ## plot
  plot_niche_signal(multi_comps[[i]], ellipses = TRUE, xlab = xlb, ylab = ylb,
                    lty = 1)
  title(main = paste("Scenario", i), line = 0.7, font.main = 1)
  
  pv <- ifelse(m_res[[i]]$`Pr(>F)`[1] < 0.05, "< 0.05", 
               paste("=", round(m_res[[i]]$`Pr(>F)`[1], 2)))
  legend("topleft", legend = substitute(paste(italic("p"), x), 
                                        env = list(x = pv)), bty = "n", 
         cex = 0.9, inset = 0.01)
}

plot.new()
legend("center", legend = c("Reconstructed host", "Reconstructed pathogen"), 
       lty = 1, lwd = 1.5, col = c("black", "red"), 
       y.intersp = 2, bty = "n")

dev.off()
# ------------------------------------------------------------------------------



# Figure with all cases and geography ------------------------------------------
# This shows virtyual niches in environmental and geographic spaces
# generating geographic predictions of virtual niches 
# host prediction 
pred_host_g <- evniche:::ell_predict(data = bio_SA, features = host_niche)

# pathogen prediction
pred_path_g <- lapply(path_niches, function(x) {
  evniche:::ell_predict(data = bio_SA, features = x)
})

# elements of figure
## titles and subtitles
mains <- c("Host niche", "Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4",
           "Scenario 5", "Scenario 6", "Scenario 7")

subs <- c("Virtual niche", "Geographic projection")

## color palette
suit_col <- bluered(255)

## limits for geographic plot
boxpam <- t(bbox(bio_SA))
boxpam <- SpatialPointsDataFrame(boxpam, data.frame(boxpam), 
                                 proj4string = bio_SA@crs)

# exporting figure
jpeg("Figures/all_cases1.jpg", width = 166, height = 195, units = "mm", res = 600)

## layout
matlay <- matrix(c(1:6, 6, 7, 8, 8, 18:22, 9, 9, 10, 11, 11, 23:27, 12, 12, 13, 
                   14, 14, 28:32, 15, 15, 16, 17, 17, 33:37, 38, 38, 39, 40, 40), 
                 nrow = 10, byrow = TRUE)

layout(matlay, widths = c(8, 8, 0.5, 8, 8), 
       heights = c(1.2, rep(c(1, 10), 4), 2.5))

## main
par(cex = 0.7)
par(mar = rep(0, 4)); plot.new(); text(0.5, 0.4, labels = subs[1])
plot.new(); text(0.5, 0.4, labels = subs[2])
plot.new(); plot.new(); text(0.5, 0.4, labels = subs[1])
plot.new(); text(0.5, 0.4, labels = subs[2])

## subtitles
par(cex = 0.65)
for (i in 1:4) {
  plot.new(); text(0.5, 0.4, labels = mains[((i * 2) - 1)])
  plot.new(); plot.new(); text(0.5, 0.4, labels = mains[(i * 2)])
}

## host niche
## virtual niche E space
par(mar = c(2.1, 2.1, 0.5, 0.5))
plot(data_T_P[, 3:4], pch = 16, col = alpha("gray55", 0.3), axes = FALSE) 
lines(host_ell, col = "black", lwd = 1.5) 
points(vd_pre_host[, 3:4], pch = 16, col = "black", cex = 0.6) 

axis(1, line = 0.005, tcl = -0.2, cex.axis = 0.6, mgp = c(0, 0.2, 0)) 
axis(2, line = 0.005, tcl = -0.2, cex.axis = 0.6, mgp = c(0, 0.2, 0))
title(ylab = "Precipitation", line = 1.1, cex.lab = 0.8)
box()

# virtual niche G prediction
par(mar = c(2.1, 0.5, 0.5, 0.5))
plot(boxpam, col = NA)
image(pred_host_g$suitability_trunc, col = suit_col, add = TRUE)
box()

par(mar = rep(0, 4)); plot.new()

## pathogen niches
for (i in 1:7) {
  ## virtual niche E space
  par(mar = c(2.1, 2.1, 0.5, 0.5))
  plot(data_T_P[, 3:4], pch = 16, col = alpha("gray55", 0.3), axes = FALSE) 
  lines(host_ell, col = "black", lwd = 1.5) 
  points(vd_pre_host[, 3:4], pch = 16, col = "black", cex = 0.6) 
  lines(path_ells[[i]], col = "red", lwd = 1) 
  points(vd_pre_path[[i]][, 3:4], pch = 3, col = "red", cex = 0.3) 
  
  axis(1, line = 0.005, tcl = -0.2, cex.axis = 0.6, mgp = c(0, 0.2, 0)) 
  axis(2, line = 0.005, tcl = -0.2, cex.axis = 0.6, mgp = c(0, 0.2, 0))
  title(ylab = "Precipitation", line = 1.1, cex.lab = 0.8)
  
  if (i %in% 6:7) {
    title(xlab = "Temperature", line = 1.1, cex.lab = 0.8)
  }
  
  box()
  
  # virtual niche G prediction
  par(mar = c(2.1, 0.5, 0.5, 0.5))
  plot(boxpam, col = NA)
  image(pred_path_g[[i]]$suitability_trunc, col = suit_col, add = TRUE)
  box()
  
  if (i %in% c(2, 4, 6)) {
    par(mar = rep(0, 4)); plot.new() 
  }
}

## legend
par(cex = 0.5)
par(mar = rep(0, 4)); plot.new() 
legend(0.1, 1, legend = c("Background", "Host records", "Pathogen records"), 
       pch = c(16, 16, 3), pt.cex = c(1, 0.8, 0.5), 
       col = c(alpha("gray5", 0.3), "black", "red"), bty = "n")
legend(0.6, 0.8, legend = c("Host niche", "Pathogen niche"), lty = 1, 
       lwd = c(1.5, 1), col = c("black", "red"), bty = "n")

plot.new(); plot.new() 
legend_bar("center", col = suit_col, width_prop = 0.5, heigh_prop = 0.15, 
           title = NULL, labels = c("Low", "High"), horizontal = TRUE, 
           cex = 1.5, labels_offset = 0.4)
text(0.5, 0.7, labels = "Suitability")

dev.off()
# ------------------------------------------------------------------------------



# Plotting univariate results --------------------------------------------------
# This shows results from univariate tests
# univariate mean comparison
jpeg("Figures/mean_univariate_all1.jpg", width = 166, height = 195, units = "mm",
     res = 600)

par(mfrow = c(7, 2), cex = 0.5, mar = c(3.5, 3.5, 0.5, 0.5))

for (i in 1:length(uni_comps)) {
  plot_niche_signal(uni_comps[[i]]$Temperature, xlab = "",
                    ylab = "", breaks = 30, lty = c(1, 2))
  title(xlab = ifelse(i == 7, "Temperature (mean)", ""), 
        ylab = paste("Scenario", i, "- Freq."), line = 2.3)
  if (i == 1) {
    legend("topright", legend = c("Observed", "CI 95%"), lty = c(1, 2), 
           col = c("blue", "black"), bty = "n", inset = -0.01, x.intersp = 0.5)
  }
  
  plot_niche_signal(uni_comps[[i]]$Precipitation, breaks = 30,
                    xlab = "", ylab = "", lty = c(1, 2))
  title(xlab = ifelse(i == 7, "Precipitation (mean)", ""), line = 2.3)
}

dev.off()


# univariate median comparison
jpeg("Figures/median_univariate_all1.jpg", width = 166, height = 195, units = "mm",
     res = 600)

par(mfrow = c(7, 2), cex = 0.5, mar = c(3.5, 3.5, 0.5, 0.5))

for (i in 1:length(uni_comps)) {
  plot_niche_signal(uni_comps[[i]]$Temperature, xlab = "", lty = c(1, 2), 
                    ylab = "", statistic = "median", breaks = 30)
  title(xlab = ifelse(i == 7, "Temperature (median)", ""), 
        ylab = paste("Scenario", i, "- Freq."), line = 2.3)
  if (i == 1) {
    legend("topright", legend = c("Observed", "CI 95%"), lty = c(1, 2), 
           col = c("blue", "black"), bty = "n", inset = -0.01, x.intersp = 0.5)
  }
  
  plot_niche_signal(uni_comps[[i]]$Precipitation, breaks = 30, lty = c(1, 2),
                    xlab = "", ylab = "", statistic = "median")
  title(xlab = ifelse(i == 7, "Precipitation (median)", ""), line = 2.3)
}

dev.off()


# univariate SD comparison
jpeg("Figures/SD_univariate_all1.jpg", width = 166, height = 195, units = "mm",
     res = 600)

par(mfrow = c(7, 2), cex = 0.5, mar = c(3.5, 3.5, 0.5, 0.5))

for (i in 1:length(uni_comps)) {
  plot_niche_signal(uni_comps[[i]]$Temperature, xlab = "", lty = c(1, 2), 
                    ylab = "", statistic = "SD", breaks = 30)
  title(xlab = ifelse(i == 7, "Temperature (SD)", ""), 
        ylab = paste("Scenario", i, "- Freq."), line = 2.3)
  
  plot_niche_signal(uni_comps[[i]]$Precipitation, breaks = 30, lty = c(1, 2), 
                    xlab = "", ylab = "", statistic = "SD")
  title(xlab = ifelse(i == 7, "Precipitation (SD)", ""), line = 2.3)
  
  if (i == 1) {
    legend("topleft", legend = c("Observed", "CI 95%"), lty = c(1, 2), 
           col = c("blue", "black"), bty = "n", inset = -0.01, x.intersp = 0.5)
  }
}

dev.off()



# univariate range comparison
jpeg("Figures/range_univariate_all1.jpg", width = 166, height = 195, units = "mm",
     res = 600)

par(mfrow = c(7, 2), cex = 0.5, mar = c(3.5, 3.5, 0.5, 0.5))

for (i in 1:length(uni_comps)) {
  plot_niche_signal(uni_comps[[i]]$Temperature, xlab = "", lty = c(1, 2),
                    ylab = "", statistic = "range", breaks = 30)
  title(xlab = ifelse(i == 7, "Temperature (range)", ""), 
        ylab = paste("Scenario", i, "- Freq."), line = 2.3)
  if (i == 1) {
    legend("topleft", legend = c("Observed", "CI 95%"), lty = c(1, 2), 
           col = c("blue", "black"), bty = "n", inset = -0.01, x.intersp = 0.5)
  }
  
  plot_niche_signal(uni_comps[[i]]$Precipitation, breaks = 30, lty = c(1, 2), 
                    xlab = "", ylab = "", statistic = "range")
  title(xlab = ifelse(i == 7, "Precipitation (range)", ""), line = 2.3)
}

dev.off()
# ------------------------------------------------------------------------------