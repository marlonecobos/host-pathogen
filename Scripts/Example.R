# ------------------------------------------------------------------------------
# Project: DETECTING SIGNALS OF SPECIES’ ECOLOGICAL NICHES IN RESULTS OF STUDIES
#          WITH DEFINED SAMPLING PROTOCOLS: EXAMPLE APPLICATION TO PATHOGEN NICHES
# Author: Marlon E. Cobos and A. Townsend Peterson
# Date: 09/02/2022
# ------------------------------------------------------------------------------

# Description ------------------------------------------------------------------
# This script contains an example of how to use the R functions developed to 
# run the analyses proposed. Plots designed to help with interpretations are
# also produced in this example.
#
# This example corresponds to Scenario 3 from the seven examples presented in 
# the manuscript. In this case, the data indicates that the signal detected for
# the pathogen ecological niche is different from random expectations 
# considering the host.
# ------------------------------------------------------------------------------


# Functions to run analyses to detect niche signal -----------------------------
# To see a complete description of all arguments that can be set in these 
# functions please see the files in the repository.
source("https://raw.githubusercontent.com/marlonecobos/host-pathogen/main/Scripts/niche_signal.R")
source("https://raw.githubusercontent.com/marlonecobos/host-pathogen/main/Scripts/plot_niche_signal.R")
# ------------------------------------------------------------------------------


# Working directory ------------------------------------------------------------
setwd("YOUR/DIRECTORY") # optional
# ------------------------------------------------------------------------------


# Example data -----------------------------------------------------------------
example_url <- "https://github.com/jsoberon/PAMs-Mexico/blob/master/PAM_IUCN_MEX.RData?raw=true"
load(url(example_url)) # the object will be scenario3
# ------------------------------------------------------------------------------


# Detecting niche signals ------------------------------------------------------
# This part performs the analyses described in the protocols

# multivariate approach (PERMANOVA)
m_comp <- niche_signal(data = scenario3, condition = "Infection", 
                       variables = c("Temperature", "Precipitation"), 
                       method = "permanova")

# univariate approach
## tenmperature
u_temp <- niche_signal(data = scenario3, condition = "Infection", 
                       variables = "Temperature", method = "univariate")

## precipitation
u_prec <- niche_signal(data = scenario3, condition = "Infection", 
                       variables = "Precipitation", method = "univariate")
# ------------------------------------------------------------------------------


# Checking and plotting results (multivariate approach) ------------------------
# results from the multivariate approach are contained in a portion of the 
# resulting list.

# checking results from PERMANOVA
m_comp$permanova_results

# Plotting results (let's see two options: points and ellipses)
plot_niche_signal(m_comp) # red color = pathogen
plot_niche_signal(m_comp, ellipses = TRUE, lty = 1)

# ellipsoids will give us a hint about point density, which is hard to see in 
# if only points are overlapped. 

## both options together
opar <- par(no.readonly = TRUE)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 0.5, 0.5)) # setting plotting area

plot_niche_signal(m_comp)
plot_niche_signal(m_comp, ellipses = TRUE, lty = 1)

par(opar) # setting plotting area back to original
# ------------------------------------------------------------------------------


# Checking and plotting results (univariate approach) --------------------------
# a summary of the results from the univariate approach that help a lot with 
# interpretations are also in the list resulting from analyses.

# checking results from univariate tests
# the summary will tell you which null hypotheses were accepted and rejected, 
# when rejected it will say whether the statistic of the positive records 
# (pathogen in this case) is lower or higher than random expectations.

## temperature
u_temp$univariate_results$hypothesis_test

## precipitation
u_prec$univariate_results$hypothesis_test


# Plotting results (four statistics are tested for each variable)
# lets plot one example
plot_niche_signal(u_temp, statistic = "mean") # blue = observed; black = 95% CI

# all results for Temperature
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 0.5, 0.5)) # setting plotting area

plot_niche_signal(u_temp, statistic = "mean")
plot_niche_signal(u_temp, statistic = "SD")
plot_niche_signal(u_temp, statistic = "median")
plot_niche_signal(u_temp, statistic = "range")

# and all results for Precipitation
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 0.5, 0.5)) # setting plotting area

plot_niche_signal(u_prec, statistic = "mean")
plot_niche_signal(u_prec, statistic = "SD")
plot_niche_signal(u_prec, statistic = "median")
plot_niche_signal(u_prec, statistic = "range")

par(opar) # setting plotting area back to original
# ------------------------------------------------------------------------------