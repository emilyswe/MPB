####Created February 11 2025 by Emily Swerdfager
####PREAMBLE: The purpose of this script is to merge the output from 4-mpb-wrangle-wildtrax-data with the covariate files for both regions


# Load libraries
library(tidyverse)

# 1. Read in Bird Data (Site x Species Matrix)
bird_data <- read.csv("/Users/Bronwyn/Documents/local-git/MPB/Output/01_MPB-bird-data-territoryA.csv")

# 2. Read in Covariate Data for Each Region
jnp_cov <- read.csv("Output/3-jnp-disturbances-tsd.csv")  # Jasper sites
sab_cov <- read.csv("Output/2_tsd-sab.csv")  # Southern Alberta sites

# 3. Join Bird Data with Covariates for JNP
jnp_bird_cov <- jnp_cov %>%
  inner_join(bird_data, by = "location")  # Merge by location

# 4. Join Bird Data with Covariates for SAB
sab_bird_cov <- sab_cov %>%
  inner_join(bird_data, by = "location")  # Merge by location

# 5. Save Final Datasets
write.csv(jnp_bird_cov, "/Users/Bronwyn/Documents/local-git/MPB/Output/05_MPB_bird_covariates_JNP.csv", row.names = FALSE)
write.csv(sab_bird_cov, "/Users/Bronwyn/Documents/local-git/MPB/Output/05_MPB_bird_covariates_SAB.csv", row.names = FALSE)

# END OF SCRIPT

jasperdat <- read.csv("Output/05_MPB_bird_covariates_JNP.csv")
sabdat <- read.csv("Output/05_MPB_bird_covariates_SAB.csv")
