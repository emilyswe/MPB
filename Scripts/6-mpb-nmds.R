#########Created February 12 2025 by Emily Swerdfager
####PREAMBLE: The purpose of this script is to perform a preliminary NMDS on the MPB bird data for jasper

# Load required packages
library(dplyr)
library(vegan)
library(ggplot2)


###########################
###########################
###########################
###########################
###########################
#JASPER DATASET NMDS:

# Read in dataset
bird_data <- read.csv("Output/05-mpb-jnp-dataset.csv")

# Define species columns (they follow the 4-letter convention)
species_columns <- grep("^[A-Z]{4}$", colnames(bird_data), value = TRUE)

# Subset NMDS data (keeping all rows, even those with zero abundance)
nmds_data <- bird_data %>%
  select(location, treatment, TSD, latitude, longitude, all_of(species_columns))

# Extract species matrix (excluding metadata)
species_matrix <- as.matrix(nmds_data[, species_columns])

# Run NMDS using Bray-Curtis distance (without removing zeros)
set.seed(999)
nmds_result <- metaMDS(species_matrix, distance = "bray", k = 2, trymax = 100)

# Extract NMDS scores and add to nmds_data
nmds_coords <- as.data.frame(scores(nmds_result, display = "sites"))
nmds_data$mdsA <- nmds_coords[, 1]
nmds_data$mdsB <- nmds_coords[, 2]

# Save NMDS results
write.csv(nmds_data[, c("location", "treatment", "TSD", "mdsA", "mdsB")], 
          "Output/6-mpb_nmds_2axis.csv", row.names = FALSE)

####### **Plot NMDS results**

# Read in NMDS results
dat3 <- read.csv("Output/6-mpb_nmds_2axis.csv")

# Plot NMDS grouped by treatment (burn, standing, harvest)
plot3 <- ggplot(dat3, aes(x = mdsA, y = mdsB)) +
  stat_ellipse(aes(colour = treatment), level = 0.68, linewidth = 1) +
  scale_color_manual(values = c("firebrick", "forestgreen", "blue"), 
                     labels = c("Burn", "Standing", "Harvest"), name = "Treatment") +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(x = "NMDS Axis 1", y = "NMDS Axis 2", title = "NMDS of Jasper Bird Communities by MPB Treatment")

# Save plot
ggsave(plot3, filename = "Output/5-mpb-jnp-nmds-mean.jpeg", width = 8, height = 6)


###########################
###########################
###########################
###########################
###########################
# SOUTHERN ALBERTA DATASET NMDS (sab):

# Load necessary libraries
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)

# Read in dataset
bird_data <- read.csv("Output/05-mpb-sab-dataset.csv")

# Define species columns (they follow the 4-letter convention)
species_columns <- grep("^[A-Z]{4}$", colnames(bird_data), value = TRUE)

# Subset NMDS data (keeping all rows, even those with zero abundance)
nmds_data <- bird_data %>%
  select(location, treatment, TSD, latitude, longitude, all_of(species_columns))

# Extract species matrix (excluding metadata)
species_matrix <- as.matrix(nmds_data[, species_columns])

# Run NMDS using Bray-Curtis distance (without removing zeros)
set.seed(999)
nmds_result <- metaMDS(species_matrix, distance = "bray", k = 2, trymax = 100)

# Extract NMDS scores and add to nmds_data
nmds_coords <- as.data.frame(scores(nmds_result, display = "sites"))
nmds_data$mdsA <- nmds_coords[, 1]
nmds_data$mdsB <- nmds_coords[, 2]

# Save NMDS results
write.csv(nmds_data[, c("location", "treatment", "TSD", "mdsA", "mdsB")], 
          "Output/6-mpb_sab_nmds_2axis.csv", row.names = FALSE)

####### **Plot NMDS results**

# Read in NMDS results
dat3 <- read.csv("Output/6-mpb_sab_nmds_2axis.csv")

# Plot NMDS grouped by treatment (burn, standing, harvest)
plot3 <- ggplot(dat3, aes(x = mdsA, y = mdsB)) +
  stat_ellipse(aes(colour = treatment), level = 0.68, linewidth = 1) +
  scale_color_manual(values = c("firebrick", "forestgreen", "blue"), 
                     labels = c("Burn", "Standing", "Harvest"), name = "Treatment") +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(x = "NMDS Axis 1", y = "NMDS Axis 2", title = "NMDS of Southern Alberta Bird Communities by MPB Treatment")

# Save plot
ggsave(plot3, filename = "Output/5-mpb-sab-nmds-mean.jpeg", width = 8, height = 6)

####### **PERMANOVA Tests**
##################### Updated April 2 2025
######################
#JASPER:

#JASPER
# PERMANOVA + delta AIC Analysis for Jasper Bird Community Data

# Load packages
library(dplyr)
library(vegan)

# Read in dataset
bird_data <- read.csv("Output/05-mpb-jnp-dataset.csv")

# Define species columns (assumes all 4-letter codes)
species_columns <- grep("^[A-Z]{4}$", colnames(bird_data), value = TRUE)

# Subset data for NMDS and PERMANOVA
nmds_data <- bird_data %>%
  select(location, treatment, TSD, latitude, longitude, all_of(species_columns))

# Extract species matrix
species_matrix <- as.matrix(nmds_data[, species_columns])

# Extract covariates for PERMANOVA
env_data <- nmds_data %>%
  select(treatment, TSD)

# --- PERMANOVA Models ---
set.seed(999)
permanova_treatment <- adonis2(species_matrix ~ treatment, data = env_data, method = "bray", permutations = 999)

set.seed(999)
permanova_tsd <- adonis2(species_matrix ~ TSD, data = env_data, method = "bray", permutations = 999)

set.seed(999)
permanova_null <- adonis2(species_matrix ~ 1, data = env_data, method = "bray", permutations = 999)

# --- AIC Calculations ---
n <- nrow(species_matrix)
k_treatment <- length(unique(env_data$treatment)) - 1  # 3 groups = 2 df
k_tsd <- 1
k_null <- 0

r2_treatment <- permanova_treatment$R2[1]
r2_tsd <- permanova_tsd$R2[1]
r2_null <- permanova_null$R2[1]

aic_treatment <- n * log(1 - r2_treatment) + 2 * k_treatment
aic_tsd <- n * log(1 - r2_tsd) + 2 * k_tsd
aic_null <- n * log(1 - r2_null) + 2 * k_null

# --- ΔAIC Calculations ---
aic_all <- c(aic_treatment, aic_tsd, aic_null)
min_aic <- min(aic_all)
delta_aic <- aic_all - min_aic

# --- Console Output ---
cat("PERMANOVA AIC Results (Jasper):\n")
cat("AIC - Disturbance Type:", round(aic_treatment, 2), "\n")
cat("AIC - TSD:", round(aic_tsd, 2), "\n")
cat("AIC - Null:", round(aic_null, 2), "\n")
cat("ΔAIC Values:", round(delta_aic, 2), "\n")

# --- Save CSV ---
aic_results <- data.frame(
  Model = c("Disturbance Type", "TSD", "Null"),
  R2 = c(r2_treatment, r2_tsd, r2_null),
  AIC = round(aic_all, 2),
  Delta_AIC = round(delta_aic, 2)
)
write.csv(aic_results, "Output/5-mpb-jnp-permanova_AIC.csv", row.names = FALSE)


######################
######################
#SOUTHERN ALBERTA:
# Load required packages
library(vegan)
library(dplyr)

# Load the Southern Rockies dataset
bird_data <- read.csv("Output/05-mpb-sab-dataset.csv")

# Define species columns (4-letter codes)
species_columns <- grep("^[A-Z]{4}$", colnames(bird_data), value = TRUE)

# Subset relevant NMDS data
nmds_data <- bird_data %>%
  select(location, treatment, TSD, all_of(species_columns))

# Create species matrix
species_matrix <- as.matrix(nmds_data[, species_columns])

# Prepare environmental data
env_data <- nmds_data %>%
  select(treatment, TSD)

# PERMANOVA - Disturbance Type
set.seed(999)
permanova_treatment <- adonis2(
  species_matrix ~ treatment,
  data = env_data,
  method = "bray",
  permutations = 999
)

# PERMANOVA - TSD
set.seed(999)
permanova_tsd <- adonis2(
  species_matrix ~ TSD,
  data = env_data,
  method = "bray",
  permutations = 999
)

# PERMANOVA - Null Model
set.seed(999)
permanova_null <- adonis2(
  species_matrix ~ 1,
  data = env_data,
  method = "bray",
  permutations = 999
)

# AIC calculation
n <- nrow(species_matrix)
k_treatment <- length(unique(env_data$treatment)) - 1  # 2 df
k_tsd <- 1
k_null <- 0

r2_treatment <- permanova_treatment$R2[1]
r2_tsd <- permanova_tsd$R2[1]
r2_null <- permanova_null$R2[1]

aic_treatment <- n * log(1 - r2_treatment) + 2 * k_treatment
aic_tsd <- n * log(1 - r2_tsd) + 2 * k_tsd
aic_null <- n * log(1 - r2_null) + 2 * k_null

# ΔAIC calculation
raw_aics <- c(aic_treatment, aic_tsd, aic_null)
delta_aic <- raw_aics - min(raw_aics)

# Print results
cat("PERMANOVA AIC Results (Southern Rockies):\n")
cat("AIC - Disturbance Type:", round(aic_treatment, 2), "\n")
cat("AIC - TSD:", round(aic_tsd, 2), "\n")
cat("AIC - Null:", round(aic_null, 2), "\n")
cat("ΔAIC Values:", round(delta_aic, 2), "\n")

# Save to CSV
aic_results <- data.frame(
  Model = c("Disturbance Type", "TSD", "Null"),
  R2 = c(r2_treatment, r2_tsd, r2_null),
  AIC = raw_aics,
  DeltaAIC = delta_aic
)
write.csv(aic_results, "Output/5-mpb-sab-permanova_AIC.csv", row.names = FALSE)


library(ggplot2)

##################################
##################################
##################################
##plot JNP NMDS results with sites

# Load NMDS data
dat_jasper <- read.csv("Output/6-mpb_nmds_2axis.csv")  

# NMDS plot with sites as points + ellipses
ggplot(dat_jasper, aes(x = mdsA, y = mdsB, color = treatment)) +
  geom_point(size = 2, alpha = 0.7) +  # Sites as points
  stat_ellipse(level = 0.68, linewidth = 1) +  # Treatment ellipses
  scale_color_manual(values = c("firebrick", "forestgreen", "blue"),
                     labels = c("Burn", "Standing", "Harvest"),
                     name = "Treatment") +
  theme_minimal() +
  labs(title = "NMDS of Jasper Bird Communities by Disturbance Type",
       x = "NMDS Axis 1", y = "NMDS Axis 2")

# Save
ggsave("Output/5-mpb-jnp-nmds-sites.jpeg", width = 8, height = 6)

##################################
##################################
##################################
##plot SAB NMDS results with sites 

# Load NMDS data
dat_sab <- read.csv("Output/6-mpb_sab_nmds_2axis.csv")

# NMDS plot with sites as points + ellipses
ggplot(dat_sab, aes(x = mdsA, y = mdsB, color = treatment)) +
  geom_point(size = 2, alpha = 0.7) +  # Sites as points
  stat_ellipse(level = 0.68, linewidth = 1) +  # Treatment ellipses
  scale_color_manual(values = c("firebrick", "forestgreen", "blue"),
                     labels = c("Burn", "Standing", "Harvest"),
                     name = "Treatment") +
  theme_minimal() +
  labs(title = "NMDS of Southern Rockies Bird Communities by Disturbance Type",
       x = "NMDS Axis 1", y = "NMDS Axis 2")

# Save
ggsave("Output/5-mpb-sab-nmds-sites.jpeg", width = 8, height = 6)


# Load required package
library(magick)

# Read in your saved NMDS plots
img_jnp <- image_read("Output/5-mpb-jnp-nmds-sites.jpeg")
img_sab <- image_read("Output/5-mpb-sab-nmds-sites.jpeg")

# Combine them side-by-side
combined_nmds <- image_append(c(img_jnp, img_sab))

# Save the combined image
image_write(combined_nmds, path = "Output/5-mpb-nmds_combined.jpeg", format = "jpeg")

