#########Created February 12 2025 by Emily Swerdfager
####PREAMBLE: The purpose of this script is to perform a preliminary NMDS on the MPB bird data for jasper

# Load required packages
packs <- c("dplyr", "vegan", "ggplot2")
for (q in 1:length(packs)) {
  if (!require(packs[q], character.only = TRUE)) {
    install.packages(packs[q])
    require(packs[q])
  }
}

install.packages("vegan", dependencies = TRUE)
library(vegan)


# Read in dataset
bird_data <- read.csv("Output/05_MPB_bird_covariates_JNP.csv")

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
  scale_color_manual(values = c("firebrick", "darkgreen", "goldenrod"), 
                     labels = c("Burn", "Standing", "Harvest"), name = "Treatment") +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(x = "NMDS Axis 1", y = "NMDS Axis 2", title = "NMDS of Bird Communities by MPB Treatment")

# Save plot
ggsave(plot3, filename = "Output/5-NMDS_MPB.jpeg", width = 8, height = 6)

####### **PERMANOVA Tests**

# Prepare environmental variables for PERMANOVA
env_data <- nmds_data %>%
  select(treatment, TSD)  # Select only relevant covariates

# Run PERMANOVA for treatment
set.seed(999)
permanova_treatment <- adonis2(
  species_matrix ~ treatment, 
  data = env_data, 
  method = "bray", 
  permutations = 999
)

# Run PERMANOVA for TSD
set.seed(999)
permanova_tsd <- adonis2(
  species_matrix ~ TSD, 
  data = env_data, 
  method = "bray", 
  permutations = 999
)

# Display results
print("PERMANOVA Results for Treatment:")
print(permanova_treatment)

print("PERMANOVA Results for Time Since Disturbance (TSD):")
print(permanova_tsd)

# Save PERMANOVA results
write.csv(as.data.frame(permanova_treatment), "Output/5-permanova_treatment.csv", row.names = TRUE)
write.csv(as.data.frame(permanova_tsd), "Output/5-permanova_tsd.csv", row.names = TRUE)
