# Load necessary packages
library(ggplot2)

###############################
###############################
###############################
###############################
###############################
#NMDS REGRESSION ON JASPER DATASET:

# Read in the NMDS output data
nmds_data <- read.csv("Output/6-mpb_nmds_2axis.csv")

# Check structure of the data
str(nmds_data)

# Load required libraries
library(ggplot2)

# Ensure treatment is a factor and set reference level to "standing"
nmds_data$treatment <- factor(nmds_data$treatment, levels = c("standing", "burn", "harvest"))

# Run linear regressions for NMDS axes
lm_mdsA <- lm(mdsA ~ treatment + TSD, data = nmds_data)
lm_mdsB <- lm(mdsB ~ treatment + TSD, data = nmds_data)

# Print model summaries
summary(lm_mdsA)
summary(lm_mdsB)

# Define custom colors
custom_colors <- c("burn" = "red", "harvest" = "blue", "standing" = "forestgreen")

# Create and save NMDS1 (mdsA) vs. TSD plot
plot_mdsA <- ggplot(nmds_data, aes(x = TSD, y = mdsA, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_color_manual(values = custom_colors) +
  theme_bw() +
  labs(title = "Regression of NMDS1 (mdsA) vs. TSD by Treatment")

ggsave("Output/NMDS1_vs_TSD.png", plot = plot_mdsA, width = 7, height = 5, dpi = 300)

# Create and save NMDS2 (mdsB) vs. TSD plot
plot_mdsB <- ggplot(nmds_data, aes(x = TSD, y = mdsB, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_color_manual(values = custom_colors) +
  theme_bw() +
  labs(title = "Regression of NMDS2 (mdsB) vs. TSD by Treatment")

ggsave("Output/NMDS2_vs_TSD.png", plot = plot_mdsB, width = 7, height = 5, dpi = 300)




###############################
###############################
###############################
###############################
###############################
#NMDS REGRESSION ON SOUTHERN ALBRERTA DATASET:


# Load necessary packages
library(ggplot2)

# Read in the NMDS output data for Southern Alberta
nmds_data_sab <- read.csv("Output/6-mpb_sab_nmds_2axis.csv")

# Check structure of the data
str(nmds_data_sab)

# Ensure treatment is a factor and set reference level to "standing"
nmds_data_sab$treatment <- factor(nmds_data_sab$treatment, levels = c("standing", "burn", "harvest"))

# Run linear regressions for NMDS axes
lm_mdsA_sab <- lm(mdsA ~ treatment + TSD, data = nmds_data_sab)
lm_mdsB_sab <- lm(mdsB ~ treatment + TSD, data = nmds_data_sab)

# Print model summaries
summary(lm_mdsA_sab)
summary(lm_mdsB_sab)

# Define custom colors
custom_colors <- c("burn" = "red", "harvest" = "blue", "standing" = "forestgreen")

# Create and save NMDS1 (mdsA) vs. TSD plot
plot_mdsA_sab <- ggplot(nmds_data_sab, aes(x = TSD, y = mdsA, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_color_manual(values = custom_colors) +
  theme_bw() +
  labs(title = "Regression of NMDS1 (mdsA) vs. TSD by Treatment (Southern Alberta)")

ggsave("Output/NMDS1_vs_TSD_SAB.png", plot = plot_mdsA_sab, width = 7, height = 5, dpi = 300)

# Create and save NMDS2 (mdsB) vs. TSD plot
plot_mdsB_sab <- ggplot(nmds_data_sab, aes(x = TSD, y = mdsB, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_color_manual(values = custom_colors) +
  theme_bw() +
  labs(title = "Regression of NMDS2 (mdsB) vs. TSD by Treatment (Southern Alberta)")

ggsave("Output/NMDS2_vs_TSD_SAB.png", plot = plot_mdsB_sab, width = 7, height = 5, dpi = 300)


