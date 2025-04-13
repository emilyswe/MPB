#Updated April 2 2025 to use GLM with poisson and log-link
# ALPHA RICHNESS ANALYSIS - Poisson GLM update
# Created April 2, 2025 by Emily Swerdfager (based on original Feb 19 script)

# 1. Load Required Packages ----
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggeffects)

# 2. Read in datasets for Jasper & Southern Alberta ----
jnp_data <- read.csv("Output/05-mpb-jnp-dataset.csv")
sab_data <- read.csv("Output/05-mpb-sab-dataset.csv")

# Ensure treatment is a factor with correct reference levels
jnp_data$treatment <- factor(jnp_data$treatment, levels = c("standing", "burn", "harvest"))
sab_data$treatment <- factor(sab_data$treatment, levels = c("standing", "burn", "harvest"))

# 3. Calculate Alpha Richness ----
species_columns_jnp <- colnames(jnp_data)[grep("^[A-Z]{4}$", colnames(jnp_data))]
species_columns_sab <- colnames(sab_data)[grep("^[A-Z]{4}$", colnames(sab_data))]

jnp_data$alpharich <- rowSums(jnp_data[, species_columns_jnp] > 0)
sab_data$alpharich <- rowSums(sab_data[, species_columns_sab] > 0)

# 4. Fit Poisson GLMs ----
model_jnp <- glm(alpharich ~ treatment + TSD, data = jnp_data, family = poisson(link = "log"))
summary(model_jnp)
saveRDS(model_jnp, "Output/8_glm_alpharich_jnp.rds")

model_sab <- glm(alpharich ~ treatment + TSD, data = sab_data, family = poisson(link = "log"))
summary(model_sab)
saveRDS(model_sab, "Output/8_glm_alpharich_sab.rds")

# 5. Generate Predictions ----
preds_jnp <- ggpredict(model_jnp, terms = c("TSD", "treatment"))
preds_sab <- ggpredict(model_sab, terms = c("TSD", "treatment"))

# 6. Plot Predictions ----
plot_preds_jnp <- ggplot(preds_jnp, aes(x = x, y = predicted, color = group, fill = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  scale_color_manual(values = c("forestgreen", "red", "blue")) +
  scale_fill_manual(values = c("forestgreen", "red", "blue")) +
  labs(title = "Predicted Alpha Richness (Jasper)", 
       x = "Time Since Disturbance", y = "Alpha Richness") +
  theme_classic()
ggsave("Output/8_predicted_alpharich_CI_jnp.png", plot = plot_preds_jnp, width = 8, height = 6)

plot_preds_sab <- ggplot(preds_sab, aes(x = x, y = predicted, color = group, fill = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  scale_color_manual(values = c("forestgreen", "red", "blue")) +
  scale_fill_manual(values = c("forestgreen", "red", "blue")) +
  labs(title = "Predicted Alpha Richness (Southern Rockies", 
       x = "Time Since Disturbance", y = "Alpha Richness") +
  theme_classic()
ggsave("Output/8_predicted_alpharich_CI_sab.png", plot = plot_preds_sab, width = 8, height = 6)

# 7. Save Predictions ----
write.csv(preds_jnp, "Output/8_predicted_alpha_richness_jnp.csv", row.names = FALSE)
write.csv(preds_sab, "Output/8_predicted_alpha_richness_sab.csv", row.names = FALSE)



#######################
#######################
#######################
#######################
#######################
#######################
#######################
#######################
# Load libraries
library(ggplot2)
library(ggeffects)
library(patchwork)

# Read predicted data
preds_jnp <- read.csv("Output/8_predicted_alpha_richness_jnp.csv")
preds_sab <- read.csv("Output/8_predicted_alpha_richness_sab.csv")

# Ensure consistent factor levels for group
preds_jnp$group <- factor(preds_jnp$group, levels = c("standing", "burn", "harvest"))
preds_sab$group <- factor(preds_sab$group, levels = c("standing", "burn", "harvest"))

# Plot for Jasper
plot_jnp <- ggplot(preds_jnp, aes(x = x, y = predicted, color = group, fill = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  scale_color_manual(values = c("forestgreen", "red", "blue")) +
  scale_fill_manual(values = c("forestgreen", "red", "blue")) +
  labs(title = "Jasper",
       x = "Time Since Disturbance (Years)", y = "Alpha Richness",
       color = "Disturbance Type", fill = "Disturbance Type") +
  theme_classic() +
  theme(legend.position = "none")

# Plot for Southern Alberta
plot_sab <- ggplot(preds_sab, aes(x = x, y = predicted, color = group, fill = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  scale_color_manual(values = c("forestgreen", "red", "blue")) +
  scale_fill_manual(values = c("forestgreen", "red", "blue")) +
  labs(title = "Southern Rockies",
       x = "Time Since Disturbance (Years)", y = "Alpha Richness",
       color = "Disturbance Type", fill = "Disturbance Type") +
  theme_classic() +
  theme(legend.position = "right")

# Combine the plots 
combined_plot <- plot_jnp + plot_sab +
  plot_annotation(
    title = "Predicted Alpha Richness by Region and Disturbance Type"
  ) &
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold")  # Increased from 16 to 18
  )

# Save the plot
ggsave("Output/8_combined_alpha_richness_panels.png", combined_plot, width = 12, height = 6)






















############################
############################
############################
############################
############################
############################
############################
############################
# 8-ALPHA RICHNESS ANALYSIS - MPB Chapter 2
# Created February 19, 2025, by Emily Swerdfager
# This script quantifies and visualizes alpha richness by treatment for both the Jasper and Southern Alberta datasets.

# 1. Load Required Packages ----
library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)
library(GGally)
library(ggeffects)

# 2. Read in datasets for Jasper & Southern Alberta ----
jnp_data <- read.csv("Output/05-mpb-jnp-dataset.csv")
sab_data <- read.csv("Output/05-mpb-sab-dataset.csv")

# Ensure treatment is a factor with correct reference levels
jnp_data$treatment <- factor(jnp_data$treatment, levels = c("standing", "burn", "harvest"))
sab_data$treatment <- factor(sab_data$treatment, levels = c("standing", "burn", "harvest"))

# 3. Calculate Alpha Richness ----
# Identify species columns using regex (assuming 4-letter bird codes)
species_columns_jnp <- colnames(jnp_data)[grep("^[A-Z]{4}$", colnames(jnp_data))]
species_columns_sab <- colnames(sab_data)[grep("^[A-Z]{4}$", colnames(sab_data))]

# Calculate alpha richness for each dataset
jnp_data$alpharich <- rowSums(jnp_data[, species_columns_jnp] > 0)
sab_data$alpharich <- rowSums(sab_data[, species_columns_sab] > 0)

# 4. Exploratory Plots ----
# Scatter plots: Alpha richness vs. TSD by treatment
plot_jnp <- ggplot(jnp_data, aes(x = TSD, y = alpharich, color = treatment)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm") +
  theme_bw() +
  scale_color_manual(values = c("forestgreen", "red", "blue")) +
  labs(title = "Alpha Richness vs. TSD (Jasper)", y = "Alpha Richness", x = "Time Since Disturbance (TSD)")
ggsave("Output/8_smoothed_scatterplot_alpharich_jnp.png", plot = plot_jnp, width = 10, height = 8)

plot_sab <- ggplot(sab_data, aes(x = TSD, y = alpharich, color = treatment)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm") +
  theme_bw() +
  scale_color_manual(values = c("forestgreen", "red", "blue")) +
  labs(title = "Alpha Richness vs. TSD (Southern Alberta)", y = "Alpha Richness", x = "Time Since Disturbance (TSD)")
ggsave("Output/8_smoothed_scatterplot_alpharich_sab.png", plot = plot_sab, width = 10, height = 8)

# 5. Linear Model ----
# Model the effect of treatment and TSD on alpha richness
model_jnp <- lm(alpharich ~ treatment + TSD, data = jnp_data, na.action = "na.fail")
summary(model_jnp)
saveRDS(model_jnp, "Output/8_lm_alpharich_jnp.rds")

model_sab <- lm(alpharich ~ treatment + TSD, data = sab_data, na.action = "na.fail")
summary(model_sab)
saveRDS(model_sab, "Output/8_lm_alpharich_sab.rds")

# 6. Predict & Visualize Model Effects with Confidence Intervals ----
# Generate predictions with confidence intervals
preds_jnp <- ggpredict(model_jnp, terms = c("TSD", "treatment"))
preds_sab <- ggpredict(model_sab, terms = c("TSD", "treatment"))

# Plot Jasper predictions with 95% CIs
plot_preds_jnp <- ggplot(preds_jnp, aes(x = x, y = predicted, color = group, fill = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  scale_color_manual(values = c("forestgreen", "red", "blue")) +  
  scale_fill_manual(values = c("forestgreen", "red", "blue")) +  
  labs(title = "Predicted Alpha Richness (Jasper) with 95% CIs", 
       x = "Time Since Disturbance", y = "Alpha Richness") +
  theme_classic()
ggsave("Output/8_predicted_alpharich_CI_jnp.png", plot = plot_preds_jnp, width = 8, height = 6)

# Plot Southern Alberta predictions with 95% CIs
plot_preds_sab <- ggplot(preds_sab, aes(x = x, y = predicted, color = group, fill = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  scale_color_manual(values = c("forestgreen", "red", "blue")) +  
  scale_fill_manual(values = c("forestgreen", "red", "blue")) +  
  labs(title = "Predicted Alpha Richness (Southern Alberta) with 95% CIs", 
       x = "Time Since Disturbance", y = "Alpha Richness") +
  theme_classic()
ggsave("Output/8_predicted_alpharich_CI_sab.png", plot = plot_preds_sab, width = 8, height = 6)

# 7. Save Model Predictions ----
write.csv(preds_jnp, "Output/8_predicted_alpha_richness_jnp.csv", row.names = FALSE)
write.csv(preds_sab, "Output/8_predicted_alpha_richness_sab.csv", row.names = FALSE)

