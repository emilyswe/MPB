# 1. Load Required Packages ----
library(dplyr)
library(ggplot2)
library(tidyr)
library(MASS)  # For Negative Binomial CI

# 2. Read in Datasets ----
jnp_data <- read.csv("Output/05-mpb-jnp-dataset.csv")
sab_data <- read.csv("Output/05-mpb-sab-dataset.csv")

# Ensure treatment is a factor with correct reference levels
jnp_data$treatment <- factor(jnp_data$treatment, levels = c("standing", "burn", "harvest"))
sab_data$treatment <- factor(sab_data$treatment, levels = c("standing", "burn", "harvest"))

# 3. Identify Species Columns ----
species_columns_jnp <- colnames(jnp_data)[grep("^[A-Z]{4}$", colnames(jnp_data))] 
species_columns_sab <- colnames(sab_data)[grep("^[A-Z]{4}$", colnames(sab_data))] 

# 4. Gamma Richness Rarefaction Setup ----
n_replicates <- 100  # Number of iterations for rarefaction
n_sampling_points <- seq(10, 100, by = 10)  # Range of sample sizes

# Function to calculate gamma richness
calculate_gamma_richness <- function(data, n) {
  if (nrow(data) < n) n <- nrow(data)  # Adjust n if too large
  sample_data <- data[sample(nrow(data), n, replace = TRUE), species_columns_jnp, drop = FALSE]
  richness_value <- sum(colSums(sample_data > 0, na.rm = TRUE) > 0)
  return(richness_value)
}

# 5. Rarefaction Function ----
run_rarefaction <- function(data) {
  gamma_results <- data.frame(treatment = character(), replicate = integer(), numpc = integer(), richness = integer())
  
  for (treat in unique(data$treatment)) {
    treat_data <- data[data$treatment == treat, ]
    for (n in n_sampling_points) {
      for (rep in 1:n_replicates) {
        if (nrow(treat_data) > 0) {
          richness_value <- calculate_gamma_richness(treat_data, n)
          gamma_results <- rbind(gamma_results, data.frame(treatment = treat, replicate = rep, numpc = n, richness = richness_value))
        }
      }
    }
  }
  return(gamma_results)
}

# 6. Compute Rarefaction Results ----
gamma_jnp <- run_rarefaction(jnp_data)
gamma_sab <- run_rarefaction(sab_data)

# 7. Compute Summary Statistics ----
compute_summary <- function(gamma_data) {
  gamma_summary <- gamma_data %>%
    group_by(treatment, numpc) %>%
    summarise(
      mean_richness = mean(richness, na.rm = TRUE),
      sd_richness = sd(richness, na.rm = TRUE),
      
      # Normal CIs
      l95_norm = mean_richness - 1.96 * (sd_richness / sqrt(n_replicates)),
      u95_norm = mean_richness + 1.96 * (sd_richness / sqrt(n_replicates)),
      
      # Poisson CIs
      l95_pois = qpois(0.025, mean_richness),
      u95_pois = qpois(0.975, mean_richness),
      
      # Negative Binomial CIs
      l95_nb = mean_richness * exp(-1.96 * (1 / sqrt(n_replicates))),
      u95_nb = mean_richness * exp(1.96 * (1 / sqrt(n_replicates))),
      
      .groups = 'drop'
    )
  return(gamma_summary)
}

gamma_jnp_summary <- compute_summary(gamma_jnp)
gamma_sab_summary <- compute_summary(gamma_sab)

# Save outputs
write.csv(gamma_jnp_summary, "Output/9-jnp_gamma_richness.csv", row.names = FALSE)
write.csv(gamma_sab_summary, "Output/9-sab_gamma_richness.csv", row.names = FALSE)

# 8. Plot Function ----
plot_gamma_richness <- function(gamma_summary, title, filename, ci_type) {
  ci_lower <- paste0("l95_", ci_type)
  ci_upper <- paste0("u95_", ci_type)
  
  # Define correct color mapping (ensures correct colors for treatments)
  color_mapping <- c("standing" = "forestgreen", "burn" = "red", "harvest" = "blue")
  fill_mapping <- c("standing" = "forestgreen", "burn" = "red", "harvest" = "blue")
  
  ggplot(gamma_summary, aes(x = numpc, y = mean_richness, color = treatment, fill = treatment)) +
    geom_ribbon(aes_string(ymin = ci_lower, ymax = ci_upper), alpha = 0.2) +
    geom_line(size = 1) +
    scale_color_manual(values = color_mapping, name = "Treatment") +  # Explicitly set correct colors
    scale_fill_manual(values = fill_mapping, name = "Treatment") +  # Explicitly set correct fills
    labs(title = title, x = "Number of Surveys", y = "Cumulative Number of Species") +
    theme_bw() +
    theme(legend.position = "bottom")
  
  ggsave(filename, width = 10, height = 5)
}

# 9. Generate Plots ----
# Normal CIs
plot_gamma_richness(gamma_jnp_summary, "Gamma Richness Rarefaction - Jasper (Normal CIs)", "Output/9-gamma_richness_jnp_norm.png", ci_type = "norm")
plot_gamma_richness(gamma_sab_summary, "Gamma Richness Rarefaction - Southern Alberta (Normal CIs)", "Output/9-gamma_richness_sab_norm.png", ci_type = "norm")

# Poisson CIs
plot_gamma_richness(gamma_jnp_summary, "Gamma Richness Rarefaction - Jasper (Poisson CIs)", "Output/9-gamma_richness_jnp_pois.png", ci_type = "pois")
plot_gamma_richness(gamma_sab_summary, "Gamma Richness Rarefaction - Southern Alberta (Poisson CIs)", "Output/9-gamma_richness_sab_pois.png", ci_type = "pois")

# Negative Binomial CIs
plot_gamma_richness(gamma_jnp_summary, "Gamma Richness Rarefaction - Jasper (Negative Binomial CIs)", "Output/9-gamma_richness_jnp_nb.png", ci_type = "nb")
plot_gamma_richness(gamma_sab_summary, "Gamma Richness Rarefaction - Southern Alberta (Negative Binomial CIs)", "Output/9-gamma_richness_sab_nb.png", ci_type = "nb")

#############
#############
############
####
#updated april 2 2025 to create combined plot for thesis doc
# Load libraries
library(ggplot2)
library(patchwork)
library(dplyr)
library(readr)


# Read in summary data
gamma_jnp_summary <- read.csv("Output/9-jnp_gamma_richness.csv")
gamma_sab_summary <- read.csv("Output/9-sab_gamma_richness.csv")

# Set consistent factor levels
gamma_jnp_summary$treatment <- factor(gamma_jnp_summary$treatment, levels = c("standing", "burn", "harvest"))
gamma_sab_summary$treatment <- factor(gamma_sab_summary$treatment, levels = c("standing", "burn", "harvest"))

# Color mapping
color_mapping <- c("standing" = "forestgreen", "burn" = "red", "harvest" = "blue")

# Jasper plot
plot_jnp <- ggplot(gamma_jnp_summary, aes(x = numpc, y = mean_richness, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = l95_norm, ymax = u95_norm), alpha = 0.2) +
  geom_line(size = 1) +
  scale_color_manual(values = color_mapping) +
  scale_fill_manual(values = color_mapping) +
  labs(title = "Jasper",
       x = "Number of Surveys", y = "Cumulative Number of Species") +
  coord_cartesian(ylim = c(20, 60)) +
  theme_classic() +
  theme(legend.position = "none")

# Southern Rockies plot
plot_sab <- ggplot(gamma_sab_summary, aes(x = numpc, y = mean_richness, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = l95_norm, ymax = u95_norm), alpha = 0.2) +
  geom_line(size = 1) +
  scale_color_manual(values = color_mapping) +
  scale_fill_manual(values = color_mapping) +
  labs(title = "Southern Rockies",
       x = "Number of Surveys", y = "Cumulative Number of Species") +
  coord_cartesian(ylim = c(20, 60)) +
  theme_classic() +
  theme(legend.position = "right")

# Combine with shared title
combined_gamma <- plot_jnp + plot_sab +
  plot_annotation(title = "Gamma Richness by Region and Disturbance Type") &
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

# Save
ggsave("Output/9_combined_gamma_richness_panels.png", combined_gamma, width = 12, height = 6)
