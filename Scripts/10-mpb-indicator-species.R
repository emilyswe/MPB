# Load required packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(indicspecies)

# Read in dataset
bird_data <- read.csv("Output/05-mpb-jnp-dataset.csv")

# Define species columns (they follow the 4-letter convention)
species_columns <- grep("^[A-Z]{4}$", colnames(bird_data), value = TRUE)

# Subset data for indicator species analysis
indval_data <- bird_data %>%
  select(location, treatment, all_of(species_columns))

# Convert species data to matrix
species_matrix <- as.matrix(indval_data[, species_columns])

# Convert treatment to factor
indval_data$treatment <- as.factor(indval_data$treatment)

# Run Indicator Species Analysis (IndVal) with 999 permutations
set.seed(123)
indval_result <- multipatt(species_matrix, indval_data$treatment, func = "IndVal.g", control = how(nperm = 999))

# Convert IndVal results into a dataframe
indval_df <- indval_result$sign
indval_df$species <- rownames(indval_result$sign)  # Add species names

# Remove rows with NA p-values and filter for significant species (p < 0.05)
indval_df <- indval_df %>%
  filter(!is.na(p.value) & p.value < 0.05)

# Identify the strongest group association for each species
numeric_columns <- c("s.burn", "s.harvest", "s.standing")  # Columns containing indicator values
group_assignments <- apply(indval_df[numeric_columns], 1, function(x) names(x)[which.max(x)]) 

# Add group assignments to the dataframe
indval_df$group <- group_assignments

# Select final columns
significant_species <- indval_df %>%
  select(species, stat, p.value, group)

# Save results to CSV
write.csv(significant_species, "Output/10-indicator_species_results.csv", row.names = FALSE)

####### **🚀 Visualization: Indicator Species Bar Plot**
indicator_plot <- ggplot(significant_species, aes(x = reorder(species, stat), y = stat, fill = group)) +
  geom_bar(stat = "identity", show.legend = TRUE) +
  coord_flip() +
  scale_fill_manual(values = c("firebrick", "forestgreen", "blue"), 
                    labels = c("Burn", "Standing", "Harvest"), name = "Treatment") +
  labs(title = "Significant Indicator Species by Treatment",
       x = "Species",
       y = "Indicator Value (IndVal)") +
  theme_minimal()

# Save plot
ggsave("Output/10-indicator_species_plot.jpeg", plot = indicator_plot, width = 8, height = 6)

#################
#################
################
#################
####### indicator species for sab (updated april 2 2025)

# Load required packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(indicspecies)
library(patchwork)

# ==== 1. Load Southern Rockies Data ====
sab_data <- read.csv("Output/05-mpb-sab-dataset.csv")

# Identify species columns (assuming 4-letter codes)
species_columns <- grep("^[A-Z]{4}$", colnames(sab_data), value = TRUE)

# Prepare input for IndVal
indval_data_sab <- sab_data %>%
  select(location, treatment, all_of(species_columns))

# Convert to matrix
species_matrix_sab <- as.matrix(indval_data_sab[, species_columns])
indval_data_sab$treatment <- as.factor(indval_data_sab$treatment)

# ==== 2. Run Indicator Species Analysis ====
set.seed(123)
indval_result_sab <- multipatt(
  species_matrix_sab,
  indval_data_sab$treatment,
  func = "IndVal.g",
  control = how(nperm = 999)
)

# Extract and filter significant species (p < 0.05)
indval_df_sab <- indval_result_sab$sign
indval_df_sab$species <- rownames(indval_df_sab)
indval_df_sab <- indval_df_sab %>%
  filter(!is.na(p.value) & p.value < 0.05)

# Identify the group with strongest association
group_assignments_sab <- apply(indval_df_sab[, c("s.burn", "s.harvest", "s.standing")],
                               1, function(x) names(x)[which.max(x)])
indval_df_sab$group <- group_assignments_sab

# Final dataframe
significant_species_sab <- indval_df_sab %>%
  select(species, stat, p.value, group)

# Save CSV
write.csv(significant_species_sab, "Output/10-indicator_species_results_sab.csv", row.names = FALSE)

# ==== 3. Create SAB Plot ====
significant_species_sab <- read.csv("Output/10-indicator_species_results_sab.csv")

plot_sab <- ggplot(significant_species_sab, aes(x = reorder(species, stat), y = stat, fill = group)) +
  geom_bar(stat = "identity", show.legend = TRUE) +
  coord_flip() +
  scale_fill_manual(values = c("firebrick", "forestgreen", "blue"),
                    labels = c("Burn", "Standing", "Harvest"),
                    name = "Disturbance Type") +
  labs(title = "Southern Rockies",
       x = "Species",
       y = "Indicator Value (IndVal)") +
  theme_bw(base_size = 12)

ggsave("Output/10-indicator_species_plot_sab.png", plot = plot_sab, width = 8, height = 6)


# ==== 4. Load Jasper Plot ====
significant_species_jnp <- read.csv("Output/10-indicator_species_results.csv")

plot_jnp <- ggplot(significant_species_jnp, aes(x = reorder(species, stat), y = stat, fill = group)) +
  geom_bar(stat = "identity", show.legend = TRUE) +
  coord_flip() +
  scale_fill_manual(values = c("firebrick", "forestgreen", "blue"), 
                    labels = c("Burn", "Standing", "Harvest"), name = "Disturbance Type") +
  labs(title = "Jasper",
       x = "Species",
       y = "Indicator Value (IndVal)") +
  theme_minimal()

# Save plot
ggsave("Output/10-indicator_species_plot.jpeg", plot = plot_jnp, width = 8, height = 6)



# ==== 5. Combine Plots ====
combined_plot <- plot_jnp + plot_sab +
  plot_annotation(title = "Indicator Species by Region and Disturbance Type") &
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

# Save combined plot
ggsave("Output/10_combined_indicator_species_panels2.png", combined_plot, width = 12, height = 6)

