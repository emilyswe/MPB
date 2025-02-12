#setwd
setwd("~/Documents/local-git/MPB")

# Load required libraries
library(sf)
library(dplyr)
library(tidyr)
library(terra)
library(raster)
library(fasterize)

#load jnp disturbance data
disturbances <- read.csv("/Users/Bronwyn/Documents/local-git/MPB/Input/jnp-disturbances.csv")
str(disturbances)

################################
################################
######Step #1. calculate disturbance year and treatment type for jasper sites

# Initialize the new columns in the dataframe to avoid errors
disturbances$disturbance_year <- NA  # Initialize the disturbance_year column with NAs
disturbances$treatment <- NA  # Initialize the treatment column with NAs

# Loop through each site in the dataframe
for (i in 1:nrow(disturbances)) {
  pointData <- disturbances[i, ]  # Get the site data
  
  # Get the disturbance years for all the relevant columns
  disturbance_years <- c(pointData$mpbyear, pointData$maxFY150m, pointData$maxFY300m, 
                         pointData$maxFY500m, pointData$maxHY150m, pointData$maxHY300m, 
                         pointData$maxHY500m)
  
  # Get the names of the disturbance columns
  disturbance_column_names <- c("mpbyear", "maxFY150m", "maxFY300m", 
                                "maxFY500m", "maxHY150m", "maxHY300m", "maxHY500m")
  
  # Remove NA values and match valid years with their column names
  valid_years <- disturbance_years[!is.na(disturbance_years)]  # Remove NAs before getting the max value
  valid_columns <- disturbance_column_names[!is.na(disturbance_years)]
  
  # If there's at least one valid year, assign the most recent year to disturbance_year
  if (length(valid_years) > 0) {
    # Get the most recent disturbance year
    max_year <- max(valid_years)
    disturbances$disturbance_year[i] <- max_year
    
    # Determine the treatment based on the column associated with the max year
    max_column <- valid_columns[which.max(valid_years)]
    
    if (max_column == "mpbyear") {
      disturbances$treatment[i] <- "standing"
    } else if (max_column %in% c("maxFY150m", "maxFY300m", "maxFY500m")) {
      disturbances$treatment[i] <- "burn"
    } else if (max_column %in% c("maxHY150m", "maxHY300m", "maxHY500m")) {
      disturbances$treatment[i] <- "harvest"
    }
  } else {
    disturbances$disturbance_year[i] <- NA  # If no valid year, set as NA
    disturbances$treatment[i] <- NA  # No treatment if no disturbance year is found
  }
}

# Check the updated dataframe
head(disturbances)

# Optional: Save the updated dataframe to CSV
write.csv(disturbances, 
          file = "/Users/Bronwyn/Documents/local-git/MPB/Output/jnp-disturbances-tsd.csv", 
          row.names = FALSE)



###############################################
###############################################
###############################################
###############################################
###############################################
###############################################

##########Step #2. Calculate TSD for jnp sites


# Step 1: Calculate TSD (Time Since Disturbance)
disturbances <- disturbances %>%
  mutate(TSD = ifelse(!is.na(disturbance_year), visityear - disturbance_year, NA))  # Calculate TSD where disturbance_year is valid

# Step 2: View the updated dataframe (optional)
View(disturbances)

# Step 3: Save the updated dataframe to a new CSV file
write.csv(disturbances, 
          file = "/Users/Bronwyn/Documents/local-git/MPB/Output/3-jnp-disturbances-tsd.csv", 
          row.names = FALSE)

###############################################
###############################################
###############################################
###############################################
###############################################
###############################################
##########Step #3. visualize the data

######make tsd by treatment plot
library(ggplot2)

# Adjust for discrete TSD values using geom_bar()
tsd_treatment_plot <- ggplot(disturbances, aes(x = factor(TSD), fill = treatment)) +
  geom_bar(
    position = "dodge",
    color = "black",
    alpha = 0.7
  ) +
  scale_fill_manual(values = c("burn" = "red", "harvest" = "blue", "standing" = "forestgreen")) +
  labs(
    title = "Distribution of TSD by Treatment (Jasper)",
    x = "Time Since Disturbance (years)",
    y = "# of Sites",
    fill = "Treatment"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

# Save the adjusted plot
ggsave(
  filename = "/Users/Bronwyn/Documents/local-git/MPB/Output/3-JNP-TSD_distribution_by_treatment.png",
  plot = tsd_treatment_plot,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)


###############################################
###############################################
# Histogram of sites by treatment

# Create the bar plot for distribution of sites by treatment
site_treatment_plot <- ggplot(disturbances, aes(x = treatment, fill = treatment)) +
  geom_bar(color = "black", alpha = 0.7) +
  scale_fill_manual(values = c("burn" = "red", "harvest" = "blue", "standing" = "forestgreen")) +
  labs(title = "Distribution of Sites by Treatment (Jasper)",
       x = "Treatment",
       y = "# of Sites",
       fill = "Treatment") +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

# Save the plot as a JPEG
ggsave("/Users/Bronwyn/Documents/local-git/MPB/Output/JNP-Sites_by_treatment.jpeg", 
       plot = site_treatment_plot, 
       width = 6, height = 4, dpi = 300)



###############################################
###############################################
###############################################
###############################################
###############################################
###############################################
########## plot severity and attack stages for jnp sites 
library(ggplot2)
library(patchwork)
library(dplyr)

# Step 1: Load the CSV file
disturbances <- read.csv("/Users/Bronwyn/Documents/local-git/MPB/Output/3-jnp-disturbances-tsd.csv")

# Step 2: Replace values of 0 with 1 in SEVERITY500 and DOMCLASS500
disturbances <- disturbances %>%
  mutate(
    SEVERITY500 = ifelse(SEVERITY500 == 0, 1, SEVERITY500),
    DOMCLASS500 = ifelse(DOMCLASS500 == 0, 1, DOMCLASS500)
  )

# Step 3: Reclassify SEVERITY500 codes (3 and 4 as "Severe")
disturbances <- disturbances %>%
  mutate(
    SEVERITY500 = case_when(
      SEVERITY500 == 1 ~ "Moderate",
      SEVERITY500 == 2 ~ "High",
      SEVERITY500 %in% c(3, 4) ~ "Severe",
      TRUE ~ NA_character_
    )
  )

# Step 4: Replace DOMCLASS500 values with descriptive labels
disturbances <- disturbances %>%
  mutate(
    DOMCLASS500 = case_when(
      DOMCLASS500 == 1 ~ "NA",
      DOMCLASS500 == 2 ~ "Red",
      DOMCLASS500 == 3 ~ "Grey",
      TRUE ~ NA_character_
    )
  )

# Step 5: Prepare data for plotting
severity_data <- data.frame(severity = na.omit(disturbances$SEVERITY500))
attack_code_data <- data.frame(attack_stage = na.omit(disturbances$DOMCLASS500))

# Step 6: Create histograms
# Histogram for SEVERITY500
p_severity <- ggplot(severity_data, aes(x = severity)) +
  geom_bar(fill = "forestgreen", color = "black", alpha = 0.7) +
  labs(
    title = "Distribution of Severity (Jasper)",
    x = "Severity",
    y = "# of Sites"
  ) +
  scale_x_discrete(limits = c("Moderate", "High", "Severe")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), # Title formatting
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

# Histogram for DOMCLASS500 (Attack Stage)
p_attack_code <- ggplot(attack_code_data, aes(x = attack_stage)) +
  geom_bar(fill = "blue", color = "black", alpha = 0.7) +
  labs(
    title = "Distribution of Attack Stage (Jasper)",
    x = "Attack Stage",
    y = "# of Sites"
  ) +
  scale_x_discrete(limits = c("NA", "Red", "Grey")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), # Title formatting
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

# Step 7: Combine the histograms into one image
combined_plot <- p_severity / p_attack_code

# Step 8: Save the combined plot
ggsave(
  "/Users/Bronwyn/Documents/local-git/MPB/Output/JNP-severity-attackcodes-histograms.jpeg",
  plot = combined_plot,
  width = 8,
  height = 10,
  dpi = 300
)



########################
########################
########################
######################## patch size calculation, code seems to be missing but found this in chat GPT convo, title MPB Covariates Jasper Plots need to check that it works 
#added Feb 11 2025

# Assign patch size based on treatment and dynamic buffer
disturbances <- disturbances %>%
  mutate(
    RelevantArea = case_when(
      treatment == "burn" & !is.na(BURNAREA150) ~ BURNAREA150,
      treatment == "burn" & !is.na(BURNAREA300) ~ BURNAREA300,
      treatment == "burn" & !is.na(BURNAREA500) ~ BURNAREA500,
      treatment == "harvest" & !is.na(HARVESTAREA150) ~ HARVESTAREA150,
      treatment == "harvest" & !is.na(HARVESTAREA300) ~ HARVESTAREA300,
      treatment == "harvest" & !is.na(HARVESTAREA500) ~ HARVESTAREA500,
      treatment == "standing" & !is.na(BEETLEAREA150) ~ BEETLEAREA150,
      treatment == "standing" & !is.na(BEETLEAREA300) ~ BEETLEAREA300,
      treatment == "standing" & !is.na(BEETLEAREA500) ~ BEETLEAREA500,
      TRUE ~ NA_real_
    )
  )

# Check for missing values and troubleshoot if needed
print(disturbances %>% filter(is.na(RelevantArea)))

# Create the histogram
library(ggplot2)

patch_size_plot <- ggplot(data = disturbances, aes(x = RelevantArea, fill = treatment)) +
  geom_histogram(binwidth = 10000, position = "stack", color = "black", alpha = 0.7) +
  scale_fill_manual(values = c("burn" = "red", "harvest" = "blue", "standing" = "forestgreen")) +
  labs(
    title = "Patch Size Distribution by Treatment (Dynamic Buffer)",
    x = "Patch Area (square meters)",
    y = "# of Sites",
    fill = "Treatment Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

# Save the plot
ggsave(
  filename = "/Users/Bronwyn/Documents/local-git/MPB/Output/patch_size_distribution_dynamic_buffer.jpeg",
  plot = patch_size_plot,
  width = 10,
  height = 8,
  dpi = 300
)

