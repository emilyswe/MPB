###Created February 11 2025 by Emily Swerdfager
###PREAMBLE: The purpose of this script is to download bird data for the MPB project and reformat it for analysis

# Install remotes if you haven't already
install.packages("remotes")

# Install wildRtrax from GitHub
remotes::install_github("ABbiodiversity/wildrtrax")

# Load libraries
library(tidyverse)
library(wildrtrax)
library(lubridate)

# Set environment variables for your WildTrax credentials
Sys.setenv(WT_USERNAME = 'emilyswe@ualberta.ca')
Sys.setenv(WT_PASSWORD = 'xeSba3-roqrob-turjoz')

# Authenticate into WildTrax (you've already done this)
wt_auth()

###############################################
###############################################
###############################################

# 1. Download Bird Data for MPB Project (ID: 2129)
mpb_raw <- wt_download_report(2129, "ARU", "main", FALSE)

# 2. Remove Non-Bird Species
mpb_filtered <- wt_tidy_species(mpb_raw, remove = c("mammal", "amphibian", "abiotic", "insect", "unknown"))

# 3. Replace TMTT Counts Using the WildTrax Function (uses prediction table)
mpb_filtered <- wt_replace_tmtt(mpb_filtered)

# 4. Convert 'individual_count' to Numeric and Handle Confidence Intervals (CI)
mpb_filtered <- mpb_filtered %>%
  mutate(individual_count = case_when(
    individual_count == "CI 1" ~ 1,  # Confidence Interval 1 → 1
    individual_count == "CI 3" ~ 3,  # Confidence Interval 3 → 3
    TRUE ~ as.numeric(individual_count)  # Convert remaining values to numbers
  ))

# 5. Filter for Breeding Window and Morning Surveys
mpb_filtered <- mpb_filtered %>%
  mutate(
    year = year(recording_date_time),
    julian = yday(recording_date_time),
    hour = hour(recording_date_time)
  ) %>%
  filter(
    julian >= 147 & julian <= 183,  # May 27 – July 2
    year %in% c(2023, 2024),        # Only 2023 and 2024
    hour >= 4 & hour < 10           # 4 AM – 10 AM surveys
  )

# 6. Get Maximum Count Per Species Per Site
#mpb_max <- mpb_filtered %>%
 # group_by(location, species_code) %>%
  #summarise(max_count = max(individual_count, na.rm = TRUE), .groups = "drop")

# 6. Get Mean Count Per Species Per Site (Rounded Up)
mpb_mean <- mpb_filtered %>%
  group_by(location, species_code) %>%
  summarise(mean_count = ceiling(mean(individual_count, na.rm = TRUE)), .groups = "drop")

# 7. Convert to Wide Format (Site by Species Matrix)
mpb_wide <- mpb_mean %>%
  pivot_wider(names_from = species_code, values_from = mean_count, values_fill = 0)

# 8. Save Final Output
write.csv(mpb_wide, "Output/1-mpb-bird-data-mean-count-wide.csv", row.names = FALSE)

#################################
################filter for territory type A species only

birds <- read.csv("Output/1-mpb-bird-data-mean-count-wide.csv")

# Explore to determine which species to remove ----



# Load libraries
library(tidyverse)

# 1. Read in Bird Data (Site x Species Matrix)
birds <- read.csv("Output/1-mpb-bird-data-mean-count-wide.csv")

# 2. Identify Species Columns (Four-Letter Uppercase Codes)
species_columns <- grep("^[A-Z]{4}$", colnames(birds), value = TRUE)

# 3. Create a List of Species Names
species_list <- data.frame(species = species_columns)

# 4. Define the Correct List of Species to Remove (Non-Territorial A & Misidentified Birds)
species_to_remove <- c(
  # Icterids
  "RWBL", "BHCO",  
  # Corvids
  "AMCR", "BLJA", "CORA", "CAJA",  
  # Aerial Insectivores
  "TRES", "VGSW", "CONI", "PUMA",  
  # Waterfowl
  "CANG", "GWTE", "GADW",  
  # Shorebirds & Marsh Birds
  "GRYE", "COLO", "SORA", "SPSA", "WISN", "SOSA",
  # Owls & Raptors
  "GHOW", "NSWO",  
  # Loons & Marsh Birds
  "SACR",  
  # Misidentified / Non-Territorial Birds
  "RUGR", "PAWA", "NAWA", "MOWA", "SWSP", 
  "DUGR", "BAWW", "BLPW", "NOPA"
)

# 5. Remove Non-Type A Species (Drop Columns)
birds_filtered <- birds %>%
  select(-all_of(species_to_remove))

# 6. Verify Remaining Species
final_species_columns <- grep("^[A-Z]{4}$", colnames(birds_filtered), value = TRUE)
final_species_list <- data.frame(species = final_species_columns)

# Print the final species list to confirm
print(final_species_list)

# 7. Save the Filtered Dataset
write.csv(birds_filtered, "/Users/Bronwyn/Documents/local-git/MPB/Output/1-MPB-bird-data-territoryA.csv", row.names = FALSE)

