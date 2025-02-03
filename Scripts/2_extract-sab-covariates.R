#######created by Lionel and Emily January 13 2025
####edited by Emily January 14 2025
###the purpose of this script is to extract covariates for sab sites, with the main initial focus being the time since disturbance
##first we turn the MPB left standing shapefile into a raster and use that to extract MPB attack year for all sab sites. We use the nearest distance value for NA values, assuming this is valid because all sites were infected by the same outbreak in the 80s.

#setwd
setwd("~/Documents/local-git/MPB")

#####Step #1. - extract disturbance year of beetle attack for sab sites
# Load required libraries
library(sf)
library(dplyr)
library(tidyr)
library(terra)
library(raster)
library(fasterize)

# Step 1: Load the Beetle Left Standing shapefile
beetles <- st_read("/Users/Bronwyn/Documents/local-git/MPB/Input/2_outputs/largerBeetleLeftStanding.shp")

# Step 2: Transform CRS to TMERC
beetlesTMERC <- beetles %>%
  st_transform(crs = "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")

# Step 3: Load Southern Alberta sites shapefile
sab_sites <- st_read("/Users/Bronwyn/Documents/local-git/MPB/Input/sab_sites.shp")

# Step 4: Transform CRS to TMERC
sab_sitesTMERC <- sab_sites %>%
  st_transform(crs = "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")

# Step 5: Define study area and crop beetles data
studyArea <- st_as_sfc(st_bbox(sab_sitesTMERC))
st_write(studyArea, "/Users/Bronwyn/Documents/local-git/MPB/Input/studyArea_sab.shp", append = FALSE)
studyArea <- st_read("/Users/Bronwyn/Documents/local-git/MPB/Input/studyArea_sab.shp") %>%
  st_transform(crs = crs(sab_sitesTMERC))
beetles_studyArea <- st_crop(beetlesTMERC, studyArea)

# Step 6: Rasterize beetle year data
byear <- beetles_studyArea["betl_yr"]
v <- vect(byear)
r <- rast(v, res = 1000)
zBeetle <- rasterize(v, r, field = "betl_yr")
plot(zBeetle)

# Step 7: Create 150m and 300m buffers
b150 <- st_buffer(sab_sitesTMERC, 150)
b300 <- st_buffer(sab_sitesTMERC, 300)

# Step 8: Extract beetle year data for each buffer
pointList <- list()

for (i in 1:nrow(sab_sitesTMERC)) {
  pointData <- sab_sitesTMERC[i, ]
  
  # Extract beetle year within 150m and 300m buffers
  try(meanBeetleYear_150m <- terra::extract(zBeetle, b150[i, ], fun = "max", na.rm = TRUE))
  try(pointData$BY150m <- as.integer(meanBeetleYear_150m$last))
  
  try(meanBeetleYear_300m <- terra::extract(zBeetle, b300[i, ], fun = "max", na.rm = TRUE))
  try(pointData$BY300m <- as.integer(meanBeetleYear_300m$last))
  
  # Handle missing values by finding the nearest valid point in pointList
  if (is.na(pointData$BY150m) & is.na(pointData$BY300m)) {
    valid_points <- do.call(rbind, pointList) %>% 
      st_as_sf(crs = st_crs(sab_sitesTMERC)) %>%
      filter(!is.na(BY150m) | !is.na(BY300m))
    
    if (nrow(valid_points) > 0) {
      distances <- st_distance(st_geometry(pointData), st_geometry(valid_points))
      nearest_idx <- which.min(distances)
      nearest_point <- valid_points[nearest_idx, ]
      
      pointData$BY150m <- nearest_point$BY150m
      pointData$BY300m <- nearest_point$BY300m
    }
  }
  
  # Append to pointList
  try(pointList[[i]] <- pointData)
  print(paste0("Line ", i, " done."))
}

# Step 9: Combine processed points into a data frame
pointBeetleYear <- do.call(rbind, pointList)
summaries <- data.frame(pointBeetleYear)

# Step 10: Save summaries to CSV
write.csv(summaries, file = "/Users/Bronwyn/Documents/local-git/MPB/Output/2_meanBeetleYear_sab.csv")



###############################################
###############################################
###############################################
###############################################
###############################################
###############################################

#########Step #2. extract disturbance year of wildfire for sab sites

# 1. Read the burn shapefile
burn_shapefile <- st_read("/Users/Bronwyn/Documents/local-git/MPB/Input/2_outputs/largerBurnedAfterBeetles.shp")

# 2. Check the structure of the burn shapefile
print(str(burn_shapefile))  # This will give us a summary of the shapefile's columns

# Function to extract the burn year for each site based on buffer distance
get_burn_year <- function(site, burn_shapefile, buffer_dist) {
  # Ensure the CRS of site and burn shapefile are the same
  if (st_crs(site) != st_crs(burn_shapefile)) {
    site <- st_transform(site, crs = st_crs(burn_shapefile))  # Transform to match burn_shapefile CRS
  }
  
  # Create a buffer around the site point (in meters)
  buffer <- st_buffer(site, dist = buffer_dist)
  
  # Check if the site falls within the burn polygon (point-in-polygon)
  burn_polygon <- st_intersects(buffer, burn_shapefile, sparse = FALSE)
  
  # If the site intersects a burn polygon, return the most recent burn year
  if (any(burn_polygon)) {
    return(max(burn_shapefile$YEAR_1[burn_polygon]))  # Most recent burn year
  } else {
    return(NA)  # If no intersection, return NA
  }
}

# Loop through each site and extract burn year for 150m, 300m, and 500m buffers
for (i in 1:nrow(pointBeetleYear)) {
  pointData <- pointBeetleYear[i,]  # Get the site
  
  # Ensure the site is an sf object for spatial operations
  pointData_sf <- st_as_sf(pointData, coords = c("longitude", "latitude"), crs = st_crs(burn_shapefile))
  
  # Extract burn year for 150m buffer
  pointData$FY150m <- get_burn_year(pointData_sf, burn_shapefile, buffer_dist = 150)
  
  # Extract burn year for 300m buffer
  pointData$FY300m <- get_burn_year(pointData_sf, burn_shapefile, buffer_dist = 300)
  
  # Extract burn year for 500m buffer
  pointData$FY500m <- get_burn_year(pointData_sf, burn_shapefile, buffer_dist = 500)
  
  # Update the main dataframe with the new columns for this row
  pointBeetleYear[i, c("FY150m", "FY300m", "FY500m")] <- pointData[c("FY150m", "FY300m", "FY500m")]
  
  print(paste0("Site ", i, " processed."))
}

# After the loop, inspect the dataframe
str(pointBeetleYear)

#save
write.csv(pointBeetleYear, "/Users/Bronwyn/Documents/local-git/MPB/Output/2_FYandBY.csv", row.names = FALSE)

###############################################
###############################################
###############################################
###############################################
###############################################
###############################################

#########Step #3. extract disturbance year of harvest for sab sites

# 1. Read the harvest shapefile
harvest_shapefile <- st_read("/Users/Bronwyn/Documents/local-git/MPB/Input/2_outputs/largerpossibleSalvageLogging.shp")

# 2. Check the structure of the harvest shapefile
print(str(harvest_shapefile))  # This will give us a summary of the shapefile's columns

# Create empty columns for harvest years (initialize with NAs)
pointBeetleYear$H150m <- NA
pointBeetleYear$H300m <- NA
pointBeetleYear$H500m <- NA

# Function to extract harvest year for each site based on buffer distance
get_harvest_year <- function(site, harvest_shapefile, buffer_dist) {
  # Ensure the CRS of site and harvest shapefile are the same
  if (st_crs(site) != st_crs(harvest_shapefile)) {
    site <- st_transform(site, crs = st_crs(harvest_shapefile))  # Transform to match harvest_shapefile CRS
  }
  
  # Create a buffer around the site point (in meters)
  buffer <- st_buffer(site, dist = buffer_dist)
  
  # Check if the site falls within the harvest polygon (point-in-polygon)
  harvest_polygon <- st_intersects(buffer, harvest_shapefile, sparse = FALSE)
  
  # If the site intersects a harvest polygon, return the most recent harvest year
  if (any(harvest_polygon)) {
    return(max(harvest_shapefile$HrvstYr[harvest_polygon]))  # Get the most recent harvest year
  } else {
    return(NA)  # If no intersection, return NA
  }
}

# 1. Loop through each site and extract harvest year for 150m, 300m, and 500m buffers
for (i in 1:nrow(pointBeetleYear)) {
  pointData <- pointBeetleYear[i,]  # Get the site
  
  # Extract harvest year for 150m buffer
  pointBeetleYear$H150m[i] <- get_harvest_year(pointData, harvest_shapefile, buffer_dist = 150)
  
  # Extract harvest year for 300m buffer
  pointBeetleYear$H300m[i] <- get_harvest_year(pointData, harvest_shapefile, buffer_dist = 300)
  
  # Extract harvest year for 500m buffer
  pointBeetleYear$H500m[i] <- get_harvest_year(pointData, harvest_shapefile, buffer_dist = 500)
  
  print(paste0("Site ", i, " processed."))
}

# Inspect the dataframe to see the new columns for harvest years
View(pointBeetleYear)

# Save the updated dataframe with beetle, burn, and harvest years to CSV
write.csv(pointBeetleYear, file = "/Users/Bronwyn/Documents/local-git/MPB/Output/2_pointBeetleYear_with_harvest.csv", row.names = FALSE)


###clean it up
# Drop the geometry column and convert to a regular dataframe
pointBeetleYear_clean <- pointBeetleYear %>%
  st_drop_geometry()  # Removes the geometry column

# Verify the cleaned dataframe
View(pointBeetleYear_clean)

# Save the cleaned dataframe to CSV
write.csv(pointBeetleYear_clean, file = "/Users/Bronwyn/Documents/local-git/MPB/Output/2_pointBeetleYear_with_harvest_clean.csv", row.names = FALSE)


###############################################
###############################################
###############################################
###############################################
###############################################
###############################################

#########Step #4. calculate disturbance year and treatment type for sab sites
# Initialize the new columns in the dataframe to avoid errors
# Initialize the new columns in the dataframe
pointBeetleYear$disturbance_year <- NA  # Initialize the disturbance_year column with NAs
pointBeetleYear$treatment <- NA  # Initialize the treatment column with NAs

# Loop through each site in the dataframe
for (i in 1:nrow(pointBeetleYear)) {
  pointData <- pointBeetleYear[i, ]  # Get the site data
  
  # Get the disturbance years for all the relevant columns
  disturbance_years <- c(pointData$BY150m, pointData$BY300m, pointData$FY150m, pointData$FY300m, 
                         pointData$FY500m, pointData$H150m, pointData$H300m, pointData$H500m)
  
  # Get the names of the disturbance columns
  disturbance_column_names <- c("BY150m", "BY300m", "FY150m", "FY300m", 
                                "FY500m", "H150m", "H300m", "H500m")
  
  # Remove NA values and match valid years with their column names
  valid_years <- disturbance_years[!is.na(disturbance_years)]  # Remove NAs before getting the max value
  valid_columns <- disturbance_column_names[!is.na(disturbance_years)]
  
  # If there's at least one valid year, assign the most recent year to disturbance_year
  if (length(valid_years) > 0) {
    # Get the most recent disturbance year
    max_year <- max(valid_years)
    pointData$disturbance_year <- max_year
    
    # Determine the treatment based on the column associated with the max year
    max_column <- valid_columns[which.max(valid_years)]
    
    if (max_column %in% c("BY150m", "BY300m")) {
      pointData$treatment <- "standing"
    } else if (max_column %in% c("FY150m", "FY300m", "FY500m")) {
      pointData$treatment <- "burn"
    } else if (max_column %in% c("H150m", "H300m", "H500m")) {
      pointData$treatment <- "harvest"
    }
  } else {
    pointData$disturbance_year <- NA  # If no valid year, set as NA
    pointData$treatment <- NA  # No treatment if no disturbance year is found
  }
  
  # Update the row in the dataframe
  pointBeetleYear[i, ] <- pointData  # Update the row
}

# Check the updated dataframe
head(pointBeetleYear)


# Step 1: Ensure `pointBeetleYear` has all the updated columns (including `disturbance_year` and `treatment`)
# (This is from your earlier loop and logic)

# Step 2: Clean the updated dataframe (drop geometry and ensure all columns are included)
pointBeetleYear_clean <- pointBeetleYear %>%
  st_drop_geometry()  # Removes the geometry column

# Step 3: Verify the cleaned dataframe (optional, for inspection)
View(pointBeetleYear_clean)

# Step 4: Save the cleaned dataframe to CSV
write.csv(pointBeetleYear_clean, 
          file = "/Users/Bronwyn/Documents/local-git/MPB/Output/2_pointBeetleYear_with_harvest_clean.csv", 
          row.names = FALSE)



###############################################
###############################################
###############################################
###############################################
###############################################
###############################################

##########Step #5. Calculate TSD for sab sites


# Step 1: Calculate TSD (Time Since Disturbance)
pointBeetleYear_clean <- pointBeetleYear_clean %>%
  mutate(TSD = ifelse(!is.na(disturbance_year), year - disturbance_year, NA))  # Calculate TSD where disturbance_year is valid

# Step 2: View the updated dataframe (optional)
View(pointBeetleYear_clean)

# Step 3: Save the updated dataframe to a new CSV file
write.csv(pointBeetleYear_clean, 
          file = "/Users/Bronwyn/Documents/local-git/MPB/Output/2_tsd-sab.csv", 
          row.names = FALSE)

###############################################
###############################################
###############################################
###############################################
###############################################
###############################################

##########Step #6. visualize the data


library(ggplot2)

# Combined histogram with fill for treatment
tsd_treatment_plot <- ggplot(pointBeetleYear_clean, aes(x = TSD, fill = treatment)) +
  geom_histogram(binwidth = 5, position = "dodge", color = "black", alpha = 0.7) +
  scale_fill_manual(values = c("burn" = "red", "harvest" = "blue", "standing" = "forestgreen")) +
  labs(
    title = "Distribution of TSD by Treatment",
    x = "Time Since Disturbance (TSD)",
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

# Save the plot with `ggsave`
ggsave(
  filename = "/Users/Bronwyn/Documents/local-git/MPB/Output/TSD_distribution_by_treatment.png",
  plot = tsd_treatment_plot,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"  # Ensure the background is white
)


###############################################
###############################################
#create TSD histogram for each treatment 

library(ggplot2)
library(patchwork)

# Create histograms for each treatment with centered titles and adjusted green color
p_standing <- ggplot(subset(pointBeetleYear_clean, treatment == "standing"), aes(x = TSD)) +
  geom_histogram(binwidth = 5, fill = "forestgreen", color = "black", alpha = 0.7) +
  labs(title = "Standing", x = "Time Since Disturbance (TSD)", y = "Count") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))  # Center title

p_burn <- ggplot(subset(pointBeetleYear_clean, treatment == "burn"), aes(x = TSD)) +
  geom_histogram(binwidth = 5, fill = "red", color = "black", alpha = 0.7) +
  labs(title = "Burn", x = "Time Since Disturbance (TSD)", y = "Count") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))  # Center title

p_harvest <- ggplot(subset(pointBeetleYear_clean, treatment == "harvest"), aes(x = TSD)) +
  geom_histogram(binwidth = 5, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Harvest", x = "Time Since Disturbance (TSD)", y = "Count") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))  # Center title

# Combine the plots into one figure
combined_plot <- p_standing / p_burn / p_harvest  # Stacks the histograms vertically

# Save the combined plot as a JPEG
ggsave("/Users/Bronwyn/Documents/local-git/MPB/Output/TSD_histograms_by_treatment.jpeg", 
       plot = combined_plot, 
       width = 8, height = 10, dpi = 300)



###############################################
###############################################
# Histogram of sites by treatment

# Create the bar plot for distribution of sites by treatment
site_treatment_plot <- ggplot(pointBeetleYear_clean, aes(x = treatment, fill = treatment)) +
  geom_bar(color = "black", alpha = 0.7) +
  scale_fill_manual(values = c("burn" = "red", "harvest" = "blue", "standing" = "forestgreen")) +
  labs(title = "Distribution of Sites by Treatment",
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
ggsave("/Users/Bronwyn/Documents/local-git/MPB/Output/Sites_by_treatment.jpeg", 
       plot = site_treatment_plot, 
       width = 6, height = 4, dpi = 300)



###############################################
###############################################
###############################################
###############################################
###############################################
###############################################

##########Step #7. Get number of attacked trees and severity of mpb attack for sab sites


library(sf)
library(dplyr)
library(terra)

# Assuming you already have the beetlesTMERC object loaded (representing the beetle shapefile)
# and the pointBeetleYear object contains your sites with buffers.

# Create buffers around the points in pointBeetleYear
b150 <- st_buffer(pointBeetleYear, dist = 150)  # 150m buffer
b300 <- st_buffer(pointBeetleYear, dist = 300)  # 300m buffer
b500 <- st_buffer(pointBeetleYear, dist = 500)  # 500m buffer

# Initialize columns for num_trs and severity
pointBeetleYear$num_trs <- NA
pointBeetleYear$severity <- NA

# Loop through each site to extract num_trs and severity
for (i in 1:nrow(pointBeetleYear)) {
  # Get the site buffer
  siteBuffer <- b500[i, ]  # Using 500m buffer as the widest net
  
  # Intersect the buffer with the beetle shapefile to get the overlapping polygons
  intersectedBeetles <- st_intersection(beetlesTMERC, siteBuffer)
  
  if (nrow(intersectedBeetles) > 0) {
    # Extract the sum of num_trs (number of attacked trees) within the buffer
    total_num_trs <- sum(intersectedBeetles$num_trs, na.rm = TRUE)
    
    # Extract the maximum severity value within the buffer
    max_severity <- max(intersectedBeetles$SEVERIT, na.rm = TRUE)
    
    # Assign the extracted values to the respective columns
    pointBeetleYear$num_trs[i] <- total_num_trs
    pointBeetleYear$severity[i] <- max_severity
  } else {
    # Assign NA if no intersected polygons
    pointBeetleYear$num_trs[i] <- NA
    pointBeetleYear$severity[i] <- NA
  }
  
  # Print progress
  print(paste0("Processed site ", i, " of ", nrow(pointBeetleYear)))
}

# View the updated dataframe
View(pointBeetleYear)
str(pointBeetleYear)

# Save the updated dataframe to a CSV file
write.csv(st_drop_geometry(pointBeetleYear), 
          file = "/Users/Bronwyn/Documents/local-git/MPB/Output/pointBeetleYear_with_num_trs_and_severity.csv", 
          row.names = FALSE)



###############################################
###############################################
###############################################
###############################################
###############################################
###############################################

##########Step #8. Visualize number of attacked trees and severity of mpb attack for sab sites

# Filter non-NA values for num_trs and severity
num_trs_data <- pointBeetleYear$num_trs[!is.na(pointBeetleYear$num_trs)]
severity_data <- pointBeetleYear$severity[!is.na(pointBeetleYear$severity)]

# Check for valid data
if (length(num_trs_data) > 0 && length(severity_data) > 0) {
  # Histogram for number of attacked trees
  p_num_trs <- ggplot(data = data.frame(num_trs = num_trs_data), aes(x = num_trs)) +
    geom_histogram(binwidth = 20, fill = "blue", color = "black", alpha = 0.7) +
    labs(
      title = "Distribution of Number of Attacked Trees",
      x = "Number of Attacked Trees",
      y = "# of Sites"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) # Center the title
  
  # Histogram for severity
  p_severity <- ggplot(data = data.frame(severity = severity_data), aes(x = severity)) +
    geom_histogram(binwidth = 1, fill = "forestgreen", color = "black", alpha = 0.7) +
    labs(
      title = "Distribution of Severity",
      x = "Severity",
      y = "# of Sites"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) # Center the title
  
  # Combine the histograms into one image
  combined_plot <- p_num_trs / p_severity
  
  # Save the combined plot
  ggsave(
    "/Users/Bronwyn/Documents/local-git/MPB/Output/num_trs_and_severity_histograms.jpeg",
    plot = combined_plot,
    width = 8, height = 10, dpi = 300
  )
} else {
  stop("No valid data in num_trs or severity columns to plot.")
}
