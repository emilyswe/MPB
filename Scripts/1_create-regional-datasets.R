#########
#title: create MPB regional datasets
#this script splits the master location shapefile into two datasets; Jasper and Southern AB, respectively. It also removes the sites north of Jasper that will be excluded from analysis.

# Load the libraries
library(sf)
library(dplyr)

# Read the locations shapefile
points <- st_read("/Users/Bronwyn/Documents/local-git/MPB/Input/All_MPB_Locations.shp")

# Filter out northern AB sites based on latitude (remove sites with latitude > 55.5)
points_filtered <- points %>%
  filter(latitude <= 55.5)  # Remove sites with latitude greater than 55.5

# Now split the filtered dataset into Jasper and Southern Alberta
# Example bounding box for Jasper (adjust these values if needed):
# Latitude for Jasper is between 52.5 and 53.5
# Longitude for Jasper is between -118.5 and -117.0
jasper_sites <- points_filtered %>%
  filter(latitude > 52.5 & latitude < 53.5 & 
           longitude > -118.5 & longitude < -117.0)

southern_ab_sites <- points_filtered %>%
  filter(!(latitude > 52.5 & latitude < 53.5 & 
             longitude > -118.5 & longitude < -117.0))

# Save the Jasper and Southern Alberta sites as shapefiles
st_write(jasper_sites, "/Users/Bronwyn/Documents/local-git/MPB/Input/jnp_sites.shp")
st_write(southern_ab_sites, "/Users/Bronwyn/Documents/local-git/MPB/Input/sab_sites.shp")

# Save the Jasper and Southern Alberta sites as CSVs
st_write(jasper_sites, "/Users/Bronwyn/Documents/local-git/MPB/Input/jnp_sites.csv")
st_write(southern_ab_sites, "/Users/Bronwyn/Documents/local-git/MPB/Input/sab_sites.csv")


