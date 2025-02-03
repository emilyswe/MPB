# ---
# title: MPB project - get some GIS
# author: Elly Knight & Emily Swerdfager
# created: Oct 25, 2024
# ---

#NOTES################################


#PREAMBLE############################
setwd("~/Documents/local-git/MPB")

#1. Load packages----
library(tidyverse) #basic data wrangling
library(sf) #work with polygons
library(terra) #for the crs function

#2. Set GD root path----
root <- "G:/.shortcut-targets-by-id/1vozK1hbdXP5HCPjjn4rkaDt7Y2MqqQbM/MPB"

#3. Get the points----
dat <- read.csv("/Users/Bronwyn/Documents/local-git/MPB/Input/All_MPB_Locations.csv")

str(dat)

#4. Turn off sci not----
options(scipen=99999)

#WRANGLING###############

#1. Check----
ggplot(dat) + 
  geom_point(aes(x=longitude, y=latitude, colour = factor(year))) +
  facet_wrap(~factor(year))

#2. Make a sf object----
pts <- st_as_sf(dat, coords=c("longitude", "latitude"), crs=4326, remove=FALSE) |> 
  st_transform(3400)

#3. Buffer our object----
buff <- st_buffer(pts, 300) |> 
  mutate(area = as.numeric(st_area(geometry)))

#save pts as shapefile
# Save the spatial object as a shapefile
st_write(pts, "/Users/Bronwyn/Documents/local-git/MPB/All_MPB_Locations.shp")

#4. Take a sample for testing----
 buff <- buff |> 
   sample_n(10)

#PINE#########

#1. Read it in----
pine <- read_sf("/Users/Bronwyn/Documents/local-git/MPB/Input/All_MPB_Locations.shp")

#2. Look at it----
str(pine)

#3. Get the pine----
pts.pine <- buff |> 
  st_transform(crs(pine)) |> #match the projections
  st_intersection() |> #get the unioned polygon
  mutate(area_purepine = as.numeric(st_area(geometry))) |> #get the area of each pine*buffer poly
  st_drop_geometry() |> 
  group_by(location, latitude, longitude, year, area) |> #pick the columsn you want to retain and group by them (location level)
  summarize(area_purepine = sum(area_purepine)) |> #sum the pine polys
  ungroup() |> 
  mutate(prop_purepine = area_purepine/area) #calculate proportion

#4. Plot----
ggplot(pts.pine) +
  geom_histogram(aes(x=prop_purepine))

#5. Make fancy plots
# Create and customize the plot as count = y-axis and prop pure pine = x-axis
pine_histogram <- ggplot(pts.pine) +
  geom_histogram(aes(x = prop_purepine), bins = 20, fill = "steelblue", color = "black") +
  labs(
    title = "Proportion of Pure Pine Across Study Sites",
    x = "Proportion of Pure Pine",
    y = "Number of Sites"
  ) +
  scale_x_continuous(
    limits = c(0.9, 1.05),  # Set range for x-axis
    breaks = seq(0.9, 1.0, by = 0.01),  # Add more granular tick marks
    labels = scales::number_format(accuracy = 1, scale = 100)  # Show numbers without '%'
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.background = element_rect(fill = "white", color = NA),  # White background for panel
    plot.background = element_rect(fill = "white", color = NA),   # White background for entire plot
    panel.grid.major = element_line(color = "gray80"),            # Lighter gridlines
    panel.grid.minor = element_blank()                            # Remove minor gridlines
  )

# Save the plot as a PNG
ggsave(
  filename = "/Users/Bronwyn/Documents/local-git/MPB/Output/pine_histogram-es.png",
  plot = pine_histogram,
  width = 8,  # Width in inches
  height = 6, # Height in inches
  dpi = 300   # Resolution
)


###################make plot showing proportion of sites in pure pine
# Create bins and calculate proportions
pts.pine <- pts.pine |> 
  mutate(bin = cut(prop_purepine, breaks = seq(0.9, 1.05, by = 0.01), include.lowest = TRUE)) |>
  group_by(bin) |>
  summarize(count = n()) |>
  mutate(proportion = count / sum(count))

# Plot with proportion on the y-axis
pine_histogram <- ggplot(pts.pine) +
  geom_bar(aes(x = bin, y = proportion), stat = "identity", fill = "steelblue", color = "black") +
  labs(
    title = "Proportion of Pure Pine Across Study Sites",
    x = "Proportion of Pure Pine",
    y = "Proportion of Sites"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_blank()
  )

# Save the plot
ggsave(
  filename = "/Users/Bronwyn/Documents/local-git/MPB/Output/pine_histogram_proportion_sites_confirmed.png",
  plot = pine_histogram,
  width = 8,
  height = 6,
  dpi = 300
)


#WRR###########

#1. Look at the layers----
st_layers(file.path(root, "WRR Areas_Emily Swerdfager_2024-01-19.gdb"))

#2. Read it in----
#get rid of the one weird polygon and then make the rest valid
wrr <- read_sf(file.path(root, "WRR Areas_Emily Swerdfager_2024-01-19.gdb"), "WRR_Complete") |> 
  dplyr::filter(!is.na(st_is_valid(read_sf(file.path(root, "WRR Areas_Emily Swerdfager_2024-01-19.gdb"), "WRR_Complete")))) |> 
  st_make_valid()

#3. Look at it----
str(wrr)

ggplot(wrr) +
  geom_sf(aes(fill=Treatment_Name))

#4. Get the patch size----
pts.wrr <- pts |> #use points not buffer
  st_transform(crs(wrr)) |> #match the projections
  st_join(wrr) |> #get the spatial join because don't care about the polygon union geometry and want to retain the columns from wrr
  st_drop_geometry() |> 
  rename(area_harvest = Shape_Area) |> 
  dplyr::select(location, latitude, longitude, area_harvest)
