# To install PrestoGP, use the command
# devtools::install_github("NIEHS/PrestoGP")
# Note that this requires the devtools package to be installed
# All other packages for this project (including devtools) are available
# on CRAN
library(sf)
library(dplyr)
library(elevatr)
library(terra)
library(tmap)
library(tigris)
library(PrestoGP)

# Read in the data
no2_data <- read.csv("https://github.com/NIEHS/LURK-Vecchia/raw/refs/heads/master/data/US_ST_NO2_Data.csv")
power_plants <- st_read("Power_Plants.shp")

# Create an sf object from the NO2 monitoring stations
# Convert the NO2 data to an sf object using latitude and longitude
no2_sf <- st_as_sf(no2_data, coords = c("Longitude", "Latitude"), crs = 4326)

# Ensure both datasets use the same coordinate reference system (CRS)
power_plants <- st_transform(power_plants, crs = st_crs(no2_sf))

# Find the mean NO2 level at each monitoring station
no2_means <- no2_data %>%
  group_by(ID, Latitude, Longitude) %>%
  summarize(Y = mean(Y, na.rm = TRUE), .groups = "drop")
no2_means_sf <- st_as_sf(no2_means, coords = c("Longitude", "Latitude"),
  crs = 4326)

# Select power plants with PrimSource as "coal", "natural gas", or "petroleum"
power_plants_co2 <- power_plants %>%
  filter(PrimSource %in% c("coal", "natural gas", "petroleum"))

# Create a 10km buffer around each point in no2_means
no2_means_buffered <- st_buffer(no2_means_sf, dist = 10000)

# Identify rows of no2_means with at least one relevant power plant within the
# buffer. Perform a spatial join to identify intersections.
intersections <- st_join(no2_means_buffered, power_plants_co2,
  join = st_intersects)

# Create a column to indicate if the monitoring station is near a power plant
no2_comparison <- intersections %>%
  group_by(ID) %>%
  summarize(Y = mean(Y),
    near_power_plant = !anyNA(PrimSource))
no2_means_sf$near_power_plant <- no2_comparison$near_power_plant

# Perform a t-test on NO2 levels
# Split the NO2 levels based on proximity to power plants
t_test_result <- t.test(no2_means_sf$Y[no2_means_sf$near_power_plant],
  no2_means_sf$Y[!no2_means_sf$near_power_plant])
print(t_test_result)

# Add a column for symbol type based on proximity
no2_means_sf <- no2_means_sf %>%
  mutate(symbol_type = ifelse(near_power_plant, "Near Power Plant",
    "Not Near Power Plant"))

# Create a base map of the United States
us_geo <- states(class = "sf", progress = FALSE)
us_geo <- st_transform(us_geo, crs = st_crs(no2_means_sf))
lower48 <- us_geo %>%
  filter(REGION != 9) %>%
  shift_geometry()
us_map <- tm_shape(lower48) + tm_borders()

# Plot the NO2 levels with different symbols for proximity
us_map + tm_shape(no2_means_sf) +
  tm_symbols(col = "Y", size = 0.05, palette = "viridis",
             style = "jenks",
             title.col = "Mean NO2 Level",
             shape = "symbol_type", title.shape = "Proximity") +
  tm_layout(legend.outside = TRUE)

# Fit a PrestoGP model to the original NO2 data:
Y <- no2_data$Y
locs <- as.matrix(cbind(no2_data$Longitude, no2_data$Latitude,
                        no2_data$YearFrac))
X <- as.matrix(no2_data[, -(1:6)])
no2_pgp <- new("VecchiaModel", n_neighbors = 25)
no2_pgp <- prestogp_fit(no2_pgp, Y, X, locs, scaling = c(1, 1, 2),
                        optim.control=list(trace = 0, reltol = 1e-2, maxit = 5000))

# Group by "ID" to get unique monitoring stations
unique_stations <- no2_sf %>% group_by(ID) %>% slice(1)

# Get elevation data for each monitoring station using the elevatr package
elevation_data <- get_elev_point(unique_stations, src = "aws", z = 10,
  override_size_check = TRUE)
unique_stations$Elevation0 <- elevation_data$elevation

# Create a 1 km buffer around each unique monitoring station
buffers <- st_buffer(unique_stations, dist = 1000)

# Retrieve a digital elevation model (DEM) for the area using elevatr
dem <- get_elev_raster(locations = buffers, z = 10,
  override_size_check = TRUE)

# Convert DEM to a terra object for easier handling
dem_raster <- rast(dem)

# Calculate the maximum elevation and average slope within each buffer
max_elevation <- rep(NA, nrow(unique_stations))
mean_slope <- rep(NA, nrow(unique_stations))
for (i in seq_len(nrow(unique_stations))) {
  # Crop the DEM to the current buffer
  cropped_dem <- crop(dem_raster, buffers[i, ])
  # Mask the DEM to the exact buffer shape
  masked_dem <- mask(cropped_dem, vect(buffers[i, ]))
  # Calculate the maximum elevation
  max_elevation[i] <- max(values(masked_dem), na.rm = TRUE)
  # Calculate the mean slope
  buffer_slope <- terrain(masked_dem, v = "slope")
  mean_slope[i] <- mean(values(buffer_slope), na.rm = TRUE)
}
buffers$Max_Elevation <- max_elevation
buffers$Mean_Slope <- mean_slope

# Use a spatial join to find power plants within the buffers
joined_data <- st_join(buffers, power_plants, join = st_intersects)

# Calculate the sum of coal/natural gas/petroleum MW capacity within each buffer
results <- joined_data %>%
  group_by(ID) %>%
  summarize(Elevation_Diff1km = mean(Max_Elevation) - mean(Elevation0),
    Mean_Slope1km = mean(Mean_Slope),
    Coal_MW1km = sum(Coal_MW, na.rm = TRUE),
    NG_MW1km = sum(NG_MW, na.rm = TRUE),
    Petrol_MW1km = sum(Crude_MW, na.rm = TRUE))

# Join the new data points to the original NO2 data set
no2_data <- left_join(no2_data, results, by = "ID")

# Drop the "geometry" column
no2_data <- no2_data[, -ncol(no2_data)]

# Now repeat the above steps for 10km and 100km buffers:
buffers <- st_buffer(unique_stations, dist = 10000)
dem <- get_elev_raster(locations = buffers, z = 10,
  override_size_check = TRUE)
dem_raster <- rast(dem)
max_elevation <- rep(NA, nrow(unique_stations))
mean_slope <- rep(NA, nrow(unique_stations))
for (i in seq_len(nrow(unique_stations))) {
  cropped_dem <- crop(dem_raster, buffers[i, ])
  masked_dem <- mask(cropped_dem, vect(buffers[i, ]))
  max_elevation[i] <- max(values(masked_dem), na.rm = TRUE)
  buffer_slope <- terrain(masked_dem, v = "slope")
  mean_slope[i] <- mean(values(buffer_slope), na.rm = TRUE)
}
buffers$Max_Elevation <- max_elevation
buffers$Mean_Slope <- mean_slope
joined_data <- st_join(buffers, power_plants, join = st_intersects)
results <- joined_data %>%
  group_by(ID) %>%
  summarize(Elevation_Diff10km = mean(Max_Elevation) - mean(Elevation0),
    Mean_Slope10km = mean(Mean_Slope),
    Coal_MW10km = sum(Coal_MW, na.rm = TRUE),
    NG_MW10km = sum(NG_MW, na.rm = TRUE),
    Petrol_MW10km = sum(Crude_MW, na.rm = TRUE))
no2_data <- left_join(no2_data, results, by = "ID")
no2_data <- no2_data[, -ncol(no2_data)]

buffers <- st_buffer(unique_stations, dist = 100000)
dem <- get_elev_raster(locations = buffers, z = 10,
  override_size_check = TRUE)
dem_raster <- rast(dem)
max_elevation <- rep(NA, nrow(unique_stations))
mean_slope <- rep(NA, nrow(unique_stations))
for (i in seq_len(nrow(unique_stations))) {
  cropped_dem <- crop(dem_raster, buffers[i, ])
  masked_dem <- mask(cropped_dem, vect(buffers[i, ]))
  max_elevation[i] <- max(values(masked_dem), na.rm = TRUE)
  buffer_slope <- terrain(masked_dem, v = "slope")
  mean_slope[i] <- mean(values(buffer_slope), na.rm = TRUE)
}
buffers$Max_Elevation <- max_elevation
buffers$Mean_Slope <- mean_slope
joined_data <- st_join(buffers, power_plants, join = st_intersects)
results <- joined_data %>%
  group_by(ID) %>%
  summarize(Elevation_Diff100km = mean(Max_Elevation) - mean(Elevation0),
    Mean_Slope100km = mean(Mean_Slope),
    Coal_MW100km = sum(Coal_MW, na.rm = TRUE),
    NG_MW100km = sum(NG_MW, na.rm = TRUE),
    Petrol_MW100km = sum(Crude_MW, na.rm = TRUE))
no2_data <- left_join(no2_data, results, by = "ID")
no2_data <- no2_data[, -ncol(no2_data)]

# Create a new set of predictors including our new variables:
X_new <- as.matrix(no2_data[, -(1:6)])

# Delete predictors that are all equal to 0:
X_new <- X_new[,-which(colSums(abs(X_new)) == 0)]

# Fit a PrestoGP model to the expanded NO2 data:
no2_new_pgp <- new("VecchiaModel", n_neighbors = 25)
no2_new_pgp <- prestogp_fit(no2_new_pgp, Y, X_new, locs, scaling = c(1, 1, 2),
  optim.control=list(trace = 0, reltol = 1e-2, maxit = 5000))

# Clean up:
rm(buffers, joined_data, no2_sf, results, unique_stations, Y, X, X_new, locs,
  elevation_data, dem, dem_raster, max_elevation, no2_means, no2_means_sf,
  power_plants_co2, no2_means_buffered, intersections, no2_comparison,
  t_test_result, us_geo, lower48, us_map, mean_slope, buffer_slope,
  cropped_dem, masked_dem, i)
