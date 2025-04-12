# To install PrestoGP, use the command
# devtools::install_github("NIEHS/PrestoGP")
# Note that this requires the devtools package to be installed
# All other packages for this project (including devtools) are available
# on CRAN
library(sf)
library(dplyr)
library(tmap)
library(spData)
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
data(us_states)
us_geo <- st_transform(us_states, crs = st_crs(no2_means_sf))
us_map <- tm_shape(us_geo) + tm_borders()

# Plot the NO2 levels with different symbols for proximity
us_map + tm_shape(no2_means_sf) +
  tm_symbols(col = "Y", size = 0.5, palette = "plasma",
             style = "jenks",
             title.col = "Mean NO2 Level", border.lwd = NA,
             shape = "symbol_type", title.shape = "Proximity",
             shapes = c(21,24)) +
  tm_layout(legend.outside = TRUE)

# Fit a PrestoGP model to the original NO2 data:
Y <- no2_data$Y
locs <- as.matrix(cbind(no2_data$Longitude, no2_data$Latitude,
                        no2_data$YearFrac))
X <- as.matrix(no2_data[, -(1:6)])
no2_pgp <- new("VecchiaModel", n_neighbors = 25)
no2_pgp <- prestogp_fit(no2_pgp, Y, X, locs, scaling = c(1, 1, 2))

# Group by "ID" to get unique monitoring stations
unique_stations <- no2_sf %>% group_by(ID) %>% slice(1)

# Create a 1 km buffer around each unique monitoring station
buffers <- st_buffer(unique_stations, dist = 1000)

# Use a spatial join to find power plants within the buffers
joined_data <- st_join(buffers, power_plants, join = st_intersects)

# Calculate the sum of coal/natural gas/petroleum MW capacity within each buffer
results <- joined_data %>%
  group_by(ID) %>%
  summarize(Coal_MW1km = sum(Coal_MW, na.rm = TRUE),
    NG_MW1km = sum(NG_MW, na.rm = TRUE),
    Petrol_MW1km = sum(Crude_MW, na.rm = TRUE))

# Join the new data points to the original NO2 data set
no2_data <- left_join(no2_data, results, by = "ID")

# Drop the "geometry" column
no2_data <- no2_data[, -ncol(no2_data)]

# Now repeat the above steps for 10km and 100km buffers:
buffers <- st_buffer(unique_stations, dist = 10000)
joined_data <- st_join(buffers, power_plants, join = st_intersects)
results <- joined_data %>%
  group_by(ID) %>%
  summarize(Coal_MW10km = sum(Coal_MW, na.rm = TRUE),
    NG_MW10km = sum(NG_MW, na.rm = TRUE),
    Petrol_MW10km = sum(Crude_MW, na.rm = TRUE))
no2_data <- left_join(no2_data, results, by = "ID")
no2_data <- no2_data[, -ncol(no2_data)]

buffers <- st_buffer(unique_stations, dist = 100000)
joined_data <- st_join(buffers, power_plants, join = st_intersects)
results <- joined_data %>%
  group_by(ID) %>%
  summarize(Coal_MW100km = sum(Coal_MW, na.rm = TRUE),
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
no2_new_pgp <- prestogp_fit(no2_new_pgp, Y, X_new, locs,
  scaling = c(1, 1, 2))
no2_new_pgp.relax <- new("VecchiaModel", n_neighbors = 25)
no2_new_pgp.relax <- prestogp_fit(no2_new_pgp, Y, X_new, locs,
  scaling = c(1, 1, 2), penalty = "relaxed")
no2_new_pgp.scad <- new("VecchiaModel", n_neighbors = 25)
no2_new_pgp.scad <- prestogp_fit(no2_new_pgp, Y, X_new, locs,
  scaling = c(1, 1, 2),  penalty = "SCAD")
no2_new_pgp.mcp <- new("VecchiaModel", n_neighbors = 25)
no2_new_pgp.mcp <- prestogp_fit(no2_new_pgp, Y, X_new, locs,
  scaling = c(1, 1, 2), penalty = "MCP")

# Clean up:
rm(buffers, joined_data, no2_sf, results, unique_stations, Y, X, X_new, locs,
  no2_means, no2_means_sf, power_plants_co2, no2_means_buffered,
  intersections, no2_comparison, t_test_result, us_geo, us_map)
