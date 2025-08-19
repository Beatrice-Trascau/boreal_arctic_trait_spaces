##----------------------------------------------------------------------------##
# PAPER 3: BOREAL AND ARCTIC PLANT SPECIES TRAIT SPACES 
# 3.2_traits_biome_boundary_distance_relationships
# This script contains code which explores ther relationships between species
# trait values and their mean distance to biome boundaries
##----------------------------------------------------------------------------##

# 1. PREPARE DATA --------------------------------------------------------------

## 1.1. Load data --------------------------------------------------------------

# Load cleaned traits
load(here("data", "derived_data", "TRY_traits_cleaned_July2025.RData"))
cleaned_traits <- cleaned_traits_July2025

# Load distances to biome boundaries
detailed_results <- readRDS(here("data", "derived_data", 
                                 "species_summaries_dist_to_biome_boundary_June25.rds"))

# Load biome boundaries
global_biomes <- st_read(here("data", "raw_data", "biomes", 
                                             "wwf_terr_ecos.shp"))

## 1.2. Prepare biomes ---------------------------------------------------------

# Load Boreal Forest (BIOME = 6)
boreal_forest <- st_union(global_biomes[global_biomes$BIOME == 6,])

# Load Tundra (BIOME = 11)
tundra <- st_union(global_biomes[global_biomes$BIOME == 11 &(global_biomes$REALM == "PA"|global_biomes$REALM == "NA"), ])

# Make sure geometries are valid
boreal_forest <- st_make_valid(boreal_forest)
tundra <- st_make_balid(tundra)

# Re-project biomes to North Pole Lamvert Azimuthal Equal Area (EPSG: 3574)
boreal_sf <- st_transform(boreal_forest, "EPSG:3574")
tundra_sf <- st_transform(tundra, "EPSG:3574")

# 2. ADD EXTRA INFORMATION TO CLEANED TRAITS -----------------------------------

## 2.1. Add species-level biome boundary distances -----------------------------

# Remove Elodea canadensis & hybrids
biome_boundaries <- detailed_results |>
  filter(!species == "Elodea canadensis",
         !grepl("×", species))

# Combine the dataframes
traits_biome_boundaries <- cleaned_traits |>
  left_join(biome_boundaries, by = c("StandardSpeciesName" = "species"))

# Rename column with mean distance to biome boundary to reflect the fact that
  # it is a species-level mean


## 2.2. Calculate distance to biome boundary for each trait record -------------

# Filter out records not within study area
traits_with_coords <- traits_biome_boundaries |>
  filter(!is.na(LON_site), !is.na(LAT_site)) |>
  filter(LON_site >= -180, LON_site <= 180,
         LAT_site >= -90, LAT_site <= 90)

# Check how many records were removed
nrow(traits_biome_boundaries)
nrow(traits_with_coords)

# Convert traits to spatial objec
traits_sf <- st_as_sf(traits_with_coords, 
                      coords = c("LON_site", "LAT_site"), 
                      crs = 4326) |>
  st_transform(crs = "EPSG:3574") 

# Get boundaries of both biomes
boreal_boundary <- st_boundary(boreal_sf)
tundra_boundary <- st_boundary(tundra_sf)

# Check which biome each point is in
in_boreal <- st_intersects(traits_sf, boreal_sf, sparse = FALSE)[,1]
in_tundra <- st_intersects(traits_sf, tundra_sf, sparse = FALSE)[,1]

# Create a distance column
traits_with_coords$record_level_distance_to_biome_boundary <- NA

# Give positive distance to points in boreal forest
if(sum(in_boreal) > 0){
  boreal_points <- traits_sf[in_boreal, ]
  boreal_distances <- st_distance(boreal_points, boreal_boundary)
  # take minimum distance to boundary
  traits_with_coords$record_level_distance_to_biome_boundary <- 
    apply(boreal_distances, 1, min) / 1000 # convert to km
}

# Give negative distance to points in tundra
if(sum(in_tundra) > 0){
  tundra_points <- traits_sf[in_tundra, ]
  tundra_distances <- st_distance(tundra_points, tundra_boundary)
  # take minimum distance to boundary
  traits_with_coords$record_level_distance_to_biome_boundary <- 
    -apply(tundra_distances, 1, min) / 1000 # convert to km
}

# Add biome classification for individual records
traits_with_coords <- traits_with_coords |>
  mutate(individual_record_biome = case_when(in_boreal ~ "boreal",
                                             in_tundra ~ "tundra",
                                             TRUE ~ "other"))

# Update the main dataframe
traits_biome_boundaries <- traits_with_coords

## 2.3. Clean growth form column -----------------------------------------------

## 2.4. Add LHS strategy for each species --------------------------------------

# Check if we have all required traits
required_traits_lhs <- c("PlantHeight", "SeedMass", "SLA")
available_traits <- unqiue(traits_biome_boundaries$TraitNameNew)

# Check if any traits are missing
missing_lhs <- setdiff(required_traits_lhs, available_traits)
available_lhs <- intersect(required_traits_lhs, available_traits)

lhs_data <- traits_biome_boundaries |>
  # filter for required traits
  fileter(TraitNameNew %in% required_traits_lhs) |>
  # remove outliers (>= 5 SDs)
  group_by(StandardSpeciesName, TraitNameNew) |>
  mutate(trait_mean = mean(StdValue, na.rm = TRUE),
         trait_sd = sd(StdValue, na.rm = TRUE)) |>
  filter(abs(StdValue - trait_mean) < 5 * trait_mean) |>
  ungroup() |>
  # select relevant columns
  select(StandardSpeciesName, TraitNameNew, StdValue) |>
  # calculate median value per species per trait
  group_by(StandardSpeciesName, TraitNameNew) |>
  summarise(MedianTraitValue = median(StdValue, na.rm = TRUE),
            .groups = "drop") |>
  # convert to wide format
  pivot_wider(names_from = TraitNameNew,
              values_from = MedianTraitValue) |>
  # rename columns to match MultiTraits requirements
  rename(Height = PlantHeight,
         SeedMass = SeedMass,
         SLA = SLA) |>
  # removed species with missing data for required traits
  filter(!is.na(PlantHeight), !is.na(SLA), !is.na(SeedMass))

# Select only the columns needed for the LHS calculations
lhs_traits_only <- lhs_data |>
  select(SLA, Height, SeedMass)

# Calculate LHS strategies
lhs_results <- LHS(lhs_traits_only)

# Add species names to results
lhs_final <- lhs_results |>
  mutate(StandardSpeciesName = lhs_data$StandardSpeciesName)

# Add LHS strategies to the main dataframe
traits_biome_boundaries_with_lhs <- traits_biome_boundaries |>
  left_join(lhs_final |>
              select(StandardSpeciesName, contains("lhs") | contains("LHS") | contains("strategy")),
            by = "StandardSpeciesName")


# 3. PLOT RELATIONSHIPS --------------------------------------------------------

# Calculate universal x-axis limits
x_min <- min(traits_biome_boundaries$mean_distance_km, na.rm = TRUE)
x_max <- max(traits_biome_boundaries$mean_distance_km, na.rm = TRUE)

## 3.1. Filter for each trait --------------------------------------------------

# Filter for plant height
plant_height <- traits_biome_boundaries |>
  filter(TraitNameNew == "PlantHeight")

# Filter for SLA
sla <- traits_biome_boundaries |>
  filter(TraitNameNew == "SLA")

# Filter for seed mass
seed_mass <- traits_biome_boundaries |>
  filter(TraitNameNew == "SeedMass")

# Filter for leaf N
leaf_N <- traits_biome_boundaries |>
  filter(TraitNameNew == "LeafN")

# Filter for leaf N:C ratio
leaf_CN_ratio <- traits_biome_boundaries |>
  filter(TraitNameNew == "LeafCN")

## 2.2 Plot all the data -------------------------------------------------------

# Plot plant height
plant_height_all_data <- ggplot(plant_height, aes(x = mean_distance_km, 
                                                  y = CleanedTraitValue)) +
  geom_point() +
  geom_smooth() +
  xlim(x_min, x_max) +
  labs(x = "Mean Distance to Biome Boundary (km)", y = "Plant height (cm)") +
  theme_classic()

# Plot SLA
sla_all_data <- ggplot(sla, aes(x = mean_distance_km,
                                y = CleanedTraitValue)) +
  geom_point() +
  geom_smooth() +
  xlim(x_min, x_max) +
  labs(x = "Mean Distance to Biome Boundary (km)", y = expression("SLA (m"^2*"/kg)")) +
  theme_classic()

# Plot seed mass
seed_mass_all_data <- ggplot(seed_mass, aes(x = mean_distance_km,
                                            y = CleanedTraitValue)) +
  geom_point() +
  geom_smooth() +
  xlim(x_min, x_max) +
  labs(x = "Mean Distance to Biome Boundary (km)", y = "Seed Mass (mg)") +
  theme_classic()

# Plot Leaf N content
leaf_N_all_data <- ggplot(leaf_N, aes(x = mean_distance_km,
                                      y = CleanedTraitValue)) +
  geom_point() +
  geom_smooth() +
  xlim(x_min, x_max) +
  labs(x = "Mean Distance to Biome Boundary (km)", y = "Leaf N (mg/g)") +
  theme_classic()

# Plot Leaf C:N ratio
leaf_CN_ratio_all_data <- ggplot(leaf_CN_ratio, aes(x = mean_distance_km,
                                                    y = CleanedTraitValue)) +
  geom_point() +
  geom_smooth() +
  xlim(x_min, x_max) +
  labs(x = "Mean Distance to Biome Boundary (km)", y = "Leaf C:N Ratio") +
  theme_classic()

# Combine plots 
trait_relationships_all_data <- plot_grid(plant_height_all_data, sla_all_data,
                                          seed_mass_all_data, leaf_N_all_data,
                                          leaf_CN_ratio_all_data,
                                          labels = c("a)", "b)", "c)", "d)", "e)"),
                                          nrow = 3,
                                          align = "hv",
                                          axis = "tblr")

## 2.3. Plot median values per species -----------------------------------------

# Calculate median values per species
traits_median <- traits_biome_boundaries |>
  group_by(StandardSpeciesName, TraitNameNew) |>
  summarise(MedianTraitValue = median(CleanedTraitValue, na.rm = TRUE),
            .groups = "drop")

# Add distance to biome boundaries
traits_median_df <- traits_median |>
  left_join(biome_boundaries, by = c("StandardSpeciesName" = "species")) |>
  filter(!is.na(mean_distance_km), !is.na(MedianTraitValue))

# Plant height
median_plant_height <- traits_median_df |>
  filter(TraitNameNew == "PlantHeight") |>
  ggplot(aes(x = mean_distance_km, y = MedianTraitValue)) +
  geom_point() +
  geom_smooth() +
  xlim(x_min, x_max) +
  labs(x = "Mean Distance to Biome Boundary (km)", y = "Median Plant height (cm)") +
  theme_classic()

# SLA
median_sla <- traits_median_df |>
  filter(TraitNameNew == "SLA") |>
  ggplot(aes(x = mean_distance_km, y = MedianTraitValue)) +
  geom_point() +
  geom_smooth() +
  xlim(x_min, x_max) +
  labs(x = "Mean Distance to Biome Boundary (km)", y = expression("Median SLA (m"^2*"/kg)")) +
  theme_classic()

# Seed mass
median_seed_mass <- traits_median_df |>
  filter(TraitNameNew == "SeedMass") |>
  ggplot(aes(x = mean_distance_km, y = MedianTraitValue)) +
  geom_point() +
  geom_smooth() +
  xlim(x_min, x_max) +
  labs(x = "Mean Distance to Biome Boundary (km)", y = "Median Seed Mass (mg)") +
  theme_classic()

# Leaf N content
median_leaf_N <- traits_median_df |>
  filter(TraitNameNew == "LeafN") |>
  ggplot(aes(x = mean_distance_km, y = MedianTraitValue)) +
  geom_point() +
  geom_smooth() +
  xlim(x_min, x_max) +
  labs(x = "Mean Distance to Biome Boundary (km)", y = "Median Leaf N (mg/g)") +
  theme_classic()

# Leaf C:N ratio
median_c_n <- traits_median_df |>
  filter(TraitNameNew == "LeafCN") |>
  ggplot(aes(x = mean_distance_km, y = MedianTraitValue)) +
  geom_point() +
  geom_smooth() +
  xlim(x_min, x_max) +
  labs(x = "Mean Distance to Biome Boundary (km)", y = "Median Leaf C:N Ratio") +
  theme_classic()

# Combine plots 
median_trait_relationships <- plot_grid(median_plant_height, median_sla,
                                        median_seed_mass, median_leaf_N,
                                        median_c_n,
                                        labels = c("a)", "b)", "c)", "d)", "e)"),
                                        nrow = 3,
                                        align = "hv",
                                        axis = "tblr")

## 2.4. Remove values >= 5 SD --------------------------------------------------

# Calculate means
plant_height_mean <- mean(plant_height$CleanedTraitValue, na.rm = TRUE)
sla_mean <- mean(sla$CleanedTraitValue, na.rm = TRUE)
seed_mass_mean <- mean(seed_mass$CleanedTraitValue, na.rm = TRUE)
leaf_N_mean <- mean(leaf_N$CleanedTraitValue, na.rm = TRUE)
leaf_CN_mean <- mean(leaf_CN_ratio$CleanedTraitValue, na.rm = TRUE)

# Calculate SD
plant_height_sd <- sd(plant_height$CleanedTraitValue, na.rm = TRUE)
sla_sd <- sd(sla$CleanedTraitValue, na.rm = TRUE)
seed_mass_sd <- sd(seed_mass$CleanedTraitValue, na.rm = TRUE)
leaf_N_sd <- sd(leaf_N$CleanedTraitValue, na.rm = TRUE)
leaf_CN_sd <- sd(leaf_CN_ratio$CleanedTraitValue, na.rm = TRUE)

# Plant height
plant_height_plot <- plant_height |>
  filter(abs(CleanedTraitValue - plant_height_mean) < 5 * plant_height_sd) |>
  ggplot(aes(x = mean_distance_km, y = CleanedTraitValue)) +
  geom_point() +
  geom_smooth() +
  xlim(x_min, x_max) +
  labs(x = "Mean Distance to Biome Boundary (km)", y = "Plant height (cm)") +
  theme_classic()

# SLA
sla_plot <- sla |>
  filter(abs(CleanedTraitValue - sla_mean) < 5 * sla_sd) |>
  ggplot(aes(x = mean_distance_km, y = CleanedTraitValue)) +
  geom_point() +
  geom_smooth() +
  xlim(x_min, x_max) +
  labs(x = "Mean Distance to Biome Boundary (km)", y = expression("SLA (m"^2*"/kg)")) +
  theme_classic()

# Seed Mass
seed_mass_plot <- seed_mass |>
  filter(abs(CleanedTraitValue - seed_mass_mean) < 5 * seed_mass_sd) |>
  ggplot(aes(x = mean_distance_km, y = CleanedTraitValue)) +
  geom_point() +
  geom_smooth() +
  xlim(x_min, x_max) +
  labs(x = "Mean Distance to Biome Boundary (km)", y = "Seed Mass (mg)") +
  theme_classic()

# Leaf N content
leaf_N_plot <- leaf_N |>
  filter(abs(CleanedTraitValue - leaf_N_mean) < 5 * leaf_N_sd) |>
  ggplot(aes(x = mean_distance_km, y = CleanedTraitValue)) +
  geom_point() +
  geom_smooth() +
  xlim(x_min, x_max) +
  labs(x = "Mean Distance to Biome Boundary (km)", y = "Leaf N (mg/g)") +
  theme_classic()

# Leaf C:N ratio
leaf_CN_plot <- leaf_CN_ratio |>
  filter(abs(CleanedTraitValue - leaf_CN_mean) < 5 * leaf_CN_sd) |>
  ggplot(aes(x = mean_distance_km, y = CleanedTraitValue)) +
  geom_point() +
  geom_smooth() +
  xlim(x_min, x_max) +
  labs(x = "Mean Distance to Biome Boundary (km)", y = "Leaf C:N Ratio") +
  theme_classic()

# Combine plots in one figure
trait_relationships_plots <- plot_grid(plant_height_plot, sla_plot,
                                       seed_mass_plot, leaf_N_plot,
                                       leaf_CN_plot,
                                       labels = c("a)", "b)", "c)", "d)", "e)"),
                                       nrow = 3,
                                       align = "hv",
                                       axis = "tblr")

# END OF SCRIPT ----------------------------------------------------------------