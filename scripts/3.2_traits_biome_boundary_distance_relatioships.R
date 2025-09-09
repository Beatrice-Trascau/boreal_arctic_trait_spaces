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
tundra <- st_make_valid(tundra)

# Re-project biomes to North Pole Lamvert Azimuthal Equal Area (EPSG: 3574)
boreal_sf <- st_transform(boreal_forest, "EPSG:3574")
tundra_sf <- st_transform(tundra, "EPSG:3574")

# 2. ADD EXTRA INFORMATION TO CLEANED TRAITS -----------------------------------

## 2.1. Add species-level biome boundary distances -----------------------------

# Remove Elodea canadensis & hybrids
biome_boundaries <- detailed_results |>
  filter(!species == "Elodea canadensis")

biome_boundaries <- biome_boundaries |>
  filter(!str_detect(species, " × "))

# Combine the dataframes
traits_biome_boundaries <- cleaned_traits |>
  left_join(biome_boundaries, by = c("StandardSpeciesName" = "species"))

# Rename column with mean distance to biome boundary to reflect the fact that
  # it is a species-level mean
traits_biome_boundaries <- traits_biome_boundaries |>
  rename(species_level_mean_distance_km = mean_distance_km)

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
  traits_with_coords$record_level_distance_to_biome_boundary[in_boreal] <- 
    apply(boreal_distances, 1, min) / 1000 # convert to km
}

# Give negative distance to points in tundra
if(sum(in_tundra) > 0){
  tundra_points <- traits_sf[in_tundra, ]
  tundra_distances <- st_distance(tundra_points, tundra_boundary)
  # take minimum distance to boundary
  traits_with_coords$record_level_distance_to_biome_boundary[in_tundra] <- 
    -apply(tundra_distances, 1, min) / 1000 # convert to km
}

# Add biome classification for individual records
traits_with_coords <- traits_with_coords |>
  mutate(individual_record_biome = case_when(in_boreal ~ "boreal",
                                             in_tundra ~ "tundra",
                                             TRUE ~ "other"))

# Update the main dataframe
traits_biome_boundaries <- traits_with_coords

## 2.3. Add LHS strategy for each species --------------------------------------

# Check if we have all required traits
required_traits_lhs <- c("PlantHeight", "SeedMass", "SLA")
available_traits <- unique(traits_biome_boundaries$TraitNameNew)

# Check if any traits are missing
missing_lhs <- setdiff(required_traits_lhs, available_traits)
available_lhs <- intersect(required_traits_lhs, available_traits)

lhs_data <- traits_biome_boundaries |>
  # filter for required traits
  filter(TraitNameNew %in% required_traits_lhs) |>
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
              values_from = MedianTraitValue)

# Rename columns to match MultiTraits requirements
lhs_data <- lhs_data |>
  rename(Height = PlantHeight,
         SeedMass = SeedMass,
         SLA = SLA) |>
  # removed species with missing data for required traits
  filter(!is.na(Height), !is.na(SLA), !is.na(SeedMass))

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
            by = "StandardSpeciesName") |>
  mutate(max_trait_value = max(StdValue, na.rm = TRUE))

# 3. CLEAN GROWTH FORMS --------------------------------------------------------

# Check how many growth forms there are
unique(traits_biome_boundaries$GrowthForm)

# Create function to search by species
species_trait_summary <- function(StandardSpeciesName) {
  traits_biome_boundaries_with_lhs |>
    filter(StandardSpeciesName == !!StandardSpeciesName) |>
    count(TraitNameNew, name = "n_records") |>
    arrange(desc(n_records)) |>
    mutate(species = StandardSpeciesName) |>
    select(species, TraitNameNew, n_records)
}

# Check species summary
species_trait_summary("Alchemilla alpina")
species_trait_summary("Pinus sylvestris")

# Function to get summary for any growth form
growthform_trait_summary <- function(growth_form) {
  traits_biome_boundaries_with_lhs |>
    filter(GrowthForm == !!growth_form) |>
    count(TraitNameNew, name = "n_records") |>
    arrange(desc(n_records)) |>
    mutate(growth_form = growth_form) |>
    select(growth_form, TraitNameNew, n_records)
}

# Usage examples:
growthform_trait_summary("Tree")

# 4. PLOT RELATIONSHIPS --------------------------------------------------------

# Calculate universal x-axis limits
x_min <- min(traits_biome_boundaries_with_lhs$record_level_distance_to_biome_boundary, na.rm = TRUE)
x_max <- max(traits_biome_boundaries_with_lhs$record_level_distance_to_biome_boundary, na.rm = TRUE)

## 4.1. Filter for each trait --------------------------------------------------

# Filter for plant height
plant_height <- traits_biome_boundaries_with_lhs |>
  filter(TraitNameNew == "PlantHeight")

# Filter for SLA
sla <- traits_biome_boundaries_with_lhs |>
  filter(TraitNameNew == "SLA")

# Filter for seed mass
seed_mass <- traits_biome_boundaries_with_lhs |>
  filter(TraitNameNew == "SeedMass")

# Filter for leaf N
leaf_N <- traits_biome_boundaries_with_lhs |>
  filter(TraitNameNew == "LeafN")

# Filter for leaf N:C ratio
leaf_CN_ratio <- traits_biome_boundaries_with_lhs |>
  filter(TraitNameNew == "LeafCN")

## 4.2 Plot all the data -------------------------------------------------------

# Plot plant height
plant_height_all_data <- ggplot(plant_height, aes(x = record_level_distance_to_biome_boundary, 
                                                  y = StdValue)) +
  geom_point() +
  geom_smooth() +
  xlim(x_min, x_max) +
  labs(x = "Distance to Biome Boundary (km)", y = "Plant height (cm)") +
  theme_classic()

# Plot SLA
sla_all_data <- ggplot(sla, aes(x = record_level_distance_to_biome_boundary,
                                y = StdValue)) +
  geom_point() +
  geom_smooth() +
  xlim(x_min, x_max) +
  labs(x = "Distance to Biome Boundary (km)", y = expression("SLA (m"^2*"/kg)")) +
  theme_classic()

# Plot seed mass
seed_mass_all_data <- ggplot(seed_mass, aes(x = record_level_distance_to_biome_boundary,
                                            y = StdValue)) +
  geom_point() +
  geom_smooth() +
  xlim(x_min, x_max) +
  labs(x = "Distance to Biome Boundary (km)", y = "Seed Mass (mg)") +
  theme_classic()

# Plot Leaf N content
leaf_N_all_data <- ggplot(leaf_N, aes(x = record_level_distance_to_biome_boundary,
                                      y = StdValue)) +
  geom_point() +
  geom_smooth() +
  xlim(x_min, x_max) +
  labs(x = "Distance to Biome Boundary (km)", y = "Leaf N (mg/g)") +
  theme_classic()

# Plot Leaf C:N ratio
leaf_CN_ratio_all_data <- ggplot(leaf_CN_ratio, aes(x = record_level_distance_to_biome_boundary,
                                                    y = StdValue)) +
  geom_point() +
  geom_smooth() +
  xlim(x_min, x_max) +
  labs(x = "Distance to Biome Boundary (km)", y = "Leaf C:N Ratio") +
  theme_classic()

# Combine plots 
trait_relationships_all_data <- plot_grid(plant_height_all_data, sla_all_data,
                                          seed_mass_all_data, leaf_N_all_data,
                                          leaf_CN_ratio_all_data,
                                          labels = c("a)", "b)", "c)", "d)", "e)"),
                                          nrow = 3,
                                          align = "hv",
                                          axis = "tblr")

## 4.3. Remove values >= 5 SD --------------------------------------------------

# Calculate means
plant_height_mean <- mean(plant_height$StdValue, na.rm = TRUE)
sla_mean <- mean(sla$StdValue, na.rm = TRUE)
seed_mass_mean <- mean(seed_mass$StdValue, na.rm = TRUE)
leaf_N_mean <- mean(leaf_N$StdValue, na.rm = TRUE)
leaf_CN_mean <- mean(leaf_CN_ratio$StdValue, na.rm = TRUE)

# Calculate SD
plant_height_sd <- sd(plant_height$StdValue, na.rm = TRUE)
sla_sd <- sd(sla$StdValue, na.rm = TRUE)
seed_mass_sd <- sd(seed_mass$StdValue, na.rm = TRUE)
leaf_N_sd <- sd(leaf_N$StdValue, na.rm = TRUE)
leaf_CN_sd <- sd(leaf_CN_ratio$StdValue, na.rm = TRUE)

# Log trait values
plant_height <- plant_height |>
  mutate(LogStdValue = log(StdValue))
sla <- sla |>
  mutate(LogStdValue = log(StdValue))
seed_mass <- seed_mass |>
  mutate(LogStdValue = log(StdValue))
leaf_N <- leaf_N |>
  mutate(LogStdValue = log(StdValue))
leaf_CN_ratio <- leaf_CN_ratio |>
  mutate(LogStdValue = log(StdValue))

# Plant height
plant_height_plot <- plant_height |>
  filter(abs(LogStdValue - plant_height_mean) < 5 * plant_height_sd) |>
  ggplot(aes(x = record_level_distance_to_biome_boundary, y = LogStdValue)) +
  geom_point() +
  geom_smooth() +
  xlim(x_min, x_max) +
  labs(x = "Distance to Biome Boundary (km)", y = "Log Plant height (cm)") +
  theme_classic()

# SLA
sla_plot <- sla |>
  filter(abs(LogStdValue - sla_mean) < 5 * sla_sd) |>
  ggplot(aes(x = record_level_distance_to_biome_boundary, y = LogStdValue)) +
  geom_point() +
  geom_smooth() +
  xlim(x_min, x_max) +
  labs(x = "Distance to Biome Boundary (km)", y = expression("Log SLA (m"^2*"/kg)")) +
  theme_classic()

# Seed Mass
seed_mass_plot <- seed_mass |>
  filter(abs(LogStdValue - seed_mass_mean) < 5 * seed_mass_sd) |>
  ggplot(aes(x = record_level_distance_to_biome_boundary, y = LogStdValue)) +
  geom_point() +
  geom_smooth() +
  xlim(x_min, x_max) +
  labs(x = "Distance to Biome Boundary (km)", y = "Log Seed Mass (mg)") +
  theme_classic()

# Leaf N content
leaf_N_plot <- leaf_N |>
  filter(abs(LogStdValue - leaf_N_mean) < 5 * leaf_N_sd) |>
  ggplot(aes(x = record_level_distance_to_biome_boundary, y = LogStdValue)) +
  geom_point() +
  geom_smooth() +
  xlim(x_min, x_max) +
  labs(x = "Distance to Biome Boundary (km)", y = "Log Leaf N (mg/g)") +
  theme_classic()

# Leaf C:N ratio
leaf_CN_plot <- leaf_CN_ratio |>
  filter(abs(LogStdValue - leaf_CN_mean) < 5 * leaf_CN_sd) |>
  ggplot(aes(x = record_level_distance_to_biome_boundary, y = LogStdValue)) +
  geom_point() +
  geom_smooth() +
  xlim(x_min, x_max) +
  labs(x = "Distance to Biome Boundary (km)", y = "Log Leaf C:N Ratio") +
  theme_classic()

# Combine plots in one figure
trait_relationships_plots <- plot_grid(plant_height_plot, sla_plot,
                                       seed_mass_plot, leaf_N_plot,
                                       leaf_CN_plot,
                                       labels = c("a)", "b)", "c)", "d)", "e)"),
                                       nrow = 3,
                                       align = "hv",
                                       axis = "tblr")

## 4.4. Plot median values per species -----------------------------------------

# Determine most common growth form per species
species_growth_forms <- traits_biome_boundaries_with_lhs |>
  filter(!is.na(GrowthForm)) |>
  group_by(StandardSpeciesName, GrowthForm) |>
  summarise(count = n(), .groups = "drop") |>
  group_by(StandardSpeciesName) |>
  slice_max(count, n = 1, with_ties = FALSE) |>
  select(StandardSpeciesName, GrowthForm)

# Calculate median values per species
traits_median <- traits_biome_boundaries_with_lhs |>
  group_by(StandardSpeciesName, TraitNameNew) |>
  summarise(max_trait_value = max(StdValue, na.rm = TRUE),
            MedianTraitValue = median(StdValue, na.rm = TRUE),
            .groups = "drop")
  
# Add distance to biome boundaries and growth form 
traits_median_df <- traits_median |>
  left_join(species_growth_forms, by = "StandardSpeciesName") |>
  left_join(biome_boundaries, by = c("StandardSpeciesName" = "species")) |>
  rename(species_level_mean_distance_km = mean_distance_km) |>
  filter(!is.na(species_level_mean_distance_km), !is.na(MedianTraitValue),
         !is.na(GrowthForm))

# Calculate universal x-axis limits
x_min_species <- min(traits_median_df$species_level_mean_distance_km, na.rm = TRUE)
x_max_species <- max(traits_median_df$species_level_mean_distance_km, na.rm = TRUE)

# Plant height
median_plant_height <- traits_median_df |>
  filter(TraitNameNew == "PlantHeight") |>
  ggplot(aes(x = species_level_mean_distance_km, y = MedianTraitValue)) +
  geom_point() +
  geom_smooth() +
  xlim(x_min_species, x_max_species) +
  labs(x = "Mean Distance to Biome Boundary (km)", y = "Median Plant height (cm)") +
  theme_classic()

# SLA
median_sla <- traits_median_df |>
  filter(TraitNameNew == "SLA") |>
  ggplot(aes(x = species_level_mean_distance_km, y = MedianTraitValue)) +
  geom_point() +
  geom_smooth() +
  xlim(x_min_species, x_max_species) +
  labs(x = "Mean Distance to Biome Boundary (km)", y = expression("Median SLA (m"^2*"/kg)")) +
  theme_classic()

# Seed mass
median_seed_mass <- traits_median_df |>
  filter(TraitNameNew == "SeedMass") |>
  ggplot(aes(x = species_level_mean_distance_km, y = MedianTraitValue)) +
  geom_point() +
  geom_smooth() +
  xlim(x_min_species, x_max_species) +
  labs(x = "Mean Distance to Biome Boundary (km)", y = "Median Seed Mass (mg)") +
  theme_classic()

# Leaf N content
median_leaf_N <- traits_median_df |>
  filter(TraitNameNew == "LeafN") |>
  ggplot(aes(x = species_level_mean_distance_km, y = MedianTraitValue)) +
  geom_point() +
  geom_smooth() +
  xlim(x_min_species, x_max_species) +
  labs(x = "Mean Distance to Biome Boundary (km)", y = "Median Leaf N (mg/g)") +
  theme_classic()

# Leaf C:N ratio
median_c_n <- traits_median_df |>
  filter(TraitNameNew == "LeafCN") |>
  ggplot(aes(x = species_level_mean_distance_km, y = MedianTraitValue)) +
  geom_point() +
  geom_smooth() +
  xlim(x_min_species, x_max_species) +
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

## 4.5. Median trait value per species per growth form -------------------------

# Get all unique growth forms (excluding NA)
all_growth_forms <- unique(traits_median_df$GrowthForm)
all_growth_forms <- all_growth_forms[!is.na(all_growth_forms)]

# Split into 4 groups (6-7 growth forms each)
panel_groups <- list(Panel1 = all_growth_forms[1:7], 
                     Panel2 = all_growth_forms[8:13],
                     Panel3 = all_growth_forms[14:19],
                     Panel4 = all_growth_forms[20:25])

# Or you could group them more logically by type:
panel_groups_logical <- list(Panel1 = c("Herb", "Graminoid", "Graminoid; Herb", 
                                        "Herb; Shrub", "Herb; Tree", "Herb; Shrub; Tree",
                                        "Hemicryptophyte; Herb; Tree"),
                             Panel2 = c("Tree", "Shrub", "Shrub; Tree", 
                                        "Woody", "Shrub; Woody", "Tree; Woody", 
                                        "Shrub; Tree; Woody"),
                             Panel3 = c("Aquatic", "Aquatic; Herb", 
                                        "Aquatic; Graminoid; Herb", "Fern", 
                                        "Fern; Herb", "Cryptogam"),
                             Panel4 = c("Shrub; Vine", "Herb; Vine", "Herb; Parasite", 
                                        "Parasite", "Bryophyte"))

# Function to create a 4-panel plot for a specific trait
create_trait_growth_form_plot <- function(trait_name, trait_label, use_logical_grouping = TRUE) {
  
  # Choose which grouping to use
  groups <- if(use_logical_grouping) panel_groups_logical else panel_groups
  
  plot_list <- list()
  
  for(i in 1:4) {
    panel_data <- traits_median_df |>
      filter(TraitNameNew == trait_name,
             GrowthForm %in% groups[[i]],
             !is.na(GrowthForm))
    
    plot_list[[i]] <- ggplot(panel_data, aes(x = species_level_mean_distance_km, 
                                             y = MedianTraitValue, 
                                             color = GrowthForm)) +
      geom_point(alpha = 0.7, size = 1.5) +
      geom_smooth(method = "loess", se = FALSE, size = 1) +
      scale_color_viridis_d(name = "Growth Form", option = "turbo") +
      labs(x = "Mean Distance to Biome Boundary (km)", 
           y = trait_label) +
      theme_classic() +
      theme(legend.position = "right",
            legend.text = element_text(size = 7),
            legend.title = element_text(size = 8)) +
      guides(color = guide_legend(ncol = 1, override.aes = list(size = 2)))
  }
  
  # Combine into 4-panel figure
  combined_plot <- plot_grid(plotlist = plot_list, 
                             nrow = 2, ncol = 2,
                             align = "hv")
  
  return(combined_plot)
}

# Create plots for each trait
plant_height_growth_forms <- create_trait_growth_form_plot("PlantHeight", "Median Plant Height (cm)")
sla_growth_forms <- create_trait_growth_form_plot("SLA", expression("Median SLA (m"^2*"/kg)"))
seed_mass_growth_forms <- create_trait_growth_form_plot("SeedMass", "Median Seed Mass (mg)")
leaf_n_growth_forms <- create_trait_growth_form_plot("LeafN", "Median Leaf N (mg/g)")
leaf_cn_growth_forms <- create_trait_growth_form_plot("LeafCN", "Median Leaf C:N Ratio")

# Display the plots
plant_height_growth_forms
sla_growth_forms
seed_mass_growth_forms
leaf_n_growth_forms
leaf_cn_growth_forms

# 5. NMDS ----------------------------------------------------------------------

## 5.1. Prepare data for NMDS --------------------------------------------------

# Categorise species as either boreal or tundra specialists
species_biome_classification <- traits_median_df |>
  distinct(StandardSpeciesName, species_level_mean_distance_km) |>
  mutate(biome_category = case_when(species_level_mean_distance_km > 0 ~ "boreal",
                                    species_level_mean_distance_km < 0 ~ "tundra",
                                    TRUE ~ "boundary"))  # exactly at boundary

# Create trait values per species in wide format for NMDS
trait_matrix <- traits_median_df |>
  filter(StandardSpeciesName %in% species_biome_classification$StandardSpeciesName) |>
  select(StandardSpeciesName, TraitNameNew, MedianTraitValue) |>
  pivot_wider(names_from = TraitNameNew, 
              values_from = MedianTraitValue) |>
  column_to_rownames("StandardSpeciesName")

# Create table of species that are missing values
missing_plant_height <- trait_matrix |>
  filter(is.na(PlantHeight))
missing_leaf_CN <- trait_matrix |>
  filter(is.na(LeafCN))
missing_leaf_N <- trait_matrix |>
  filter(is.na(LeafN))
missing_SLA <- trait_matrix |>
  filter(is.na(SLA))
missing_seed_mass <- trait_matrix |>
  filter(is.na(SeedMass))

# Remove species with missing data for any trait
complete_trait_matrix <- trait_matrix[complete.cases(trait_matrix), ]

# Check which traits are available
colnames(complete_trait_matrix)

## 5.2. Run NMDS ---------------------------------------------------------------

# Set seed for NMDS
set.seed(53135)

# Run NMDS
nmds_result <- metaMDS(complete_trait_matrix, 
                       distance = "euclidean",
                       k = 3,
                       trymax = 100)

# Check stress value (should be < 0.2, ideally < 0.1)
nmds_result$stress # 0.04216575

# Check stress value per dimension
dimcheck_out <- dimcheckMDS(complete_trait_matrix,
                            distance = "bray",
                            k = 5)

## 5.3. Plot output ------------------------------------------------------------

# Extract NMDS scores
nmds_scores <- as.data.frame(nmds_result$points)
nmds_scores$StandardSpeciesName <- rownames(nmds_scores)

# Add biome classification
nmds_plot_data <- nmds_scores |>
  left_join(species_biome_classification, by = "StandardSpeciesName") |>
  filter(!is.na(biome_category))

# Create polygon hull
boreal_scores <- nmds_plot_data[nmds_plot_data$biome_category == "boreal", ][chull(nmds_plot_data[nmds_plot_data$biome_category == 
                                                                           "boreal", c("MDS1", "MDS2")]), ]  # hull values for boreal
tundra_scores <- nmds_plot_data[nmds_plot_data$biome_category == "tundra", ][chull(nmds_plot_data[nmds_plot_data$biome_category == 
                                                                   "tundra", c("MDS1", "MDS2")]), ]  # hull values for grp B

hull_data <- rbind(boreal_scores, tundra_scores)  #combine grp.a and grp.b
hull.data

# Create NMDS plot
nmds_plot <- ggplot(nmds_plot_data, aes(x = MDS1, y = MDS2, color = biome_category)) +
  geom_polygon(data=hull_data,aes(x=MDS1,y=MDS2,fill=biome_category,group=biome_category),alpha=0.30) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = c("boreal" = "darkgreen", "tundra" = "skyblue"),
                     name = "Biome") +
  labs(x = "NMDS1", y = "NMDS2") +
  theme_classic() +
  theme(legend.position = "bottom") +
  # add stress value as subtitle
  labs(subtitle = paste("Based on", ncol(trait_matrix), "traits,", nrow(nmds_plot_data), "species"))

# Add ellipses around each group
nmds_plot_with_hulls <- nmds_plot +
  stat_ellipse(aes(fill = biome_category), alpha = 0.2, geom = "polygon",
               level = 0.95) +
  scale_fill_manual(values = c("boreal" = "darkgreen", "tundra" = "skyblue"),
                    guide = "none")

## 5.4. Plot output with trait vectors -----------------------------------------

# Fit trait vectors to the NMDS ordination
trait_fit <- envfit(nmds_result, trait_matrix, permutations = 999)

# Extract the vectors
trait_vectors <- as.data.frame(scores(trait_fit, "vectors"))
trait_vectors$trait <- rownames(trait_vectors)

# Check significance of vectors
print(trait_fit)

# Plot with proper trait vectors (no arbitrary scaling)
nmds_plot_with_proper_vectors <- nmds_plot_with_hulls +
  geom_segment(data = trait_vectors, 
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "black", 
               size = 1,
               inherit.aes = FALSE) +
  geom_text(data = trait_vectors,
            aes(x = NMDS1 * 1.1, y = NMDS2 * 1.1, label = trait),
            color = "black", 
            size = 3, 
            fontface = "bold",
            inherit.aes = FALSE)

## 5.5. Statistical tests ------------------------------------------------------

# Test for significant differences between groups using PERMANOVA
permanova_result <- adonis2(trait_matrix ~ biome_category, 
                            data = species_biome_classification[species_biome_classification$StandardSpeciesName %in% rownames(trait_matrix), ])

print(permanova_result)

# Test for homogeneity of dispersions
dispersion_test <- betadisper(vegdist(trait_matrix), 
                              species_biome_classification$biome_category[species_biome_classification$StandardSpeciesName %in% rownames(trait_matrix)])
permutest(dispersion_test)

# 6. GLLVM? --------------------------------------------------------------------

# Fit gllvm with trait_matrix and biome classifications
gllvm_model <- gllvm(trait_matrix, 
                     X = data.frame(biome = nmds_plot_data$biome_category),
                     num.lv = 2,
                     family = "gaussian",
                     seed = 52164)

# Compare to NMDS stress
gllvm_model$logL

# Extract and plot results
ordiplot(gllvm_model, biplot = TRUE)
coefplot(gllvm_model)  

# END OF SCRIPT ----------------------------------------------------------------