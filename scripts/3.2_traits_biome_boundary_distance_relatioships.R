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


## 1.2. Combine traits and biome boundary distances ----------------------------

# Remove Elodea canadensis & hybrids
biome_boundaries <- detailed_results |>
  filter(!species == "Elodea canadensis",
         !grepl("×", species))

# Combine the dataframes
traits_biome_boundaries <- cleaned_traits |>
  left_join(biome_boundaries, by = c("StandardSpeciesName" = "species"))

# 2. PLOT RELATIONSHIPS --------------------------------------------------------  

## 2.1. Filter for each trait --------------------------------------------------

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
  labs(x = "Mean Distance to Biome Boundary (km)", y = "Plant height (cm)") +
  theme_classic()

# Plot SLA
sla_all_data <- ggplot(sla, aes(x = mean_distance_km,
                                y = CleanedTraitValue)) +
  geom_point() +
  geom_smooth() +
  labs(x = "Mean Distance to Biome Boundary (km)", y = expression("SLA (m"^2*"/kg)")) +
  theme_classic()

# Plot seed mass
seed_mass_all_data <- ggplot(seed_mass, aes(x = mean_distance_km,
                                            y = CleanedTraitValue)) +
  geom_point() +
  geom_smooth() +
  labs(x = "Mean Distance to Biome Boundary (km)", y = "Seed Mass (mg)") +
  theme_classic()

# Plot Leaf N content
leaf_N_all_data <- ggplot(leaf_N, aes(x = mean_distance_km,
                                      y = CleanedTraitValue)) +
  geom_point() +
  geom_smooth() +
  labs(x = "Distance to Biome Boundary (km)", y = "Leaf N (mg/g)") +
  theme_classic()

# Plot Leaf C:N ratio
leaf_CN_ratio_all_data <- ggplot(leaf_CN_ratio, aes(x = mean_distance_km,
                                                    y = CleanedTraitValue)) +
  geom_point() +
  geom_smooth() +
  labs(x = "Distance to Biome Boundary (km)", y = "Leaf C:N Ratio") +
  theme_classic()

# Combine plots 
trait_relationships_all_data <- plot_grid(plant_height_all_data, sla_all_data,
                                          seed_mass_all_data, leaf_N_all_data,
                                          leaf_CN_ratio_all_data,
                                          labels = c("a)", "b)", "c)", "d)", "e)"),
                                          ncols = 2)

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
  labs(x = "Mean Distance to Biome Boundary (km)", y = "Median Plant height (cm)") +
  theme_classic()

# SLA
median_sla <- traits_median_df |>
  filter(TraitNameNew == "SLA") |>
  ggplot(aes(x = mean_distance_km, y = MedianTraitValue)) +
  geom_point() +
  geom_smooth() +
  labs(x = "Mean Distance to Biome Boundary (km)", y = expression("Median SLA (m"^2*"/kg)")) +
  theme_classic()

# Seed mass
median_seed_mass <- traits_median_df |>
  filter(TraitNameNew == "SeedMass") |>
  ggplot(aes(x = mean_distance_km, y = MedianTraitValue)) +
  geom_point() +
  geom_smooth() +
  labs(x = "Mean Distance to Biome Boundary (km)", y = "Median Seed Mass (mg)") +
  theme_classic()

# Leaf N content
median_leaf_N <- traits_median_df |>
  filter(TraitNameNew == "LeafN") |>
  ggplot(aes(x = mean_distance_km, y = MedianTraitValue)) +
  geom_point() +
  geom_smooth() +
  labs(x = "Mean Distance to Biome Boundary (km)", y = "Median Leaf N (mg/g)") +
  theme_classic()

# Leaf C:N ratio
median_c_n <- traits_median_df |>
  filter(TraitNameNew == "LeafCN") |>
  ggplot(aes(x = mean_distance_km, y = MedianTraitValue)) +
  geom_point() +
  geom_smooth() +
  labs(x = "Mean Distance to Biome Boundary (km)", y = "Median Leaf C:N Ratio") +
  theme_classic()

# Combine plots 
median_trait_relationships <- plot_grid(median_plant_height, median_sla,
                                        median_seed_mass, median_leaf_N,
                                        median_c_n,
                                        labels = c("a)", "b)", "c)", "d)", "e)"),
                                        ncols = 2)

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
  labs(x = "Mean Distance to Biome Boundary (km)", y = "Plant height (cm)") +
  theme_classic()

# SLA
sla_plot <- sla |>
  filter(abs(CleanedTraitValue - sla_mean) < 5 * sla_sd) |>
  ggplot(aes(x = mean_distance_km, y = CleanedTraitValue)) +
  geom_point() +
  geom_smooth() +
  labs(x = "Mean Distance to Biome Boundary (km)", y = expression("SLA (m"^2*"/kg)")) +
  theme_classic()

# Seed Mass
seed_mass_plot <- seed_mass |>
  filter(abs(CleanedTraitValue - seed_mass_mean) < 5 * seed_mass_sd) |>
  ggplot(aes(x = mean_distance_km, y = CleanedTraitValue)) +
  geom_point() +
  geom_smooth() +
  labs(x = "Mean Distance to Biome Boundary (km)", y = "Seed Mass (mg)") +
  theme_classic()

# Leaf N content
leaf_N_plot <- leaf_N |>
  filter(abs(CleanedTraitValue - leaf_N_mean) < 5 * leaf_N_sd) |>
  ggplot(aes(x = mean_distance_km, y = CleanedTraitValue)) +
  geom_point() +
  geom_smooth() +
  labs(x = "Mean Distance to Biome Boundary (km)", y = "Leaf N (mg/g)") +
  theme_classic()

# Leaf C:N ratio
leaf_CN_plot <- leaf_CN_ratio |>
  filter(abs(CleanedTraitValue - leaf_CN_mean) < 5 * leaf_CN_sd) |>
  ggplot(aes(x = mean_distance_km, y = CleanedTraitValue)) +
  geom_point() +
  geom_smooth() +
  labs(x = "Mean Distance to Biome Boundary (km)", y = "Leaf C:N Ratio") +
  theme_classic()

# Combine plots in one figure
trait_relationships_plots <- plot_grid(plant_height_plot, sla_plot,
                                          seed_mass_plot, leaf_N_plot,
                                          leaf_CN_plot,
                                          labels = c("a)", "b)", "c)", "d)", "e)"),
                                          ncols = 2)

# END OF SCRIPT ----------------------------------------------------------------