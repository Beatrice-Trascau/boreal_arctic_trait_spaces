##----------------------------------------------------------------------------##
# PAPER 3: BOREAL AND ARCTIC PLANT SPECIES TRAIT SPACES 
# 3.2_traits_biome_boundary_distance_relationships
# This script contains code which explores ther relationships between species
# trait values and their mean distance to biome boundaries
##----------------------------------------------------------------------------##

# 1. PREPARE DATA --------------------------------------------------------------

## 1.1. Load data --------------------------------------------------------------

# Load cleaned traits
load(here("data", "derived_data", "TRY_traits_cleaned_species_names.RData"))
cleaned_traits <- traits_cleaned_species_names

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

## 2.1. Clean separate df for Plant Height

# Filter for plant height
plant_height <- traits_biome_boundaries |>
  filter(TraitNameNew == "PlantHeight")

# Check df
glimpse(plant_height) # trait values considered characters
unique(plant_height$OrigUnitStr) # and 4 different measurements for plant height....... (m, cm, meter, mm)

# Convert everything to cm
plant_height_fixed <- plant_height |>
  mutate(OrigValueStr = as.numeric(OrigValueStr),
         FixedTraitValue = case_when(OrigUnitStr %in% c("m", "meter") ~ OrigValueStr * 100,
                                     OrigUnitStr == "cm" ~ OrigValueStr,
                                     OrigUnitStr == "mm" ~ OrigValueStr / 10,
                                     .default = OrigValueStr))

b <- plant_height |>
  filter(ObservationID %in% c("1322294", "1322298", "1322299"))


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
leaf_NC_ratio <- traits_biome_boundaries |>
  filter(TraitNameNew == "LeafCN")

# Plot plant height
plant_height_all_data <- ggplot(plant_height, aes(x = median_distance_km, y = TraitValue)) +
  geom_point() +
  geom_smooth() +
  labs(x = "Distance to Biome Boundary (km)", y = "Plant height ")