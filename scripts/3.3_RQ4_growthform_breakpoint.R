##----------------------------------------------------------------------------##
# PAPER 3: BOREAL AND ARCTIC PLANT SPECIES TRAIT SPACES 
# 3.3_RQ4_growthform_breakpoint
# This script contains code for breakpoint regression analysis testing whether
# trait-distance relationships differ on either side of the biome boundary
# broken down by growth form
##----------------------------------------------------------------------------##

# 1. SETUP ---------------------------------------------------------------------

## 1.1. Load data --------------------------------------------------------------

# Load packages and functions
library(here)
source(here("scripts", "0_setup.R"))

# Load cleaned traits
load(here("data", "derived_data", "TRY_traits_cleaned_July2025.RData"))
cleaned_traits <- cleaned_traits_July2025

# Load distances to biome boundaries
detailed_results <- readRDS(here("data", "derived_data", 
                                 "species_summaries_dist_to_biome_boundary_June25.rds"))

# Load biome boundaries
global_biomes <- st_read(here("data", "raw_data", "biomes", 
                              "wwf_terr_ecos.shp"))

# Load CAFF quality check
caff_check <- read.xlsx(here("data", "derived_data", "caff_quality_check.xlsx"),
                        sheet = 1, skipEmptyRows = TRUE)

# Cleaned growth forms 
# Validated at workshop by Mariana & James
validated_growthforms <- read.xlsx(here("data", "derived_data",
                                        "BIEN_growth_forms_corrected_at_workshop.xlsx"),
                                   sheet = 1, skipEmptyRows = TRUE)

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

# Get boundaries for distance calculations
boreal_boundary <- st_boundary(boreal_sf)
tundra_boundary <- st_boundary(tundra_sf)

# 2. ADD EXTRA INFORMATION TO CLEANED TRAITS -----------------------------------

## 2.1. Add species-level biome boundary distances -----------------------------

# Remove Elodea canadensis & hybrids
biome_boundaries <- detailed_results |>
  filter(!species == "Elodea canadensis") |>
  filter(!str_detect(species, " × "))

# REVERSE THE SIGN: multiply by -1 to follow new convention
# Original: positive = boreal, negative = tundra
# New: negative = boreal, positive = tundra
biome_boundaries <- biome_boundaries |>
  mutate(mean_distance_km = -mean_distance_km,
         median_distance_km = -median_distance_km,
         min_distance_km = -min_distance_km,
         max_distance_km = -max_distance_km,
         mean_distance_boreal_km = -mean_distance_boreal_km,
         median_distance_boreal_km = -median_distance_boreal_km,
         mean_distance_tundra_km = -mean_distance_tundra_km,
         median_distance_tundra_km = -median_distance_tundra_km,
         weighted_mean_distance_km = -weighted_mean_distance_km)

# Combine the dataframes
traits_biome_boundaries <- cleaned_traits |>
  left_join(biome_boundaries, by = c("StandardSpeciesName" = "species"))

# Rename column with mean distance to biome boundary to reflect the fact that
# it is a species-level mean
traits_biome_boundaries <- traits_biome_boundaries |>
  rename(species_level_mean_distance_km = mean_distance_km)

## 2.2. Calculate distance to biome boundary for each trait record -------------
## WITH REVERSED SIGN CONVENTION

# Filter to records with valid coordinates
traits_with_coords <- traits_biome_boundaries |>
  filter(!is.na(LON_site), !is.na(LAT_site)) |>
  filter(LON_site >= -180, LON_site <= 180,
         LAT_site >= -90, LAT_site <= 90)

# Check how many records were removed
nrow(traits_biome_boundaries)
nrow(traits_with_coords)

# Convert traits to spatial objec5
traits_sf <- st_as_sf(traits_with_coords, 
                      coords = c("LON_site", "LAT_site"), 
                      crs = 4326) |>
  st_transform(crs = "EPSG:3574")

# Check which biome each point is in
in_boreal <- st_intersects(traits_sf, boreal_sf, sparse = FALSE)[,1]
in_tundra <- st_intersects(traits_sf, tundra_sf, sparse = FALSE)[,1]

# Initialize distance column
traits_with_coords$distance_to_boundary_km <- NA

# REVERSED CONVENTION:
# Calculate distances for boreal points (NEGATIVE)
if(sum(in_boreal) > 0) {
  boreal_points <- traits_sf[in_boreal, ]
  boreal_distances <- st_distance(boreal_points, boreal_boundary)
  traits_with_coords$distance_to_boundary_km[in_boreal] <- 
    -apply(boreal_distances, 1, min) / 1000  # NEGATIVE for boreal
}

# Calculate distances for tundra points (POSITIVE)
if(sum(in_tundra) > 0) {
  tundra_points <- traits_sf[in_tundra, ]
  tundra_distances <- st_distance(tundra_points, tundra_boundary)
  traits_with_coords$distance_to_boundary_km[in_tundra] <- 
    apply(tundra_distances, 1, min) / 1000  # POSITIVE for tundra
}

# Add biome classification
traits_with_coords$biome <- NA
traits_with_coords$biome[in_boreal] <- "boreal"
traits_with_coords$biome[in_tundra] <- "tundra"

# Create site identifier
traits_with_coords <- traits_with_coords |>
  mutate(site_name = paste0(round(LON_site, 3), "_", round(LAT_site, 3)))

# Check how many records are in boreal and tundra biome
sum(traits_with_coords$biome == "boreal", na.rm = TRUE)
sum(traits_with_coords$biome == "tundra", na.rm = TRUE)

## 2.3. CAFF quality control ---------------------------------------------------

# Rename species name column in CAFF check
caff_check <- caff_check |>
  rename(StandardSpeciesName = SPECIES_CLEAN) |>
  dplyr::select(StandardSpeciesName, final.category)

# Add CAFF classification (keep all columns from traits_with_coords)
traits_with_caff <- traits_with_coords |>
  left_join(caff_check, by = "StandardSpeciesName")

# Filter to valid species only
cleaned_traits_final <- traits_with_caff |>
  mutate(caff_biome_category = final.category) |>
  filter(!is.na(caff_biome_category)) |>
  filter(caff_biome_category != "remove") |>
  filter(!is.na(distance_to_boundary_km)) |>
  filter(!is.na(biome))

# Check how many species there were originally
length(unique(traits_with_coords$StandardSpeciesName))

# Check how many species are left after CAFF filtering
length(unique(cleaned_traits_final$StandardSpeciesName))

# Check the total number of records after filtering
nrow(cleaned_traits_final)

# Create piecewise distance variables with REVERSED convention
# Now: boreal_side uses negative distances, tundra_side uses positive distances
cleaned_traits_final <- cleaned_traits_final |>
  mutate(distance_boreal_side = pmin(distance_to_boundary_km, 0),  # negative values (boreal)
         distance_tundra_side = pmax(distance_to_boundary_km, 0))  # positive values (tundra)

## 2.4. Add cleaned growth forms -----------------------------------------------

# Check how many growth forms there are (and check that you don't have different spellings for the same thing)
unique(validated_growthforms$`James'.confirmation`)

# Fix column names before merging
growthforms <- validated_growthforms |>
  clean_names() |>
  rename(StandardSpeciesName = standard_species_name) |>
  dplyr::select(-x1)

# Check which values are left

# Add to the dataframe
