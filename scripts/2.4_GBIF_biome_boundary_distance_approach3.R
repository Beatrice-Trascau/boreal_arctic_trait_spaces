##----------------------------------------------------------------------------##
# PAPER 3: BOREAL AND ARCTIC PLANT SPECIES TRAIT SPACES 
# 2.4_GBIF_biome_boundary_distance_approach3
# This script contains code which calcualtes the distance to biome boundary for
# each species in our list with a third approach
##----------------------------------------------------------------------------##

# 1. SETUP ---------------------------------------------------------------------

## 1.1. Libraries --------------------------------------------------------------

# Load packages and functions
library(here)
source(here("scripts", "0_setup.R"))

## 1.2. GBIF credentials -------------------------------------------------------

# Add your details
# options(gbif_user = "your_username")
# options(gbif_pwd = "your_password") 
# options(gbif_email = "your_email")

## 1.3. Species data -----------------------------------------------------------

# Load species list 
load(here("data", "derived_data", "all_filtered_standardised_species.RData"))

# Extract species list 
species_list <- unique(corrected_species_list$CheckedSpeciesName)

## 1.4. Load and prepare biomes ------------------------------------------------

# Load WWF biomes
# Citation: Olson, D. M., Dinerstein, E., Wikramanayake, E. D., Burgess, N. D., Powell, G. V. N., Underwood, E. C., D'Amico, J. A., Itoua, I., Strand, H. E., Morrison, J. C., Loucks, C. J., Allnutt, T. F., Ricketts, T. H., Kura, Y., Lamoreux, J. F., Wettengel, W. W., Hedao, P., Kassem, K. R. 2001. Terrestrial ecoregions of the world: a new map of life on Earth. Bioscience 51(11):933-938.
global_biomes <- st_read(here("data", "raw_data", "biomes", "wwf_terr_ecos.shp"))

# Load Boreal Forest (BIOME = 6)
boreal_forest <- st_union(global_biomes[global_biomes$BIOME == 6,])

# Load Tundra (BIOME = 11) - Palearctic and Nearctic only
tundra <- st_union(global_biomes[global_biomes$BIOME == 11 & 
                                   (global_biomes$REALM == "PA" | global_biomes$REALM == "NA"), ])

# Make sure geometries are valid
boreal_forest <- st_make_valid(boreal_forest)
tundra <- st_make_valid(tundra)

# Transform to WGS84 for GBIF compatibility
boreal_forest <- st_transform(boreal_forest, "EPSG:4326")
tundra <- st_transform(tundra, "EPSG:4326")

## 1.4. Create analysis grid ---------------------------------------------------

# Combine biomes to use for the grid extent
combined_biomes <- st_union(boreal_forest, tundra)

# Create bounding box from the combined biomes
combined_extent <- st_bbox(combined_biomes)

# Create grid with 0.45 degrees resolution (~50km at high latitudes)
grid <- rast(extent = combined_extent, 
             resolution = c(0.45, 0.45),
             crs = "EPSG:4326")

# Convert to polygons
polygrid <- as.polygons(grid)

# Give each cell and ID
polygrid$cell_id <- 1:nrow(polygrid)

# Convert to sf object
polygrid_sf <- st_as_sf(polygrid)

# Filter cells within boreal biome
cells_in_boreal <- st_intersects(polygrid_sf, boreal_forest, sparse = FALSE)[,1]

# Filter cells within tundra biome
cells_in_tundra <- st_intersects(polygrid_sf, tundra, sparse = FALSE)[,1]

# Get cells in both biomes
cells_in_biomes <- cells_in_boreal | cells_in_tundra

# Keep only cells within biome
polygrid_filtered <- polygrid_sf[cells_in_biomes, ]

# Get boreal cells
polygrid_filtered$in_boreal <- cells_in_boreal[cells_in_biomes]

# Get tundra cells
polygrid_filtered$in_tundra <- cells_in_tundra[cells_in_biomes]

# Summary
cat("Grid created:", nrow(polygrid_filtered), "cells within biomes\n")
cat("  Boreal cells:", sum(polygrid_filtered$in_boreal), "\n")
cat("  Tundra cells:", sum(polygrid_filtered$in_tundra), "\n")









