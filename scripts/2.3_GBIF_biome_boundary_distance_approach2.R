##----------------------------------------------------------------------------##
# PAPER 3: BOREAL AND ARCTIC PLANT SPECIES TRAIT SPACES 
# 2.3_GBIF_biome_boundary_distance_approach2
# This script contains code which calcualtes the distance to biome boundary for
# each species in our list with a second approach
##----------------------------------------------------------------------------##

# 1. SETUP ---------------------------------------------------------------------

## 1.1. Species data -----------------------------------------------------------

# Load packages and functions
library(here)
source(here("scripts", "0_setup.R"))

# Load species list 
load(here("data", "derived_data", "all_filtered_standardised_species.RData"))

# Extract species list 
species_list <- unique(corrected_species_list$CheckedSpeciesName)

## 1.2. WWF biomes -------------------------------------------------------------

# Load WWF biomes
# Citation: Olson, D. M., Dinerstein, E., Wikramanayake, E. D., Burgess, N. D., Powell, G. V. N., Underwood, E. C., D'Amico, J. A., Itoua, I., Strand, H. E., Morrison, J. C., Loucks, C. J., Allnutt, T. F., Ricketts, T. H., Kura, Y., Lamoreux, J. F., Wettengel, W. W., Hedao, P., Kassem, K. R. 2001. Terrestrial ecoregions of the world: a new map of life on Earth. Bioscience 51(11):933-938.
global_biomes <- st_read (here("data", "raw_data", "biomes", "wwf_terr_ecos.shp"))

#Load Boreal Forest (BIOME = 6)
boreal_forest <- st_union(global_biomes[global_biomes$BIOME == 6,])

# Load Tundra (BIOME = 11)
tundra <- st_union(global_biomes[global_biomes$BIOME == 11 &(global_biomes$REALM == "PA"|global_biomes$REALM == "NA"), ])

# Make sure boreal_forest and tundra are valid geometries
boreal_forest <- st_make_valid(boreal_forest)
tundra <- st_make_valid(tundra)

# Check CRS
biome_crs <- st_crs(boreal_forest)

## 1.3. Create 1km x 1km grid --------------------------------------------------

# Combine biomes for grid extent
combined_biomes <- st_union(boreal_forest, tundra)

# Create bounding box for grid creation
bbox <- st_bbox(combined_biomes)

# Create grid with 1km resolution
# Note: grid resolution depends on CRS units. Assuming meters for projected CRS
grid_resoluton <- 1000

# Create raster template
grid_raster <- rast(xmin = bbox["xmin"], xmax = bbox["xmax"],
                    ymin = bbox["ymin"], ymax = bbox["ymax"],
                    resolution = grid_resolution,
                    crs = st_crs(combined_biomes)$wkt)

# Convert to polygons and then to sf
grid_polygons <- as.polygons(grid_raster)
grid_sf <- st_as_sf(grid_polygons)

# Add cell ID
grid_sf$cell_id <- 1:nrow(grid_sf)








