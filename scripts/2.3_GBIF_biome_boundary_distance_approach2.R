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
grid_resoluton <- 10000

# Create raster template
grid_raster <- rast(xmin = bbox["xmin"], xmax = bbox["xmax"],
                    ymin = bbox["ymin"], ymax = bbox["ymax"],
                    resolution = grid_resoluton,
                    crs = st_crs(combined_biomes)$wkt)

# Convert to polygons and then to sf
grid_polygons <- as.polygons(grid_raster)
grid_sf <- st_as_sf(grid_polygons)

# Add cell ID
grid_sf$cell_id <- 1:nrow(grid_sf)

## 1.4. Filter grid to biome areas and calculate distances ---------------------

# Keep only cells that intersect with either biome
grid_intersects <- st_intersects(grid_sf, combined_biomes)
grid_in_biomes <- grid_sf[lengths(grid_intersects) > 0, ]

# Check how many cells there are within biomes
cat("Grid cells within biomes:", nrow(grid_in_biomes), "\n")

# Get centroids for the biome boundary distance calculations
grid_centroids <- st_centroid(grid_in_biomes)

# Determine which biome each cell belongs to
boreal_cells <- st_intersects(grid_centroids, boreal_forest)
tundra_cells <- st_intersects(grid_centroids, tundra)

# Create biome classification for each cell
grid_in_biomes$biome <- "mixed" # default
grid_in_biomes$biome[lengths(boreal_cells) > 0 & lengths(tundra_cells) == 0] <- "boreal"
grid_in_biomes$biome[lengths(tundra_cells) > 0 & lengths(boreal_cells) == 0] <- "tundra"

# Calculate distances to biome boundaries
grid_in_biomes$distance_to_boundary <- NA

# Boreal cells: +ve distance to biome boundary
boreal_mask <- grid_in_biomes$biome == "boreal"
if(sum(boreal_mask) > 0){
  boreal_centroids <- st_centroid(grid_in_biomes[boreal_mask, ])
  boreal_distances <- as.numeric(st_distance(boreal_centroids, st_boundary(tundra)))
  grid_in_biomes$distance_to_boundary[boreal_mask] <- boreal_distances
}

# Tundra cells: -ve distance to boreal boundary
tundra_mask <- grid_in_biomes$biome == "tundra"
if(sum(tundra_mask) > 0){
  tundra_centroids <- st_centroid(grid_in_biomes[tundra_mask, ])
  tundra_distances <- as.numeric(st_distance(tundra_centroids, st_boundary(boreal_forest)))
  grid_in_biomes$distance_to_boundary[tundra_mask] <- -tundra_distances
}

# Remove mixed cells
grid_final <- grid_in_biomes[grid_in_biomes$biome %in% c("boreal", "tundra"), ]

# Quick summary of cells
cat("Final grid cells after removing mixed cells:", nrow(grid_final), "\n")
cat("Boreal cells:", sum(grid_final$biome == "boreal"), "\n")
cat("Tundra cells:", sum(grid_final$biome == "tundra"), "\n")

# Convert to lat/lon for GBIF queries
grid_wgs84 <- st_transform(grid_final, crs = 4326)
grid_centroids_wgs84 <- st_centroid(grid_wgs84)
grid_coords <- st_coordinates(grid_centroids_wgs84)

# Handle case where coordinates don't have column names
if(is.null(colnames(grid_coords))) {
  # Assume standard longitude (X) and latitude (Y)
  colnames(grid_coords) <- c("X", "Y")
}

# Extract lat & long for final grid
grid_final$longitude <- grid_coords[, "X"]
grid_final$latitude  <- grid_coords[, "Y"]

## 1.5. Prepare df to store results --------------------------------------------

# Create a dataframe to store results in
results_df <- data.frame(species_name = character(),
                         avg_distance_to_boundary = numeric(),
                         classification = character(),
                         n_records_used = integer(),
                         stringsAsFactors = FALSE)

# Create tracking variables
succesful_species <- character()
failed_species <- character()
n_total_species <- length(species_list)

# 2. PROCESS EACH SPECIES ------------------------------------------------------

for(i in seq_along(species_list)){
  species_name <- species_list[i]
  cat("Processing species", i, "of", n_total_species, ":", species_name, "\n")
  
  tryCatch({
    ## 2.1. Get occurrence counts for each grid cell ---------------------------
    
    # Initialise presence vector
    cell_presence <- rep(0, nrow(grid_final))
    
    # Query each cell for occurrences
    for(j in 1: nrow(grid_final)){
      
      # Ensure coordinates are correctly formatted as numbers
      lat_center <- as.numeric(grid_final$latitude[j])
      lon_center <- as.numeric(grid_final$longitude[j])
      
      # Create coordinate bounds 
      lat_min <- lat_center - 0.005
      lat_max <- lat_center + 0.005
      lon_min <- lon_center - 0.005
      lon_max <- lon_center + 0.005
      
      # Create proper WKT polygon format for the bounding box
      wkt_polygon <- paste0("POLYGON((", 
                            lon_min, " ", lat_min, ",",
                            lon_max, " ", lat_min, ",", 
                            lon_max, " ", lat_max, ",",
                            lon_min, " ", lat_max, ",",
                            lon_min, " ", lat_min, "))")
      
      # Count occurrences
      cell_count <- occ_count(scientificName = species_name,
                              hasCoordinate = TRUE,
                              coordinateUncertaintyInMeters = "0,1000",
                              geometry = wkt_polygon)
      
      # Mark as present if any occurrences are found
      if(cell_count > 0) {
        cell_presence[j] <- 1
      }
      
      # Progress indicator for large grids
      if(j %% 1000 == 0) {
        cat("    Processed", j, "of", nrow(grid_final), "cells\n")
      }
    }
    
    # Check if species found in any cells
    n_occupied_cells <- sum(cell_presence)
    if(n_occupied_cells == 0) {
      cat("  No occurrences found for", species_name, "\n")
      failed_species <- c(failed_species, paste(species_name, "(no occurrences)"))
      next
    }
    
    ## 2.2. Calculate median distance across occupied cells --------------------
    
    # Get distances for occupied cells only
    occupied_distances <- grid_final$distance_to_boundary[cell_presence == 1]
    
    # Remove NA distances
    occupied_distances <- occupied_distances[!is.na(occupied_distances)]
    
    # Check if there are only NA distances
    if(length(occupied_distances) == 0) {
      cat("  No valid distances calculated for", species_name, "\n")
      failed_species <- c(failed_species, paste(species_name, "(no valid distances)"))
      next
    }
    
    # Calculate median distance
    median_distance <- median(occupied_distances, na.rm = TRUE)
    
    # Classify species
    classification <- ifelse(median_distance > 0, "boreal", "arctic")
    
    # Store results
    new_row <- data.frame(species_name = species_name,
                          avg_distance_to_boundary = median_distance,
                          classification = classification,
                          n_records_used = length(occupied_distances),
                          stringsAsFactors = FALSE)
    
    # Add to results df
    results_df <- rbind(results_df, new_row)
    
    # Store successful species
    succesful_species <- c(succesful_species, species_name)
    
  }, error = function(e){
    cat("  Error processing", species_name, ":", e$message, "\n")
    failed_species <- c(failed_species, paste(species_name, "(error:", e$message, ")"))
  })
  
  cat("\n")
}

# 3. FINAL SUMMARY -------------------------------------------------------------

# Quick summary of successful and failed species
cat("=== PROCESSING COMPLETE ===\n")
cat("Total species processed:", n_total_species, "\n")
cat("Successful species:", length(succesful_species), "\n")
cat("Failed species:", length(failed_species), "\n")

# Inspect reason for failure in failed species
if(length(failed_species) > 0) {
  cat("\nFailed species: ")
  for(failed in failed_species){
    cat("  -", failed, "\n")
  }
}

# Save results
save(results_df, file = here("data", "derived_data", 
                             "species_biome_classification_approach2.RData"))

# END OF SCRIPT ----------------------------------------------------------------

