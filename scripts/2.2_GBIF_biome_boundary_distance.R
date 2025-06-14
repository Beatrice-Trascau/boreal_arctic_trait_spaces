##----------------------------------------------------------------------------##
# PAPER 3: BOREAL AND ARCTIC PLANT SPECIES TRAIT SPACES 
# 2.2_GBIF_biome_boundary_distance
# This script contains code which calcualtes the distance to biome boundary for
# each species in our list
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

## 1.3. Create geographic filter for GBIF search -------------------------------

# Get the souther boundary (minimum latitude) of boreal forest
boreal_bbox <- st_bbox(st_transform(boreal_forest, crs = 4326))
min_latitude <- boreal_bbox["ymin"]

# Create simple boundy box WKT that covers northern regions (from boreal souther boundary to 90°N across all longitudes)
northern_region_wkt <- sprintf("POLYGON((-180 %.6f, 180 %.6f, 180 90, -180 90, -180 %.6f))", 
                               min_latitude, min_latitude, min_latitude)

## 1.4. Prepare df to store results --------------------------------------------

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
    ## 2.1. Get occurrence data from GBIF --------------------------------------
    cat("  Downloading GBIF data...\n")
    
    # Set up search in GBIF
    gbif_data <- occ_search(scientificName = species_name,
                            hasCoordinate = TRUE,
                            coordinateUncertaintyInMeters = "0,1000",
                            geometry = northern_region_wkt,
                            limit = 200000)
    
    # Check if data was returned
    if(is.null(gbif_data$data) || nrow(gbif_data$data) == 0){
      cat(" No occurrence data foud for", species_name, "\n")
      failed_species <- c(failed_species, paste(species_name, "(no data)"))
      next
    }
    occurrence_data <- gbif_data$data
    cat("  Found", nrow(occurrence_data), "occurrence records\n")
    
    ## 2.2. Convert occurrences to spatial points ------------------------------
    
    # Remove records without coordinates
    occurrence_data <- occurrence_data[!is.na(occurrence_data$decimalLongitude) &
                                         !is.na(occurrence_data$decimalLatitude), ]
    
    if(nrow(occurrence_data) == 0) {
      cat("  No valid coordinates found for", species_name, "\n")
      failed_species <- c(failed_species, paste(species_name, "(no coordinates)"))
      next
    }
    
    # Create sf points object
    occurrence_points <- st_as_sf(occurrence_data,
                                  coords = c("decimalLongitude", "decimalLatitude"),
                                  crs = 4326) #WGS84
    
    # Transform to same CRS as biome polygons
    occurrence_points <- st_transform(occurrence_points, crs = biome_crs)
    
    ## 2.3. Filter points to those within the biomes ---------------------------
    
    # Find points within each biome
    boreal_points <- occurrence_points[lengths(st_within(occurrence_points, boreal_forest)) > 0,]
    tundra_points <- occurrence_points[lengths(st_within(occurrence_points, tundra)) > 0,]
    
    # Display the number of points found
    cat("  Points in boreal:", nrow(boreal_points), "\n")
    cat("  Points in tundra:", nrow(tundra_points), "\n")
    
    # Check if there are points in at least one biome
    if(nrow(boreal_points) == 0 && nrow(tundra_points) == 0){
      cat("  No points found in either biome for", species_name, "\n")
      failed_species <- c(failed_species, paste(species_name, "(no biome overlap)"))
      next
    }
    
    ## 2.4. Calculate distance to biome boundaries -----------------------------
    
    # Create variable to store distance
    distances <- numeric()
    
    # Boreal points: positive distance to tundra boundary
    if(nrow(boreal_points) > 0) {
      boreal_distances <- as.numeric(st_distance(boreal_points, st_boundary(tundra)))
      distances <- c(distances, boreal_distances)
    }
    
    # Tundra points: negative distance to boreal boundary
    if(nrow(tundra_points) > 0){
      tundra_distances <- -as.numeric(st_distance(tundra_points, st_boundary(boreal_forest)))
      distances <- c(distances, tundra_distances)
    }
    
    ## 2.5. Calculate species-level metrics ------------------------------------
    
    # Get average and sample size
    avg_distance <- mean(distances, na.rm = TRUE)
    n_records <- length(distances) 
    
    # Classify species based on average distance
    classification <- ifelse(avg_distance > 0, "boreal", "arctic")
    
    # Add classification to results dataframe
    new_row <- data.frame(species_name = species_name,
                          avg_distance_to_boundary = avg_distance,
                          classification = classification,
                          n_records_used = n_records,
                          stringsAsFactors = FALSE)
    
    # Add to results dataframe
    results_df <- rbind(results_df, new_row)
    succesful_species <- c(succesful_species, species_name)
    
    ## 2.6. Clean up temporary data --------------------------------------------
    rm(gbif_data, occurrence_data, occurrence_points, boreal_points, tundra_points, 
       boreal_distances, tundra_distances, distances)
    gc()
    
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
                             "species_biome_classification.RData"))


# END OF SCRIPT ----------------------------------------------------------------