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

## 1.3. Prepare df to store results --------------------------------------------

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
                            coordinateUncertaintyInMeters = "0,1000")
    
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
    
  })
}