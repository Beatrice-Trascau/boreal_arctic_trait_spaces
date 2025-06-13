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

# 2. CALCULATE DISTANCE TO BIOME BOUNDARY --------------------------------------

## 2.1. Process species --------------------------------------------------------

process_species_chunk_with_distances <- function(species_list, grid_polygons, boreal_sf, tundra_sf, chunk_number = 1){
  cat("=== Processing Chunk", chunk_number, "===\n")
  cat("Species:", paste(species_list, collapse = ", "), "\n")
  
  # Get study area WKT
  study_area_wkt <- st_as_text(st_as_sfc(st_bbox(st_union(boreal_sf, tundra_sf))))
  
  ### 2.1.1. Get taxon keys for species ----------------------------------------
  
  # Create data frame to store taxon information in
  taxon_info <- data.frame(species = species_list,
                           taxon_key = NA,
                           stringsAsFactors = FALSE)
  
  # Loop to look up taxon keys for each item on the species list
  for(i in 1:length(species_list)) {
    tryCatch({
      backbone <- name_backbone(species_list[i])
      if(!is.null(backbone$usageKey) && !is.na(backbone$usageKey)) {
        taxon_info$taxon_key[i] <- backbone$usageKey
        cat("  ✓", species_list[i], ":", backbone$usageKey, "\n")
      } else {
        cat("  ✗ No taxon key for", species_list[i], "\n")
      }
    }, error = function(e) {
      cat("  ✗ Error with", species_list[i], ":", e$message, "\n")
    })
    Sys.sleep(0.5)
  }
  
  # Store valid keys in the dataframe
  valid_keys <- taxon_info$taxon_key[!is.na(taxon_info$taxon_key)]
  
  # Check if there were any valid taxon keys found
  if(length(valid_keys) == 0) {
    cat("No valid taxon keys found for this chunk!\n")
    return(list(results = data.frame(), download_metadata = NULL))
  }
  
  # Get summary of how many taxon keys were found
  cat("Found", length(valid_keys), "valid taxon keys\n")
  
  ### 2.1.2. Create and process GBIF download ----------------------------------
  tryCatch({
    download_key <- occ_download(pred_in("taxonKey", valid_keys),
                                 pred("geometry", study_area_wkt),
                                 pred("hasCoordinate", TRUE),
                                 pred("hasGeospatialIssue", FALSE),
                                 pred_lte("coordinateUncertaintyInMeters", 50000)
                                 format = "SIMPLE_CSV")
    
    # Wait for download
    occ_download_wait(download_key, status_ping = 60)
    
    # Get download information and store complete metadata
    download_info <- occ_download_meta(download_key)
    cat("Download completed. Total records:", download_info$totalRecords, "\n")
    
    # Store metadata for citation
    download_metadata <- list(chunk_number = chunk_number,
                              download_key = download_key,
                              doi = download_info$doi,
                              created = download_info$created,
                              modified = download_info$modified,
                              total_records = download_info$totalRecords,
                              species_in_chunk = species_list,
                              taxon_keys = valid_keys,
                              download_link = paste0("https://doi.org/", download_info$doi),
                              citation = paste0("GBIF Occurrence Download ", download_info$doi, 
                                                " accessed via GBIF.org on ", Sys.Date()))
    # Check if any records were found
    if(download_info$totalRecords == 0) {
      cat("No records found for these species\n")
      return(list(results = data.frame(), download_metadata = download_metadata))
    }
    
    # Download file
    cat("Downloading file...\n")
    zip_file <- occ_download_get(download_key, path = tempdir(), overwrite = TRUE)
    
    # Process the donwloaded data with distances (using function defined above)
    results <- process_occurrence_data_with_distances(zip_file, grid_polygons, boreal_sf, tundra_sf)
    
    # Clean up
    file.remove(zip_file)
    
    # Add chunk info and download key to results
    if(nrow(results) > 0) {
      results$chunk <- chunk_number
      results$download_key <- download_key
    }
    
    return(list(results = results, download_metadata = download_metadata))
}, error = function(e) {
  cat("Error in download process:", e$message, "\n")
  return(list(results = data.frame(), download_metadata = NULL))
})
}

## 2.2. Function to calculate distances ----------------------------------------
process_occurrence_data_with_distances <- function(zip_file, grid_polygons, boreal_sf, tundra_sf) {
  
  ### 2.2.1. Extract and read occurrence data ----------------------------------
  
  # Extract and read CSV
  temp_dir <- tempdir()
  unzip(zip_file, exdir = temp_dir, overwrite = TRUE)
  csv_files <- list.files(temp_dir, pattern = "\\.csv$", full.names = TRUE)
  
  if(length(csv_files) == 0) {
    return(data.frame())
  }
  
  csv_file <- csv_files[1]
  
  # Read occurrence data
  occ_data <- fread(csv_file, 
                    select = c("taxonKey", "species", "decimalLongitude", "decimalLatitude"),
                    quote = "")
  
  # Check how many occurrences there are
  cat("  Read", nrow(occ_data), "occurrence records\n")
  cat("  Species:", paste(unique(occ_data$species), collapse = ", "), "\n")
  
  if(nrow(occ_data) == 0) {
    return(data.frame())
  }
  
  # Remove records with missing coordinates
  occ_data <- occ_data[!is.na(decimalLongitude) & !is.na(decimalLatitude)]
  
  ### 2.2.2. Spatial join with grid --------------------------------------------
  
  # Convert to spatial data
  occ_sf <- st_as_sf(occ_data, 
                     coords = c("decimalLongitude", "decimalLatitude"),
                     crs = 4326)
  
  # Spatial join with grid
  joined <- st_join(occ_sf, grid_polygons)
  joined_clean <- joined[!is.na(joined$cell_id), ]
  
  # Check how many occurrences there were
  cat("  Found", nrow(joined_clean), "occurrences within grid cells\n")
  
  if(nrow(joined_clean) == 0) {
    return(data.frame())
  }
  
  ### 2.2.3. Count occurrences per cell per species ----------------------------
  
  # Get a dataframe of the joined grid & occurrences
  joined_df <- st_drop_geometry(joined_clean)
  
  # Get a count of the occurrences per cell
  counts <- joined_df |>
    group_by(cell_id, species, in_boreal, in_tundra) |>
    summarise(count = n(), .groups = 'drop')
  
  # Filter out cells with too few occurrences
  counts_filtered <- counts[counts$count > 1, ]
  
  if(nrow(counts_filtered) == 0) {
    return(data.frame())
  }
  
  ### 2.2.4. Calculate distances to biome boundaries ---------------------------
  
  # Get boundaries of both biomes
  boreal_boundary <- st_boundary(boreal_sf)
  tundra_boundary <- st_boundary(tundra_sf)
  
  # Get centroids of cells with occurrences
  unique_cells <- unique(counts_filtered$cell_id)
  cell_geometries <- grid_polygons[grid_polygons$cell_id %in% unique_cells, ]
  cell_centroids <- st_centroid(cell_geometries)
  
  # Calculate distances for each cell
  distance_lookup <- data.frame(cell_id = cell_centroids$cell_id,
                                distance_to_boundary = NA,
                                stringsAsFactors = FALSE)
  
  # Add cell information 
  for(i in 1:nrow(distance_lookup)) {
    cell_id <- distance_lookup$cell_id[i]
    centroid <- cell_centroids[cell_centroids$cell_id == cell_id, ]
    cell_info <- grid_polygons[grid_polygons$cell_id == cell_id, ]
    
    if(cell_info$in_boreal) {
      # Cell is in boreal forest - distance to boreal boundary (positive)
      dist_to_boundary <- as.numeric(st_distance(centroid, boreal_boundary))
      distance_lookup$distance_to_boundary[i] <- mean(dist_to_boundary)
      
    } else if(cell_info$in_tundra) {
      # Cell is in tundra - distance to tundra boundary (negative)
      dist_to_boundary <- as.numeric(st_distance(centroid, tundra_boundary))
      distance_lookup$distance_to_boundary[i] <- -mean(dist_to_boundary)
    }
  }
  
  # Convert distances from degrees to kilometers
  distance_lookup$distance_to_boundary_km <- distance_lookup$distance_to_boundary * 111
  
  
  ### 2.2.5. Merge and finalize results ----------------------------------------
  
  # Merge distances with counts
  final_results <- merge(counts_filtered, 
                         distance_lookup[, c("cell_id", "distance_to_boundary", 
                                             "distance_to_boundary_km")], 
                         by = "cell_id")
  
  # Add biome classification
  final_results <- merge(counts_filtered, 
                         distance_lookup[, c("cell_id", "distance_to_boundary", 
                                             "distance_to_boundary_km")], 
                         by = "cell_id")
  
  # Check how many distance caluclations were completed
  cat("  Completed distance calculations for", nrow(final_results), "records\n")
  return(as.data.frame(final_results))
}

## 2.3. Function to calculate species-level summaries --------------------------
calculate_species_summaries <- function(results_df){
  
  # Create df with summaries
  species_summaries <- results_df |>
    group_by(species) |>
    summarise(total_cells = n(),
              total_occurrences = sum(count),
              # biome distribution
              cells_in_boreal = sum(biome == "boreal"),
              cells_in_tundra = sum(biome == "tundra"),
              pct_cells_boreal = round(100 * sum(biome == "boreal") / n(), 1),
              pct_cells_tundra = round(100 * sum(biome == "tundra") / n(), 1),
              # distance statistics (all cells)
              mean_distance_km = round(mean(distance_to_boundary_km, na.rm = TRUE), 2),
              median_distance_km = round(median(distance_to_boundary_km, na.rm = TRUE), 2),
              sd_distance_km = round(sd(distance_to_boundary_km, na.rm = TRUE), 2),
              min_distance_km = round(min(distance_to_boundary_km, na.rm = TRUE), 2),
              max_distance_km = round(max(distance_to_boundary_km, na.rm = TRUE), 2),
              # distance statistics for boreal cells only
              mean_distance_boreal_km = round(mean(distance_to_boundary_km[biome == "boreal"], 
                                                   na.rm = TRUE), 2),
              median_distance_boreal_km = round(median(distance_to_boundary_km[biome == "boreal"], 
                                                       na.rm = TRUE), 2),
              # distance statistics for tundra cells only  
              mean_distance_tundra_km = round(mean(distance_to_boundary_km[biome == "tundra"], 
                                                   na.rm = TRUE), 2),
              median_distance_tundra_km = round(median(distance_to_boundary_km[biome == "tundra"], 
                                                       na.rm = TRUE), 2),
              # occurrence-weighted distances (cells with more occurrences get more weight)
              weighted_mean_distance_km = round(weighted.mean(distance_to_boundary_km, count, 
                                                              na.rm = TRUE), 2),
              # range metrics
              latitudinal_range_km = round(max_distance_km - min_distance_km, 2),
              
              .groups = 'drop')
  
  # Replace NaN with NA for cleaner output
  species_summaries[sapply(species_summaries, is.nan)] <- NA
  
  # Check how many species had summaries calculated
  cat("Species summaries calculated for", nrow(species_summaries), "species\n")

  return(as.data.frame(species_summaries))
}



