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
source(here("scripts", "0_GBIF_creds.R"))

## 1.2. Species data -----------------------------------------------------------

# Load species list 
load(here("data", "derived_data", "corrected_species_list_27June2025.RData"))

# Extract species list 
species_list <- unique(corrected_species_list$SpeciesName)

# Initialise global species-taxon tracking
global_species_taxon_mapping <- data.frame(species = character(),
                                           taxon_key = character(),
                                           match_type = character(),
                                           confidence = character(),
                                           match_status = character(),
                                           chunk_number = integer(),
                                           stringsAsFactors = FALSE)

## 1.3. Load and prepare biomes ------------------------------------------------

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
boreal_forest_wgs84 <- st_transform(boreal_forest, "EPSG:4326")
tundra_wgs84 <- st_transform(tundra, "EPSG:4326")

# Transform biomes to North Pole Lambert Azimuthal Equal Area (EPSG: 3574) for spatial analysis
boreal_sf <- st_transform(boreal_forest, "EPSG:3574")
tundra_sf <- st_transform(tundra, "EPSG:3574")

## 1.4. Create analysis grid ---------------------------------------------------

# Combine biomes to use for the grid extent (in EPSG:3574)
combined_biomes <- st_union(boreal_sf, tundra_sf)

# Create bounding box from the combined biomes
combined_extent <- st_bbox(combined_biomes)

# Create grid with 0.45 degrees resolution (50km at high latitudes)
grid <- rast(extent = combined_extent, 
             resolution = c(50000, 50000),  # 50km 
             crs = "EPSG:3574")

# Convert to polygons
polygrid <- as.polygons(grid)

# Give each cell and ID
polygrid$cell_id <- 1:nrow(polygrid)

# Convert to sf object
polygrid_sf <- st_as_sf(polygrid)

# Filter cells within boreal biome
cells_in_boreal <- st_intersects(polygrid_sf, boreal_sf, sparse = FALSE)[,1]

# Filter cells within tundra biome
cells_in_tundra <- st_intersects(polygrid_sf, tundra_sf, sparse = FALSE)[,1]

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
  
  # Get study area WKT (in WGS84 for GBIF compatibility)
  study_area_wkt <- st_as_text(st_as_sfc(st_bbox(st_union(boreal_forest_wgs84, tundra_wgs84))))
  
  ### 2.1.1. Get taxon keys for species ----------------------------------------
  
  # Create data frame to store taxon information in
  taxon_info <- data.frame(species = species_list,
                           taxon_key = NA,
                           match_type = NA,
                           confidence = NA,
                           match_status = "FAILED",
                           chunk_number = chunk_number,
                           stringsAsFactors = FALSE)
  
  # Loop to look up taxon keys for each item on the species list
  for(i in 1:length(species_list)) {
    tryCatch({
      backbone <- name_backbone(species_list[i])
      if(!is.null(backbone$usageKey) && !is.na(backbone$usageKey)) {
        taxon_info$taxon_key[i] <- as.character(backbone$usageKey)
        taxon_info$match_type[i] <- ifelse(is.null(backbone$matchType), "UNKNOWN", backbone$matchType)
        taxon_info$confidence[i] <- ifelse(is.null(backbone$confidence), "UNKNOWN", 
                                           as.character(backbone$confidence))
        taxon_info$match_status[i] <- "SUCCESS"
        cat("  ✓", species_list[i], ":", backbone$usageKey, 
            "(", backbone$matchType, ",", backbone$confidence, ")\n")
      } else {
        taxon_info$match_type[i] <- ifelse(is.null(backbone$matchType), "NO_MATCH", backbone$matchType)
        taxon_info$confidence[i] <- "N/A"
        cat("  ✗ No taxon key for", species_list[i], 
            "- Match type:", taxon_info$match_type[i], "\n")
      }
    }, error = function(e) {
      taxon_info$match_type[i] <- "ERROR"
      taxon_info$confidence[i] <- "N/A"
      cat("  ✗ Error with", species_list[i], ":", e$message, "\n")
    })
    Sys.sleep(0.5)
  }
  
  # Update global tracking
  global_species_taxon_mapping <<- rbind(global_species_taxon_mapping, taxon_info)
  
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
                                 pred_lte("coordinateUncertaintyInMeters", 50000),
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
    
    # Process the donwloaded data with distances (using function defined below)
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
  
  # Create UNIQUE temp directory for each download
  temp_dir <- file.path(tempdir(), paste0("gbif_", basename(zip_file), "_", Sys.time() |> as.numeric()))
  dir.create(temp_dir, recursive = TRUE)
  
  # Extract and read CSV
  unzip(zip_file, exdir = temp_dir, overwrite = TRUE)
  csv_files <- list.files(temp_dir, pattern = "\\.csv$", full.names = TRUE)
  
  # DEBUG: Show what we're actually reading
  cat("  DEBUG: Temp dir:", temp_dir, "\n")
  cat("  DEBUG: CSV files found:", length(csv_files), "\n")
  cat("  DEBUG: Reading file:", basename(csv_files[1]), "\n")
  
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
  
  # Clean up the unique temp directory
  unlink(temp_dir, recursive = TRUE)
  
  # Remove records with missing coordinates
  occ_data <- occ_data[!is.na(decimalLongitude) & !is.na(decimalLatitude)]
  
  ### 2.2.2. Spatial join with grid --------------------------------------------
  
  # Convert to spatial data
  occ_sf <- st_as_sf(occ_data, 
                     coords = c("decimalLongitude", "decimalLatitude"),
                     crs = 4326)
  
  # Transform occurrences to EPSG:3574 to match the grid
  occ_sf_proj <- st_transform(occ_sf, "EPSG:3574")
  
  # Spatial join with grid
  joined <- st_join(occ_sf_proj, grid_polygons)
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
                                distance_to_boundary_m = NA,
                                stringsAsFactors = FALSE)
  
  # Add cell information 
  for(i in 1:nrow(distance_lookup)) {
    cell_id <- distance_lookup$cell_id[i]
    centroid <- cell_centroids[cell_centroids$cell_id == cell_id, ]
    cell_info <- grid_polygons[grid_polygons$cell_id == cell_id, ]
    
    if(cell_info$in_boreal) {
      # Cell is in boreal forest - distance to boreal boundary (positive)
      dist_to_boundary <- as.numeric(st_distance(centroid, boreal_boundary))
      distance_lookup$distance_to_boundary_m[i] <- mean(dist_to_boundary)
      
    } else if(cell_info$in_tundra) {
      # Cell is in tundra - distance to tundra boundary (negative)
      dist_to_boundary <- as.numeric(st_distance(centroid, tundra_boundary))
      distance_lookup$distance_to_boundary_m[i] <- -mean(dist_to_boundary)
    }
  }
  
  # Convert distances from meters to kilometers
  distance_lookup$distance_to_boundary_km <- distance_lookup$distance_to_boundary_m / 1000
  
  
  ### 2.2.5. Merge and finalize results ----------------------------------------
  
  # Merge distances with counts
  final_results <- merge(counts_filtered, 
                         distance_lookup[, c("cell_id", "distance_to_boundary_m", 
                                             "distance_to_boundary_km")], 
                         by = "cell_id")
  
  # Add biome classification
  final_results$biome <- ifelse(final_results$in_boreal, "boreal", 
                                ifelse(final_results$in_tundra, "tundra", "other"))
  
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

## 2.4. Function to manage GBIF citations --------------------------------------

save_gbif_citations <- function(download_metadata_list, filename = "gbif_download_citations.txt"){
 
  # Check if there is any metadata to save
  if(length(download_metadata_list) == 0) {
    cat("No download metadata to save\n")
    return()
  }
  
  ### 2.4.1. Create file to store citations ------------------------------------
  citation_lines <- c("GBIF Download Citations",
                      "=" %>% rep(50) %>% paste(collapse = ""),
                      paste("Generated on:", Sys.time()),"",
                      "Individual Downloads:",
                      "-" %>% rep(30) %>% paste(collapse = ""),"")
  
  # Add metadata
  for(i in 1:length(download_metadata_list)) {
    metadata <- download_metadata_list[[i]]
    if(!is.null(metadata)) {
      citation_lines <- c(citation_lines,
                          paste("Chunk", metadata$chunk_number, ":"),
                          paste("  Download Key:", metadata$download_key),
                          paste("  DOI:", metadata$doi),
                          paste("  Download Link:", metadata$download_link),
                          paste("  Total Records:", metadata$total_records),
                          paste("  Species (", length(metadata$species_in_chunk), "):", paste(metadata$species_in_chunk, collapse = ", ")),
                          paste("  Created:", metadata$created),
                          paste("  Citation:", metadata$citation),
                          ""
      )
    }
  }
  
  ### 2.4.2. Add combined citation ---------------------------------------------
  
  all_dois <- sapply(download_metadata_list, function(x) if(!is.null(x)) x$doi else NULL)
  all_dois <- all_dois[!sapply(all_dois, is.null)]
  
  if(length(all_dois) > 0) {
    citation_lines <- c(citation_lines,
                        "",
                        "Combined Citation:",
                        "-" %>% rep(30) %>% paste(collapse = ""),
                        "",
                        "GBIF.org (", format(Sys.Date(), "%Y"), ") GBIF Occurrence Downloads: ",
                        paste(all_dois, collapse = "; "),
                        " accessed via GBIF.org on ", format(Sys.Date(), "%d %B %Y"), ".",
                        "",
                        "BibTeX format:",
                        "@misc{gbif_downloads,",
                        "  author = {GBIF.org},",
                        paste0("  title = {GBIF Occurrence Downloads},"),
                        paste0("  year = {", format(Sys.Date(), "%Y"), "},"),
                        paste0("  note = {", paste(all_dois, collapse = "; "), "},"),
                        paste0("  url = {https://www.gbif.org},"),
                        paste0("  urldate = {", format(Sys.Date(), "%Y-%m-%d"), "}"),
                        "}"
    )
  }
  
  ### 2.4.3. Write files -------------------------------------------------------
  
  # Write to file
  writeLines(citation_lines, filename)
  cat("Citations saved to:", filename, "\n")
  
  # Also save as RDS for programmatic access
  saveRDS(download_metadata_list, gsub("\\.txt$", ".rds", filename))
  cat("Metadata saved to:", gsub("\\.txt$", ".rds", filename), "\n")
}

## 2.5. Function to save species-taxon key mappings ----------------------------
save_species_taxon_mapping <- function(mapping_df, filename_base = "species_taxon_key_mapping"){
  
  # Create summary statistics
  total_species <- nrow(mapping_df)
  successful_matches <- sum(mapping_df$match_status == "SUCCESS")
  failed_matches <- sum(mapping_df$match_status == "FAILED")
  error_matches <- sum(mapping_df$match_type == "ERROR")
  
  # Save as CSV
  csv_filename <- paste0(filename_base, "_", format(Sys.Date(), "%B%d"), ".csv")
  write.csv(mapping_df, here("data", "derived_data", csv_filename), row.names = FALSE)
  
  # Create a detailed text report
  txt_filename <- paste0(filename_base, "_", format(Sys.Date(), "%B%d"), ".txt")
  
  report_lines <- c(
    "Species-Taxon Key Mapping Report",
    "=" %>% rep(50) %>% paste(collapse = ""),
    paste("Generated on:", Sys.time()),
    paste("Analysis completed using EPSG:3574 projection"),
    "",
    "SUMMARY STATISTICS:",
    "-" %>% rep(30) %>% paste(collapse = ""),
    paste("Total species processed:", total_species),
    paste("Successful matches:", successful_matches, paste0("(", round(100*successful_matches/total_species, 1), "%)")),
    paste("Failed matches:", failed_matches, paste0("(", round(100*failed_matches/total_species, 1), "%)")),
    paste("Errors:", error_matches, paste0("(", round(100*error_matches/total_species, 1), "%)")),
    "",
    "MATCH TYPE BREAKDOWN:",
    "-" %>% rep(30) %>% paste(collapse = ""),
    ""
  )
  
  # Add match type breakdown
  if(nrow(mapping_df) > 0) {
    match_type_summary <- table(mapping_df$match_type)
    for(match_type in names(match_type_summary)) {
      count <- match_type_summary[match_type]
      pct <- round(100 * count / total_species, 1)
      report_lines <- c(report_lines, paste("  ", match_type, ":", count, paste0("(", pct, "%)")))
    }
  }
  
  report_lines <- c(report_lines,
                    "",
                    "CONFIDENCE BREAKDOWN (for successful matches):",
                    "-" %>% rep(30) %>% paste(collapse = ""),
                    ""
  )
  
  # Add confidence breakdown for successful matches
  successful_mapping <- mapping_df[mapping_df$match_status == "SUCCESS" & mapping_df$confidence != "N/A", ]
  if(nrow(successful_mapping) > 0) {
    confidence_summary <- table(successful_mapping$confidence)
    for(confidence in names(confidence_summary)) {
      count <- confidence_summary[confidence]
      pct <- round(100 * count / successful_matches, 1)
      report_lines <- c(report_lines, paste("  ", confidence, ":", count, paste0("(", pct, "% of successful)")))
    }
  }
  
  report_lines <- c(report_lines,
                    "",
                    "DETAILED MAPPINGS:",
                    "-" %>% rep(30) %>% paste(collapse = ""),
                    "See the CSV file for complete species-by-species details.",
                    "",
                    paste("CSV file saved as:", csv_filename),
                    paste("RDS file saved as:", gsub("\\.csv$", ".rds", csv_filename))
  )
  
  # Write text report
  writeLines(report_lines, here("data", "derived_data", txt_filename))
  
  # Save as RDS for programmatic access
  rds_filename <- gsub("\\.csv$", ".rds", csv_filename)
  saveRDS(mapping_df, here("data", "derived_data", rds_filename))
  
  cat("Species-taxon key mapping saved:\n")
  cat("  Text report:", txt_filename, "\n")
  cat("  CSV file:", csv_filename, "\n")
  cat("  RDS file:", rds_filename, "\n")
  cat("Summary: ", successful_matches, "/", total_species, " species successfully matched\n")
}

# 3. MAIN ANALYSIS FUNCTIONS ---------------------------------------------------

## 3.1. Main function to process all species -----------------------------------

analyze_species_list <- function(species_list, chunk_size = 5, start_chunk = 1){
  
  # Add some way to track things
  cat("=== Starting Analysis of", length(species_list), "Species ===\n")
  cat("Chunk size:", chunk_size, "\n")
  cat("Starting from chunk:", start_chunk, "\n")
  
  ### 3.1.1. Use the biomes and grid created before ----------------------------
  cat("Using pre-loaded biomes and grid...\n")
  cat("Grid cells available:", nrow(polygrid_filtered), "\n")
  
  ### 3.1.2. Process species in chunks -----------------------------------------
  
  # Split species into chunks
  species_chunks <- split(species_list, ceiling(seq_along(species_list)/chunk_size))
  total_chunks <- length(species_chunks)
  
  # Initialize results storage
  all_results <- list()
  all_download_metadata <- list()
  
  # Process chunks starting from start_chunk
  for(chunk_i in start_chunk:total_chunks) {
    cat("\n" , rep("=", 50), "\n")
    cat("CHUNK", chunk_i, "of", total_chunks, "\n")
    cat(rep("=", 50), "\n")
    
    current_species <- species_chunks[[chunk_i]]
    
    tryCatch({
      chunk_output <- process_species_chunk_with_distances(
        current_species, 
        polygrid_filtered, 
        boreal_sf, 
        tundra_sf,
        chunk_i)
      
      # Extract results and metadata
      chunk_results <- chunk_output$results
      chunk_metadata <- chunk_output$download_metadata
      
      if(nrow(chunk_results) > 0) {
        all_results[[chunk_i]] <- chunk_results
        cat("✓ Chunk", chunk_i, "completed:", nrow(chunk_results), "results\n")
        
        # Save intermediate results
        saveRDS(chunk_results, here("data", "derived_data", paste0("chunk_", chunk_i, "_results.rds")))
      } else {
        cat("✗ Chunk", chunk_i, "produced no results\n")
      }
      
      # Store metadata regardless of results
      if(!is.null(chunk_metadata)) {
        all_download_metadata[[chunk_i]] <- chunk_metadata
        saveRDS(chunk_metadata, here("data", "derived_data", paste0("chunk_", chunk_i, "_metadata.rds")))
      }
      
      # Save progress
      if(length(all_results) > 0) {
        combined_results <- do.call(rbind, all_results)
        saveRDS(combined_results, here("data", "derived_data", "species_analysis_progress.rds"))
        cat("Progress saved. Total results so far:", nrow(combined_results), "\n")
      }
      
      # Save citation progress
      if(length(all_download_metadata) > 0) {
        save_gbif_citations(all_download_metadata, here("data", "derived_data", "gbif_citations_progress.txt"))
      }
      
      # Be nice to GBIF between chunks
      if(chunk_i < total_chunks) {
        cat("Waiting 2 minutes before next chunk...\n")
        Sys.sleep(120)
      }
      
    }, error = function(e) {
      cat("ERROR in chunk", chunk_i, ":", e$message, "\n")
      cat("Continuing to next chunk...\n")
    })
  }
  
  ### 3.1.3. Combine results and create final outputs --------------------------
  if(length(all_results) > 0) {
    final_results <- do.call(rbind, all_results)
    
    # Calculate species summaries
    species_summaries <- calculate_species_summaries(final_results)
    
    # Save final results
    saveRDS(final_results, here("data", "derived_data", 
                                "dist_to_biome_boundary_June25.rds"))
    write.csv(final_results, here("data", "derived_data", 
                                  "dist_to_biome_boundary_June25.csv"), 
              row.names = FALSE)
    
    # Save species summaries
    saveRDS(species_summaries, here("data", "derived_data", 
                                    "species_summaries_dist_to_biome_boundary_June25.rds"))
    write.csv(species_summaries, here("data", "derived_data", 
                                      "species_summaries_dist_to_biome_boundary_June25.csv"), 
              row.names = FALSE)
    
    # Save final citations
    save_gbif_citations(all_download_metadata, here("data", "derived_data", 
                                                    "biome_boundary_gbif_download_citations_final_June25.txt"))
    
    # Save final species-taxon key mapping
    save_species_taxon_mapping(global_species_taxon_mapping, "species_taxon_key_mapping_final")
    
    cat("\n=== ANALYSIS COMPLETE ===\n")
    cat("Total cell-species results:", nrow(final_results), "\n")
    cat("Species analyzed:", length(unique(final_results$species)), "\n")
    cat("Species-taxon mappings tracked:", nrow(global_species_taxon_mapping), "\n")
    cat("Files saved:\n")
    cat("  - dist_to_biome_boundary_EPSG3574_June27.rds/.csv (all cell-level data)\n")
    cat("  - species_summaries_dist_to_biome_boundary_EPSG3574_June27.rds/.csv (species-level statistics)\n")
    cat("  - biome_boundary_gbif_download_citations_final_EPSG3574_June27.txt (citation information)\n")
    cat("  - species_taxon_key_mapping_final_June27.csv/.txt/.rds (species-taxon mappings with confidence)\n")
    
    # Summary statistics
    cat("\nOverall Summary:\n")
    cat("Projection used: EPSG:3574 (North Pole Lambert Azimuthal Equal Area)\n")
    cat("Boreal cells:", sum(final_results$biome == "boreal"), "\n")
    cat("Tundra cells:", sum(final_results$biome == "tundra"), "\n")
    cat("Distance range (km):", range(final_results$distance_to_boundary_km, na.rm = TRUE), "\n")
    cat("GBIF downloads used:", length(all_download_metadata), "\n")
    
    # Species mapping summary
    successful_matches <- sum(global_species_taxon_mapping$match_status == "SUCCESS")
    total_species <- nrow(global_species_taxon_mapping)
    cat("Species successfully matched to GBIF:", successful_matches, "/", total_species, 
        "(", round(100*successful_matches/total_species, 1), "%)\n")
    
    return(list(
      results = final_results,
      species_summaries = species_summaries,
      download_metadata = all_download_metadata,
      species_taxon_mapping = global_species_taxon_mapping
    ))
  } else {
    cat("No results produced\n")
    # Still save the species-taxon mapping even if no occurrence data was found
    if(nrow(global_species_taxon_mapping) > 0) {
      save_species_taxon_mapping(global_species_taxon_mapping, "species_taxon_key_mapping_final")
    }
    return(list(results = data.frame(), 
                species_summaries = data.frame(), 
                download_metadata = list(),
                species_taxon_mapping = global_species_taxon_mapping))
  }
}

# 4. RUN ANALYSIS --------------------------------------------------------------

cat("Script loaded successfully. Starting analysis...\n")

# Test with a small subset first 
#test_results <- analyze_species_list(species_list[1:10], chunk_size = 2, start_chunk = 1)

# Run full analysis
results <- analyze_species_list(species_list, chunk_size = 5, start_chunk = 1)

# END OF SCRIPT ----------------------------------------------------------------