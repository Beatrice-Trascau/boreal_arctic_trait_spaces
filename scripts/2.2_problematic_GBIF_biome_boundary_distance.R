##----------------------------------------------------------------------------##
# PAPER 3: BOREAL AND ARCTIC PLANT SPECIES TRAIT SPACES 
# 2.2_problematic_GBIF_biome_boundary_distance
# This script contains code which calcualtes the distance to biome boundary for
# Picea sitchensis and Ilex aquifolium - two of the problematic species from
# the previous calculations
##----------------------------------------------------------------------------##

# 1. SETUP ---------------------------------------------------------------------

## 1.1. Libraries --------------------------------------------------------------

# Load packages and functions
library(here)
source(here("scripts", "0_setup.R"))
source(here("scripts", "0_GBIF_creds.R"))

## 1.2. Species list -----------------------------------------------------------

# Extract species list 
target_species <- c("Picea sitchensis", "Ilex aquifolium")

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

# Get boundaries for distance calculations
boreal_boundary <- st_boundary(boreal_sf)
tundra_boundary <- st_boundary(tundra_sf)

# Get study area for GBIF query
study_area_wkt <- st_as_text(st_as_sfc(st_bbox(st_union(boreal_forest_wgs84, tundra_wgs84))))

# 3. PROCESS EACH SPECIES ------------------------------------------------------

# Initialise list to store results
all_occurrences <- list()
download_info <- list()

# Process each of the two species
for(species_name in target_species) {
  
  # Get taxon key
  backbone <- name_backbone(species_name)
  if(is.null(backbone$usageKey) || is.na(backbone$usageKey)) {
    cat("No taxon key found for", species_name, "\n")
    next
  }
  
  taxon_key <- backbone$usageKey
  cat("Found taxon key:", taxon_key, "\n")
  
  # Download GBIF data
  download_key <- occ_download(pred_in("taxonKey", taxon_key),
                               pred("geometry", study_area_wkt),
                               pred("hasCoordinate", TRUE),
                               pred("hasGeospatialIssue", FALSE),
                               pred_lte("coordinateUncertaintyInMeters", 50000),
                               format = "SIMPLE_CSV")
  
  # Wait for download
  occ_download_wait(download_key, status_ping = 60)
  
  # Get download metadata
  download_meta <- occ_download_meta(download_key)
  cat("Download completed. Total records:", download_meta$totalRecords, "\n")
  
  # Store download info for citations
  download_info[[species_name]] <- list(species = species_name,
                                        download_key = download_key,
                                        doi = download_meta$doi,
                                        total_records = download_meta$totalRecords,
                                        citation = paste0("GBIF Occurrence Download ", download_meta$doi, 
                                                          " accessed via GBIF.org on ", Sys.Date()))
  
  if(download_meta$totalRecords == 0) {
    cat("No records found for", species_name, "\n")
    next
  }
  
  # Download and save file permanently
  safe_species_name <- gsub(" ", "_", species_name)
  zip_filename <- paste0("GBIF_", safe_species_name, "_", format(Sys.Date(), "%Y%m%d"), ".zip")
  zip_file <- occ_download_get(download_key, 
                               path = here("data", "raw_data"), 
                               file = zip_filename,
                               overwrite = TRUE)
  
  # Extract directly to permanent location
  extract_dir <- here("data", "derived_data", "gbif_extracts")
  if(!dir.exists(extract_dir)) dir.create(extract_dir, recursive = TRUE)
  
  species_extract_dir <- file.path(extract_dir, safe_species_name)
  unzip(zip_file, exdir = species_extract_dir, overwrite = TRUE)
  csv_files <- list.files(species_extract_dir, pattern = "\\.csv$", full.names = TRUE)
  
  # Read occurrence data directly from permanent location
  occ_data <- fread(csv_files[1], 
                    select = c("gbifID", "taxonKey", "species", "decimalLongitude", 
                               "decimalLatitude", "coordinateUncertaintyInMeters",
                               "country", "stateProvince", "year"),
                    quote = "")
  
  # Remove records with missing coordinates
  occ_data_clean <- occ_data[!is.na(decimalLongitude) & !is.na(decimalLatitude)]
  cat("Processing", nrow(occ_data_clean), "records with coordinates\n")
  
  if(nrow(occ_data_clean) == 0) next
  
  # Convert to spatial data and transform to EPSG:3574 for distance calculations
  occ_sf <- st_as_sf(occ_data_clean, 
                     coords = c("decimalLongitude", "decimalLatitude"),
                     crs = 4326) |>
    st_transform("EPSG:3574")
  
  # Determine which biome each occurrence is in
  in_boreal <- st_intersects(occ_sf, boreal_sf, sparse = FALSE)[,1]
  in_tundra <- st_intersects(occ_sf, tundra_sf, sparse = FALSE)[,1]
  
  # Initialize columns
  occ_data_clean$distance_to_boundary_m <- NA
  occ_data_clean$distance_to_boundary_km <- NA
  occ_data_clean$biome <- NA
  
  # Calculate distances for boreal occurrences (positive distances)
  if(sum(in_boreal) > 0) {
    boreal_points <- occ_sf[in_boreal, ]
    boreal_distances <- st_distance(boreal_points, boreal_boundary)
    min_distances <- apply(boreal_distances, 1, min)
    occ_data_clean$distance_to_boundary_m[in_boreal] <- as.numeric(min_distances)
    occ_data_clean$biome[in_boreal] <- "boreal"
    cat("  Calculated distances for", sum(in_boreal), "boreal occurrences\n")
  }
  
  # Calculate distances for tundra occurrences (negative distances)
  if(sum(in_tundra) > 0) {
    tundra_points <- occ_sf[in_tundra, ]
    tundra_distances <- st_distance(tundra_points, tundra_boundary)
    min_distances <- apply(tundra_distances, 1, min)
    occ_data_clean$distance_to_boundary_m[in_tundra] <- -as.numeric(min_distances)
    occ_data_clean$biome[in_tundra] <- "tundra"
    cat("  Calculated distances for", sum(in_tundra), "tundra occurrences\n")
  }
  
  # Convert to kilometers
  occ_data_clean$distance_to_boundary_km <- occ_data_clean$distance_to_boundary_m / 1000
  
  # Keep only occurrences within study biomes
  occ_final <- occ_data_clean[!is.na(occ_data_clean$biome), ]
  cat("Final dataset:", nrow(occ_final), "occurrences within study biomes\n")
  
  # Store results
  all_occurrences[[species_name]] <- occ_final
  
  cat("Waiting 30 seconds before next species...\n\n")
  Sys.sleep(30)
}

# 4. COMBINE RESULTS AND CALCULATE SUMMARIES -----------------------------------

if(length(all_occurrences) > 0) {
  
  # Combine all occurrence data
  final_occurrences <- do.call(rbind, all_occurrences)
  
  # Calculate species-level summary statistics
  species_summaries <- final_occurrences |>
    group_by(species) |>
    # add count statistics
    summarise(total_occurrences = n(),
              # biome distribution
              occurrences_in_boreal = sum(biome == "boreal", na.rm = TRUE),
              occurrences_in_tundra = sum(biome == "tundra", na.rm = TRUE),
              pct_occurrences_boreal = round(100 * sum(biome == "boreal", na.rm = TRUE) / n(), 1),
              pct_occurrences_tundra = round(100 * sum(biome == "tundra", na.rm = TRUE) / n(), 1),
              # distance statistics (all occurrences)
              mean_distance_km = round(mean(distance_to_boundary_km, na.rm = TRUE), 2),
              median_distance_km = round(median(distance_to_boundary_km, na.rm = TRUE), 2),
              sd_distance_km = round(sd(distance_to_boundary_km, na.rm = TRUE), 2),
              min_distance_km = round(min(distance_to_boundary_km, na.rm = TRUE), 2),
              max_distance_km = round(max(distance_to_boundary_km, na.rm = TRUE), 2),
              # distance statistics for boreal occurrences only
              mean_distance_boreal_km = round(mean(distance_to_boundary_km[biome == "boreal"], na.rm = TRUE), 2),
              median_distance_boreal_km = round(median(distance_to_boundary_km[biome == "boreal"], na.rm = TRUE), 2),
              # distance statistics for tundra occurrences only  
              mean_distance_tundra_km = round(mean(distance_to_boundary_km[biome == "tundra"], na.rm = TRUE), 2),
              median_distance_tundra_km = round(median(distance_to_boundary_km[biome == "tundra"], na.rm = TRUE), 2),
              # weighted mean (same as regular mean for individual occurrences)
              weighted_mean_distance_km = round(mean(distance_to_boundary_km, na.rm = TRUE), 2),
              # range metrics
              distance_range_km = round(max_distance_km - min_distance_km, 2),
              .groups = 'drop')
  
  # Replace NaN with NA
  species_summaries[sapply(species_summaries, is.nan)] <- NA
  
  # Save final results
  date_suffix <- format(Sys.Date(), "%Y%m%d")
  
  # Save occurrence data
  saveRDS(final_occurrences, here("data", "derived_data", paste0("target_species_occurrences_", date_suffix, ".rds")))
  write.csv(final_occurrences, here("data", "derived_data", paste0("target_species_occurrences_", date_suffix, ".csv")), row.names = FALSE)
  
  # Save species summaries  
  saveRDS(species_summaries, here("data", "derived_data", paste0("target_species_summaries_", date_suffix, ".rds")))
  write.csv(species_summaries, here("data", "derived_data", paste0("target_species_summaries_", date_suffix, ".csv")), row.names = FALSE)
  
  # Save download information for citations
  citation_lines <- c("GBIF Download Citations for Target Species Analysis",
                      paste("Generated on:", Sys.time()), "")
  
  for(species in names(download_info)) {
    info <- download_info[[species]]
    citation_lines <- c(citation_lines,
                        paste("Species:", info$species),
                        paste("  DOI:", info$doi),
                        paste("  Total Records:", info$total_records),
                        paste("  Citation:", info$citation), "")
  }
  
  writeLines(citation_lines, here("data", "derived_data", paste0("target_species_citations_", date_suffix, ".txt")))
  
  
} else {
  cat("No results produced\n")
}

# END OF SCRIPT ----------------------------------------------------------------