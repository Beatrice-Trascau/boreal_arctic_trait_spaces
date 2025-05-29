##----------------------------------------------------------------------------##
# PAPER 3: BOREAL AND ARCTIC PLANT SPECIES TRAIT SPACES 
# 2.1_GBIF_download
# This script contains code which downloads the GBIF Occurrence records used in
# the analysis
##----------------------------------------------------------------------------##

# 1. SETUP ---------------------------------------------------------------------

# Load data & packages
library(here)
source(here("scripts", "0_setup.R"))
load(here("data", "derived_data", "all_filtered_standardised_species.RData"))

# Setting username, password and email 
Sys.setenv(GBIF_USER = "my_username") # change with your details
Sys.setenv(GBIF_PWD = "my_password")
Sys.setenv(GBIF_EMAIL = "my_email")

# 2. PREPARE SPECIES LIST ------------.-----------------------------------------

# Get species list
species_list <- unique(corrected_species_list$CheckedSpeciesName)

# Convert species names to GBIF taxon keys
taxon_keys <- c()
failed_species <- c()

for(i in 1:length(species_list)){
  # get taxon keys for each species
  key_lookup <- name_backbone(species_list[i])
  
  # store in either taxon_keys or failed_species based on success
  if(!is.null(key_lookup$usageKey) && key_lookup$matchType != "NONE"){
    taxon_keys <- c(taxon_keys, key_lookup$usageKey)
  } else {
    failed_species <- c(failed_species, species_list[i])
  }
}

# Print failed species
if(length(failed_species) > 0){
  print(failed_species)
}

# 3. CREATE DOWNLOAD REUQEST ---------------------------------------------------

# Create download request
download_key <- occ_download(pred_in("taxonKey", taxon_keys),
                             pred_gte("decimalLatitude", 40),
                             pred_lte("decimalLatitude", 90),
                             pred("hasCoordinate", TRUE),
                             format = "DWCA")

# Save download key
cat("Download key:", download_key, "\n") #0006611-250525065834625

# Check download status
occ_download_wait(download_key)

# Get download metadata
download_meta <- occ_download_meta(download_key)

# 4. IMPORT AND SAVE DATA ------------------------------------------------------

# Import the download data
occurrence_species <- occ_download_get(download_key) |>
  occ_download_import()

# Save the data
save(occurrence_species, here::here("data", "raw_data", 
                                    "occurrence_species_northern_hemisphere.rda"))

# Save download key metadata
download_info <- list(download_key = download_key,
                      species_list = species_list,
                      taxon_keys = taxon_keys,
                      failed_species = failed_species,
                      download_date = Sys.Date(),
                      total_records = nrow(occurrence_species))

# Save download info
save(download_info, file = here::here("data", "raw_data", "download_info_species.rda"))

# END OF SCRIPT ----------------------------------------------------------------
