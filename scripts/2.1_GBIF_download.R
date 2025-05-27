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

# 2. DOWNLOAD GBIF OCCURRENCES -------------------------------------------------

# Define polar projection
polarproj <- '+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'

# List of northern countries
northern_countries <- c('NO', 'SE', 'FI','RU','CA','IS','GL','SJ','MN','JP')

# Extract list of species from df
species_list <- unique(corrected_species_list$CheckedSpeciesName)

# Function to download GBIF occurrences for each species
download_species_occs <- function(species_name, verbose = TRUE){
  if(verbose) cat("Downloading:", species_name, "\n")
  # search for occurrences
  gbif_result <- occ_search(scientificName = species_name,
                            hasCoordinate = TRUE,
                            hasGeospatialIssue = FALSE,
                            decimalLatitude = "40,90",
                            decimalLongitude = "-180,180",
                            basisOfRecord = c("HUMAN_OBSERVATION", "OBSERVATION", 
                                              "MACHINE_OBSERVATION", "PRESERVED_SPECIMEN", 
                                              "LIVING_SPECIMEN"),
                            fields = c("key", "scientificName", "decimalLatitude",
                                       "decimalLongitude", "year", "month", "country", 
                                       "basisOfRecord", "coordinateUncertaintyInMeters"))
  
  # check if data was returned
  if(gbif_result$meta$count == 0) {
    if(verbose) cat("  No occurrences found for:", species_name, "\n")
    return(NULL)
  }
  
  # extract data
  occurrences <- gbif_result$data
  
  # quick data clean
  occurrences_clean <- occurrences |>
    # remove any records without coordinates
    filer(!is.na(decimalLatitude), !is.na(decimalLongitude)) |>
    # remove records with high coordinate uncertainty 
    filter(is.na(coordinateUncertaintyInMeters) | coordinateUncertaintyInMeters <= 10000) |>
    # add clean species name
    mutate(species_clean = gsub(" ", ".", species_name)) |>
    # select relevant columns
    select(species = scientificName, species_clean, decimalLongitude, decimalLatitude, 
           year, month, country, basisOfRecord, coordinateUncertaintyInMeters)
  
  # get information about how many records were cleaned out
  if(verbose) {
    cat("  Found:", nrow(gbif_result$data), "total records\n")
    cat("  After cleaning:", nrow(occurrences_clean), "records\n")
  }
  
  return(occurrences_clean)
}

# Download GBIF data
gbif_results <- download_species_occs(species_list)

