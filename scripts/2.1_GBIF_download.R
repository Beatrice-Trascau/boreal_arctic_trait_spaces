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
species_name <- list(corrected_species_list$CheckedSpeciesName)

# Function to download GBIF occurrences for each species