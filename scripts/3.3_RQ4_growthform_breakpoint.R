##----------------------------------------------------------------------------##
# PAPER 3: BOREAL AND ARCTIC PLANT SPECIES TRAIT SPACES 
# 3.3_RQ4_growthform_breakpoint
# This script contains code for breakpoint regression analysis testing whether
# trait-distance relationships differ on either side of the biome boundary
# broken down by growth form
##----------------------------------------------------------------------------##

# 1. SETUP ---------------------------------------------------------------------

## 1.1. Load data --------------------------------------------------------------

# Load packages and functions
library(here)
source(here("scripts", "0_setup.R"))

# Load cleaned traits
load(here("data", "derived_data", "TRY_traits_cleaned_July2025.RData"))
cleaned_traits <- cleaned_traits_July2025

# Load distances to biome boundaries
detailed_results <- readRDS(here("data", "derived_data", 
                                 "species_summaries_dist_to_biome_boundary_June25.rds"))

# Load biome boundaries
global_biomes <- st_read(here("data", "raw_data", "biomes", 
                              "wwf_terr_ecos.shp"))

# Load CAFF quality check
caff_check <- read.xlsx(here("data", "derived_data", "caff_quality_check.xlsx"),
                        sheet = 1, skipEmptyRows = TRUE)

# Cleaned growth forms 
# Validated at workshop by Mariana & James
validated_growthforms <- read.xlsx(here("data", "derived_data",
                                        "BIEN_growth_forms_corrected_at_workshop.xlsx"),
                                   sheet = 1, skipEmptyRows = TRUE)