##----------------------------------------------------------------------------##
# PAPER 3: BOREAL AND ARCTIC PLANT SPECIES TRAIT SPACES 
# 5.1_authorship_check
##----------------------------------------------------------------------------##

# 1. SETUP ---------------------------------------------------------------------

# Load packages
library(here)
source(here("scripts", "0_setup.R"))

# Load trait data
load(here("data", "raw_data", "try_beatrice.RData"))
try_raw <- try.final.control

# Load authorship request spreadsheet
authroship_request <- read.xlsx(here("data", "raw_data", "TRY_authorship_req.xlsx"),
                                sheet = 1, skipEmptyRows = TRUE)

# 2. EXTRACT TRAIT IDS ---------------------------------------------------------

# Check out the authorship request file
View(authroship_request)

# Rename column to match global traits dataframe
authroship_request <- authroship_request |>
  rename(TraitID = Traits.ID)

# Extract the unique trait ID from collapsed column
cleaned_authorship <- authroship_request |>
  separate_rows(TraitID, sep = ",") |>
  filter(TraitID != " ")
  
# Check TraitID values
glimpse(cleaned_authorship) #TraitID = character
glimpse(try_raw) # TraitID = integer

# Convert TraitID column to integer in cleaned_authorsphi df
cleaned_authorship <- cleaned_authorship |>
  mutate(TraitID = as.numeric(TraitID))

# Extract the TraitIDs present in global_traits dataframe
cleaned_authorship$TraitID_present_in_global_traits <- cleaned_authorship$TraitID %in% try_raw$TraitID
cleaned_authorship$DatasetID_present_in_global_traits <- cleaned_authorship$DatasetID %in% try_raw$DatasetID

#Keep only 