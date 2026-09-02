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

# Rename column to match
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

#Keep only records with TraitID_present_in_global_traits == TRUE
to_request <- cleaned_authorship |>
  filter(TraitID_present_in_global_traits == TRUE)

# Extract unique list of emails
request_emails <- unique(to_request$Email)

# 3. CLEAN SPECIES NAMES -------------------------------------------------------

## 3.1. Remove morphospecies and suspect names ---------------------------------

global_traits1 <- try_raw |>
  filter(!grepl("\\bsp\\.$|\\bsp\\b", AccSpeciesName) &
           !AccSpeciesName %in% c("Grass", "Fern", "Unknown", "Graminoid") &
           !AccSpeciesName %in% c("Hieracium sect.")) |>
  mutate(RawSpeciesName = AccSpeciesName,
         AccSpeciesName = if_else(AccSpeciesName == "Eri sch", 
                                  "Eriophorum scheuchzeri", 
                                  AccSpeciesName))

# Compare number of species brfore and after removing the morphospecies
length(unique(try_raw$AccSpeciesName))
length(unique(global_traits1$AccSpeciesName))

## 3.2. Standardise subspecies/varieties to species level ----------------------

global_traits2 <- global_traits1 |>
  mutate(StandardSpeciesName = str_replace(AccSpeciesName, 
                                           " (subsp\\.|var\\.|sect\\.).*$", ""))

# Remove 'x' at the end (hybrids)
global_traits3 <- global_traits2 |>
  mutate(StandardSpeciesName = str_replace(StandardSpeciesName, "\\sx$", ""))

# Remove Hieracium (generic name)
global_traits4 <- global_traits3 |>
  filter(StandardSpeciesName != "Hieracium")

## 3.3. Clean species names based on taxonomic check ---------------------------

global_traits5 <- global_traits4 |>
  mutate(StandardSpeciesName = case_when(
    StandardSpeciesName == "Casteleja occidens" ~ "Castilleja occidentalis",
    StandardSpeciesName == "Sausarrea angustifolium" ~ "Saussurea angustifolia",
    StandardSpeciesName == "Spirodela polyrrhiza" ~ "Spirodela polyrhiza",
    StandardSpeciesName == "Salix doniana" ~ "Salix purpurea",
    StandardSpeciesName == "Silene samojedora" ~ "Silene samojedorum",
    StandardSpeciesName == "Salix myrtifolia" ~ "Salix myrtillifolia",
    StandardSpeciesName == "Calamagrostis purpuras" ~ "Calamagrostis purpurea",
    StandardSpeciesName == "Pedicularis vertisilata" ~ "Pedicularis verticillata",
    StandardSpeciesName == "Peticites frigidus" ~ "Petasites frigidus",
    StandardSpeciesName == "Senecio atropurpuris" ~ "Senecio atropurpureus",
    StandardSpeciesName == "Gentia glauca" ~ "Gentiana glauca",
    StandardSpeciesName == "Polemonium acutifolium" ~ "Polemonium acutiflorum",
    StandardSpeciesName == "Rumex lapponum" ~ "Rumex lapponicus",
    StandardSpeciesName == "Pedicularis vertillis" ~ "Pedicularis verticillata",
    StandardSpeciesName == "Sabulina rossii" ~ "Sabulina rosei",
    StandardSpeciesName == "Echinops crispus" ~ "Echinops ritro",
    StandardSpeciesName == "Salix fuscenses" ~ "Salix fuscescens",
    StandardSpeciesName == "Salix laponicum" ~ "Salix lapponum",
    StandardSpeciesName == "Senecio atropupuris" ~ "Senecio atropurpureus",
    StandardSpeciesName == "Salix argyocarpon" ~ "Salix argyrocarpa",
    StandardSpeciesName == "Salix herbaceae-polaris" ~ "Salix herbacea",
    .default = StandardSpeciesName))

# Check how many species are left
length(unique(global_traits5$StandardSpeciesName))

# 4. EXTRACT DATASETS AND REFERENCES PROVIDING THE MOST RECORDS ----------------

## 4.1. Highest number of records ----------------------------------------------

# Get a list of the dataset names ordere by the number of records
dataset_counts <- global_traits5 |>
  count(Dataset, name = "n_records") |>
  arrange(desc(n_records))

# View the top datasets
print(dataset_counts)

# Get a list of references ordered by the number of records
reference_counts <- global_traits5 |>
  count(Reference, name = "n_records") |>
  arrange(desc(n_records))

# View the top references
print(reference_counts)

## 4.2. Datasets and references with 50+ records per trait or species ----------

# Get dataset names for datasets giving >= 50 records per trait
datasets_50plus_per_trait <- global_traits5 |>
  group_by(Dataset, TraitNameNew) |>
  summarise(n_records = n(), .groups = "drop") |>
  filter(n_records >= 200) |>
  arrange(desc(n_records))

# Get dataset names for datasets giving >= 50 records per species
datasets_50plus_per_species <- global_traits5 |>
  group_by(Dataset, StandardSpeciesName) |>
  summarise(n_records = n(), .groups = "drop") |>
  filter(n_records >= 200) |>
  arrange(desc(n_records))

# Get reference names for datasets giving >= 50 records per trait
references_50plus_per_trait <- global_traits5 |>
  group_by(Reference, TraitNameNew) |>
  summarise(n_records = n(), .groups = "drop") |>
  filter(n_records >= 200) |>
  arrange(desc(n_records))

# Get reference names for datasets giving >= 50 records per species
references_50plus_per_species <- global_traits5 |>
  group_by(Reference, StandardSpeciesName) |>
  summarise(n_records = n(), .groups = "drop") |>
  filter(n_records >= 200) |>
  arrange(desc(n_records))

# View results
print("Datasets with 50+ records per trait:")
print(datasets_50plus_per_trait)

print("Datasets with 50+ records per species:")

print(datasets_50plus_per_species)

print("References with 50+ records per trait:")
print(references_50plus_per_trait)
unique(references_50plus_per_trait$Reference)
print("References with 50+ records per species:")
print(references_50plus_per_species)
unique(references_50plus_per_species$Reference)

# Extract unique references from each dataframe
unique_refs_per_trait <- references_50plus_per_trait |>
  distinct(Reference) |>
  pull(Reference)

unique_refs_per_species <- references_50plus_per_species |>
  distinct(Reference) |>
  pull(Reference)

# Combine the unique values for >= 50 records per species or trait
unique_refs_combined <- unique(c(unique_refs_per_trait, unique_refs_per_species))

# Save to file
write_csv(tibble(Reference = unique_refs_combined), 
          here("data", "derived_data","unique_references_200plus_combined.csv"))

# END OF SCRIPT ----------------------------------------------------------------