##----------------------------------------------------------------------------##
# PAPER 3: BOREAL AND ARCTIC PLANT SPECIES TRAIT SPACES 
# 1.4_global_TRY_cleaned_traits
# This script contains code which cleans raw version of the traits dataset,
# which had not been filtered to only contain records from within the boreal
# forest or tundra - this will be used for RQ1 and RQ2 and the NMDS
##----------------------------------------------------------------------------##

# 1. LOAD DATA -----------------------------------------------------------------

# Setup
library(here)
source(here("scripts", "0_setup.R"))

# Load raw traits data
load(here("data", "raw_data", "try_beatrice.RData"))
try_raw <- try.final.control

# Load distances to biome boundaries
detailed_results <- readRDS(here("data", "derived_data", 
                                 "species_summaries_dist_to_biome_boundary_June25.rds"))

# 2. CLEAN SPECIES NAMES -------------------------------------------------------

## 2.1. Remove morphospecies and subspecies ------------------------------------

# Check for morphospecies
sp_endings <- sum(grepl("\\bsp\\.$|\\bsp\\b", try_raw$AccSpeciesName)) #1478 records
generic_names <- sum(try_raw$AccSpeciesName %in% c("Grass", "Fern", "Unknown")) #43 records
suspect_names <- sum(try_raw$AccSpeciesName %in% c("Hieracium sect.")) #59 records
cf_names <- sum(grepl("\\bcf\\b", try_raw$AccSpeciesName)) #0
aff_names <- sum(grepl("\\baff\\b", try_raw$AccSpeciesName)) #0
single_word_names <- sum(!grepl("\\s", try_raw$AccSpeciesName) & 
                           nchar(try_raw$AccSpeciesName) > 2) #76 records

# Check the single word names
single_examples <- try_raw$AccSpeciesName[!grepl("\\s", try_raw$AccSpeciesName) & 
                                               nchar(try_raw$AccSpeciesName) > 2]
print(single_examples) #"Unknown" "Grass"   "Graminoid", "Fern

# Filter out morphospecies from list
global_traits1 <- try_raw |>
  filter(!grepl("\\bsp\\.$|\\bsp\\b", AccSpeciesName) &
           # remove generic species names
           !AccSpeciesName %in% c("Grass", "Fern", "Unknown", "Graminoid") &
           # remove suspect species names
           !AccSpeciesName %in% c("Hieracium sect.")) |>
  mutate(RawSpeciesName = AccSpeciesName,
         AccSpeciesName = if_else(AccSpeciesName == "Eri sch", "Eriophorum scheuchzeri", AccSpeciesName))

# Check how many records were removed
cat("Original species count:", length(unique(try_raw$AccSpeciesName)), "\n") # 1799
cat("Filtered species count:", length(unique(global_traits1$AccSpeciesName)), "\n") # 1754
cat("Number of morphospecies/generic entries removed:", 
    length(unique(try_raw$AccSpeciesName)) - length(unique(global_traits1$AccSpeciesName)), "\n") # 45

# Check which records have been removed
removed_species <- anti_join(try_raw, global_traits1, 
                             by = "AccSpeciesName")
unique(removed_species$AccSpeciesName)

# Check if any 'sp.' names remain
remaining_sp <- sum(grepl("\\bsp\\.$|\\bsp\\b", global_traits1$AccSpeciesName))
cat("'sp.' names remaining after filtering:", remaining_sp, "\n") #0

# Check if any generic names remain
remaining_generic <- sum(global_traits1$AccSpeciesName %in% c("Grass", "Fern", "Unknown", "Graminoid"))
cat("Generic names remaining after filtering:", remaining_generic, "\n") #0

# Check if any suspect names remain
remaining_suspect <- sum(global_traits1$AccSpeciesName %in% c("Hieracium sect."))
cat("Suspect names remaining after filtering:", remaining_suspect, "\n") #0

# Check for single word names in filtered data
remaining_single_words <- sum(!grepl("\\s", global_traits1$AccSpeciesName) & 
                                nchar(global_traits1$AccSpeciesName) > 2)
cat("Single word names remaining after filtering:", remaining_single_words, "\n") #0

## 2.2. Standardise subspecies/varieties to species level ----------------------

# Extract list of subspecies/varieties
subspecies_varieties <- global_traits1 |>
  filter(grepl("subsp\\.", AccSpeciesName) | 
           grepl("var\\.", AccSpeciesName) |
           grepl("sect\\.", AccSpeciesName))

# Count how many records there were in subspecies group
subspecies_count <- global_traits1 |>
  filter(grepl("subsp\\.", AccSpeciesName))

# Count how many records there were in varieties group
varieties_count <- global_traits1 |>
  filter(grepl("var\\.", AccSpeciesName))

# Count how many records there were in sections group
sections_count <- global_traits1 |>
  filter(grepl("sect\\.", AccSpeciesName))

# Print results
cat("Total records with subspecies/varieties/sections:", length(unique(subspecies_varieties$AccSpeciesName)), "\n") #44
cat("Subspecies (subsp.):", length(unique(subspecies_count$AccSpeciesName)), "\n") #37
cat("Varieties (var.):", length(unique(varieties_count$AccSpeciesName)), "\n") #3
cat("Sections (sect.):", length(unique(sections_count$AccSpeciesName)), "\n") #4

# Standardise records
global_traits2 <- global_traits1 |>
  # replace everything after and including "subsp.", "var.", or "sect." with nothing
  mutate(StandardSpeciesName = str_replace(AccSpeciesName, " (subsp\\.|var\\.|sect\\.).*$", ""))

# Get species names that do not have subspecies/varieties etc
regular_species <- global_traits1 |>
  filter(!grepl("subsp\\.", AccSpeciesName) &
           !grepl("var\\.", AccSpeciesName) &
           !grepl("sect\\.", AccSpeciesName))

# Check the counts
cat("Original filtered species count:", length(unique(global_traits1$AccSpeciesName)), "\n") #1754
cat("Number of subspecies/varieties found:", length(unique(subspecies_varieties$AccSpeciesName)), "\n") #44
cat("Number of regular species:", length(unique(regular_species$AccSpeciesName)), "\n") #1710
cat("Final standardized species count:", length(unique(global_traits2$StandardSpeciesName)), "\n") #1722

# Check if there are any subspecies left
remaining_subspecies <- global_traits2 |>
  filter(grepl("subsp\\.", StandardSpeciesName))
cat("Subspecies (subsp.) remaining after standardization:", length(unique(remaining_subspecies$StandardSpeciesName)), "\n") #0

# Check if any varieties remain in the final list
remaining_varieties <- global_traits2 |>
  filter(grepl("var\\.", StandardSpeciesName))
cat("Varieties (var.) remaining after standardization:",  length(unique(remaining_varieties$StandardSpeciesName)), "\n") #0

# Check if any sections remain in the final list
remaining_sections <- global_traits2 |>
  filter(grepl("sect\\.", StandardSpeciesName))
cat("Sections (sect.) remaining after standardization:", length(unique(remaining_sections$StandardSpeciesName)), "\n") #0

# Check for names with more than 2 words (potential missed subspecies/varieties)
multi_word_names <- global_traits2$StandardSpeciesName[lengths(strsplit(global_traits2$StandardSpeciesName, "\\s+")) > 2]
remaining_multi_word <- length(multi_word_names)
cat("Names with more than 2 words remaining:", remaining_multi_word, "\n") #681
multi_word_names # all have 'x' at the end

# Remove 'x' at the end of the species names
global_traits3  <- global_traits2 |>
  mutate(StandardSpeciesName = str_replace(StandardSpeciesName, "\\sx$", ""))

# Verify there are no multi-word records left
remaining_multi_word_fixed <- global_traits3$StandardSpeciesName[
  lengths(strsplit(global_traits3$StandardSpeciesName, "\\s+")) > 2] #0

# Remove Hieracium (again)
global_traits4 <- global_traits3 |>
  filter(StandardSpeciesName != "Hieracium")
