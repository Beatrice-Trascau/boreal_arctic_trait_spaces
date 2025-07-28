##----------------------------------------------------------------------------##
# PAPER 3: BOREAL AND ARCTIC PLANT SPECIES TRAIT SPACES 
# 1.1_TRY_species_list
# This script contains code which prepares the species list for the distance to
# biome boundary calculation
##----------------------------------------------------------------------------##

# 1. SETUP ---------------------------------------------------------------------

# Load filtered TRY data
library(here)
source(here("scripts", "0_setup.R"))
load(here("data", "derived_data", "try_filtered_all_data_27July2025.RData"))

# Inspect data
glimpse(try_filtered_all_data_27July2025)

# 2. CLEAN SPECIES LIST --------------------------------------------------------

## 2.1. Extract species list ---------------------------------------------------

# Check for morphospecies
sp_endings <- sum(grepl("\\bsp\\.$|\\bsp\\b", try_filtered_all_data_27July2025$AccSpeciesName)) #949 records
generic_names <- sum(try_filtered_all_data_27July2025$AccSpeciesName %in% c("Grass", "Fern", "Unknown")) #43 records
suspect_names <- sum(try_filtered_all_data_27July2025$AccSpeciesName %in% c("Hieracium sect.")) #59 records
cf_names <- sum(grepl("\\bcf\\b", try_filtered_all_data_27July2025$AccSpeciesName)) #0
aff_names <- sum(grepl("\\baff\\b", try_filtered_all_data_27July2025$AccSpeciesName)) #0
single_word_names <- sum(!grepl("\\s", try_filtered_all_data_27July2025$AccSpeciesName) & 
                           nchar(try_filtered_all_data_27July2025$AccSpeciesName) > 2) #43 records

# Check the single word names
single_examples <- try_filtered_all_data_27July2025$AccSpeciesName[!grepl("\\s", try_filtered_all_data_27July2025$AccSpeciesName) & 
                                                                     nchar(try_filtered_all_data_27July2025$AccSpeciesName) > 2]
print(single_examples) #"Unknown" "Grass" "Fern"

# Filter out morphospecies from list
filtered_species_list_1 <- try_filtered_all_data_27July2025 |>
  filter(!grepl("\\bsp\\.$|\\bsp\\b", AccSpeciesName) &
           # remove generic species names
           !AccSpeciesName %in% c("Grass", "Fern", "Unknown") &
           # remove suspect species names
           !AccSpeciesName %in% c("Hieracium sect.")) |>
  mutate(AccSpeciesName = if_else(AccSpeciesName == "Eri sch", "Eriophorum scheuchzeri", AccSpeciesName))

# Check how many records were removed
cat("Original species count:", length(unique(try_filtered_all_data_27July2025$AccSpeciesName)), "\n") # 800
cat("Filtered species count:", length(unique(filtered_species_list_1$AccSpeciesName)), "\n") # 770
cat("Number of morphospecies/generic entries removed:", 
    length(unique(try_filtered_all_data_27July2025$AccSpeciesName)) - length(unique(filtered_species_list_1$AccSpeciesName)), "\n") # 30

# Check which records have been removed
removed_species <- anti_join(try_filtered_all_data_27July2025, filtered_species_list_1, 
                             by = "AccSpeciesName")
unique(removed_species$AccSpeciesName)

# Check if any 'sp.' names remain
remaining_sp <- sum(grepl("\\bsp\\.$|\\bsp\\b", filtered_species_list_1$AccSpeciesName))
cat("'sp.' names remaining after filtering:", remaining_sp, "\n") #0

# Check if any generic names remain
remaining_generic <- sum(filtered_species_list_1$SpeciesName %in% c("Grass", "Fern", "Unknown"))
cat("Generic names remaining after filtering:", remaining_generic, "\n") #0

# Check if any suspect names remain
remaining_suspect <- sum(filtered_species_list_1$AccSpeciesName %in% c("Hieracium sect."))
cat("Suspect names remaining after filtering:", remaining_suspect, "\n") #0

# Check for single word names in filtered data
remaining_single_words <- sum(!grepl("\\s", filtered_species_list_1$AccSpeciesName) & 
                                nchar(filtered_species_list_1$AccSpeciesName) > 2)
cat("Single word names remaining after filtering:", remaining_single_words, "\n") #0

## 2.2. Standardise subspecies/varieties to species level ----------------------

# Extract list of subspecies/varieties
subspecies_varieties <- filtered_species_list_1 |>
  filter(grepl("subsp\\.", AccSpeciesName) | 
           grepl("var\\.", AccSpeciesName) |
           grepl("sect\\.", AccSpeciesName))

# Count how many records there were in subspecies group
subspecies_count <- filtered_species_list_1 |>
  filter(grepl("subsp\\.", AccSpeciesName))

# Count how many records there were in varieties group
varieties_count <- filtered_species_list_1 |>
  filter(grepl("var\\.", AccSpeciesName))

# Count how many records there were in sections group
sections_count <- filtered_species_list_1 |>
  filter(grepl("sect\\.", AccSpeciesName))

# Print results
cat("Total records with subspecies/varieties/sections:", length(unique(subspecies_varieties$AccSpeciesName)), "\n") #22
cat("Subspecies (subsp.):", length(unique(subspecies_count$AccSpeciesName)), "\n") #19
cat("Varieties (var.):", length(unique(varieties_count$AccSpeciesName)), "\n") #2
cat("Sections (sect.):", length(unique(sections_count$AccSpeciesName)), "\n") #1

# Standardise records
filtered_species_list_2 <- filtered_species_list_1 |>
  # replace everything after and including "subsp.", "var.", or "sect." with nothing
  mutate(StandardSpeciesName = str_replace(AccSpeciesName, " (subsp\\.|var\\.|sect\\.).*$", ""))

# Get species names that do not have subspecies/varieties etc
regular_species <- filtered_species_list_1 |>
  filter(!grepl("subsp\\.", AccSpeciesName) &
           !grepl("var\\.", AccSpeciesName) &
           !grepl("sect\\.", AccSpeciesName))

# Check the counts
cat("Original filtered species count:", length(unique(filtered_species_list_1$AccSpeciesName)), "\n") #770
cat("Number of subspecies/varieties found:", length(unique(subspecies_varieties$AccSpeciesName)), "\n") #22
cat("Number of regular species:", length(unique(regular_species$AccSpeciesName)), "\n") #748
cat("Final standardized species count:", length(unique(filtered_species_list_2$StandardSpeciesName)), "\n") #754

# Check if there are any subspecies left
remaining_subspecies <- filtered_species_list_2 |>
  filter(grepl("subsp\\.", StandardSpeciesName))
cat("Subspecies (subsp.) remaining after standardization:", length(unique(remaining_subspecies$StandardSpeciesName)), "\n") #0

# Check if any varieties remain in the final list
remaining_varieties <- filtered_species_list_2 |>
  filter(grepl("var\\.", StandardSpeciesName))
cat("Varieties (var.) remaining after standardization:",  length(unique(remaining_varieties$StandardSpeciesName)), "\n") #0

# Check if any sections remain in the final list
remaining_sections <- filtered_species_list_2 |>
  filter(grepl("sect\\.", StandardSpeciesName))
cat("Sections (sect.) remaining after standardization:", length(unique(remaining_sections$StandardSpeciesName)), "\n") #0

# Check for names with more than 2 words (potential missed subspecies/varieties)
multi_word_names <- filtered_species_list_2$StandardSpeciesName[lengths(strsplit(filtered_species_list_2$StandardSpeciesName, "\\s+")) > 2]
remaining_multi_word <- length(multi_word_names)
cat("Names with more than 2 words remaining:", remaining_multi_word, "\n") #22
multi_word_names # all have 'x' at the end

# Remove 'x' at the end of the species names
filtered_species_list_3  <- filtered_species_list_2 |>
  mutate(StandardSpeciesName = str_replace(StandardSpeciesName, "\\sx$", ""))

# Verify there are no multi-word records left
remaining_multi_word_fixed <- filtered_species_list_3$StandardSpeciesName[
  lengths(strsplit(filtered_species_list_3$StandardSpeciesName, "\\s+")) > 2] #0

# Remove Hieracium (again)
filtered_species_list_4 <- filtered_species_list_3 |>
  filter(StandardSpeciesName != "Hieracium")

# Save cleaned species list
# save(filtered_species_list_4,
#      file = here("data", "derived_data", "cleaned_species_list_27July2025.RData"))

# 3. TAXON CHECK ---------------------------------------------------------------

## 3.1. Run taxon check --------------------------------------------------------

# Load cleaned species list
load(here("data", "derived_data", "cleaned_species_list_27July2025.RData"))

# Get list of species
spp <- unique(filtered_species_list_4$StandardSpeciesName)

# Identify empty species names
empty <- filtered_species_list_4 |> 
  filter(StandardSpeciesName == " ") # no empty cells

# Load WFO data
library(WorldFlora)
WFO.remember('data/WFO_Backbone/classification.csv')

# Create dataframe with unique species names only
sp_names_only <- filtered_species_list_4 |>
  distinct(StandardSpeciesName) # 752 records 

# Run taxon check
taxon_check <- WFO.match(spec.data = sp_names_only,
                         spec.name = "StandardSpeciesName",
                         WFO.file = 'data/WFO_Backbone/classification.csv',
                         no.dates = TRUE)

# Save taxon check to file
# write.csv(taxon_check, here("data", "derived_data",
#                             "WFO_taxon_check_26June2025.csv"))

## 3.2. Check record count mismatch --------------------------------------------

# Load taxon check
taxon_check <- read.csv(here("data", "derived_data",
                             "WFO_taxon_check_26June2025.csv"))

# Spp list = 752 but taxon_check = 1072

# Check for duplicates in taxon_check
duplicates_in_taxon_check <- taxon_check |>
  group_by(SpeciesName.ORIG) |>
  summarise(count = n(), .groups = "drop") |>
  filter(count > 1) |>
  arrange(desc(count))
nrow(duplicates_in_taxon_check) # 191 species appear multiple times in taxon check

# Check why there are multiple matches
multi_match_analysis <- taxon_check |>
  group_by(SpeciesName.ORIG) |>
  summarise(match_count = n(),
            unique_scientific_names = n_distinct(scientificName),
            fuzzy_matches = sum(Fuzzy, na.rm = TRUE),
            accepted_matches = sum(taxonomicStatus == "Accepted", na.rm = TRUE),
            .groups = "drop") |>
  filter(match_count > 1)
print(head(multi_match_analysis, 10))

# Check sequence information
sequence_info <- taxon_check |>
  select(SpeciesName.ORIG, scientificName, OriSeq, Subseq, Matched, Fuzzy) |>
  arrange(OriSeq, Subseq) |>
  head(20)

print(sequence_info)

## 3.3. Check taxon match issues -----------------------------------------------

# Check for unmatched records
unmatched <- taxon_check |> 
  filter(Matched == FALSE | is.na(Matched)) # 2 records: Casteleja occidens & Sausarrea angustifolium
if(nrow(unmatched) > 0) {
  cat("These species names were not found in WFO:\n")
  print(unmatched |> select(SpeciesName, SpeciesName.ORIG) |> head(10))
  if(nrow(unmatched) > 10) cat("... and", nrow(unmatched) - 10, "more\n")
}

# Check for fuzzy matches
fuzzy_matches <- taxon_check |> 
  filter(Fuzzy == TRUE) # 38 records with fuzzy match
if(nrow(fuzzy_matches) > 0) {
  fuzzy_comparison <- fuzzy_matches |> 
    select(Original = SpeciesName.ORIG, Matched = scientificName, Distance = Fuzzy.dist) |>
    head(10)
  print(fuzzy_comparison)
  if(nrow(fuzzy_matches) > 10) cat("... and", nrow(fuzzy_matches) - 10, "more\n")
}

# Check for hybrids
hybrids <- taxon_check |> 
  filter(Hybrid == TRUE)
if(nrow(hybrids) > 0) {
  cat("These are identified as hybrids:\n")
  print(hybrids |> select(SpeciesName.ORIG, scientificName) |> head(10))
} # no hybrids

# Check for records with brackets
brackets <- taxon_check |> filter(Brackets.detected == TRUE)
if(nrow(brackets) > 0) {
  cat("Bracket examples:\n")
  print(brackets |> select(SpeciesName.ORIG) |> head(5))
} # no records with brackets

# Check for records with numbers
numbers <- taxon_check |> filter(Number.detected == TRUE)
if(nrow(numbers) > 0) {
  cat("Number examples:\n")
  print(numbers |> select(SpeciesName.ORIG) |> head(5)) 
} # no records with numbers

## 3.4. Fix taxon check issues -------------------------------------------------

# Fix unmatched records (i.e. misspelled records)
corrected_species_list_1 <- filtered_species_list_5 |>
  mutate(SpeciesName_original = SpeciesName, # keep original for reference
         SpeciesName = case_when(SpeciesName == "Casteleja occidens" ~ "Castilleja occidentalis",
                                 SpeciesName == "Sausarrea angustifolium" ~ "Saussurea angustifolia",
                                 .default = SpeciesName))

# Review fuzzy matches
fuzzy_review <- fuzzy_matches |>
  select(Original = SpeciesName.ORIG, Matched = scientificName, 
         Distance = Fuzzy.dist, Status = taxonomicStatus) |>
  arrange(Distance)
# View(fuzzy_review)

# Manual fix for fuzzy matches
corrected_species_list_2 <- corrected_species_list_1 |>
  mutate(SpeciesName = case_when(
    SpeciesName == "Spirodela polyrrhiza" ~ "Spirodela polyrhiza",
    SpeciesName == "Salix doniana" ~ "Salix purpurea",
    SpeciesName == "Silene samojedora" ~ "Silene samojedorum",
    SpeciesName == "Salix myrtifolia" ~ "Salix myrtillifolia",
    SpeciesName == "Calamagrostis purpuras" ~ "Calamagrostis purpurea",
    SpeciesName == "Pedicularis vertisilata" ~ "Pedicularis verticillata",
    SpeciesName == "Peticites frigidus" ~ "Petasites frigidus",
    SpeciesName == "Senecio atropurpuris" ~ "Senecio atropurpureus",
    SpeciesName == "Gentia glauca" ~ "Gentiana glauca",
    SpeciesName == "Polemonium acutifolium" ~ "Polemonium acutiflorum",
    SpeciesName == "Rumex lapponum" ~ "Rumex lapponicus",
    SpeciesName == "Pedicularis vertillis" ~ "Pedicularis verticillata",
    SpeciesName == "Sabulina rossii" ~ "Sabulina rosei",
    SpeciesName == "Echinops crispus" ~ "Echinops ritro",
    SpeciesName == "Salix fuscenses" ~ "Salix fuscescens",
    SpeciesName == "Salix laponicum" ~ "Salix lapponum",
    SpeciesName == "Senecio atropupuris" ~ "Senecio atropurpureus",
    SpeciesName == "Salix argyocarpon" ~ "Salix argyrocarpa",
    SpeciesName == "Salix herbaceae-polaris" ~ "Salix herbacea",
    .default = SpeciesName))

# Check for duplicate records
duplicated_species <- corrected_species_list_2 |>
  group_by(SpeciesName) |>
  summarise(Count = n(),
            .groups = "drop") |>
  filter(Count > 1) |>
  arrange(desc(Count)) # all duplicates look normal

# Remove duplicates
corrected_species_list <- corrected_species_list_2 |>
  distinct(SpeciesName, .keep_all = TRUE) # 739 species

# Write corrected species list to file
# save(corrected_species_list, file = here("data", "derived_data", 
#                                          "corrected_species_list_27June2025.RData"))

# END OF SCRIPT ----------------------------------------------------------------  