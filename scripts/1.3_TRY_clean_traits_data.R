##----------------------------------------------------------------------------##
# PAPER 3: BOREAL AND ARCTIC PLANT SPECIES TRAIT SPACES 
# 1.3_TRY_clean_traits_data
# This script contains code which cleans the species names in the traits data
##----------------------------------------------------------------------------##

# 1. SETUP ---------------------------------------------------------------------

# Load filtered TRY data
library(here)
source(here("scripts", "0_setup.R"))
load(here("data", "derived_data", "try_filtered_all_data_27July2025.RData"))
raw_traits <- try_filtered_all_data_27July2025

# Inspect data
glimpse(raw_traits)

# 2. CLEAN SPECIES LIST --------------------------------------------------------

## 2.1. Remove morphospecies and suspect names ---------------------------------

# Check for morphospecies
sp_endings <- sum(grepl("\\bsp\\.$|\\bsp\\b", raw_traits$AccSpeciesName)) #949 records
generic_names <- sum(raw_traits$AccSpeciesName %in% c("Grass", "Fern", "Unknown")) #43 records
suspect_names <- sum(raw_traits$AccSpeciesName %in% c("Hieracium sect.")) #59 records
cf_names <- sum(grepl("\\bcf\\b", raw_traits$AccSpeciesName)) #0
aff_names <- sum(grepl("\\baff\\b", raw_traits$AccSpeciesName)) #0
single_word_names <- sum(!grepl("\\s", raw_traits$AccSpeciesName) & 
                           nchar(raw_traits$AccSpeciesName) > 2) #43 records

# Check the single word names
single_examples <- raw_traits$AccSpeciesName[!grepl("\\s", raw_traits$AccSpeciesName) & 
                                                      nchar(raw_traits$AccSpeciesName) > 2]
print(single_examples) #"Unknown" "Grass"   "Fern"

# Filter out morphospecies from list
filtered_traits_1 <- raw_traits |>
  filter(!grepl("\\bsp\\.$|\\bsp\\b", AccSpeciesName) &
           # remove generic species names
           !AccSpeciesName %in% c("Grass", "Fern", "Unknown") &
           # remove suspect species names
           !AccSpeciesName %in% c("Hieracium sect.")) |>
  mutate(RawSpeciesName = AccSpeciesName,
         AccSpeciesName = if_else(AccSpeciesName == "Eri sch", "Eriophorum scheuchzeri", AccSpeciesName))

# Check how many records were removed
cat("Original species count:", length(unique(raw_traits$AccSpeciesName)), "\n") # 800
cat("Filtered species count:", length(unique(filtered_traits_1$AccSpeciesName)), "\n") # 770
cat("Number of morphospecies/generic entries removed:", 
    length(unique(raw_traits$AccSpeciesName)) - length(unique(filtered_traits_1$AccSpeciesName)), "\n") # 30

# Check which records have been removed
removed_species <- anti_join(raw_traits, filtered_traits_1, 
                             by = "AccSpeciesName")
unique(removed_species$AccSpeciesName)

# Check if any 'sp.' names remain
remaining_sp <- sum(grepl("\\bsp\\.$|\\bsp\\b", filtered_traits_1$AccSpeciesName))
cat("'sp.' names remaining after filtering:", remaining_sp, "\n") #0

# Check if any generic names remain
remaining_generic <- sum(filtered_traits_1$AccSpeciesName %in% c("Grass", "Fern", "Unknown"))
cat("Generic names remaining after filtering:", remaining_generic, "\n") #0

# Check if any suspect names remain
remaining_suspect <- sum(filtered_traits_1$AccSpeciesName %in% c("Hieracium sect."))
cat("Suspect names remaining after filtering:", remaining_suspect, "\n") #0

# Check for single word names in filtered data
remaining_single_words <- sum(!grepl("\\s", filtered_traits_1$AccSpeciesName) & 
                                nchar(filtered_traits_1$AccSpeciesName) > 2)
cat("Single word names remaining after filtering:", remaining_single_words, "\n") #0

## 2.2. Standardise subspecies/varieties to species level ----------------------

# Extract list of subspecies/varieties
subspecies_varieties <- filtered_traits_1 |>
  filter(grepl("subsp\\.", AccSpeciesName) | 
           grepl("var\\.", AccSpeciesName) |
           grepl("sect\\.", AccSpeciesName))

# Count how many records there were in subspecies group
subspecies_count <- filtered_traits_1 |>
  filter(grepl("subsp\\.", AccSpeciesName))

# Count how many records there were in varieties group
varieties_count <- filtered_traits_1 |>
  filter(grepl("var\\.", AccSpeciesName))

# Count how many records there were in sections group
sections_count <- filtered_traits_1 |>
  filter(grepl("sect\\.", AccSpeciesName))

# Print results
cat("Total records with subspecies/varieties/sections:", length(unique(subspecies_varieties$AccSpeciesName)), "\n") #22
cat("Subspecies (subsp.):", length(unique(subspecies_count$AccSpeciesName)), "\n") #19
cat("Varieties (var.):", length(unique(varieties_count$AccSpeciesName)), "\n") #2
cat("Sections (sect.):", length(unique(sections_count$AccSpeciesName)), "\n") #1

# Standardise records
filtered_traits_2 <- filtered_traits_1 |>
  # replace everything after and including "subsp.", "var.", or "sect." with nothing
  mutate(StandardSpeciesName = str_replace(AccSpeciesName, " (subsp\\.|var\\.|sect\\.).*$", ""))

# Get species names that do not have subspecies/varieties etc
regular_species <- filtered_traits_1 |>
  filter(!grepl("subsp\\.", AccSpeciesName) &
           !grepl("var\\.", AccSpeciesName) &
           !grepl("sect\\.", AccSpeciesName))

# Check the counts
cat("Original filtered species count:", length(unique(filtered_traits_1$AccSpeciesName)), "\n") #770
cat("Number of subspecies/varieties found:", length(unique(subspecies_varieties$AccSpeciesName)), "\n") #22
cat("Number of regular species:", length(unique(regular_species$AccSpeciesName)), "\n") #748
cat("Final standardized species count:", length(unique(filtered_traits_2$StandardSpeciesName)), "\n") #754

# Check if there are any subspecies left
remaining_subspecies <- filtered_traits_2 |>
  filter(grepl("subsp\\.", StandardSpeciesName))
cat("Subspecies (subsp.) remaining after standardization:", length(unique(remaining_subspecies$StandardSpeciesName)), "\n") #0

# Check if any varieties remain in the final list
remaining_varieties <- filtered_traits_2 |>
  filter(grepl("var\\.", StandardSpeciesName))
cat("Varieties (var.) remaining after standardization:",  length(unique(remaining_varieties$StandardSpeciesName)), "\n") #0

# Check if any sections remain in the final list
remaining_sections <- filtered_traits_2 |>
  filter(grepl("sect\\.", StandardSpeciesName))
cat("Sections (sect.) remaining after standardization:", length(unique(remaining_sections$StandardSpeciesName)), "\n") #0

# Check for names with more than 2 words (potential missed subspecies/varieties)
multi_word_names <- filtered_traits_2$StandardSpeciesName[lengths(strsplit(filtered_traits_2$StandardSpeciesName, "\\s+")) > 2]
remaining_multi_word <- length(multi_word_names)
cat("Names with more than 2 words remaining:", remaining_multi_word, "\n") #22
multi_word_names # all have 'x' at the end

# Remove 'x' at the end of the species names
filtered_traits_3  <- filtered_traits_2 |>
  mutate(StandardSpeciesName = str_replace(StandardSpeciesName, "\\sx$", ""))

# Verify there are no multi-word records left
remaining_multi_word_fixed <- filtered_traits_3$StandardSpeciesName[
  lengths(strsplit(filtered_traits_3$StandardSpeciesName, "\\s+")) > 2] #0

# Remove Hieracium (again)
filtered_traits_4 <- filtered_traits_3 |>
  filter(StandardSpeciesName != "Hieracium")

## 2.3. Clean species names ----------------------------------------------------

# This is done in accordance with the taxonomic check from script 1.2 which was
# run on the list of species from this dataframe
filtered_traits_5 <- filtered_traits_4 |>
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

# 3. CLEAN TRAITS VALUES -------------------------------------------------------

## 3.1. Solve character issue --------------------------------------------------

# Check structure of traits df
glimpse(filtered_traits_5) # trait values are characters - why?


# Find out why the column OrigValueStr was read in as a character
filtered_traits_5 |>
  # get the number of rows that have letters = 1015
  summarise(has_letters = sum(grepl("[a-zA-Z]", OrigValueStr)),
            # number of rows with ranges = 1012
            has_ranges = sum(grepl("-", OrigValueStr)),
            # number of rows with special characters = 1015
            has_special_characters = sum(grepl("[^0-9.-]", OrigValueStr)),
            total_rows = n())

# Inspect rows with letters
rows_with_letters <- filtered_traits_5 |>
  filter(grepl("[a-zA-Z]", OrigValueStr))

# Inspect rows with ranges
rows_with_ranges <- filtered_traits_5 |>
  filter(grepl("-", OrigValueStr))

# Inspect rows with special characters
rows_with_characters <- filtered_traits_5 |>
  filter(grepl("[^0-9.-]", OrigValueStr))

# Rows are the same across rows_with_letters, rows_with_ranges, rows_with_characters
# 2 types of issues: *cm at the end of values and scientific notations
# remove trailing text after values for the rows with *cm
# the scientific notations will be converted correctly to numerical by as.numeric()
filtered_traits_6 <- filtered_traits_5 |>
  mutate(OrigValueStr = str_remove(OrigValueStr, " cm\\*$"),
         CleanedTraitValue = as.numeric(OrigValueStr))

## 3.2. Check consistency of units of measurement ------------------------------

# Check if there are any rows where OrigUnitStr and UnitName are different
mismatched_units <- filtered_traits_6 |>
  filter(OrigUnitStr != UnitName) |>
  select(StandardSpeciesName, TraitName, OrigValueStr, OrigUnitStr,
         UnitName, CleanedTraitValue, Dataset)

# Create an interactive table to check mismatches 
unit_explorer_table <- datatable(
  mismatched_units,
  
  # Basic options
  options = list(pageLength = 25,
                 scrollX = TRUE,
                 scrollY = "500px",
                 dom = 'Bfrtip',
                 buttons = c('copy', 'csv', 'excel'),
                 search = list(regex = TRUE, caseInsensitive = TRUE),
    
    # Column definitions for better formatting
    columnDefs = list(list(width = '150px', targets = c(0, 1)),  # Species and Trait names
      list(width = '100px', targets = c(2, 3, 4, 5)),  # Value columns
      list(width = '200px', targets = c(6)),  # Dataset and 
      list(className = 'dt-center', targets = c(2, 3, 4, 5)))),  # Center align values
   
  # Add filter widgets
  filter = 'top',
  
  # Extensions for additional functionality
  extensions = 'Buttons',
  
  # Column names
  colnames = c('Species' = 'StandardSpeciesName',
               'Trait' = 'TraitName', 
               'Original Value' = 'OrigValueStr',
               'Original Unit' = 'OrigUnitStr',
               'Standard Unit' = 'UnitName',
               'Cleaned Value' = 'CleanedTraitValue',
               'Dataset' = 'Dataset')) |>
  
  # Format specific columns with colors
  formatStyle('Original Unit', backgroundColor = '#ffebee',
              fontWeight = 'bold', color = '#c62828') |>
  
  formatStyle('Standard Unit', backgroundColor = '#e8f5e8',
              fontWeight = 'bold', color = '#2e7d32')

# Display the table
unit_explorer_table

# Check what to filter for
unique(mismatched_units$OrigUnitStr)

# Save cleaned traits dataframe
save(traits_cleaned_species_names, file = here("data", "derived_data",
                                               "TRY_traits_cleaned_species_names.RData"))
# END OF SCRIPT ----------------------------------------------------------------