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
load(here("data", "derived_data", "try_filtered.RData"))

# Inspect data
glimpse(all_filtered_species)

# 2. CLEAN SPECIES LIST --------------------------------------------------------

## 2.1. Extract species list ---------------------------------------------------

# Check for morphospecies
sp_endings <- sum(grepl("\\bsp\\.$|\\bsp\\b", all_filtered_species$SpeciesName)) #25
generic_names <- sum(all_filtered_species$SpeciesName %in% c("Grass", "Fern", "Unknown")) #3
suspect_names <- sum(all_filtered_species$SpeciesName %in% c("Hieracium sect.")) #1
cf_names <- sum(grepl("\\bcf\\b", all_filtered_species$SpeciesName)) #0
aff_names <- sum(grepl("\\baff\\b", all_filtered_species$SpeciesName)) #0
single_word_names <- sum(!grepl("\\s", all_filtered_species$SpeciesName) & 
                           nchar(all_filtered_species$SpeciesName) > 2) #3

# Check the single word names
single_examples <- all_filtered_species$SpeciesName[!grepl("\\s", all_filtered_species$SpeciesName) & 
                                                      nchar(all_filtered_species$SpeciesName) > 2]
print(single_examples) #"Unknown" "Grass"   "Fern"

# Filter out morphospecies from list
filtered_species_list_1 <- all_filtered_species |>
  filter(!grepl("\\bsp\\.$|\\bsp\\b", SpeciesName) &
           # remove generic species names
           !SpeciesName %in% c("Grass", "Fern", "Unknown") &
           # remove suspect species names
           !SpeciesName %in% c("Hieracium sect.")) |>
  mutate(SpeciesName = if_else(SpeciesName == "Eri sch", "Eriophorum scheuchzeri", SpeciesName))

# Check how many records were removed
cat("Original species count:", nrow(all_filtered_species), "\n") # 800
cat("Filtered species count:", nrow(filtered_species_list_1), "\n") # 771
cat("Number of morphospecies/generic entries removed:", 
    nrow(all_filtered_species) - nrow(filtered_species_list_1), "\n") # 29

# Check which records have been removed
removed_species <- anti_join(all_filtered_species, filtered_species_list_1, 
                             by = "SpeciesName")
print(removed_species)

# Check if any 'sp.' names remain
remaining_sp <- sum(grepl("\\bsp\\.$|\\bsp\\b", filtered_species_list_1$SpeciesName))
cat("'sp.' names remaining after filtering:", remaining_sp, "\n") #0

# Check if any generic names remain
remaining_generic <- sum(filtered_species_list_1$SpeciesName %in% c("Grass", "Fern", "Unknown"))
cat("Generic names remaining after filtering:", remaining_generic, "\n") #0

# Check if any suspect names remain
remaining_suspect <- sum(filtered_species_list_1$SpeciesName %in% c("Hieracium sect."))
cat("Suspect names remaining after filtering:", remaining_suspect, "\n") #0

# Check for single word names in filtered data
remaining_single_words <- sum(!grepl("\\s", filtered_species_list_1$SpeciesName) & 
                                nchar(filtered_species_list_1$SpeciesName) > 2)
cat("Single word names remaining after filtering:", remaining_single_words, "\n") #0

## 2.2. Standardise subspecies/varieties to species level ----------------------

# Extract list of subspecies/varieties
subspecies_varieties <- filtered_species_list_1 |>
  filter(grepl("subsp\\.", SpeciesName) | 
           grepl("var\\.", SpeciesName) |
           grepl("sect\\.", SpeciesName))

# Count how many records there were in each group
subspecies_count <- sum(grepl("subsp\\.", filtered_species_list_1$SpeciesName))
varieties_count <- sum(grepl("var\\.", filtered_species_list_1$SpeciesName))
sections_count <- sum(grepl("sect\\.", filtered_species_list_1$SpeciesName))

# Print results
cat("Total records with subspecies/varieties/sections:", nrow(subspecies_varieties), "\n") #22
cat("Subspecies (subsp.):", subspecies_count, "\n") #19
cat("Varieties (var.):", varieties_count, "\n") #2
cat("Sections (sect.):", sections_count, "\n") #1

# Standardise records
standardised_subspecies <- subspecies_varieties |>
  # replace everything after and including "subsp.", "var.", or "sect." with nothing
  mutate(StandardSpeciesName = str_replace(SpeciesName, " (subsp\\.|var\\.|sect\\.).*$", "")) |>
  # use the standardized name instead of the original
  select(SpeciesName = StandardSpeciesName, Biome, SharedAcrossBiomes)

# Get species names that do not have subspecies/varieties etc
regular_species <- filtered_species_list_1 |>
  filter(!grepl("subsp\\.", SpeciesName) &
           !grepl("var\\.", SpeciesName) &
           !grepl("sect\\.", SpeciesName))

# Combine the datasets and remove duplicates
filtered_species_list_2 <- bind_rows(regular_species,
                                     standardised_subspecies) |>
  distinct() |>
  mutate(BiomeCategory = case_when(
    SharedAcrossBiomes ~ "BorealArctic",
    Biome == "Boreal" ~ "Boreal specialist",
    Biome == "Tundra" ~ "Arctic specialist"))

# Check the counts
cat("Original filtered species count:", nrow(filtered_species_list_1), "\n") #771
cat("Number of subspecies/varieties found:", nrow(subspecies_varieties), "\n") #22
cat("Number of regular species:", nrow(regular_species), "\n") #749
cat("Final standardized species count:", nrow(filtered_species_list_2), "\n") #763

# Check if there are any subspecies left
remaining_subspecies <- sum(grepl("subsp\\.", filtered_species_list_2$SpeciesName))
cat("Subspecies (subsp.) remaining after standardization:", remaining_subspecies, "\n") #0

# Check if any varieties remain in the final list
remaining_varieties <- sum(grepl("var\\.", filtered_species_list_2$SpeciesName))
cat("Varieties (var.) remaining after standardization:", remaining_varieties, "\n") #0

# Check if any sections remain in the final list
remaining_sections <- sum(grepl("sect\\.", filtered_species_list_2$SpeciesName))
cat("Sections (sect.) remaining after standardization:", remaining_sections, "\n") #0

# Check for names with more than 2 words (potential missed subspecies/varieties)
multi_word_names <- filtered_species_list_2$SpeciesName[lengths(strsplit(filtered_species_list_2$SpeciesName, "\\s+")) > 2]
remaining_multi_word <- length(multi_word_names)
cat("Names with more than 2 words remaining:", remaining_multi_word, "\n") #5!
multi_word_names # all have 'x' at the end


# Remove 'x' at the end of the species names
filtered_species_list_3  <- filtered_species_list_2 |>
  mutate(SpeciesName = str_replace(SpeciesName, "\\sx$", ""))

# Verify there are no multi-word records left
remaining_multi_word_fixed <- filtered_species_list_3$SpeciesName[
  lengths(strsplit(filtered_species_list_3$SpeciesName, "\\s+")) > 2] #0

## 2.3. Visually inspect cleaned species list ----------------------------------

# Create summary statistics for the table
species_summary <- filtered_species_list_3 |>
  group_by(SpeciesName, BiomeCategory) |>
  summarise(Biomes = paste(unique(Biome), collapse = ", "),
            SharedAcrossBiomes = first(SharedAcrossBiomes),
            .groups = "drop") |>
  arrange(SpeciesName)

# Add some additional useful columns
species_summary <- species_summary |>
  mutate(Genus = sapply(strsplit(SpeciesName, "\\s+"), `[`, 1),
         WordCount = lengths(strsplit(SpeciesName, "\\s+"))) |>
  relocate(SpeciesName, Genus, BiomeCategory, Biomes, SharedAcrossBiomes, WordCount)

# Create interactive table
interactive_table <- datatable(species_summary,
  # table options
  options = list(pageLength = 25,
                 lengthMenu = c(10, 25, 50, 100, -1),
                 scrollX = TRUE,
                 search = list(regex = TRUE, caseInsensitive = TRUE),
                 columnDefs = list(
                   list(width = '300px', targets = 0),
                   list(width = '150px', targets = 1), 
                   list(className = 'dt-center', targets = c(2, 4, 5))),
                 dom = 'Blfrtip',
                 buttons = c('copy', 'csv', 'excel', 'pdf', 'print')),
  # add extensions for buttons
  extensions = 'Buttons',
  # column names
  colnames = c('Species Name' = 'SpeciesName',
               'Genus' = 'Genus',
               'Biome Category' = 'BiomeCategory',
               'Biomes Found' = 'Biomes',
               'Shared Across Biomes' = 'SharedAcrossBiomes',
               'Word Count' = 'WordCount'),
  # add filtering
  filter = 'top',
  # row names
  rownames = FALSE) |>
  # format columns (use the display names, not original column names)
  formatStyle('Biome Category',
              backgroundColor = styleEqual(c('Boreal specialist', 'Arctic specialist', 'BorealArctic'),
                                           c('#d4edda', '#cce7ff', '#fff3cd'))) |>
  formatStyle('Shared Across Biomes',
              backgroundColor = styleEqual(c(TRUE, FALSE),c('#fff3cd', '#f8f9fa')))

# Display the table
interactive_table

# Save cleaned species list
save(filtered_species_list_3, 
     file = here("data", "derived_data", "cleaned_species_list_26June2025.RData"))

# 3. TAXON CHECK ---------------------------------------------------------------

# Load cleaned species list
load(here("data", "derived_data", "cleaned_species_list_26June2025.RData"))

# Get list of species
spp <- unique(filtered_species_list_3$SpeciesName)

# Identify empty species names
empty <- filtered_species_list_3 |> 
  filter(SpeciesName == " ") # no empty cells

# Load WFO data
library(WorldFlora)
WFO.remember('data/WFO_Backbone/classification.csv')

# Create dataframe with unique species names only
sp_names_only <- filtered_species_list_3 |>
  distinct(SpeciesName) # 753 records => 10 duplicate records in the filtered_species_list

# Check which records are duplicated
duplicated_species <- filtered_species_list_3 |>
  group_by(SpeciesName) |>
  summarise(Count = n(),
            .groups = "drop") |>
  filter(Count > 1) |>
  arrange(desc(Count))
print(duplicated_species)


# Run taxon check
taxon_check <- WFO.match(spec.data = sp_names_only,
                         spec.name = "SpeciesName",
                         WFO.file = 'data/WFO_Backbone/classification.csv',
                         no.dates = TRUE)

# Save taxon check to file
write.csv(taxon_check, here("data", "derived_data",
                            "WFO_taxon_check_26June2025.csv"))
