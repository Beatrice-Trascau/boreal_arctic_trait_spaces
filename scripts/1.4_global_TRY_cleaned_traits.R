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

# Load CAFF quality check
caff_check <- read.xlsx(here("data", "derived_data", "caff_quality_check.xlsx"),
                        sheet = 1, skipEmptyRows = TRUE)

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

## 2.3. Clean species names ----------------------------------------------------

# This is done in accordance with the taxonomic check from script 1.2 which was
# run on the list of species from this dataframe
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

# 3. ADD DISTANCE TO BIOME BOUNDARY  -------------------------------------------

# Remove Elodea canadensis & hybrids
biome_boundaries <- detailed_results |>
  filter(!species == "Elodea canadensis")

# Remove hybrids from the distance to biome boundary df
biome_boundaries <- biome_boundaries |>
  filter(!str_detect(species, " × "))

# Combine the dataframes
global_cleaned_traits1 <- global_traits5 |>
  left_join(biome_boundaries, by = c("StandardSpeciesName" = "species"))

# Rename column with mean distance to biome boundary to reflect the fact that
# it is a species-level mean
global_cleaned_traits2 <- global_cleaned_traits1 |>
  rename(species_level_mean_distance_km = mean_distance_km)

# Remove species that do not have a species level mean distance to biome boundary
global_cleaned_traits3 <- global_cleaned_traits2 |>
  filter(!is.na(species_level_mean_distance_km))

# Save dataframe
save(global_cleaned_traits3, file = here("data", "derived_data", "global_cleaned_traits_August2025.RData"))
write.csv(global_cleaned_traits3, file = here("data", "derived_data", "global_cleaned_traits_August2025.csv"))

# Extract species list from global_cleaned_traits3
species_list_boreal_tundra <- global_cleaned_traits3 |>
  distinct(StandardSpeciesName, .keep_all = TRUE) |>
  mutate(biome_category = case_when(species_level_mean_distance_km > 0 ~ "boreal",
                                    species_level_mean_distance_km < 0 ~ "tundra",
                                    TRUE ~ "boundary")) |>
  select(c(Dataset, StandardSpeciesName, biome_category))

# Write to file
write.csv(species_list_boreal_tundra, here("data", "derived_data",
                                           "species_list_boreal_tundra_Sept2025.csv"))

# 4. CLASSIFICATION QUALITY CONTROL --------------------------------------------

# Check columns of dfs
colnames(global_cleaned_traits3)
colnames(caff_check)

# Change species name column to match the cleaned traits
caff_check <- caff_check |>
  rename(StandardSpeciesName = SPECIES_CLEAN)

# Add CAFF classification to the dataframe of cleaned traits
caff_global_cleaned_traits <- global_cleaned_traits3 |>
  left_join(caff_check |> select(StandardSpeciesName, final.category), 
            by = "StandardSpeciesName")

# Create colum with the corrected classification
caff_global_cleaned_traits2 <- caff_global_cleaned_traits |>
  # create new column with the caff category and log-transform values
  mutate(caff_biome_category = final.category) |>
  filter(!caff_biome_category == "remove")

# 4. NMDS ----------------------------------------------------------------------

# Check filtering
trial <- caff_global_cleaned_traits2 |>
  group_by(StandardSpeciesName, TraitNameNew) |>
  mutate(number_of_records = length(StdValue)) 

# Calculate species-level means for the traits
traits_median <- caff_global_cleaned_traits2 |>
  group_by(StandardSpeciesName, TraitNameNew) |>
  mutate(number_of_records = length(StdValue)) |>
  filter(number_of_records > 4) |>
  summarise(max_trait_value = max(StdValue, na.rm = TRUE),
            MedianTraitValue = median(StdValue, na.rm = TRUE),
            .groups = "drop")

# Add distance to biome boundaries and growth form 
traits_median_df <- traits_median |>
  left_join(biome_boundaries, by = c("StandardSpeciesName" = "species")) |>
  rename(species_level_mean_distance_km = mean_distance_km) |>
  filter(!is.na(species_level_mean_distance_km), !is.na(MedianTraitValue)) |>
  mutate(log_median_trait_value = log(MedianTraitValue))

# Categorise species as either boreal or tundra specialists
species_biome_classification <- traits_median_df |>
  distinct(StandardSpeciesName, species_level_mean_distance_km) |>
  mutate(pipeline_biome_category = case_when(species_level_mean_distance_km > 0 ~ "boreal",
                                    species_level_mean_distance_km < 0 ~ "tundra",
                                    TRUE ~ "boundary"))  # exactly at boundary

# Save list of tundra and boreal species 
# boreal_classification <- species_biome_classification |>
#   filter(biome_category == "boreal")
# write.csv(boreal_classification, here("data", "derived_data", 
#                                       "boreal_species_classification.csv"))

# Create trait values per species in wide format for NMDS
trait_matrix <- traits_median_df |>
  select(StandardSpeciesName, TraitNameNew, log_median_trait_value) |>
  pivot_wider(names_from = TraitNameNew, 
              values_from = log_median_trait_value) |>
  column_to_rownames("StandardSpeciesName")

# Remove Leaf C: N ratio
trait_matrix <- trait_matrix |>
  select(-LeafCN)

# Remove species with missing data for any trait
complete_trait_matrix <- trait_matrix[complete.cases(trait_matrix), ]

# Check which traits are available
View(complete_trait_matrix)

## 5.2. Run NMDS ---------------------------------------------------------------

# Set seed for NMDS
set.seed(84145)

# Run NMDS
nmds_result <- metaMDS(complete_trait_matrix, 
                       distance = "euclidean",
                       k = 3,
                       trymax = 100)

# Check stress value per dimension
dimcheck_out <- dimcheckMDS(complete_trait_matrix,
                            distance = "euclidean",
                            k = 4)

## 5.3. Plot output ------------------------------------------------------------

# Extract NMDS scores
nmds_scores <- as.data.frame(nmds_result$points)
nmds_scores$StandardSpeciesName <- rownames(nmds_scores)

# Extract biome classification
caff_biomes <- caff_global_cleaned_traits2 |>
  select(StandardSpeciesName, caff_biome_category) |>
  distinct(StandardSpeciesName, .keep_all = TRUE)

# Add biome classification
nmds_plot_data <- nmds_scores |>
  left_join(caff_biomes, by = "StandardSpeciesName") |>
  filter(!is.na(caff_biome_category))


# Get polygon hull from borel and tundra
boreal_scores <- nmds_plot_data[nmds_plot_data$caff_biome_category == "boreal", ][chull(nmds_plot_data[nmds_plot_data$caff_biome_category == 
                                                                                                    "boreal", c("MDS1", "MDS2")]), ]  # hull values for boreal
tundra_scores <- nmds_plot_data[nmds_plot_data$caff_biome_category == "tundra", ][chull(nmds_plot_data[nmds_plot_data$caff_biome_category == 
                                                                                                    "tundra", c("MDS1", "MDS2")]), ]  # hull values for grp B
# Combine into one hull
hull_data <- rbind(boreal_scores, tundra_scores)  #combine grp.a and grp.b

# Create NMDS plot
(nmds_plot <- ggplot(nmds_plot_data, aes(x = MDS1, y = MDS2, color = caff_biome_category)) +
  geom_polygon(data=hull_data,aes(x=MDS1,y=MDS2,fill=caff_biome_category,group=caff_biome_category),alpha=0.30) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = c("boreal" = "darkgreen", "tundra" = "black"),
                     name = "Biome") +
  labs(x = "NMDS1", y = "NMDS2") +
  theme_classic() +
  theme(legend.position = "bottom") +
  # add stress value as subtitle
  labs(subtitle = paste("Based on", ncol(trait_matrix), "traits,", nrow(nmds_plot_data), "species")))
  
# Fit trait vectors to the NMDS ordination
trait_fit <- envfit(nmds_result, complete_trait_matrix, permutations = 999, na.rm = TRUE)

# Extract the vectors
trait_vectors <- as.data.frame(scores(trait_fit, "vectors"))
trait_vectors$trait <- rownames(trait_vectors)

# Check significance of vectors
print(trait_fit)

# Plot NMDS with hull and vectors
(nmds_plot_with_proper_vectors <- nmds_plot +
  geom_segment(data = trait_vectors, 
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "black", 
               size = 1,
               inherit.aes = FALSE) +
  geom_text(data = trait_vectors,
            aes(x = NMDS1 * 1.1, y = NMDS2 * 1.1, label = trait),
            color = "black", 
            size = 3, 
            fontface = "bold",
            inherit.aes = FALSE))

# 5. NMDS WITHOUT SEED MASS ----------------------------------------------------

# Remove seed mass from list of traits
seedless_traits <- trait_matrix |>
  select(-SeedMass)

# Remove species with missing data for any trait
seedless_complete_trait_matrix <- seedless_traits[complete.cases(seedless_traits), ]

# Set seed for NMDS
set.seed(12345)

# Run NMDS
seedless_nmds_result <- metaMDS(seedless_complete_trait_matrix, 
                       distance = "euclidean",
                       k = 3,
                       trymax = 100)

# Check stress value per dimension
seedless_dimcheck_out <- dimcheckMDS(seedless_complete_trait_matrix,
                            distance = "euclidean",
                            k = 3)

# Extract NMDS scores
seedless_nmds_scores <- as.data.frame(seedless_nmds_result$points)
seedless_nmds_scores$StandardSpeciesName <- rownames(seedless_nmds_scores)

# Add biome classification
seedless_nmds_plot_data <- seedless_nmds_scores |>
  left_join(caff_biomes, by = "StandardSpeciesName") |>
  filter(!is.na(caff_biome_category))

# Get polygon hull from borel and tundra
seedless_boreal_scores <- seedless_nmds_plot_data[seedless_nmds_plot_data$caff_biome_category == "boreal", ][chull(seedless_nmds_plot_data[seedless_nmds_plot_data$caff_biome_category == 
                                                                                                         "boreal", c("MDS1", "MDS2")]), ]  # hull values for boreal
seedless_tundra_scores <- seedless_nmds_plot_data[seedless_nmds_plot_data$caff_biome_category == "tundra", ][chull(seedless_nmds_plot_data[seedless_nmds_plot_data$caff_biome_category == 
                                                                                                         "tundra", c("MDS1", "MDS2")]), ]  # hull values for grp B
# Combine into one hull
seedless_hull_data <- rbind(seedless_boreal_scores, seedless_tundra_scores) 

# Create NMDS plot
(seedless_nmds_plot <- ggplot(seedless_nmds_plot_data, aes(x = MDS1, y = MDS2, color = caff_biome_category)) +
    geom_polygon(data=seedless_hull_data,aes(x=MDS1,y=MDS2,fill=caff_biome_category,group=caff_biome_category),alpha=0.30) +
    geom_point(size = 2, alpha = 0.7) +
    scale_color_manual(values = c("boreal" = "darkgreen", "tundra" = "black"),
                       name = "Biome") +
    labs(x = "NMDS1", y = "NMDS2") +
    theme_classic() +
    theme(legend.position = "bottom") +
    # add stress value as subtitle
    labs(subtitle = paste("Based on", ncol(seedless_complete_trait_matrix), "traits,", nrow(seedless_nmds_plot_data), "species")))


# 6. PLOT TRAIT RELATIONSHIPS --------------------------------------------------

# Add biome classification
traits_biomes <- traits_median_df |>
  mutate(biome_category = case_when(species_level_mean_distance_km > 0 ~ "boreal",
                                    species_level_mean_distance_km < 0 ~ "tundra",
                                    TRUE ~ "boundary"))

# Plot sla & plant height
sla_height <- traits_biomes |>
  select(c("StandardSpeciesName", "TraitNameNew", "biome_category", "MedianTraitValue")) |>
  filter(TraitNameNew %in% c("PlantHeight", "SLA")) |>
  filter(!StandardSpeciesName %in% c("Picea sitchensis", "Ilex aquifolium"))

# Convert to wide format
sla_height_wide <- sla_height |>
  pivot_wider(names_from = TraitNameNew, 
              values_from = MedianTraitValue) |>
  column_to_rownames("StandardSpeciesName") |>
  filter(!is.na(PlantHeight) & !is.na(SLA)) |>
  mutate(sqrt_plant_height = sqrt(PlantHeight),
         sqrt_sla = sqrt(SLA))

# Plot 
p1<- ggplot(sla_height_wide, aes(x=sqrt_plant_height, y=sqrt_sla, 
                                 color=biome_category)) +
  geom_point() +
  theme_classic() +
  theme(legend.position="right")

ggMarginal(p1, groupColour = TRUE, groupFill = TRUE)




