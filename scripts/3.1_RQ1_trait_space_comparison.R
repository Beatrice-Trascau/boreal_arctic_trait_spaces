##----------------------------------------------------------------------------##
# PAPER 3: BOREAL AND ARCTIC PLANT SPECIES TRAIT SPACES 
# 3.1_RQ1_trait_space_comparison
# This script addresses RQ1: Are boreal species functionally distinct from 
# Arctic species?
##----------------------------------------------------------------------------##

# 1. SETUP ---------------------------------------------------------------------

# Load packages
library(here)
source(here("scripts", "0_setup.R"))

# Load data
load(here("data", "raw_data", "try_beatrice.RData"))
try_raw <- try.final.control

# Load distances to biome boundaries
detailed_results <- readRDS(here("data", "derived_data", 
                                 "species_summaries_dist_to_biome_boundary_June25.rds"))

# Load CAFF quality check
caff_check <- read.xlsx(here("data", "derived_data", "caff_quality_check.xlsx"),
                        sheet = 1, skipEmptyRows = TRUE)

# 2. DATA QUALITY CHECK --------------------------------------------------------

## 2.1. Check trait distribution and skewness ----------------------------------

# Get summary statistics
trait_summary_supp <- try_raw |>
  filter(!is.na(StdValue)) |>
  group_by(TraitNameNew) |>
  summarise(n_species = n_distinct(AccSpeciesName),
            n_records = n(),
            min = min(StdValue, na.rm = TRUE),
            Q1 = quantile(StdValue, 0.25, na.rm = TRUE),
            median = median(StdValue, na.rm = TRUE),
            mean = mean(StdValue, na.rm = TRUE),
            Q3 = quantile(StdValue, 0.75, na.rm = TRUE),
            max = max(StdValue, na.rm = TRUE),
            sd = sd(StdValue, na.rm = TRUE),
            cv = sd / mean,
            skewness = skewness(StdValue, na.rm = TRUE),
            range_ratio = max / min,
            .groups = "drop") |>
  arrange(TraitNameNew)

# Print summary statistics
print(trait_summary_supp)

# Calculate skewness
trait_skewness <- try_raw |>
  filter(!is.na(StdValue)) |>
  group_by(TraitNameNew) |>
  summarise(skewness = skewness(StdValue, na.rm = TRUE),
            .groups = "drop") |>
  arrange(desc(abs(skewness)))

# Print skewness summary
print(trait_skewness)

## 2.2. Visual inspection of distributions -------------------------------------

# Create histograms - RAW values
p_raw <- try_raw |>
  filter(!is.na(StdValue)) |>
  ggplot(aes(x = StdValue)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  facet_wrap(~ TraitNameNew, scales = "free") +
  theme_bw() +
  labs(x = "Trait value",
       y = "Count")

# Save plot
ggsave(here("figures", "FigureS3_distribution_raw_traits.png"), 
       plot = p_raw, width = 12, height = 8, dpi = 300)

# Create histograms - LOG-TRANSFORMED values
p_log <- try_raw |>
  filter(!is.na(StdValue), StdValue > 0) |>
  ggplot(aes(x = log(StdValue))) +
  geom_histogram(bins = 50, fill = "darkgreen", alpha = 0.7) +
  facet_wrap(~ TraitNameNew, scales = "free") +
  theme_bw() +
  labs( x = "log(Trait value)",
       y = "Count")

# Save plot
ggsave(here("figures", "FigureS4_distribution_log_traits.png"), 
       plot = p_log, width = 12, height = 8, dpi = 300)

## 2.3. Check species record counts --------------------------------------------

# Count records per species per trait
species_trait_counts <- try_raw |>
  filter(!is.na(StdValue)) |>
  group_by(AccSpeciesName, TraitNameNew) |>
  summarise(n_records = n(), .groups = "drop")

# Key traits for NMDS
key_traits <- c("PlantHeight", "SLA", "LeafN", "SeedMass")

# Check how many species have >=3 records for all 4 traits
species_complete_data_threshold4 <- species_trait_counts |>
  filter(TraitNameNew %in% key_traits) |>
  filter(n_records >= 3) |> 
  group_by(AccSpeciesName) |>
  summarise(n_traits_with_data = n(), .groups = "drop") |>
  filter(n_traits_with_data == 4)

# Check how many species have complete data
nrow(species_complete_data_threshold4)

# Also check breakdown by trait
trait_availability <- species_trait_counts |>
  filter(TraitNameNew %in% key_traits) |>
  filter(n_records >= 3) |>
  group_by(TraitNameNew) |>
  summarise(n_species_available = n(), .groups = "drop") |>
  arrange(desc(n_species_available))

# Check availability broken down by trait
print(trait_availability)


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

# 4. ADD DISTANCE TO BIOME BOUNDARY --------------------------------------------

# Remove Elodea canadensis & hybrids from biome boundaries
biome_boundaries <- detailed_results |>
  filter(!species == "Elodea canadensis") |>
  filter(!str_detect(species, " × "))

# Combine with trait data
global_cleaned_traits1 <- global_traits5 |>
  left_join(biome_boundaries, by = c("StandardSpeciesName" = "species"))

# Rename distance column
global_cleaned_traits2 <- global_cleaned_traits1 |>
  rename(species_level_mean_distance_km = mean_distance_km)

# Remove species without distance data
global_cleaned_traits3 <- global_cleaned_traits2 |>
  filter(!is.na(species_level_mean_distance_km))

# Check the number of unique species left after cleaning
length(unique(global_cleaned_traits3$StandardSpeciesName))

# 5. CLASSIFICATION QUALITY CONTROL --------------------------------------------

# Prepare CAFF classification
caff_check <- caff_check |>
  rename(StandardSpeciesName = SPECIES_CLEAN)

# Add CAFF classification
caff_global_cleaned_traits <- global_cleaned_traits3 |>
  left_join(caff_check |> select(StandardSpeciesName, final.category), 
            by = "StandardSpeciesName")

# Filter out species marked for removal
caff_global_cleaned_traits2 <- caff_global_cleaned_traits |>
  mutate(caff_biome_category = final.category) |>
  filter(!caff_biome_category == "remove")

# Check how many unique species there are in the newly cleaned df
length(unique(caff_global_cleaned_traits2$StandardSpeciesName))

# Check classification breakdown
classification_summary <- caff_global_cleaned_traits2 |>
  distinct(StandardSpeciesName, caff_biome_category) |>
  count(caff_biome_category)

# Look at the classification summary
print(classification_summary)

# 6. CALCULATE SPECIES-LEVEL TRAIT MEDIANS ------------------------------------

# Calculate medians for species with >=3 records per trait
traits_median <- caff_global_cleaned_traits2 |>
  group_by(StandardSpeciesName, TraitNameNew) |>
  mutate(number_of_records = length(StdValue)) |>
  filter(number_of_records > 2) |>  # >=3 records threshold
  summarise(max_trait_value = max(StdValue, na.rm = TRUE),
            MedianTraitValue = median(StdValue, na.rm = TRUE),
            n_records = n(),
            .groups = "drop")

# Check how many species-trait combinations there are left after filtering
nrow(traits_median)

# Add distance to biome boundaries
traits_median_df <- traits_median |>
  left_join(biome_boundaries, by = c("StandardSpeciesName" = "species")) |>
  rename(species_level_mean_distance_km = mean_distance_km) |>
  filter(!is.na(species_level_mean_distance_km), !is.na(MedianTraitValue)) |>
  mutate(log_median_trait_value = log(MedianTraitValue))

# Check how many species-trait combinations with complete data there are
nrow(traits_median_df)

# Categorise species by biome (based on distance)
species_biome_classification <- traits_median_df |>
  distinct(StandardSpeciesName, species_level_mean_distance_km) |>
  mutate(pipeline_biome_category = case_when(species_level_mean_distance_km > 0 ~ "boreal",
                                             species_level_mean_distance_km < 0 ~ "tundra",
                                             TRUE ~ "boundary"))

# 7. NMDS WITH 4 TRAITS  -------------------------------------------------------

## 7.1. Create trait matrix (PlantHeight, SLA, LeafN, SeedMass) ----------------

# Create wide format trait matrix
trait_matrix_4trait <- traits_median_df |>
  filter(TraitNameNew %in% key_traits) |>
  select(StandardSpeciesName, TraitNameNew, log_median_trait_value) |>
  pivot_wider(names_from = TraitNameNew, 
              values_from = log_median_trait_value) |>
  column_to_rownames("StandardSpeciesName")

# Remove species with missing data for any trait
complete_trait_matrix_4trait <- trait_matrix_4trait[complete.cases(trait_matrix_4trait), ]

# Check how many species will be included in the NMDS
nrow(complete_trait_matrix_4trait) #101

# Double check the traits included
paste(colnames(complete_trait_matrix_4trait), collapse = ", ")

## 7.2. Run NMDS ---------------------------------------------------------------

# Set seed
set.seed(532826)

# Run NMDS
nmds_4trait <- metaMDS(complete_trait_matrix_4trait, 
                       distance = "euclidean",
                       k = 3,
                       trymax = 100)

# Get the stress value
round(nmds_4trait$stress, 3) # 0.019

# Check stress across dimensions
dimcheck_4trait <- dimcheckMDS(complete_trait_matrix_4trait,
                               distance = "euclidean",
                               k = 4) # not the happiest

## 7.3. Extract NMDS scores and add classification ----------------------------

# Get NMDS scores
nmds_scores_4trait <- as.data.frame(nmds_4trait$points)
nmds_scores_4trait$StandardSpeciesName <- rownames(nmds_scores_4trait)

# Get CAFF biome classification
caff_biomes <- caff_global_cleaned_traits2 |>
  select(StandardSpeciesName, caff_biome_category) |>
  distinct(StandardSpeciesName, .keep_all = TRUE)

# Add classification to NMDS scores
nmds_plot_data_4trait <- nmds_scores_4trait |>
  left_join(caff_biomes, by = "StandardSpeciesName") |>
  filter(!is.na(caff_biome_category))

# Check how many species are in the final NMDS plot
nrow(nmds_plot_data_4trait) #101 (correct)


## 8.3. Fit trait vectors ------------------------------------------------------

# Fit trait vectors to ordination
trait_fit_4trait <- envfit(nmds_4trait, complete_trait_matrix_4trait, 
                           permutations = 999, na.rm = TRUE)

# Check values
print(trait_fit_4trait)

# Extract vectors
trait_vectors_4trait <- as.data.frame(scores(trait_fit_4trait, "vectors"))
trait_vectors_4trait$trait <- rownames(trait_vectors_4trait)

## 8.4. Statistical tests ------------------------------------------------------

# PERMANOVA with all species
permanova_4trait_full <- adonis2(complete_trait_matrix_4trait ~ caff_biome_category, 
                                 data = nmds_plot_data_4trait,
                                 method = "euclidean")

# Check PERMANOVA results
print(permanova_4trait_full)

# PERMANOVA without Silene acaulis (sensitivity check)
if("Silene acaulis" %in% rownames(complete_trait_matrix_4trait)) {
  cat("\nSensitivity analysis: Removing Silene acaulis...\n")
  
  trait_matrix_no_silene <- complete_trait_matrix_4trait[
    !rownames(complete_trait_matrix_4trait) %in% c("Silene acaulis"), ]
  
  permanova_4trait_no_silene <- adonis2(trait_matrix_no_silene ~ caff_biome_category, 
                                        data = nmds_plot_data_4trait[
                                          nmds_plot_data_4trait$StandardSpeciesName != "Silene acaulis", ],
                                        method = "euclidean")
  
  cat("\nPERMANOVA results (without Silene acaulis):\n")
  print(permanova_4trait_no_silene)
} else {
  cat("\nNote: Silene acaulis not in dataset\n")
}

# Test for homogeneity of dispersions
dispersion_4trait <- betadisper(vegdist(complete_trait_matrix_4trait, method = "euclidean"), 
                                nmds_plot_data_4trait$caff_biome_category)
dispersion_test_4trait <- permutest(dispersion_4trait)

print(dispersion_test_4trait)

## 8.5. Create NMDS plot -------------------------------------------------------

# Get convex hulls
boreal_scores_4trait <- nmds_plot_data_4trait[nmds_plot_data_4trait$caff_biome_category == "boreal", ][
  chull(nmds_plot_data_4trait[nmds_plot_data_4trait$caff_biome_category == "boreal", c("MDS1", "MDS2")]), ]

tundra_scores_4trait <- nmds_plot_data_4trait[nmds_plot_data_4trait$caff_biome_category == "tundra", ][
  chull(nmds_plot_data_4trait[nmds_plot_data_4trait$caff_biome_category == "tundra", c("MDS1", "MDS2")]), ]

hull_data_4trait <- rbind(boreal_scores_4trait, tundra_scores_4trait)

# Create base plot
nmds_plot_4trait <- ggplot(nmds_plot_data_4trait, 
                           aes(x = MDS1, y = MDS2, color = caff_biome_category)) +
  geom_polygon(data = hull_data_4trait,
               aes(x = MDS1, y = MDS2, fill = caff_biome_category, group = caff_biome_category),
               alpha = 0.30) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = c("boreal" = "darkgreen", "tundra" = "black"),
                     name = "Biome") +
  scale_fill_manual(values = c("boreal" = "darkgreen", "tundra" = "black"),
                    name = "Biome") +
  labs(x = "NMDS1", 
       y = "NMDS2",
       subtitle = paste0("Based on ", ncol(complete_trait_matrix_4trait), " traits, ", 
                         nrow(nmds_plot_data_4trait), " species | Stress = ", 
                         round(nmds_4trait$stress, 3))) +
  theme_classic() +
  theme(legend.position = "bottom")

# Add trait vectors
nmds_plot_4trait_vectors <- nmds_plot_4trait +
  geom_segment(data = trait_vectors_4trait, 
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "black", 
               size = 1,
               inherit.aes = FALSE) +
  geom_text(data = trait_vectors_4trait,
            aes(x = NMDS1 * 1.1, y = NMDS2 * 1.1, label = trait),
            color = "black", 
            size = 3, 
            fontface = "bold",
            inherit.aes = FALSE)

# Check the plot
print(nmds_plot_4trait_vectors)

# Save plot
ggsave(here("figures", "Figure2_RQ1_NMDS_4traits.png"), 
       plot = nmds_plot_4trait_vectors, width = 10, height = 8, dpi = 300)

# 10. SAVE RESULTS -------------------------------------------------------------

# Save NMDS objects
save(nmds_4trait, nmds_plot_data_4trait,
     trait_vectors_4trait,
     permanova_4trait_full,
     file = here("data", "derived_data", "RQ1_NMDS_results.RData"))

# Save species lists
species_4trait <- data.frame(StandardSpeciesName = rownames(complete_trait_matrix_4trait))
write.csv(species_4trait, 
          here("data", "derived_data", "RQ1_species_4trait_NMDS.csv"),
          row.names = FALSE)

# END OF SCRIPT ----------------------------------------------------------------