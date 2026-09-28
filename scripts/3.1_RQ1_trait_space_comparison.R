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

# Load the corrected distances to biome boundaries (corrected version from script 2.3)
detailed_results <- readRDS(here("data", "derived_data", 
                                 "species_summaries_dist_to_biome_boundary_corrected_Sept2026.rds"))

# Load CAFF quality check
caff_check <- read.xlsx(here("data", "derived_data", "caff_quality_check.xlsx"),
                        sheet = 1, skipEmptyRows = TRUE)

# Set biome colours and legend labels
biome_colours <- c("boreal" = "darkgreen", "tundra" = "darkblue")
biome_labels  <- c("boreal" = "Boreal", "tundra" = "Tundra")

# Set text sizes in point 
fig_axis_text    <- 16 
fig_axis_title   <- 18
fig_legend_text  <- 16
fig_legend_title <- 18
fig_annot_text   <- 16
fig_panel_label  <- 22

# Set a theme to be used by all trait-space figures
theme_traitspace <- function(legend_position = "bottom") {
  theme_classic() +
    theme(axis.text = element_text(size = fig_axis_text),
          axis.title = element_text(size = fig_axis_title, face = "bold"),
          legend.title = element_text(size = fig_legend_title, face = "bold"),
          legend.text = element_text(size = fig_legend_text),
          legend.key.size = unit(1, "cm"),
          plot.subtitle = element_text(size = fig_annot_text),
          legend.position = legend_position)
}

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
# TraitNameNew n_species n_records    min    Q1 median   mean    Q3    max    sd    cv skewness
# <chr>            <int>     <int>  <dbl> <dbl>  <dbl>  <dbl> <dbl>  <dbl> <dbl> <dbl>    <dbl>
#   1 LeafCN             374      3820 1      15.8  20.3   23.6   27.7    83.4 11.7  0.496    1.73 
# 2 LeafN              875     15856 0.2    14.3  20.1   21.2   26.5    69.3  9.24 0.435    0.771
# 3 PlantHeight       1004     55949 0.0001  0.05  0.12   0.450  0.28   55    1.87 4.16    10.8  
# 4 SLA               1347     51791 0.118  10.3  15.2   17.1   21.6  2000   16.5  0.963   58.7  
# 5 SeedMass           381      2797 0.001   0.24  0.857  2.45   2.35 1076.  24.3  9.92    38.3 

# Calculate skewness
trait_skewness <- try_raw |>
  filter(!is.na(StdValue)) |>
  group_by(TraitNameNew) |>
  summarise(skewness = skewness(StdValue, na.rm = TRUE),
            .groups = "drop") |>
  arrange(desc(abs(skewness)))

# Print skewness summary
print(trait_skewness)
# TraitNameNew skewness
# <chr>           <dbl>
#   1 SLA            58.7  
# 2 SeedMass       38.3  
# 3 PlantHeight    10.8  
# 4 LeafCN          1.73 
# 5 LeafN           0.771

## 2.2. Visual inspection of distributions -------------------------------------

# Create histograms - RAW values
(p_raw <- try_raw |>
  filter(!is.na(StdValue)) |>
  ggplot(aes(x = StdValue)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  facet_wrap(~ TraitNameNew, scales = "free") +
  theme_bw() +
  labs(x = "Trait value",
       y = "Count"))

# Save plot
ggsave(here("figures", "FigureS3_distribution_raw_traits.png"), 
       plot = p_raw, width = 12, height = 8, dpi = 300)

# Create histograms - LOG-TRANSFORMED values
(p_log <- try_raw |>
    filter(!is.na(StdValue), StdValue > 0) |>
    ggplot(aes(x = log(StdValue))) +
    geom_histogram(bins = 50, fill = "darkgreen", alpha = 0.7) +
    facet_wrap(~ TraitNameNew, scales = "free") +
    theme_bw() +
    labs( x = "log(Trait value)",
          y = "Count"))

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
nrow(species_complete_data_threshold4) #107

# Also check breakdown by trait
trait_availability <- species_trait_counts |>
  filter(TraitNameNew %in% key_traits) |>
  filter(n_records >= 3) |>
  group_by(TraitNameNew) |>
  summarise(n_species_available = n(), .groups = "drop") |>
  arrange(desc(n_species_available))

# Check availability broken down by trait
print(trait_availability)
# TraitNameNew n_species_available
# <chr>                      <int>
# 1 SLA                          992
# 2 PlantHeight                  571
# 3 LeafN                        528
# 4 SeedMass                     176

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
length(unique(try_raw$AccSpeciesName)) #1799
length(unique(global_traits1$AccSpeciesName)) #1754

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
length(unique(global_traits5$StandardSpeciesName)) #1706

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
length(unique(global_cleaned_traits3$StandardSpeciesName)) #1386

# 5. CLASSIFICATION QUALITY CONTROL --------------------------------------------

# Prepare CAFF classification
caff_check <- caff_check |>
  rename(StandardSpeciesName = SPECIES_CLEAN)

# Add CAFF classification
caff_global_cleaned_traits <- global_cleaned_traits3 |>
  left_join(caff_check |> dplyr::select(StandardSpeciesName, final.category), 
            by = "StandardSpeciesName")

# Filter out species marked for removal
caff_global_cleaned_traits2 <- caff_global_cleaned_traits |>
  mutate(caff_biome_category = final.category) |>
  filter(!caff_biome_category == "remove")

# Check how many unique species there are in the newly cleaned df
length(unique(caff_global_cleaned_traits2$StandardSpeciesName)) #1345

# Check classification breakdown
classification_summary <- caff_global_cleaned_traits2 |>
  distinct(StandardSpeciesName, caff_biome_category) |>
  count(caff_biome_category)

# Look at the classification summary
print(classification_summary)
# caff_biome_category     n
# <chr>               <int>
# 1 boreal               1214
# 2 tundra                131

# One row per species with its CAFF classification (used throughout)
caff_biomes <- caff_global_cleaned_traits2 |>
  dplyr::select(StandardSpeciesName, caff_biome_category) |>
  distinct(StandardSpeciesName, .keep_all = TRUE)

# 6. CALCULATE SPECIES-LEVEL TRAIT MEDIANS -------------------------------------

## 6.1. Species-level medians --------------------------------------------------

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
nrow(traits_median) # 2150

# Add distance to biome boundaries
traits_median_df <- traits_median |>
  left_join(biome_boundaries, by = c("StandardSpeciesName" = "species")) |>
  rename(species_level_mean_distance_km = mean_distance_km) |>
  filter(!is.na(species_level_mean_distance_km), !is.na(MedianTraitValue)) |>
  mutate(log_median_trait_value = log(MedianTraitValue))

# Check how many species-trait combinations with complete data there are
nrow(traits_median_df) #2150

## 6.2. Compare distance-based classification with CAFF ------------------------

# Categorise species by biome (based on distance)
species_biome_classification <- traits_median_df |>
  distinct(StandardSpeciesName, species_level_mean_distance_km) |>
  mutate(pipeline_biome_category = case_when(species_level_mean_distance_km > 0 ~ "boreal",
                                             species_level_mean_distance_km < 0 ~ "tundra",
                                             TRUE ~ "boundary"))

# Compare with the CAFF classification used in the analyses
classification_agreement <- species_biome_classification |>
  left_join(caff_biomes, by = "StandardSpeciesName")

# Cross-table: rows = distance-based, columns = CAFF
print(table(pipeline = classification_agreement$pipeline_biome_category,
            caff = classification_agreement$caff_biome_category))
#          caff
# pipeline boreal tundra
# boreal    850      0
# tundra     21     92

# Check what % of records are classified in agreement with CAFF
cat("Agreement between distance-based and CAFF classification:",
    round(100 * mean(classification_agreement$pipeline_biome_category ==
                       classification_agreement$caff_biome_category, na.rm = TRUE), 1), "%\n") # 97.8 %

# Save species where the two classifications disagree
write.csv(classification_agreement |>
            filter(pipeline_biome_category != caff_biome_category),
          here("data", "derived_data", "RQ1_classification_disagreements_pipeline_vs_CAFF.csv"),
          row.names = FALSE)

# 7. NMDS WITH 4 TRAITS  -------------------------------------------------------

## 7.1. Create trait matrix (PlantHeight, SLA, LeafN, SeedMass) ----------------

#seedless_traits <- c("PlantHeight", "SLA", "LeafN", "SeedMass")

# Create wide format trait matrix
trait_matrix_4trait <- traits_median_df |>
  filter(TraitNameNew %in% key_traits) |>
  dplyr::select(StandardSpeciesName, TraitNameNew, log_median_trait_value) |>
  pivot_wider(names_from = TraitNameNew, 
              values_from = log_median_trait_value) |>
  column_to_rownames("StandardSpeciesName")

# Remove species with missing data for any trait
complete_trait_matrix_4trait <- trait_matrix_4trait[complete.cases(trait_matrix_4trait), ]

# Check how many species will be included in the NMDS
nrow(complete_trait_matrix_4trait) #101

# Double check the traits included
paste(colnames(complete_trait_matrix_4trait), collapse = ", ")

# Standardise log-traits (mean 0, SD 1) across the species in the analysis,
# so that each trait contributes equally to the Euclidean distances
trait_scaling <- data.frame(Trait = colnames(complete_trait_matrix_4trait),
                            mean_log = colMeans(complete_trait_matrix_4trait),
                            sd_log = apply(complete_trait_matrix_4trait, 2, sd),
                            row.names = NULL)

# Check the spread of each trait before standardising
print(trait_scaling)
#         Trait   mean_log    sd_log
# 1       LeafN  3.0787352 0.3404609
# 2         SLA  2.7663334 0.3990198
# 3 PlantHeight -1.9799131 1.0788766
# 4    SeedMass -0.7691451 1.4297480

# Convert the trait matrix to data frame
complete_trait_matrix_4trait_std <- as.data.frame(scale(complete_trait_matrix_4trait))

# Check: every trait should now have mean 0 and SD 1
round(colMeans(complete_trait_matrix_4trait_std), 10)
apply(complete_trait_matrix_4trait_std, 2, sd) # All good!

## 7.2. Run NMDS ---------------------------------------------------------------

# Set seed
set.seed(532826)

# Run NMDS
nmds_4trait <- metaMDS(complete_trait_matrix_4trait_std, 
                       distance = "euclidean",
                       k = 2,
                       trymax = 100)

# Get the stress value
round(nmds_4trait$stress, 3) # 0.147

# Check stress across dimensions
dimcheck_4trait <- dimcheckMDS(complete_trait_matrix_4trait_std,
                               distance = "euclidean",
                               k = 4) # 1 = 0.29, 2 = 0.15, 3 = 0.07, 4 = 0

## 7.3. Extract NMDS scores and add classification ----------------------------

# Get NMDS scores
nmds_scores_4trait <- as.data.frame(nmds_4trait$points)
nmds_scores_4trait$StandardSpeciesName <- rownames(nmds_scores_4trait)

# Add classification to NMDS scores
nmds_plot_data_4trait <- nmds_scores_4trait |>
  left_join(caff_biomes, by = "StandardSpeciesName") |>
  filter(!is.na(caff_biome_category))

# Check how many species are in the final NMDS plot
nrow(nmds_plot_data_4trait) #101 (correct)

# Check that species order matches the trait matrix (needed for PERMANOVA/betadisper)
if (!identical(rownames(complete_trait_matrix_4trait_std), 
               nmds_plot_data_4trait$StandardSpeciesName)) {
  stop("ERROR: Species order mismatch between trait matrix and NMDS classification!")
} # all good!

## 7.4. Fit trait vectors ------------------------------------------------------

# Set seed
set.seed(532826)

# Fit trait vectors to ordination
trait_fit_4trait <- envfit(nmds_4trait, complete_trait_matrix_4trait_std, 
                           permutations = 999, na.rm = TRUE)

# Check values
print(trait_fit_4trait)
# NMDS1    NMDS2     r2 Pr(>r)    
# LeafN        0.82834  0.56022 0.6906  0.001 ***
# SLA          0.50681  0.86206 0.7791  0.001 ***
# PlantHeight  0.73091 -0.68247 0.6733  0.001 ***
# SeedMass     0.86245 -0.50614 0.7135  0.001 ***

# Extract vectors
trait_vectors_4trait <- as.data.frame(scores(trait_fit_4trait, "vectors"))
trait_vectors_4trait$trait <- rownames(trait_vectors_4trait)

## 7.5. Statistical tests ------------------------------------------------------

# Set seed
set.seed(532826)

# PERMANOVA with all species
permanova_4trait_full <- adonis2(complete_trait_matrix_4trait_std ~ caff_biome_category, 
                                 data = nmds_plot_data_4trait,
                                 method = "euclidean",
                                 permutations = 999)

# Check PERMANOVA results
print(permanova_4trait_full)
# adonis2(formula = complete_trait_matrix_4trait_std ~ caff_biome_category, data = nmds_plot_data_4trait, permutations = 999, method = "euclidean")
# Df SumOfSqs      R2      F Pr(>F)   
# Model      1     17.3 0.04324 4.4741  0.007 **
# Residual  99    382.7 0.95676                 
# Total    100    400.0 1.00000 

# Test for homogeneity of dispersions
dispersion_4trait <- betadisper(vegdist(complete_trait_matrix_4trait_std, method = "euclidean"), 
                                nmds_plot_data_4trait$caff_biome_category)
set.seed(532826)
dispersion_test_4trait <- permutest(dispersion_4trait, permutations = 999)

# Check the output of the dispersion test
print(dispersion_test_4trait)
# Response: Distances
#            Df Sum Sq Mean Sq      F N.Perm Pr(>F)
# Groups     1  1.377 1.37663 1.4461    999  0.229
# Residuals 99 94.243 0.95195 

## 7.6. Create NMDS plot -------------------------------------------------------

# Get convex hulls
boreal_scores_4trait <- nmds_plot_data_4trait[nmds_plot_data_4trait$caff_biome_category == "boreal", ][
  chull(nmds_plot_data_4trait[nmds_plot_data_4trait$caff_biome_category == "boreal", c("MDS1", "MDS2")]), ]

tundra_scores_4trait <- nmds_plot_data_4trait[nmds_plot_data_4trait$caff_biome_category == "tundra", ][
  chull(nmds_plot_data_4trait[nmds_plot_data_4trait$caff_biome_category == "tundra", c("MDS1", "MDS2")]), ]

hull_data_4trait <- rbind(boreal_scores_4trait, tundra_scores_4trait)

# Create base plot
(nmds_plot_4trait <- ggplot(nmds_plot_data_4trait, 
                            aes(x = MDS1, y = MDS2, color = caff_biome_category)) +
    geom_polygon(data = hull_data_4trait,
                 aes(x = MDS1, y = MDS2, fill = caff_biome_category, group = caff_biome_category),
                 alpha = 0.30) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_manual(values = biome_colours, labels = biome_labels, name = "Biome") +
    scale_fill_manual(values = biome_colours, labels = biome_labels, name = "Biome") +
    labs(x = "NMDS1", 
         y = "NMDS2",
         subtitle = paste0("Based on ", ncol(complete_trait_matrix_4trait_std), 
                           " standardised log-traits, ", 
                           nrow(nmds_plot_data_4trait), " species | Stress = ", 
                           round(nmds_4trait$stress, 3))) +
    theme_traitspace(legend_position = "bottom"))

# Add trait vectors
(nmds_plot_4trait_vectors <- nmds_plot_4trait +
  geom_segment(data = trait_vectors_4trait, 
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "black", 
               linewidth = 1,
               inherit.aes = FALSE) +
  geom_text(data = trait_vectors_4trait,
            aes(x = NMDS1 * 1.1, y = NMDS2 * 1.1, label = trait),
            color = "black", 
            size = fig_annot_text / .pt,   # convert points to ggplot's mm units
            fontface = "bold",
            inherit.aes = FALSE))

# Check the plot
print(nmds_plot_4trait_vectors)

# Save plot
ggsave(here("figures", "Figure2_RQ1_NMDS_4traits.png"), 
       plot = nmds_plot_4trait_vectors, width = 10, height = 8, dpi = 600)

## 7.7. Save NMDS results ------------------------------------------------------

# Save NMDS objects
save(nmds_4trait, nmds_plot_data_4trait,
     trait_vectors_4trait,
     permanova_4trait_full,
     dispersion_test_4trait,
     complete_trait_matrix_4trait_std,
     trait_scaling,
     file = here("data", "derived_data", "RQ1_NMDS_results.RData"))

# Save the means and SDs used for standardising (for the methods/supplement)
write.csv(trait_scaling,
          here("data", "derived_data", "RQ1_trait_standardisation_values.csv"),
          row.names = FALSE)

# Save species lists
species_4trait <- data.frame(StandardSpeciesName = rownames(complete_trait_matrix_4trait))
write.csv(species_4trait, 
          here("data", "derived_data", "RQ1_species_4trait_NMDS.csv"),
          row.names = FALSE)

# 8. GLLVM --------------------------------------------------------------------

## 8.1. Prepare data for GLLVM -------------------------------------------------

# Use the same standardised trait matrix as the NMDS
gllvm_data <- complete_trait_matrix_4trait_std

# Get biome classification for the SAME species that are in gllvm_data
gllvm_species_names <- rownames(gllvm_data)

# Create biome data frame with exact species match
gllvm_biome_full <- caff_biomes |>
  filter(StandardSpeciesName %in% gllvm_species_names)

# Reorder to match gllvm_data exactly
gllvm_biome_full <- gllvm_biome_full[match(gllvm_species_names, gllvm_biome_full$StandardSpeciesName), ]

# Extract just the biome column
gllvm_biome <- data.frame(biome = gllvm_biome_full$caff_biome_category)

# Check that the matching is correct
if(!all(rownames(gllvm_data) == gllvm_biome_full$StandardSpeciesName)) {
  stop("ERROR: Species order mismatch between trait data and biome classification!")
} # all good!

# Check how many species are in each biome
print(table(gllvm_biome$biome))
# boreal tundra 
# 88     13 

## 8.2. Fit GLLVM with biome predictor -----------------------------------------

# Set seed (same one as for the NMDS)
set.seed(532826)

# Fit GLLVM with 1 latent variable (4 traits limits us to num.lv = 1)
gllvm_model <- gllvm(y = gllvm_data,
                     X = gllvm_biome,
                     family = "gaussian",
                     num.lv = 1,
                     formula = ~ biome,
                     seed = 532826)

# Store and check model summary
gllvm_summary <- summary(gllvm_model)
print(gllvm_summary)
# Coefficients predictors:
# Estimate Std. Error z value Pr(>|z|)   
# biometundra:LeafN        0.08236    0.29554   0.279   0.7805   
# biometundra:SLA         -0.49535    0.29152  -1.699   0.0893 . 
# biometundra:PlantHeight -0.89484    0.28193  -3.174   0.0015 **
# biometundra:SeedMass    -0.68852    0.28761  -2.394   0.0167 * 

## 9.3. Fit null model and test biome effect -----------------------------------

# Set the same seed
set.seed(532826)

# Fit null model (without biome predictor)
gllvm_null <- gllvm(y = gllvm_data,
                    family = "gaussian",
                    num.lv = 1,
                    seed = 532826)

# Likelihood ratio test
gllvm_anova <- anova(gllvm_null, gllvm_model)
print(gllvm_anova)
#   Resid.Df        D Df.diff     P.value
# 1      392  0.00000       0            
# 2      388 20.16205       4 0.000463924

# Extract p-value
biome_p_value <- gllvm_anova$`Pr(>Chi)`[2]

# Manual LRT calculation (in case anova() gives wrong p-value)
ll_null <- as.numeric(logLik(gllvm_null))
ll_full <- as.numeric(logLik(gllvm_model))
lr_stat <- -2 * (ll_null - ll_full)
df_diff <- attr(logLik(gllvm_model), "df") - attr(logLik(gllvm_null), "df")
p_value_manual <- pchisq(lr_stat, df = df_diff, lower.tail = FALSE)

# Calculate pseudo R-squared
pseudo_r2 <- 1 - (ll_full / ll_null) # 0.01824379

## 9.4. Extract trait-specific biome effects -----------------------------------

# Get coefficients with standard errors from summary
trait_effects <- as.data.frame(gllvm_summary$Coef.tableX)

# Clean up row names to get trait names
trait_effects$Trait <- gsub("biometundra:", "", rownames(trait_effects))

# Coefficients are in SD units of the log-trait (tundra relative to boreal).
# Back-transform: SD units -> log units -> multiplicative ratio (tundra / boreal)
trait_effects <- trait_effects |>
  dplyr::select(Trait, Coefficient = Estimate, SE = `Std. Error`,
                Z_value = `z value`, P_value = `Pr(>|z|)`) |>
  left_join(trait_scaling |> dplyr::select(Trait, sd_log), by = "Trait") |>
  mutate(Lower95 = Coefficient - 1.96 * SE,
         Upper95 = Coefficient + 1.96 * SE,
         Coefficient_log = Coefficient * sd_log,
         Ratio_tundra_boreal = exp(Coefficient_log),
         Ratio_lower95 = exp(Lower95 * sd_log),
         Ratio_upper95 = exp(Upper95 * sd_log),
         Significance = case_when(P_value < 0.001 ~ "***",
                                  P_value < 0.01 ~ "**",
                                  P_value < 0.05 ~ "*",
                                  P_value < 0.1 ~ ".",
                                  TRUE ~ "ns"),
         Sig_binary = ifelse(P_value < 0.05, "Significant (p<0.05)", "Not significant"))

# Check that all traits were matched to their SD
if (any(is.na(trait_effects$sd_log))) {
  stop("ERROR: Could not match GLLVM trait names to the standardisation values!")
} # all ok!

# Display
print(trait_effects, row.names = FALSE)

# Print the important information for each trait 
for(i in 1:nrow(trait_effects)) {
  cat("  ", trait_effects$Trait[i], ":\n")
  cat("      β =", sprintf("%.3f", trait_effects$Coefficient[i]), 
      "±", sprintf("%.3f", trait_effects$SE[i]), "SD",
      trait_effects$Significance[i], "\n")
  cat("      Z =", sprintf("%.2f", trait_effects$Z_value[i]), 
      ", p =", format.pval(trait_effects$P_value[i], digits = 4), "\n")
  cat("      Tundra/boreal ratio =", sprintf("%.2f", trait_effects$Ratio_tundra_boreal[i]),
      "(95% CI", sprintf("%.2f", trait_effects$Ratio_lower95[i]), "-",
      sprintf("%.2f", trait_effects$Ratio_upper95[i]), ")\n\n")
}
# LeafN :
#   β = 0.082 ± 0.296 SD ns 
# Z = 0.28 , p = 0.7805 
# Tundra/boreal ratio = 1.03 (95% CI 0.84 - 1.25 )
# 
# SLA :
#   β = -0.495 ± 0.292 SD . 
# Z = -1.70 , p = 0.08928 
# Tundra/boreal ratio = 0.82 (95% CI 0.65 - 1.03 )
# 
# PlantHeight :
#   β = -0.895 ± 0.282 SD ** 
#   Z = -3.17 , p = 0.001504 
# Tundra/boreal ratio = 0.38 (95% CI 0.21 - 0.69 )
# 
# SeedMass :
#   β = -0.689 ± 0.288 SD * 
#   Z = -2.39 , p = 0.01667 
# Tundra/boreal ratio = 0.37 (95% CI 0.17 - 0.84 )

# Save model output
write.csv(trait_effects,
          here("data", "derived_data", "RQ1_GLLVM_trait_effects.csv"),
          row.names = FALSE)

## 9.5. Extract latent variable loadings ---------------------------------------

# Get loadings from model parameters
lv_loadings <- gllvm_model$params$theta

# Create loadings data frame
loadings_df <- data.frame(Trait = rownames(lv_loadings),
                          Loading_LV1 = lv_loadings[, 1]) |>
  arrange(desc(abs(Loading_LV1)))

# Get loadings
print(loadings_df) # latent variable captures the variation not explained by biome
# higher absolute loading value = trait contributes more to the latent variable
# Trait Loading_LV1
# LeafN             LeafN   1.0000000
# SLA                 SLA   0.4642256
# SeedMass       SeedMass   0.3504976
# PlantHeight PlantHeight   0.1883649

## 9.6. Visualize trait-specific effects --------------------------------------

# Create forest plot
gllvm_forest_plot <- ggplot(trait_effects, 
                            aes(x = reorder(Trait, Coefficient), 
                                y = Coefficient,
                                color = Sig_binary)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_errorbar(aes(ymin = Lower95, ymax = Upper95), 
                width = 0.2, linewidth = 1) +
  geom_point(size = 4) +
  scale_color_manual(values = c("Significant (p<0.05)" = "steelblue", 
                                "Not significant" = "gray60"),
                     name = "") +
  coord_flip() +
  labs(x = "",
       y = "Coefficient (tundra relative to boreal, SD units)") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        legend.position = "bottom",
        legend.text = element_text(size = 14))

# Check the plots
print(gllvm_forest_plot)

# Save plot
ggsave(here("figures", "RQ1_GLLVM_trait_effects.png"),
       plot = gllvm_forest_plot, width = 8, height = 6, dpi = 300)

## 9.7. Model diagnostics ------------------------------------------------------

# Check convergence
gllvm_model$convergence # TRUE

# Create diagnostic plots
png(here("figures", "RQ1_GLLVM_diagnostics.png"), 
    width = 10, height = 10, units = "in", res = 300)
par(mfrow = c(2, 2))
plot(gllvm_model, which = 1:4)
par(mfrow = c(1, 1))
dev.off()

# 10. PAIRWISE TRAIT COMPARISON ------------------------------------------------

# N.B: uses log values (not standardised) so that we can interpret the axes
# Each biome is shpwn as a convex hull to match the NMDS wiht points for all species and a large diamond at the biome mean of the two traits

## 10.1. Prepare data for pairwise plots ---------------------------------------

# Use the trait matrix BEFORE removing incomplete cases
pairwise_df_all <- as.data.frame(trait_matrix_4trait) |>
  mutate(Species = rownames(trait_matrix_4trait)) |>
  left_join(caff_biomes |> dplyr::select(StandardSpeciesName, caff_biome_category),
            by = c("Species" = "StandardSpeciesName")) |>
  rename(Biome = caff_biome_category) |>
  filter(!is.na(Biome))

# Check how many species there are in total
nrow(pairwise_df_all) #963

# Check how they break down by biome
print(table(pairwise_df_all$Biome))
# boreal tundra 
# 871     92 

# Define list of trait pair
trait_pairs_list <- list(c("PlantHeight", "SLA"),
                         c("PlantHeight", "LeafN"),
                         c("PlantHeight", "SeedMass"),
                         c("SLA", "LeafN"),
                         c("SLA", "SeedMass"),
                         c("LeafN", "SeedMass"))

# Count the number of species available for each trait pair
for(pair in trait_pairs_list) {
  trait1 <- pair[1]
  trait2 <- pair[2]
  
  # get total number of species
  n_species <- sum(!is.na(pairwise_df_all[[trait1]]) & !is.na(pairwise_df_all[[trait2]]))
  
  # get number of boreal species
  n_boreal <- sum(!is.na(pairwise_df_all[[trait1]]) & 
                    !is.na(pairwise_df_all[[trait2]]) & 
                    pairwise_df_all$Biome == "boreal")
  
  # get number of tundra species
  n_tundra <- sum(!is.na(pairwise_df_all[[trait1]]) & 
                    !is.na(pairwise_df_all[[trait2]]) & 
                    pairwise_df_all$Biome == "tundra")
  
  # display number of species per each trait pair broken down by biome
  cat("  ", trait1, "×", trait2, ": n =", n_species, 
      "(", n_boreal, "boreal,", n_tundra, "tundra )\n")
}
# PlantHeight × SLA : n = 411 ( 354 boreal, 57 tundra )
# PlantHeight × LeafN : n = 301 ( 276 boreal, 25 tundra )
# PlantHeight × SeedMass : n = 141 ( 109 boreal, 32 tundra )
# SLA × LeafN : n = 396 ( 365 boreal, 31 tundra )
# SLA × SeedMass : n = 155 ( 124 boreal, 31 tundra )
# LeafN × SeedMass : n = 110 ( 97 boreal, 13 tundra )

## 10.2. Function to build pairwise plot ---------------------------------------

# Define the function
plot_trait_pair <- function(df, xvar, yvar, xlab, ylab,
                            show_points = TRUE, show_legend = FALSE,
                            marginal_type = "density", marginal_alpha = 0.4) {
  
  # Species with both traits (columns renamed to x and y so ggMarginal works)
  plot_data <- df |>
    filter(!is.na(.data[[xvar]]), !is.na(.data[[yvar]])) |>
    transmute(Species, Biome, x = .data[[xvar]], y = .data[[yvar]])
  
  # Convex hull per biome
  hull_data <- plot_data |>
    group_by(Biome) |>
    slice(chull(x, y)) |>
    ungroup()
  
  # Biome means of the two traits
  biome_means <- plot_data |>
    group_by(Biome) |>
    summarise(x = mean(x), y = mean(y), .groups = "drop")
  
  # Sample sizes
  n_total <- nrow(plot_data)
  n_boreal <- sum(plot_data$Biome == "boreal")
  n_tundra <- sum(plot_data$Biome == "tundra")
  
  # Base plot: hulls first, then species points (the first point layer is what
  # ggMarginal uses for the marginal plots), then biome means
  p <- ggplot(plot_data, aes(x = x, y = y, colour = Biome, fill = Biome)) +
    geom_polygon(data = hull_data, alpha = 0.25, linewidth = 0.8)
  
  if (show_points) {
    p <- p + geom_point(alpha = 0.6, size = 2.5)
  } else {
    # invisible points so ggMarginal still has a point layer to work from
    p <- p + geom_point(alpha = 0, size = 2.5)
  }
  
  p <- p +
    geom_point(data = biome_means, size = 5, shape = 18) +
    annotate("text", x = -Inf, y = Inf, 
             label = paste0("n = ", n_total, " (", n_boreal, " boreal, ", n_tundra, " tundra)"), 
             hjust = -0.1, vjust = 1.5, size = fig_annot_text / .pt, fontface = "bold") +
    scale_colour_manual(values = biome_colours, labels = biome_labels, name = "Biome") +
    scale_fill_manual(values = biome_colours, labels = biome_labels, name = "Biome") +
    labs(x = xlab, y = ylab) +
    theme_traitspace(legend_position = if (show_legend) "right" else "none")
  
  # Add marginal plots: one semi-transparent density curve per biome
  # (set marginal_type = "boxplot" or "histogram" to change the style)
  ggMarginal(p, type = marginal_type, groupColour = TRUE, groupFill = TRUE,
             alpha = marginal_alpha)
}

## 10.3. Create the six pairwise plots -----------------------------------------

# Plant Height vs SLA
(plot1 <- plot_trait_pair(pairwise_df_all, "PlantHeight", "SLA",
                          "Plant Height (log)", "SLA (log)"))

# Plant Height vs Leaf N
(plot2 <- plot_trait_pair(pairwise_df_all, "PlantHeight", "LeafN",
                          "Plant Height (log)", "Leaf N (log)"))

# Plant Height vs Seed Mass
(plot3 <- plot_trait_pair(pairwise_df_all, "PlantHeight", "SeedMass",
                          "Plant Height (log)", "Seed Mass (log)"))

# SLA vs Leaf N
(plot4 <- plot_trait_pair(pairwise_df_all, "SLA", "LeafN",
                          "SLA (log)", "Leaf N (log)"))

# SLA vs Seed Mass
(plot5 <- plot_trait_pair(pairwise_df_all, "SLA", "SeedMass",
                          "SLA (log)", "Seed Mass (log)"))

# Lead N vs Seed Mass
(plot6 <- plot_trait_pair(pairwise_df_all, "LeafN", "SeedMass",
                          "Leaf N (log)", "Seed Mass (log)",
                          show_legend = TRUE))

# Combine into a single figure
(trait_bagplots <- plot_grid(plot1, plot2, plot3, plot4, plot5, plot6,
                            labels = c('a)', 'b)', 'c)', 'd)', 'e)', 'f)'),
                            nrow = 2))

# Save combined figure as png and pdf
ggsave(here("figures", "Figure3_pairwise_trait_comparisons_all_data.png"),
       plot = trait_bagplots, width = 20, height = 15, dpi = 600)
ggsave(here("figures", "Figure3_pairwise_trait_comparisons_all_data.pdf"),
       plot = trait_bagplots, width = 20, height = 15, dpi = 600)

# Save individual plots
ggsave(here("figures", "Figure3a_pairwise_PlantHeight_SLA.png"), plot = plot1, 
       width = 6, height = 5, dpi = 600)
ggsave(here("figures", "Figure3b_pairwise_PlantHeight_LeafN.png"), plot = plot2, 
       width = 6, height = 5, dpi = 600)
ggsave(here("figures", "Figure3c_pairwise_PlantHeight_SeedMass.png"), plot = plot3, 
       width = 6, height = 5, dpi = 600)
ggsave(here("figures", "Figure3d_pairwise_SLA_LeafN.png"), plot = plot4, 
       width = 6, height = 5, dpi = 600)
ggsave(here("figures", "Figure3e_pairwise_SLA_SeedMass.png"), plot = plot5, 
       width = 6, height = 5, dpi = 600)
ggsave(here("figures", "figure3f_pairwise_LeafN_SeedMass.png"), plot = plot6, 
       width = 6, height = 5, dpi = 600)

## 10.5 PERMANOVA per trait pair -----------------------------------------------

# Test whether the biome centroids differ 
# The traits are standardised within each pairs
# Also testing for homogeneity because PERMANOVA is sensitive to differences in spread

# Set the number of permutations
n_permutations_pairs <- 9999

# Define a function to run the PERMANOVAs
run_pair_permanova <- function(df, xvar, yvar) {
  
  # species with both traits
  d <- df |>
    filter(!is.na(.data[[xvar]]), !is.na(.data[[yvar]]))
  
  # standardise the two log-traits (mean 0, SD 1) across these species
  m <- scale(as.matrix(d[, c(xvar, yvar)]))
  
  # PERMANOVA: do the biome centroids differ?
  set.seed(532826)
  perm <- adonis2(m ~ Biome, data = d, method = "euclidean",
                  permutations = n_permutations_pairs)
  
  # dispersion: do the biomes differ in spread around their centroid?
  disp <- betadisper(dist(m, method = "euclidean"), d$Biome)
  set.seed(532826)
  disp_test <- permutest(disp, permutations = n_permutations_pairs)
  disp_means <- tapply(disp$distances, d$Biome, mean)
  
  # biome means on the log scale (the diamonds in the plots)
  means <- d |>
    group_by(Biome) |>
    summarise(mx = mean(.data[[xvar]]), my = mean(.data[[yvar]]), .groups = "drop")
  
  data.frame(Trait_pair = paste(xvar, "×", yvar),
             n_total = nrow(d),
             n_boreal = sum(d$Biome == "boreal"),
             n_tundra = sum(d$Biome == "tundra"),
             mean_x_boreal = means$mx[means$Biome == "boreal"],
             mean_x_tundra = means$mx[means$Biome == "tundra"],
             mean_y_boreal = means$my[means$Biome == "boreal"],
             mean_y_tundra = means$my[means$Biome == "tundra"],
             df_biome = perm$Df[1],
             df_residual = perm$Df[2],
             SS_biome = perm$SumOfSqs[1],
             pseudo_F = perm$F[1],
             R2 = perm$R2[1],
             p_permanova = perm$`Pr(>F)`[1],
             dispersion_boreal = disp_means[["boreal"]],
             dispersion_tundra = disp_means[["tundra"]],
             dispersion_F = disp_test$tab$F[1],
             p_dispersion = disp_test$tab$`Pr(>F)`[1])
}

# Run for all six trait pairs
pairwise_permanova <- bind_rows(lapply(trait_pairs_list, function(pair) {
  run_pair_permanova(pairwise_df_all, pair[1], pair[2])
}))

# Correct for running six tests (Holm)
pairwise_permanova <- pairwise_permanova |>
  mutate(p_permanova_holm = p.adjust(p_permanova, method = "holm"),
         p_dispersion_holm = p.adjust(p_dispersion, method = "holm"),
         .after = p_permanova) |>
  relocate(p_dispersion_holm, .after = p_dispersion)

# Look at the results
print(pairwise_permanova |>
        dplyr::select(Trait_pair, n_boreal, n_tundra, pseudo_F, R2,
                      p_permanova, p_permanova_holm, dispersion_F,
                      p_dispersion, p_dispersion_holm),
      digits = 3)
# Trait_pair n_boreal n_tundra pseudo_F     R2 p_permanova p_permanova_holm dispersion_F
# 1      PlantHeight × SLA      354       57    18.65 0.0436      0.0001           0.0006       10.380
# 2    PlantHeight × LeafN      276       25    11.19 0.0361      0.0003           0.0015        1.833
# 3 PlantHeight × SeedMass      109       32     5.35 0.0371      0.0066           0.0220        0.445
# 4            SLA × LeafN      365       31     4.02 0.0101      0.0263           0.0526        1.125
# 5         SLA × SeedMass      124       31     5.20 0.0329      0.0055           0.0220        0.174
# 6       LeafN × SeedMass       97       13     2.92 0.0263      0.0571           0.0571        1.588
# p_dispersion p_dispersion_holm
# 1       0.0013            0.0078
# 2       0.1722            0.8610
# 3       0.5055            1.0000
# 4       0.2851            0.8610
# 5       0.6761            1.0000
# 6       0.2031            0.8610

# Formatted version for the Supplementary table
pairwise_permanova_table <- pairwise_permanova |>
  transmute(`Trait pair` = Trait_pair,
            `n (boreal / tundra)` = paste0(n_total, " (", n_boreal, " / ", n_tundra, ")"),
            `Boreal mean (x, y)` = paste0(sprintf("%.2f", mean_x_boreal), ", ",
                                          sprintf("%.2f", mean_y_boreal)),
            `Tundra mean (x, y)` = paste0(sprintf("%.2f", mean_x_tundra), ", ",
                                          sprintf("%.2f", mean_y_tundra)),
            `df` = paste0(df_biome, ", ", df_residual),
            `Pseudo-F` = sprintf("%.2f", pseudo_F),
            `R²` = sprintf("%.3f", R2),
            `p` = format.pval(p_permanova, digits = 3, eps = 1e-4),
            `p (Holm)` = format.pval(p_permanova_holm, digits = 3, eps = 1e-4),
            `Dispersion F` = sprintf("%.2f", dispersion_F),
            `Dispersion p` = format.pval(p_dispersion, digits = 3, eps = 1e-4),
            `Dispersion p (Holm)` = format.pval(p_dispersion_holm, digits = 3, eps = 1e-4))

print(pairwise_permanova_table)

# Save: full numeric results and the formatted table
write.csv(pairwise_permanova,
          here("data", "derived_data", "RQ1_pairwise_PERMANOVA_full_results.csv"),
          row.names = FALSE)
write.csv(pairwise_permanova_table,
          here("data", "derived_data", "TableS_RQ1_pairwise_PERMANOVA.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# 11. VIOLIN PLOTS FOR INDIVIDUAL TRAITS ---------------------------------------

# N.B. uses non-standardised log values so that axes are still interpretable

## 11.1. Remove outliers (>5 SD from mean) per biome ---------------------------

# Plant Height - remove outliers per biome
plant_height_clean <- pairwise_df_all |>
  filter(!is.na(PlantHeight), !is.na(Biome)) |>
  group_by(Biome) |>
  mutate(mean_val = mean(PlantHeight, na.rm = TRUE),
         sd_val = sd(PlantHeight, na.rm = TRUE),
         lower_bound = mean_val - 5 * sd_val,
         upper_bound = mean_val + 5 * sd_val) |>
  filter(PlantHeight >= lower_bound & PlantHeight <= upper_bound) |>
  ungroup() |>
  dplyr::select(-mean_val, -sd_val, -lower_bound, -upper_bound)

# SLA - remove outliers per biome
sla_clean <- pairwise_df_all |>
  filter(!is.na(SLA), !is.na(Biome)) |>
  group_by(Biome) |>
  mutate(mean_val = mean(SLA, na.rm = TRUE),
         sd_val = sd(SLA, na.rm = TRUE),
         lower_bound = mean_val - 5 * sd_val,
         upper_bound = mean_val + 5 * sd_val) |>
  filter(SLA >= lower_bound & SLA <= upper_bound) |>
  ungroup() |>
  dplyr::select(-mean_val, -sd_val, -lower_bound, -upper_bound)

# Leaf N - remove outliers per biome
leafn_clean <- pairwise_df_all |>
  filter(!is.na(LeafN), !is.na(Biome)) |>
  group_by(Biome) |>
  mutate(mean_val = mean(LeafN, na.rm = TRUE),
         sd_val = sd(LeafN, na.rm = TRUE),
         lower_bound = mean_val - 5 * sd_val,
         upper_bound = mean_val + 5 * sd_val) |>
  filter(LeafN >= lower_bound & LeafN <= upper_bound) |>
  ungroup() |>
  dplyr::select(-mean_val, -sd_val, -lower_bound, -upper_bound)

# Seed Mass - remove outliers per biome
seedmass_clean <- pairwise_df_all |>
  filter(!is.na(SeedMass), !is.na(Biome)) |>
  group_by(Biome) |>
  mutate(mean_val = mean(SeedMass, na.rm = TRUE),
         sd_val = sd(SeedMass, na.rm = TRUE),
         lower_bound = mean_val - 5 * sd_val,
         upper_bound = mean_val + 5 * sd_val) |>
  filter(SeedMass >= lower_bound & SeedMass <= upper_bound) |>
  ungroup() |>
  dplyr::select(-mean_val, -sd_val, -lower_bound, -upper_bound)

# Check how many data points were removed for each trait
cat("\nOutliers removed (>5 SD from mean, calculated per biome):\n")
cat("  PlantHeight:", nrow(pairwise_df_all |> filter(!is.na(PlantHeight), !is.na(Biome))) - nrow(plant_height_clean), 
    "removed,", nrow(plant_height_clean), "remaining\n") # 0 removed, 495 remaining
cat("  SLA:", nrow(pairwise_df_all |> filter(!is.na(SLA), !is.na(Biome))) - nrow(sla_clean), 
    "removed,", nrow(sla_clean), "remaining\n") # 0 removed, 839 remaining
cat("  LeafN:", nrow(pairwise_df_all |> filter(!is.na(LeafN), !is.na(Biome))) - nrow(leafn_clean), 
    "removed,", nrow(leafn_clean), "remaining\n") # 2 removed, 437 remaining
cat("  SeedMass:", nrow(pairwise_df_all |> filter(!is.na(SeedMass), !is.na(Biome))) - nrow(seedmass_clean), 
    "removed,", nrow(seedmass_clean), "remaining\n") # 0 removed, 166 remaining

## 11.2. Create violin plots ---------------------------------------------------

# Plant Height
(violin_height <- ggplot(plant_height_clean, 
                        aes(x = Biome, y = PlantHeight, fill = Biome)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.2, size = 1.5) +
  scale_fill_manual(values = c("boreal" = "darkgreen", "tundra" = "darkblue"),
                    labels = c("boreal" = "Boreal", "tundra" = "Tundra"),
                    name = "Biome") +
  scale_x_discrete(labels = c("boreal" = "Boreal", "tundra" = "Tundra")) +
  labs(x = "", y = "Plant Height (log)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 14, face = "bold")))

# SLA
(violin_sla <- ggplot(sla_clean, 
                     aes(x = Biome, y = SLA, fill = Biome)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.2, size = 1.5) +
  scale_fill_manual(values = c("boreal" = "darkgreen", "tundra" = "darkblue"),
                    labels = c("boreal" = "Boreal", "tundra" = "Tundra"),
                    name = "Biome") +
  scale_x_discrete(labels = c("boreal" = "Boreal", "tundra" = "Tundra")) +
  labs(x = "", y = "SLA (log)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 14, face = "bold")))

# Leaf N
(violin_leafn <- ggplot(leafn_clean, 
                       aes(x = Biome, y = LeafN, fill = Biome)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.2, size = 1.5) +
  scale_fill_manual(values = c("boreal" = "darkgreen", "tundra" = "darkblue"),
                    labels = c("boreal" = "Boreal", "tundra" = "Tundra"),
                    name = "Biome") +
  scale_x_discrete(labels = c("boreal" = "Boreal", "tundra" = "Tundra")) +
  labs(x = "", y = "Leaf N (log)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 14, face = "bold")))

# Seed Mass (with legend)
(violin_seedmass <- ggplot(seedmass_clean, 
                          aes(x = Biome, y = SeedMass, fill = Biome)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.2, size = 1.5) +
  scale_fill_manual(values = c("boreal" = "darkgreen", "tundra" = "darkblue"),
                    labels = c("boreal" = "Boreal", "tundra" = "Tundra"),
                    name = "Biome") +
  scale_x_discrete(labels = c("boreal" = "Boreal", "tundra" = "Tundra")) +
  labs(x = "", y = "Seed Mass (log)") +
  theme_classic() +
  theme(legend.position = "right",
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 13),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 14, face = "bold")))

## 11.3. Combine violin plots --------------------------------------------------

# Combine into single figure
(trait_violins <- plot_grid(violin_height, violin_sla, violin_leafn, violin_seedmass,
                           labels = c('a)', 'b)', 'c)', 'd)'),
                           nrow = 2, ncol = 2))

# Save combined figure
ggsave(here("figures", "FigureS5_violin_plots_4traits.png"),
       plot = trait_violins, width = 15, height = 10, dpi = 600)

ggsave(here("figures", "FigureS5_violin_plots_4traits.pdf"),
       plot = trait_violins, width = 12, height = 10, dpi = 600)

# Save individual plots
ggsave(here("figures", "FigureS5a_violin_PlantHeight.png"), plot = violin_height, 
       width = 6, height = 5, dpi = 600)
ggsave(here("figures", "FigureS5b_violin_SLA.png"), plot = violin_sla, 
       width = 6, height = 5, dpi = 600)
ggsave(here("figures", "FigureS5c_violin_LeafN.png"), plot = violin_leafn, 
       width = 6, height = 5, dpi = 600)
ggsave(here("figures", "FigureS5d_violin_SeedMass.png"), plot = violin_seedmass, 
       width = 6, height = 5, dpi = 600)

## 11.4. Linear models to compare traits between biomes ------------------------

# Plant Height
plant_height_biome_lm <- lm(PlantHeight ~ Biome, data = plant_height_clean)
summary(plant_height_biome_lm)

# SLA
sla_biome_lm <- lm(SLA ~ Biome, data = sla_clean)
summary(sla_biome_lm)

# Leaf N
leafn_biome_lm <- lm(LeafN ~ Biome, data = leafn_clean)
summary(leafn_biome_lm)

# Seed Mass
seed_mass_lm <- lm(SeedMass ~ Biome, data = seedmass_clean)
summary(seed_mass_lm)

## 11.5 Figure and table of the linear model results ---------------------------

# The biome coefficient is the difference in mean log-trait (tundra - boreal).
# exp(coefficient) = ratio of geometric means (tundra / boreal), which is
# easier to interpret. A Welch t-test p-value is included as a check, because
# lm() assumes equal variances in both biomes and group sizes are unequal.

# Create a list of the trait labels
trait_labels <- c(PlantHeight = "Plant Height", SLA = "SLA",
                  LeafN = "Leaf N", SeedMass = "Seed Mass")

# Create a list of the lm output
trait_lms <- list(PlantHeight = plant_height_biome_lm, SLA = sla_biome_lm,
                  LeafN = leafn_biome_lm, SeedMass = seed_mass_lm)

# Create a list of the data used in the lm
trait_lm_data <- list(PlantHeight = plant_height_clean, SLA = sla_clean,
                      LeafN = leafn_clean, SeedMass = seedmass_clean)

# Crete a table of all the results
trait_lm_results <- bind_rows(lapply(names(trait_lms), function(tr) {
  mod <- trait_lms[[tr]]
  d <- trait_lm_data[[tr]]
  co <- summary(mod)$coefficients["Biometundra", ]
  ci <- confint(mod)["Biometundra", ]
  welch <- t.test(d[[tr]] ~ d$Biome, var.equal = FALSE)
  data.frame(Trait = tr,
             n_boreal = sum(d$Biome == "boreal"),
             n_tundra = sum(d$Biome == "tundra"),
             estimate_log = co[["Estimate"]],
             SE = co[["Std. Error"]],
             t_value = co[["t value"]],
             p = co[["Pr(>|t|)"]],
             lower95_log = ci[[1]],
             upper95_log = ci[[2]],
             R2 = summary(mod)$r.squared,
             p_welch = welch$p.value)
})) |>
  mutate(ratio = exp(estimate_log),
         ratio_lower95 = exp(lower95_log),
         ratio_upper95 = exp(upper95_log),
         p_holm = p.adjust(p, method = "holm"),
         p_welch_holm = p.adjust(p_welch, method = "holm"))

# Check the table
print(trait_lm_results, digits = 3)

# Forest plot of tundra / boreal ratios
(lm_ratio_plot <- trait_lm_results |>
  mutate(Trait_label = factor(trait_labels[Trait], levels = rev(trait_labels[key_traits])),
         Significance = ifelse(p_holm < 0.05, "Significant (Holm p < 0.05)", "Not significant")) |>
  ggplot(aes(x = ratio, y = Trait_label, colour = Significance)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.8) +
  geom_errorbar(aes(xmin = ratio_lower95, xmax = ratio_upper95),
                width = 0.2, linewidth = 1, orientation = "y") +
  geom_point(size = 4) +
  geom_text(aes(label = paste0("n = ", n_boreal, " / ", n_tundra)),
            vjust = -1.3, size = 0.8 * fig_annot_text / .pt, colour = "black") +
  scale_x_log10() +
  scale_colour_manual(values = c("Significant (Holm p < 0.05)" = "steelblue",
                                 "Not significant" = "grey60"),
                      name = NULL) +
  labs(x = "Tundra / boreal ratio (geometric means, 95% CI)", y = NULL) +
  theme_traitspace(legend_position = "bottom"))

# Save figure to file
ggsave(here("figures", "FigureS6_RQ1_trait_lm_ratios.png"),
       plot = lm_ratio_plot, width = 9, height = 6, dpi = 600)

# Create a formatted table for the supplement
trait_lm_table <- trait_lm_results |>
  transmute(Trait = trait_labels[Trait],
            `n (boreal / tundra)` = paste0(n_boreal, " / ", n_tundra),
            `Estimate (log)` = sprintf("%.2f", estimate_log),
            SE = sprintf("%.2f", SE),
            t = sprintf("%.2f", t_value),
            `R²` = sprintf("%.3f", R2),
            `Ratio tundra/boreal (95% CI)` = paste0(sprintf("%.2f", ratio), " (",
                                                    sprintf("%.2f", ratio_lower95), "–",
                                                    sprintf("%.2f", ratio_upper95), ")"),
            p = format.pval(p, digits = 3, eps = 1e-4),
            `p (Holm)` = format.pval(p_holm, digits = 3, eps = 1e-4),
            `Welch p (Holm)` = format.pval(p_welch_holm, digits = 3, eps = 1e-4))

# Quick check of the table
print(trait_lm_table)

# Save the table to file
write.csv(trait_lm_results,
          here("data", "derived_data", "RQ1_trait_lm_full_results.csv"), row.names = FALSE)
write.csv(trait_lm_table,
          here("data", "derived_data", "TableS_RQ1_trait_lm.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

## 11.6. Check representativeness of NMDS subset -------------------------------

# For the NMDS, we only used species with data for all 4 traits
# Here, we check that the species used are a representative subsample of all species for each trait and biome

# To do this, we draw many random subset of the same size from the pool and compare the NMDS subset with them
# (permutation / randomisation test), for both the mean and the SD of the trait.

# Notes:
#  - The NMDS subset is PART of the pool, so the two cannot be compared with a
#    two-group test (lm or t-test): the groups are not independent. The
#    permutation test handles this correctly.
#  - Uses pairwise_df_all (no outlier removal), i.e. the values the NMDS used.
#  - Small subsets give low power, so non-significant does not prove the subset
#    is representative. Look at the standardised difference as well:
#    |std diff| < 0.2 negligible, 0.2-0.5 small, 0.5-0.8 moderate, > 0.8 large.
#  - p-values are not corrected for multiple tests on purpose: here the concern
#    is missing a real difference, and correction would make that more likely.

# Get the data from the trait matrix
nmds_species <- rownames(complete_trait_matrix_4trait)

# Set number of permutations
n_perm_repr <- 9999

# Create a function to check representativeness
check_representativeness <- function(values, in_subset, n_perm) {
  pool_mean <- mean(values)
  pool_sd <- sd(values)
  sub <- values[in_subset]
  n_sub <- length(sub)
  obs_mean <- mean(sub)
  obs_sd <- sd(sub)
  
  # random subsets of the same size from the pool (without replacement)
  perm <- replicate(n_perm, {
    s <- sample(values, n_sub, replace = FALSE)
    c(mean(s), sd(s))
  })
  perm_means <- perm[1, ]
  perm_sds <- perm[2, ]
  
  # two-sided p-values (with +1 correction so p is never exactly 0)
  p_mean <- (sum(abs(perm_means - pool_mean) >= abs(obs_mean - pool_mean)) + 1) / (n_perm + 1)
  p_sd <- (sum(abs(perm_sds - mean(perm_sds)) >= abs(obs_sd - mean(perm_sds))) + 1) / (n_perm + 1)
  
  list(summary = data.frame(n_pool = length(values),
                            n_subset = n_sub,
                            pct_of_pool = 100 * n_sub / length(values),
                            pool_mean = pool_mean,
                            subset_mean = obs_mean,
                            std_mean_diff = (obs_mean - pool_mean) / pool_sd,
                            ratio_subset_pool = exp(obs_mean - pool_mean),
                            percentile_of_mean = 100 * mean(perm_means <= obs_mean),
                            p_mean = p_mean,
                            pool_sd = pool_sd,
                            subset_sd = obs_sd,
                            p_sd = p_sd),
       perm_means = perm_means)
}

# Set the seed
set.seed(532826)

# Create empty lists to sacve the results
repr_results <- list()
repr_perm <- list()

# Check representativeness
for (tr in key_traits) {
  for (b in c("boreal", "tundra")) {
    d <- pairwise_df_all |>
      filter(Biome == b, !is.na(.data[[tr]]))
    in_sub <- d$Species %in% nmds_species
    if (sum(in_sub) < 2) next
    res <- check_representativeness(d[[tr]], in_sub, n_perm_repr)
    repr_results[[paste(tr, b)]] <- cbind(Trait = tr, Biome = b, res$summary)
    repr_perm[[paste(tr, b)]] <- data.frame(Trait = tr, Biome = b,
                                            perm_mean = res$perm_means,
                                            subset_mean = res$summary$subset_mean,
                                            pool_mean = res$summary$pool_mean)
  }
}

# Combine everything into single table and df
repr_table <- bind_rows(repr_results)
repr_perm_df <- bind_rows(repr_perm)

# Check: the subset size should equal the number of NMDS species of that biome
n_nmds_per_biome <- table(nmds_plot_data_4trait$caff_biome_category)
if (any(repr_table$n_subset != as.numeric(n_nmds_per_biome[repr_table$Biome]))) {
  stop("ERROR: subset sizes do not match the number of NMDS species per biome!")
} # all ok!

# Check the table
print(repr_table, digits = 3)

# Formatted table for the supplement
repr_table_formatted <- repr_table |>
  transmute(Trait = trait_labels[Trait],
            Biome = biome_labels[Biome],
            `NMDS species / all species` = paste0(n_subset, " / ", n_pool,
                                                  " (", round(pct_of_pool), "%)"),
            `Mean all (log)` = sprintf("%.2f", pool_mean),
            `Mean NMDS (log)` = sprintf("%.2f", subset_mean),
            `Standardised difference` = sprintf("%.2f", std_mean_diff),
            `Ratio NMDS/all` = sprintf("%.2f", ratio_subset_pool),
            `Percentile` = sprintf("%.0f", percentile_of_mean),
            `p (mean)` = format.pval(p_mean, digits = 3, eps = 1e-4),
            `SD all` = sprintf("%.2f", pool_sd),
            `SD NMDS` = sprintf("%.2f", subset_sd),
            `p (SD)` = format.pval(p_sd, digits = 3, eps = 1e-4))

# Check the table
print(repr_table_formatted)

# Save table to file to use in the supplementary
write.csv(repr_table,
          here("data", "derived_data", "RQ1_NMDS_subset_representativeness_full.csv"),
          row.names = FALSE)
write.csv(repr_table_formatted,
          here("data", "derived_data", "TableS_RQ1_NMDS_subset_representativeness.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

## 11.7. Representativeness check figure ---------------------------------------

# Create labels
facet_labels <- labeller(Trait = as_labeller(trait_labels),
                         Biome = as_labeller(biome_labels))

# (a) Distribution of each trait: all species vs NMDS species
repr_long <- pairwise_df_all |>
  pivot_longer(all_of(key_traits), names_to = "Trait", values_to = "value") |>
  filter(!is.na(value))

# Combine into a single df
repr_density_data <- bind_rows(
  repr_long |> mutate(Group = "All species with data"),
  repr_long |> filter(Species %in% nmds_species) |> mutate(Group = "NMDS species")) |>
  mutate(Trait = factor(Trait, levels = key_traits))

# Group by trait and summarise the mean value
repr_density_means <- repr_density_data |>
  group_by(Trait, Biome, Group) |>
  summarise(mean_value = mean(value), .groups = "drop")

# Plot the representativeness plot
(repr_density_plot <- ggplot(repr_density_data, aes(x = value, fill = Group, colour = Group)) +
  geom_density(alpha = 0.35, linewidth = 0.7) +
  geom_vline(data = repr_density_means,
             aes(xintercept = mean_value, colour = Group),
             linetype = "dashed", linewidth = 0.8) +
  facet_grid(Biome ~ Trait, scales = "free", labeller = facet_labels) +
  scale_fill_manual(values = c("All species with data" = "grey50", "NMDS species" = "#D55E00"),
                    name = NULL) +
  scale_colour_manual(values = c("All species with data" = "grey30", "NMDS species" = "#D55E00"),
                      name = NULL) +
  labs(x = "Trait value (log)", y = "Density") +
  theme_traitspace(legend_position = "bottom") +
  theme(strip.text = element_text(size = fig_axis_title, face = "bold"),
        strip.background = element_blank()))

# Save to file
ggsave(here("figures", "FigureS7_RQ1_NMDS_subset_distributions.png"),
       plot = repr_density_plot, width = 16, height = 8, dpi = 600)

# (b) Permutation distributions: where does the NMDS subset mean fall among
#     the means of random subsets of the same size?
repr_perm_df <- repr_perm_df |>
  mutate(Trait = factor(Trait, levels = key_traits))

# Keep only unique values
repr_perm_lines <- repr_perm_df |>
  distinct(Trait, Biome, subset_mean, pool_mean)

# Factor traits and create nice labels
repr_perm_labels <- repr_table |>
  mutate(Trait = factor(Trait, levels = key_traits),
         label = paste0(round(percentile_of_mean), "th percentile\np = ",
                        format.pval(p_mean, digits = 2, eps = 1e-4)))

# Plot figure
(repr_perm_plot <- ggplot(repr_perm_df, aes(x = perm_mean)) +
  geom_histogram(bins = 40, fill = "grey75", colour = "grey55") +
  geom_vline(data = repr_perm_lines, aes(xintercept = pool_mean),
             colour = "grey30", linetype = "dashed", linewidth = 0.8) +
  geom_vline(data = repr_perm_lines, aes(xintercept = subset_mean),
             colour = "#D55E00", linewidth = 1.2) +
  geom_text(data = repr_perm_labels, aes(x = -Inf, y = Inf, label = label),
            hjust = -0.05, vjust = 1.3, size = 0.8 * fig_annot_text / .pt) +
  facet_grid(Biome ~ Trait, scales = "free", labeller = facet_labels) +
  labs(x = "Mean of random subsets (log trait value)", y = "Number of random subsets") +
  theme_traitspace(legend_position = "none") +
  theme(strip.text = element_text(size = fig_axis_title, face = "bold"),
        strip.background = element_blank()))

# Save figure to file
ggsave(here("figures", "FigureS8_RQ1_NMDS_subset_permutation.png"),
       plot = repr_perm_plot, width = 16, height = 8, dpi = 600)

# END OF SCRIPT ----------------------------------------------------------------