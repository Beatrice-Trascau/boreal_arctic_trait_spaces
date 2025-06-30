##----------------------------------------------------------------------------##
# PAPER 3: BOREAL AND ARCTIC PLANT SPECIES TRAIT SPACES 
# 3.1_exploration_biome_boundary_distance
# This script contains code which calcualtes the distance to biome boundary for
# each species in our list with a third approach
##----------------------------------------------------------------------------##

# 1. LOAD DATA -----------------------------------------------------------------

# Load detailed results
detailed_results <- readRDS(here("data", "derived_data", 
                                 "species_summaries_dist_to_biome_boundary_June25.rds"))
# Load species list 
load(here("data", "derived_data", "corrected_species_list_27June2025.RData"))

# Extract species list 
species_list <- unique(corrected_species_list$CheckedSpeciesName)


# 2. EXPLORE DATA --------------------------------------------------------------

# Check overall structure
glimpse(detailed_results)

# Check the number of species
length(unique(detailed_results$species)) #8390 - how are there even more???

# Check how many rows with NA values for distances there are 
sum(is.na(detailed_results$mean_distance_boreal_km)) #261
sum(is.na(detailed_results$mean_distance_tundra_km)) #5639

# Plot violins of distances to biome and tundra
hist(detailed_results$mean_distance_tundra_km)
hist(detailed_results$mean_distance_boreal_km)
hist(detailed_results$mean_distance_km)

# 3. INTERACTIVE SPECIES TABLE -------------------------------------------------

# Filter the summary to only include the species in the original list
filtered_summaries <- detailed_results |>
  filter(species %in% species_list)

# Check which species from original list are missing from summaries
missing_species <- setdiff(species_list, detailed_results$species)
if(length(missing_species) > 0) {
  cat("Species from original list not found in summaries:", length(missing_species), "\n")
  cat("Examples:", paste(head(missing_species, 5), collapse = ", "), "\n")
} # No missing species
 
sort(missing_species)
# Check which species in summaries are not in original list
extra_species <- setdiff(detailed_results$species, species_list)
if(length(extra_species) > 0) {
  cat("Extra species in summaries (filtered out):", length(extra_species), "\n")
  cat("Examples:", paste(head(extra_species, 5), collapse = ", "), "\n")
}

# Create the interactive table with search functionality
interactive_table <- datatable(detailed_results,
  # define table options
  options = list(
    # enable search box
    searching = TRUE,
    # Set page length options
    pageLength = 25,
    lengthMenu = c(10, 25, 50, 100, -1),
    # enable column filtering
    searchCols = list(list(search = ""),  # Species column - empty search to start
      rep(list(NULL), ncol(filtered_summaries) - 1)),  # Other columns
    # set scroll options for wide tables
    scrollX = TRUE,
    scrollY = "600px",
    # enable column reordering
    colReorder = TRUE,
    # enable buttons for export
    dom = 'Bfrtip',
    buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
  
  # enable column filters
  columnDefs = list(list(targets = "_all", searchable = TRUE))),
  # add extensions for extra functionality
  extensions = c('Buttons', 'ColReorder'),
  # table styling
  class = 'cell-border stripe hover',
  # row names
  rownames = FALSE,
  # column names
  colnames = c('Species' = 'species',
                'Total Cells' = 'total_cells',
                'Total Occurrences' = 'total_occurrences',
                'Cells in Boreal' = 'cells_in_boreal',
                'Cells in Tundra' = 'cells_in_tundra',
                '% Cells Boreal' = 'pct_cells_boreal',
                '% Cells Tundra' = 'pct_cells_tundra',
                'Mean Distance (km)' = 'mean_distance_km',
                'Median Distance (km)' = 'median_distance_km',
                'SD Distance (km)' = 'sd_distance_km',
                'Min Distance (km)' = 'min_distance_km',
                'Max Distance (km)' = 'max_distance_km',
                'Mean Distance Boreal (km)' = 'mean_distance_boreal_km',
                'Median Distance Boreal (km)' = 'median_distance_boreal_km',
                'Mean Distance Tundra (km)' = 'mean_distance_tundra_km',
                'Median Distance Tundra (km)' = 'median_distance_tundra_km',
                'Weighted Mean Distance (km)' = 'weighted_mean_distance_km',
                'Latitudinal Range (km)' = 'latitudinal_range_km'),
  # add column filter
  filter = 'top') |>
  
  # format numeric columns using the NEW column names (after renaming)
  formatRound(columns = c('Mean Distance (km)', 'Median Distance (km)', 'SD Distance (km)',
                          'Min Distance (km)', 'Max Distance (km)', 'Mean Distance Boreal (km)',
                          'Median Distance Boreal (km)', 'Mean Distance Tundra (km)',
                          'Median Distance Tundra (km)', 'Weighted Mean Distance (km)',
                          'Latitudinal Range (km)'), digits = 1) |>
  
  # format percentage columns using NEW column names
  formatRound(columns = c('% Cells Boreal', '% Cells Tundra'), digits = 1) |>
  # add some styling for better readability
  formatStyle(columns = 'Species',
              fontWeight = 'bold') |>
  # highlight cells based on biome preference
  formatStyle(columns = '% Cells Boreal',
    backgroundColor = styleInterval(cuts = 50, values = c('lightblue', 'lightgreen'))) |>
  formatStyle(columns = '% Cells Tundra',
    backgroundColor = styleInterval(cuts = 50, values = c('lightcoral', 'lightgray')))


# Display the table
#interactive_table

# END OF SCRIPT ----------------------------------------------------------------