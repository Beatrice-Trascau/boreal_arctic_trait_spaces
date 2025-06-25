##----------------------------------------------------------------------------##
# PAPER 3: BOREAL AND ARCTIC PLANT SPECIES TRAIT SPACES 
# 3.1_exploration_biome_boundary_distance
# This script contains code which calcualtes the distance to biome boundary for
# each species in our list with a third approach
##----------------------------------------------------------------------------##

# 1. LOAD DATA -----------------------------------------------------------------

species_distances <- readRDS(here("data", "derived_data",
                                  "species_summaries_dist_to_biome_boundary_June25.rds"))

# Load species list 
load(here("data", "derived_data", "all_filtered_standardised_species.RData"))

# Extract species list 
species_list <- unique(corrected_species_list$CheckedSpeciesName)


# 2. EXPLORE DATA --------------------------------------------------------------

# Check overall structure
glimpse(species_distances)

# Check the number of species
length(unique(species_distances$species)) #7573 - many more than originally intended

# Check how many rows with NA values for distances there are 
sum(is.na(species_distances$mean_distance_boreal_km)) #273
sum(is.na(species_distances$mean_distance_tundra_km)) # 4995

# Plot violins of distances to biome and tundra
hist(species_distances$mean_distance_tundra_km)

# 3. INTERACTIVE SPECIES TABLE -------------------------------------------------

# Filter the summary to only include the species in the original list
filtered_summaries <- species_distances |>
  filter(species %in% species_list)

# Check which species from original list are missing from summaries
missing_species <- setdiff(species_list, species_distances$species)
if(length(missing_species) > 0) {
  cat("Species from original list not found in summaries:", length(missing_species), "\n")
  cat("Examples:", paste(head(missing_species, 5), collapse = ", "), "\n")
}

# Check which species in summaries are not in original list
extra_species <- setdiff(species_distances$species, species_list)
if(length(extra_species) > 0) {
  cat("Extra species in summaries (filtered out):", length(extra_species), "\n")
  cat("Examples:", paste(head(extra_species, 5), collapse = ", "), "\n")
}

# Create the interactive table with search functionality
interactive_table <- datatable(filtered_summaries,
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
    buttons = c('copy', 'csv', 'excel', 'pdf', 'print')),
  
  # add extensions for extra functionality
  extensions = c('Buttons', 'ColReorder'),
  # table styling
  class = 'cell-border stripe hover',
  # row names
  rownames = FALSE,
  # column names
  colnames = c( 'Species' = 'species',
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
                'Latitudinal Range (km)' = 'latitudinal_range_km')) |>
  
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
interactive_table
