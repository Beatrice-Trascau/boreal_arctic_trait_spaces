##----------------------------------------------------------------------------##
# PAPER 3: BOREAL AND ARCTIC PLANT SPECIES TRAIT SPACES 
# 1.1_TRY_filtering
# This script contains code which loads/installs necessary packages and defines
# functions used in the analysis
##----------------------------------------------------------------------------##

# 1. TRY DATA -----------------------------------------------------------------

# Setup
library(here)
source(here("scripts", "0_setup.R"))

# Load data
load(here("data", "raw_data", "try_beatrice.RData"))
try_raw <- try.final.control

# Inspect data
glimpse(try_raw)

# 2. DOWNLOAD BIOMES -----------------------------------------------------------

# Load WWF biomes
# Citation: Olson, D. M., Dinerstein, E., Wikramanayake, E. D., Burgess, N. D., Powell, G. V. N., Underwood, E. C., D'Amico, J. A., Itoua, I., Strand, H. E., Morrison, J. C., Loucks, C. J., Allnutt, T. F., Ricketts, T. H., Kura, Y., Lamoreux, J. F., Wettengel, W. W., Hedao, P., Kassem, K. R. 2001. Terrestrial ecoregions of the world: a new map of life on Earth. Bioscience 51(11):933-938.
global_biomes <- st_read (here("data", "raw_data", "biomes", "wwf_terr_ecos.shp"))

#Load Boreal Forest (BIOME = 6)
boreal_forest <- st_union(global_biomes[global_biomes$BIOME == 6,])

# Load Tundra (BIOME = 11)
tundra <- st_union(global_biomes[global_biomes$BIOME == 11 &(global_biomes$REALM == "PA"|global_biomes$REALM == "NA"), ])

# Make sure boreal_forest and tundra are valid geometries
boreal_forest <- st_make_valid(boreal_forest)
tundra <- st_make_valid(tundra)

# Check CRS
biome_crs <- st_crs(boreal_forest)

# 3. FILTER TRY RECORDS --------------------------------------------------------

## 3.1. Convert TRY database to spatial object ---------------------------------

# Convert TRY to sf object
try_sf <- st_as_sf(try_raw, coords = c("LON_site", "LAT_site"), crs = 4326) |>
  st_transform(biome_crs)  # make sure CRS is the same as the biomes

# Make sure that the CRS of the biomes and TRY data match
cat("CRS of TRY data:", st_crs(try_sf)$input, "\n")
cat("CRS of boreal forest:", st_crs(boreal_forest)$input, "\n")
cat("CRS of tundra:", st_crs(tundra)$input, "\n") # looks good

# Check geometries are valid
boreal_forest <- st_make_valid(boreal_forest)
tundra <- st_make_valid(tundra)
try_sf <- st_make_valid(try_sf)

## 3.2. Filter out records falling outside of the biome boundaries -------------

# Filter records within biomes
boreal_try <- st_intersects(try_sf, boreal_forest)
tundra_try <- st_intersects(try_sf, tundra)

# Convert the sparse geometry index to regular indices
boreal_indices <- which(lengths(boreal_try) > 0)
tundra_indices <- which(lengths(tundra_try) > 0)

# Extract original data for both biomes
boreal_try_data <- try_raw[boreal_indices, ]
tundra_try_data <- try_raw[tundra_indices, ]

# Add biome classification
boreal_try_data$biome <- "boreal"
tundra_try_data$biome <- "tundra"

# Combine the two
try_filtered <- rbind(boreal_try_data, tundra_try_data)
try_filtered_all_data_27July2025 <- rbind(boreal_try_data, tundra_try_data)

# Compare number of records in the "raw" and "filtered" trait datasets
nrow(try_raw) # 130 213
nrow(try_filtered_all_data_27July2025) # 63 094

# 4. CREATE SPECIES LIST -------------------------------------------------------

# Extract unique species for each biome
boreal_species <- unique(boreal_try_data$AccSpeciesName)
tundra_species <- unique(tundra_try_data$AccSpeciesName)

# Get species present in both biomes
common_species <- intersect(boreal_species, tundra_species)

# Find unique species to each biome
unique_boreal_species <- setdiff(boreal_species, tundra_species)
unique_tundra_species <- setdiff(tundra_species, boreal_species)

# Create summary dataframe
all_filtered_species <- data.frame(SpeciesName = c(boreal_species, setdiff(tundra_species, 
                                                                  boreal_species)),
                          Biome = c(rep("Boreal", length(boreal_species)), 
                                    rep("Tundra", length(setdiff(tundra_species, boreal_species)))),
                          SharedAcrossBiomes = c(boreal_species %in% common_species, 
                         setdiff(tundra_species, boreal_species) %in% common_species))


# Save filtered data to file
#save(try_filtered_all_data_27July2025, file = here("data", "derived_data", "try_filtered_all_data_27July2025.RData"))
#save(all_filtered_species, file = here("data", "derived_data", "try_filtered.RData"))

# 5. CHECK CLASSIFICATION OF SPECIES -------------------------------------------

# Load Mariana's dataset
classification <- read.csv(here("data", "raw_data", "species_ab_class_jan2025.csv"))

# Check types of classification
unique(classification$ClassNew)

# Standardise classifications for ease of comparison
classification <- classification |>
  mutate(ClassNew = ifelse(ClassNew %in% c("Boreal-tundra boundary", "Ubiquitous"),
                           "BorealArctic", ClassNew))

# Ensure species names are in the same formats in both dfs
comparison_results <- all_filtered_species |>
  # select only relevant columns
  select(SpeciesName, Biome) |>
  # join with Mariana's df
  left_join(classification |>
              select(SpeciesName = SPECIES_CLEAN, PaperClassification = ClassNew), 
            by = "SpeciesName") |>
  # create a match indicator for species that are in both datasets
  mutate(FoundInBoth = !is.na(PaperClassification),
         ClassificationMatch = case_when(is.na(PaperClassification) ~ NA,
                                         Biome == PaperClassification ~ TRUE,
                                         TRUE ~ FALSE))

# Extract list of discrepancies
discrepancies <- comparison_results |>
  filter(FoundInBoth & !ClassificationMatch) |>
  select(SpeciesName, MyClassification = Biome, PaperClassification)

# 6. PLOT MAP WITH BIOME AND DATAPOITNS ----------------------------------------

## 6.1. Map with biomes and only datapoints within the biomes ------------------

# Get world basemap
world <- ne_countries(scale = "medium", returnclass = "sf")

# Define polar projection
proj_choice<-"+proj=laea +lat_0=90 +lon_0=0 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs "

# Transform all spatial objects to the polar projection
world_polar <- st_transform(world, crs = proj_choice)
boreal_forest_polar <- st_transform(boreal_forest, crs = proj_choice)
tundra_polar <- st_transform(tundra, crs = proj_choice)
try_sf_polar <- st_transform(try_sf, crs = proj_choice)

# Define custom colors for biomes
mycols <- c("#6EC90D", "skyblue")

# Combine the boreal and tundra data for the legend
biome_data <- rbind(cbind(st_drop_geometry(data.frame(BIOME = 6)), 
                          st_geometry(boreal_forest_polar)) |> st_sf(),
                    cbind(st_drop_geometry(data.frame(BIOME = 11)), 
                          st_geometry(tundra_polar)) |> st_sf())

# Define dataframe for the point legend
point_legend_data <- data.frame(category = c("Within Boreal Forest Biome",
                                             " Within Tundra Biome",
                                             "Outside of Biomes"))

# Create map with polar projection
(biomes_polar_proj <- ggplot() +
  # add biomes
  geom_sf(data = biome_data, aes(fill = factor(BIOME)), color = NA) +
  # add country outlines
  geom_sf(data = world_polar, fill = NA, color = "darkgray", size = 0.2) +
  # add points
  geom_sf(data = try_sf_polar[boreal_indices,], color = "black", size = 1, alpha = 0.7) +
  geom_sf(data = try_sf_polar[tundra_indices,], color = "#2121C7", size = 1, alpha = 0.7) +
  # set extent
  coord_sf(crs = proj_choice, 
           ylim = c(-703086, 7071423), 
           xlim = c(-505347.4, 8526158)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = c(0.2, 0.1),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        plot.title = element_text(face = "bold", size = 12),
        panel.grid.major = element_line(color = "gray90", linetype = "dashed"),
        panel.grid.minor = element_blank()) +
  # add custom colours and legend values
  scale_fill_manual(values = mycols,
                    labels = c("Boreal Forest Biome", "Tundra Biome"),
                    name = "Biome") +
  # add point legend
  guides(fill = guide_legend(title = "Biome"),
         color = guide_legend(title = "Records")))

# Save figure
ggsave(filename = here("figures", "Figure1_TRY_datapoints_within_biomes.png"),
       plot = biomes_polar_proj, width = 16, height = 12, dpi = 300)

## 6.2. Map with biomes and all datapoints -------------------------------------

# Recalculate intersections using the polar-projected data to ensure consistency
boreal_try_polar <- st_intersects(try_sf_polar, boreal_forest_polar)
tundra_try_polar <- st_intersects(try_sf_polar, tundra_polar)

# Convert to regular indices
boreal_indices_polar <- which(lengths(boreal_try_polar) > 0)
tundra_indices_polar <- which(lengths(tundra_try_polar) > 0)

# Get indices for points outside both biomes
all_indices_polar <- 1:nrow(try_sf_polar)
outside_indices_polar <- setdiff(all_indices_polar, c(boreal_indices_polar, tundra_indices_polar))

# Create map with polar projection
(all_points_biomes_polar_proj <- ggplot() +
  # add biomes
  geom_sf(data = biome_data, aes(fill = factor(BIOME)), color = NA) +
  # add country outlines
  geom_sf(data = world_polar, fill = NA, color = "darkgray", size = 0.2) +
  # add points
  geom_sf(data = try_sf_polar[outside_indices_polar,], aes(color = "Outside Target Biomes"), size = 0.8, alpha = 0.5) +
  geom_sf(data = try_sf_polar[boreal_indices_polar,], aes(color = "Inside Boreal Forest Biome"), size = 0.8, alpha = 0.7) +
  geom_sf(data = try_sf_polar[tundra_indices_polar,], aes(color = "Inside Tundra Biome"), size = 0.8, alpha = 0.7) +
  # set extent
  coord_sf(crs = proj_choice, 
           ylim = c(-703086, 7071423), 
           xlim = c(-505347.4, 8526158)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = c(0.15, 0.125),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        plot.title = element_text(face = "bold"),
        panel.grid.major = element_line(color = "gray90", linetype = "dashed"),
        panel.grid.minor = element_blank()) +
  # add custome legend for points
  scale_color_manual(values = c("Inside Boreal Forest Biome" = "black", 
                                "Inside Tundra Biome" = "#2121C7", 
                                "Outside Target Biomes" = "#8B5F65")) +
  # add custom colours and legend values
  scale_fill_manual(values = mycols,
                    labels = c("Boreal Forest Biome", "Tundra Biome"),
                    name = "Biome") +
  # change legends size
  guides(fill = guide_legend(title = "Biome"),
         color = guide_legend(title = "Records",
                              override.aes = list(size = 2.5))))

# Save figure
# ggsave(filename = here("figures", "FigureS1_biomes_and_all_points.png"),
#        plot = biomes_polar_proj, width = 16, height = 12, dpi = 300)

## 6.3. Maps with biomes and Plant Height datapoints ---------------------------

# Keep only plant height records
plant_height_records <- try_sf_polar |>
  filter(TraitNameNew == "PlantHeight")

# Recalculate intersections using the polar-projected data to ensure consistency
ph_boreal_try_polar <- st_intersects(plant_height_records, boreal_forest_polar)
ph_tundra_try_polar <- st_intersects(plant_height_records, tundra_polar)

# Convert to regular indices
ph_boreal_indices_polar <- which(lengths(ph_boreal_try_polar) > 0)
ph_tundra_indices_polar <- which(lengths(ph_tundra_try_polar) > 0)

# Get indices for points outside both biomes
ph_all_indices_polar <- 1:nrow(plant_height_records)
ph_outside_indices_polar <- setdiff(ph_all_indices_polar, 
                                    c(ph_boreal_indices_polar, ph_tundra_indices_polar))

# Create map
(plant_height_points <- ggplot() +
    # add biomes
    geom_sf(data = biome_data, aes(fill = factor(BIOME)), color = NA) +
    # add country outlines
    geom_sf(data = world_polar, fill = NA, color = "darkgray", size = 0.2) +
    # add points
    geom_sf(data = plant_height_records[ph_outside_indices_polar,], aes(color = "Outside Target Biomes"), size = 0.8, alpha = 0.5) +
    geom_sf(data = plant_height_records[ph_boreal_indices_polar,], aes(color = "Inside Boreal Forest Biome"), size = 0.8, alpha = 0.7) +
    geom_sf(data = plant_height_records[ph_tundra_indices_polar,], aes(color = "Inside Tundra Biome"), size = 0.8, alpha = 0.7) +
    # set extent
    coord_sf(crs = proj_choice, 
             ylim = c(-703086, 7071423), 
             xlim = c(-505347.4, 8526158)) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "none",
          plot.title = element_text(face = "bold"),
          panel.grid.major = element_line(color = "gray90", linetype = "dashed"),
          panel.grid.minor = element_blank()) +
    # add custome legend for points
    scale_color_manual(values = c("Inside Boreal Forest Biome" = "black", 
                                  "Inside Tundra Biome" = "#2121C7", 
                                  "Outside Target Biomes" = "#8B5F65")) +
    # add custom colours and legend values
    scale_fill_manual(values = mycols,
                      labels = c("Boreal Forest Biome", "Tundra Biome"),
                      name = "Biome") +
    # change legends size
    guides(fill = guide_legend(title = "Biome"),
           color = guide_legend(title = "Records",
                                override.aes = list(size = 2.5))))

## 6.4. Map with biomes and SLA points -----------------------------------------

# Keep only sla records
sla_records <- try_sf_polar |>
  filter(TraitNameNew == "SLA")

# Recalculate intersections using the polar-projected data to ensure consistency
sla_boreal_try_polar <- st_intersects(sla_records, boreal_forest_polar)
sla_tundra_try_polar <- st_intersects(sla_records, tundra_polar)

# Convert to regular indices
sla_boreal_indices_polar <- which(lengths(sla_boreal_try_polar) > 0)
sla_tundra_indices_polar <- which(lengths(sla_tundra_try_polar) > 0)

# Get indices for points outside both biomes
sla_all_indices_polar <- 1:nrow(sla_records)
sla_outside_indices_polar <- setdiff(sla_all_indices_polar, 
                                    c(sla_boreal_indices_polar, sla_tundra_indices_polar))

# Create map
(sla_points <- ggplot() +
    # add biomes
    geom_sf(data = biome_data, aes(fill = factor(BIOME)), color = NA) +
    # add country outlines
    geom_sf(data = world_polar, fill = NA, color = "darkgray", size = 0.2) +
    # add points
    geom_sf(data = sla_records[sla_outside_indices_polar,], aes(color = "Outside Target Biomes"), size = 0.8, alpha = 0.5) +
    geom_sf(data = sla_records[sla_boreal_indices_polar,], aes(color = "Inside Boreal Forest Biome"), size = 0.8, alpha = 0.7) +
    geom_sf(data = sla_records[sla_tundra_indices_polar,], aes(color = "Inside Tundra Biome"), size = 0.8, alpha = 0.7) +
    # set extent
    coord_sf(crs = proj_choice, 
             ylim = c(-703086, 7071423), 
             xlim = c(-505347.4, 8526158)) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "none",
          plot.title = element_text(face = "bold"),
          panel.grid.major = element_line(color = "gray90", linetype = "dashed"),
          panel.grid.minor = element_blank()) +
    # add custome legend for points
    scale_color_manual(values = c("Inside Boreal Forest Biome" = "black", 
                                  "Inside Tundra Biome" = "#2121C7", 
                                  "Outside Target Biomes" = "#8B5F65")) +
    # add custom colours and legend values
    scale_fill_manual(values = mycols,
                      labels = c("Boreal Forest Biome", "Tundra Biome"),
                      name = "Biome") +
    # change legends size
    guides(fill = guide_legend(title = "Biome"),
           color = guide_legend(title = "Records",
                                override.aes = list(size = 2.5))))

## 6.5. Map with biomes and Seed Mass points -----------------------------------------

# Keep only sla records
seed_mass_records <- try_sf_polar |>
  filter(TraitNameNew == "SeedMass")

# Recalculate intersections using the polar-projected data to ensure consistency
seed_mass_boreal_try_polar <- st_intersects(seed_mass_records, boreal_forest_polar)
seed_mass_tundra_try_polar <- st_intersects(seed_mass_records, tundra_polar)

# Convert to regular indices
seed_mass_boreal_indices_polar <- which(lengths(seed_mass_boreal_try_polar) > 0)
seed_mass_tundra_indices_polar <- which(lengths(seed_mass_tundra_try_polar) > 0)

# Get indices for points outside both biomes
seed_mass_all_indices_polar <- 1:nrow(seed_mass_records)
seed_mass_outside_indices_polar <- setdiff(seed_mass_all_indices_polar, 
                                     c(seed_mass_boreal_indices_polar, 
                                       seed_mass_tundra_indices_polar))

# Create map
(seed_mass_points <- ggplot() +
    # add biomes
    geom_sf(data = biome_data, aes(fill = factor(BIOME)), color = NA) +
    # add country outlines
    geom_sf(data = world_polar, fill = NA, color = "darkgray", size = 0.2) +
    # add points
    geom_sf(data = seed_mass_records[seed_mass_outside_indices_polar,], aes(color = "Outside Target Biomes"), size = 0.8, alpha = 0.5) +
    geom_sf(data = seed_mass_records[seed_mass_boreal_indices_polar,], aes(color = "Inside Boreal Forest Biome"), size = 0.8, alpha = 0.7) +
    geom_sf(data = seed_mass_records[seed_mass_tundra_indices_polar,], aes(color = "Inside Tundra Biome"), size = 0.8, alpha = 0.7) +
    # set extent
    coord_sf(crs = proj_choice, 
             ylim = c(-703086, 7071423), 
             xlim = c(-505347.4, 8526158)) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "none",
          plot.title = element_text(face = "bold"),
          panel.grid.major = element_line(color = "gray90", linetype = "dashed"),
          panel.grid.minor = element_blank()) +
    # add custome legend for points
    scale_color_manual(values = c("Inside Boreal Forest Biome" = "black", 
                                  "Inside Tundra Biome" = "#2121C7", 
                                  "Outside Target Biomes" = "#8B5F65")) +
    # add custom colours and legend values
    scale_fill_manual(values = mycols,
                      labels = c("Boreal Forest Biome", "Tundra Biome"),
                      name = "Biome") +
    # change legends size
    guides(fill = guide_legend(title = "Biome"),
           color = guide_legend(title = "Records",
                                override.aes = list(size = 2.5))))

## 6.6. Map with biomes and Leaf N points -----------------------------------------

# Keep only sla records
leafN_records <- try_sf_polar |>
  filter(TraitNameNew == "LeafN")

# Recalculate intersections using the polar-projected data to ensure consistency
leafN_boreal_try_polar <- st_intersects(leafN_records, boreal_forest_polar)
leafN_tundra_try_polar <- st_intersects(leafN_records, tundra_polar)

# Convert to regular indices
leafN_boreal_indices_polar <- which(lengths(leafN_boreal_try_polar) > 0)
leafN_tundra_indices_polar <- which(lengths(leafN_tundra_try_polar) > 0)

# Get indices for points outside both biomes
leafN_all_indices_polar <- 1:nrow(leafN_records)
leafN_outside_indices_polar <- setdiff(leafN_all_indices_polar,
                                       c(leafN_boreal_indices_polar,
                                         leafN_tundra_indices_polar))

# Create map
(leafN_points <- ggplot() +
    # add biomes
    geom_sf(data = biome_data, aes(fill = factor(BIOME)), color = NA) +
    # add country outlines
    geom_sf(data = world_polar, fill = NA, color = "darkgray", size = 0.2) +
    # add points
    geom_sf(data = leafN_records[leafN_outside_indices_polar,], aes(color = "Outside Target Biomes"), size = 0.8, alpha = 0.5) +
    geom_sf(data = leafN_records[leafN_boreal_indices_polar,], aes(color = "Inside Boreal Forest Biome"), size = 0.8, alpha = 0.7) +
    geom_sf(data = leafN_records[leafN_tundra_indices_polar,], aes(color = "Inside Tundra Biome"), size = 0.8, alpha = 0.7) +
    # set extent
    coord_sf(crs = proj_choice, 
             ylim = c(-703086, 7071423), 
             xlim = c(-505347.4, 8526158)) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "none",
          plot.title = element_text(face = "bold"),
          panel.grid.major = element_line(color = "gray90", linetype = "dashed"),
          panel.grid.minor = element_blank()) +
    # add custome legend for points
    scale_color_manual(values = c("Inside Boreal Forest Biome" = "black", 
                                  "Inside Tundra Biome" = "#2121C7", 
                                  "Outside Target Biomes" = "#8B5F65")) +
    # add custom colours and legend values
    scale_fill_manual(values = mycols,
                      labels = c("Boreal Forest Biome", "Tundra Biome"),
                      name = "Biome") +
    # change legends size
    guides(fill = guide_legend(title = "Biome"),
           color = guide_legend(title = "Records",
                                override.aes = list(size = 2.5))))

## 6.7. Combine all maps into single figure ------------------------------------

# Create dummy plot to extract legend
dummy_plot <- ggplot() +
  # add biomes
  geom_sf(data = biome_data, aes(fill = factor(BIOME)), color = NA) +
  geom_point(data = data.frame(x = c(1, 2, 3), y = c(1, 2, 3), 
                               type = c("Inside Boreal Forest Biome", "Inside Tundra Biome", "Outside Target Biomes")),
             aes(x = x, y = y, color = type), size = 3, alpha = 0.8) +
  scale_color_manual(values = c("Inside Boreal Forest Biome" = "black", 
                                "Inside Tundra Biome" = "#2121C7", 
                                "Outside Target Biomes" = "#8B5F65")) +
  # add the same scale_fill_manual as your maps
  scale_fill_manual(values = mycols, 
                    labels = c("Boreal Forest Biome", "Tundra Biome"),
                    name = "Biomes") +
  theme_void() +
  theme(legend.text = element_text(size = 11),
        legend.title = element_text(size = 12, face = "bold"),
        legend.spacing.y = unit(0.3, "cm")) +
  guides(fill = guide_legend(title = "Biomes", 
                             order = 1,
                             override.aes = list(alpha = 0.8, size = 0)),
         color = guide_legend(title = "Records", order = 2,
                              override.aes = list(size = 3, alpha = 0.8)))

# Extract legend from dummy plot
combined_legend <- get_legend(dummy_plot)

# Combine all maps into single panel
all_maps_with_legend <- plot_grid(plot_grid(plant_height_points, sla_points, 
                             seed_mass_points, leafN_points,
                             labels = c("Plant Height", "SLA", "Seed Mass", "Leaf N"),
                             ncol = 2, nrow = 2),
                             # add legend
                             combined_legend,
                             # adjust relative widths - 80% for plots, 20% for legend
                             rel_widths = c(0.8, 0.2),
                             ncol = 2)

# Save figure
ggsave(filename = here("figures", "FigureS2_biomes_points_trati_breakdown.png"),
       plot = all_maps_with_legend, width = 16, height = 12, dpi = 300)

# END OF SCRIPT ----------------------------------------------------------------