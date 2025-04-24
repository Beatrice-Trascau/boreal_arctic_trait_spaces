##----------------------------------------------------------------------------##
# PAPER 3: BOREAL AND ARCTIC PLANT SPECIES TRAIT SPACES 
# 1.1_TRY_filtering
# This script contains code which loads/installs necessary packages and defines
# functions used in the analysis
##----------------------------------------------------------------------------##

# 1. TRY DATA -----------------------------------------------------------------

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

# Save filtered data to file
#save(try_filtered, file = here("data", "derived_data", "try_filtered.RData"))

## 3.3 Extract species list ----------------------------------------------------

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

## 3.4. Remove morphospecies from list -----------------------------------------

# Filter out morphospecies from list
filtered_species_list <- all_filtered_species |>
  filter(!grepl("\\bsp\\.$|\\bsp\\b", SpeciesName) &
           # remove generic species names
           !SpeciesName %in% c("Grass", "Fern", "Unknown") &
           # remove suspect species names
           !SpeciesName %in% c("Eri sch", "Hieracium sect."))

# Check how many records were removed
cat("Original species count:", nrow(all_filtered_species), "\n") # 800
cat("Filtered species count:", nrow(filtered_species_list), "\n") # 770
cat("Number of morphospecies/generic entries removed:", 
    nrow(all_filtered_species) - nrow(filtered_species_list), "\n") # 30

# Check which records have been removed
removed_species <- anti_join(all_filtered_species, filtered_species_list, 
                             by = "SpeciesName")
print(removed_species)

# Save species list
#save(all_filtered_species, file = here("data", "derived_data", "all_filtered_species.RData"))

# 4. PLOT MAP WITH BIOME AND DATAPOITNS ----------------------------------------

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
mycols <- c("forestgreen", "skyblue")

# Combine the boreal and tundra data for the legend
biome_data <- rbind(
  cbind(st_drop_geometry(data.frame(BIOME = 6)), st_geometry(boreal_forest_polar)) |> st_sf(),
  cbind(st_drop_geometry(data.frame(BIOME = 11)), st_geometry(tundra_polar)) |> st_sf())

# Create map with polar projection
biomes_polar_proj <- ggplot() +
  # Add biomes
  geom_sf(data = biome_data, aes(fill = factor(BIOME)), color = NA) +
  # Add country outlines
  geom_sf(data = world_polar, fill = NA, color = "darkgray", size = 0.2) +
  # Add points
  geom_sf(data = try_sf_polar[boreal_indices,], color = "darkgreen", size = 0.8, alpha = 0.7) +
  geom_sf(data = try_sf_polar[tundra_indices,], color = "blue", size = 0.8, alpha = 0.7) +
  # Set extent
  coord_sf(crs = proj_choice, 
           ylim = c(-703086, 7071423), 
           xlim = c(-505347.4, 8526158)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = c(0.2, 0.9),
        plot.title = element_text(face = "bold"),
        panel.grid.major = element_line(color = "gray90", linetype = "dashed"),
        panel.grid.minor = element_blank()) +
  # Add custom colours and legend values
  scale_fill_manual(values = mycols,
                    labels = c("Boreal forest", "Arctic tundra"),
                    name = "Biome") +
  # Add point legend
  guides(fill = guide_legend(title = "Biome"),
         color = guide_legend(title = "Records"))

# Save figure
# ggsave(filename = here("figures", "Figure1_biomes_and_points.png"),
#        plot = biomes_polar_proj, width = 16, height = 12, dpi = 300)

# END OF SCRIPT ----------------------------------------------------------------