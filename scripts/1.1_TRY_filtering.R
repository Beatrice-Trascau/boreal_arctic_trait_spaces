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
global_biomes <-st_read (here("data", "raw_data", "biomes", "wwf_terr_ecos.shp"))

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
