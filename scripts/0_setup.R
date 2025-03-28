##----------------------------------------------------------------------------##
# PAPER 3: BOREAL AND ARCTIC PLANT SPECIES TRAIT SPACES 
# 0_setup
# This script contains code which loads/installs necessary packages and defines
# functions used in the analysis
##----------------------------------------------------------------------------##

# 1. LOAD/INSTALL PACKAGES NEEDED FOR ANALYIS ----------------------------------

# Function to check to install/load packages

# Define function
install_load_package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, repos = "http://cran.us.r-project.org")
  }
  require(x, character.only = TRUE)
}

# Define list of packages
package_vec <- c("here", "terra", "sf", "geodata", "mapview",
                 "tidyverse", "dplyr", "ggplot2","gt", "cowplot",
                 "data.table","patchwork", "styler", "scales",
                 "plotly", "tidyterra", "ggspatial", "htmlwidgets",
                 "htmltools", "patchwork", "webshot2", "CoordinateCleaner",
                 "car", "kableExtra", "readr", "rnaturalearth", "rnaturalearthdata")

# Execute the function
sapply(package_vec, install_load_package)

# 2. CREATE NECESSARY FILE STRUCTURE -------------------------------------------

# Function to create the file structure needed to run the analysis smoothly
create_project_structure <- function(base_path = "boreal_arctic_trait_space") {
  # Define the directory structure
  dirs <- c(
    file.path(base_path),
    file.path(base_path, "data"),
    file.path(base_path, "scripts"),
    file.path(base_path, "figure"),
    file.path(base_path, "data", "raw_data"),
    file.path(base_path, "data", "derived_data"),
    file.path(base_path, "data", "raw_data", "biomes")
  )
  
  # Create directories if they don't exist
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      cat("Created directory:", dir, "\n")
    } else {
      cat("Directory already exists:", dir, "\n")
    }
  }
  
  cat("\nProject structure setup complete!\n")
}

# Run function
create_project_structure

# END OF SCRIPT ----------------------------------------------------------------