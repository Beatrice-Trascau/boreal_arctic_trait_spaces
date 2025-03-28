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

# END OF SCRIPT ----------------------------------------------------------------