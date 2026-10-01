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
                 "car", "kableExtra", "readr", "rnaturalearth", "rnaturalearthdata",
                 "rgbif", "purr", "DT", "MultiTraits", "BIEN", "vegan",
                 "openxlsx", "goeveg", "moments", "gllvm", "ggExtra", "nlme",
                 "segmented", "bbmle", "mgcv", "DHARMa", "janitor",
                 "googledrive")

# Execute the function
sapply(package_vec, install_load_package)

# 2. CREATE NECESSARY FILE STRUCTURE -------------------------------------------

# Function to create the file structure needed to run the analysis smoothly
create_project_structure <- function(base_path = "boreal_arctic_trait_space") {
  # Define the directory structure
  dirs <- c(file.path(base_path),
            file.path(base_path, "data"),
            file.path(base_path, "scripts"),
            file.path(base_path, "figure"),
            file.path(base_path, "data", "raw_data"),
            file.path(base_path, "data", "WFO_Backbone"),
            file.path(base_path, "data", "derived_data"),
            file.path(base_path, "data", "raw_data", "biomes"))
  
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

# 3. DOWLOAD DATA FROM DRIVE ---------------------------------------------------

# Authenticate with Google - will open a new browser window
#drive_auth()
# When running this for the first time:
# 1. New browser window will open
# 2. You will be asked to sign in to your Google account (you will need one)
# 3. You will be asked to give permission to the googledrive package
# 4. You can close the window after you approve
# 5. A success message should appear in R

# Check that authentication worked
#drive_user() 

# Give file ID for derived_data
#file_id <- "1pHDlFHwG1KPwGuundvh-G-qbOy9H8EKS"

# Download derived_data file from drive
# drive_download(file = as_id(file_id),
#                path = here("data", "derived_data"),  
#                overwrite = FALSE)

# Give file ID for raw_data
#raw_data_file_id <- "1-aF8IVkStnY0qTNoUQSG-RO049y74fH8"

# Download raw_data file from drive
# drive_download(file = as_id(raw_data_file_id),
#                path = here("data", "raw_data"),  
#                overwrite = FALSE)

# END OF SCRIPT ----------------------------------------------------------------