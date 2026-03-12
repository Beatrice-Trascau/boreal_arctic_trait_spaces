##----------------------------------------------------------------------------##
# PAPER 3: BOREAL AND ARCTIC PLANT SPECIES TRAIT SPACES 
# 3.3_RQ3_breakpoint_analysis
# This script contains code for breakpoint regression analysis testing whether
# trait-distance relationships differ on either side of the biome boundary
##----------------------------------------------------------------------------##

# 1. SETUP ---------------------------------------------------------------------

## 1.1. Load data --------------------------------------------------------------

# Load packages and functions
library(here)
source(here("scripts", "0_setup.R"))

# Load cleaned traits
load(here("data", "derived_data", "TRY_traits_cleaned_July2025.RData"))
cleaned_traits <- cleaned_traits_July2025

# Load distances to biome boundaries
detailed_results <- readRDS(here("data", "derived_data", 
                                 "species_summaries_dist_to_biome_boundary_June25.rds"))

# Load biome boundaries
global_biomes <- st_read(here("data", "raw_data", "biomes", 
                              "wwf_terr_ecos.shp"))

# Load CAFF quality check
caff_check <- read.xlsx(here("data", "derived_data", "caff_quality_check.xlsx"),
                        sheet = 1, skipEmptyRows = TRUE)

## 1.2. Prepare biomes ---------------------------------------------------------

# Load Boreal Forest (BIOME = 6)
boreal_forest <- st_union(global_biomes[global_biomes$BIOME == 6,])

# Load Tundra (BIOME = 11)
tundra <- st_union(global_biomes[global_biomes$BIOME == 11 &(global_biomes$REALM == "PA"|global_biomes$REALM == "NA"), ])

# Make sure geometries are valid
boreal_forest <- st_make_valid(boreal_forest)
tundra <- st_make_valid(tundra)

# Re-project biomes to North Pole Lamvert Azimuthal Equal Area (EPSG: 3574)
boreal_sf <- st_transform(boreal_forest, "EPSG:3574")
tundra_sf <- st_transform(tundra, "EPSG:3574")

# Get boundaries for distance calculations
boreal_boundary <- st_boundary(boreal_sf)
tundra_boundary <- st_boundary(tundra_sf)

# 2. ADD EXTRA INFORMATION TO CLEANED TRAITS -----------------------------------

## 2.1. Add species-level biome boundary distances -----------------------------

# Remove Elodea canadensis & hybrids
biome_boundaries <- detailed_results |>
  filter(!species == "Elodea canadensis") |>
  filter(!str_detect(species, " × "))

# REVERSE THE SIGN: multiply by -1 to follow new convention
# Original: positive = boreal, negative = tundra
# New: negative = boreal, positive = tundra
biome_boundaries <- biome_boundaries |>
  mutate(mean_distance_km = -mean_distance_km,
         median_distance_km = -median_distance_km,
         min_distance_km = -min_distance_km,
         max_distance_km = -max_distance_km,
         mean_distance_boreal_km = -mean_distance_boreal_km,
         median_distance_boreal_km = -median_distance_boreal_km,
         mean_distance_tundra_km = -mean_distance_tundra_km,
         median_distance_tundra_km = -median_distance_tundra_km,
         weighted_mean_distance_km = -weighted_mean_distance_km)

# Combine the dataframes
traits_biome_boundaries <- cleaned_traits |>
  left_join(biome_boundaries, by = c("StandardSpeciesName" = "species"))

# Rename column with mean distance to biome boundary to reflect the fact that
# it is a species-level mean
traits_biome_boundaries <- traits_biome_boundaries |>
  rename(species_level_mean_distance_km = mean_distance_km)

## 2.2. Calculate distance to biome boundary for each trait record -------------
## WITH REVERSED SIGN CONVENTION

# Filter to records with valid coordinates
traits_with_coords <- traits_biome_boundaries |>
  filter(!is.na(LON_site), !is.na(LAT_site)) |>
  filter(LON_site >= -180, LON_site <= 180,
         LAT_site >= -90, LAT_site <= 90)

# Check how many records were removed
nrow(traits_biome_boundaries)
nrow(traits_with_coords)

# Convert traits to spatial objec5
traits_sf <- st_as_sf(traits_with_coords, 
                      coords = c("LON_site", "LAT_site"), 
                      crs = 4326) |>
  st_transform(crs = "EPSG:3574")

# Check which biome each point is in
in_boreal <- st_intersects(traits_sf, boreal_sf, sparse = FALSE)[,1]
in_tundra <- st_intersects(traits_sf, tundra_sf, sparse = FALSE)[,1]

# Initialize distance column
traits_with_coords$distance_to_boundary_km <- NA

# REVERSED CONVENTION:
# Calculate distances for boreal points (NEGATIVE)
if(sum(in_boreal) > 0) {
  boreal_points <- traits_sf[in_boreal, ]
  boreal_distances <- st_distance(boreal_points, boreal_boundary)
  traits_with_coords$distance_to_boundary_km[in_boreal] <- 
    -apply(boreal_distances, 1, min) / 1000  # NEGATIVE for boreal
}

# Calculate distances for tundra points (POSITIVE)
if(sum(in_tundra) > 0) {
  tundra_points <- traits_sf[in_tundra, ]
  tundra_distances <- st_distance(tundra_points, tundra_boundary)
  traits_with_coords$distance_to_boundary_km[in_tundra] <- 
    apply(tundra_distances, 1, min) / 1000  # POSITIVE for tundra
}

# Add biome classification
traits_with_coords$biome <- NA
traits_with_coords$biome[in_boreal] <- "boreal"
traits_with_coords$biome[in_tundra] <- "tundra"

# Create site identifier
traits_with_coords <- traits_with_coords |>
  mutate(site_name = paste0(round(LON_site, 3), "_", round(LAT_site, 3)))

# Check how many records are in boreal and tundra biome
sum(traits_with_coords$biome == "boreal", na.rm = TRUE)
sum(traits_with_coords$biome == "tundra", na.rm = TRUE)

## 2.3. CAFF quality control ---------------------------------------------------

# Rename species name column in CAFF check
caff_check <- caff_check |>
  rename(StandardSpeciesName = SPECIES_CLEAN) |>
  dplyr::select(StandardSpeciesName, final.category)

# Add CAFF classification (keep all columns from traits_with_coords)
traits_with_caff <- traits_with_coords |>
  left_join(caff_check, by = "StandardSpeciesName")

# Filter to valid species only
cleaned_traits_final <- traits_with_caff |>
  mutate(caff_biome_category = final.category) |>
  filter(!is.na(caff_biome_category)) |>
  filter(caff_biome_category != "remove") |>
  filter(!is.na(distance_to_boundary_km)) |>
  filter(!is.na(biome))

# Check how many species there were originally
length(unique(traits_with_coords$StandardSpeciesName))

# Check how many species are left after CAFF filtering
length(unique(cleaned_traits_final$StandardSpeciesName))

# Check the total number of records after filtering
nrow(cleaned_traits_final)

# Create piecewise distance variables with REVERSED convention
# Now: boreal_side uses negative distances, tundra_side uses positive distances
cleaned_traits_final <- cleaned_traits_final |>
  mutate(distance_boreal_side = pmin(distance_to_boundary_km, 0),  # negative values (boreal)
         distance_tundra_side = pmax(distance_to_boundary_km, 0))  # positive values (tundra)
# 3. PLANT HEIGHT --------------------------------------------------------------

## 3.1. Prepare Plant Height data ----------------------------------------------

# Filter to PlantHeight records
plant_height_final <- cleaned_traits_final |>
  filter(TraitNameNew == "PlantHeight") |>
  filter(!is.na(CleanedTraitValue))

# Check number of records, species and sited
nrow(plant_height_final) # 29601
sum(plant_height_final$biome == "boreal") # 14072
sum(plant_height_final$biome == "tundra") # 15529
length(unique(plant_height_final$StandardSpeciesName)) # 441
length(unique(plant_height_final$site_name)) # 300

## 3.2. Fit breakpoint model ---------------------------------------------------

# Fit model
model_plant_height <- lme(log(CleanedTraitValue) ~ distance_tundra_side + distance_boreal_side,
                          random = list(StandardSpeciesName = ~ 1, site_name = ~ 1),
                          data = plant_height_final,
                          method = "REML")

# Check diagnostics
round(AIC(model_plant_height), 2)
round(BIC(model_plant_height), 2)
round(logLik(model_plant_height), 2)

# Extract coefficients
coefs_ph <- fixef(model_plant_height)
coefs_ph

# Get random effects variance
print(VarCorr(model_plant_height))

## 3.3. Get model diagnostics --------------------------------------------------

# Extract residuals
residuals_ph <- residuals(model_plant_height, type = "normalized")

# Create diagnostic plots
png(here("figures", "RQ3_PlantHeight_breakpoint_diagnostics.png"),
    width = 10, height = 8, units = "in", res = 300)

par(mfrow = c(2, 2))

plot(fitted(model_plant_height), residuals_ph,
     xlab = "Fitted values", ylab = "Normalized residuals",
     main = "Plant Height - Residuals vs Fitted")
abline(h = 0, col = "red", lty = 2)

qqnorm(residuals_ph, main = "Plant Height - Q-Q Plot")
qqline(residuals_ph, col = "red")

plot(fitted(model_plant_height), sqrt(abs(residuals_ph)),
     xlab = "Fitted values", ylab = "√|Normalized residuals|",
     main = "Plant Height - Scale-Location")

plot(plant_height_final$distance_to_boundary_km, residuals_ph,
     xlab = "Distance to boundary (km)", ylab = "Normalized residuals",
     main = "Plant Height - Residuals vs Distance")
abline(h = 0, col = "red", lty = 2)
abline(v = 0, col = "blue", lty = 2)

par(mfrow = c(1, 1))
dev.off()

## 3.4. Model visualisation ----------------------------------------------------

# Create prediction data
pred_data_ph <- data.frame(distance_to_boundary_km = seq(min(plant_height_final$distance_to_boundary_km),
                                                         max(plant_height_final$distance_to_boundary_km),
                                                         length.out = 200)) |>
  mutate(distance_tundra_side = pmin(distance_to_boundary_km, 0),
         distance_boreal_side = pmax(distance_to_boundary_km, 0))

# Get predictions
pred_data_ph$predicted <- predict(model_plant_height, newdata = pred_data_ph, level = 0)
pred_data_ph$predicted_original <- exp(pred_data_ph$predicted)

# Create plot
plot_ph <- ggplot() +
  geom_point(data = plant_height_final,
             aes(x = distance_to_boundary_km, y = CleanedTraitValue, color = biome),
             alpha = 0.2, size = 1) +
  geom_line(data = pred_data_ph,
            aes(x = distance_to_boundary_km, y = predicted_original),
            color = "#dcd0ff", linewidth = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_color_manual(values = c("boreal" = "darkgreen", "tundra" = "darkblue"),
                     name = "Biome") +
  scale_y_log10() +
  labs(x = "Distance to Biome Boundary (km)",
       y = "Plant Height (m, log scale)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"))

# Check the plot
print(plot_ph)

# Save output
ggsave(here("figures", "Figure4a_PlantHeight_breakpoint_fitted.png"),
       plot = plot_ph, width = 10, height = 6, dpi = 600)

## 3.5. Run breakpoint model with estimated breakpoint -------------------------

# Fit simple lm to find breakpoint
simple_lm_ph <- lm(log(CleanedTraitValue) ~ distance_to_boundary_km,
                   data = plant_height_final)

# Fit segmented model
tryCatch({
  seg_model_ph <- segmented(simple_lm_ph, 
                            seg.Z = ~ distance_to_boundary_km,
                            psi = 0)
  
  breakpoint_estimate_ph <- seg_model_ph$psi[, "Est."]
  breakpoint_se_ph <- seg_model_ph$psi[, "St.Err"]
  
  cat("  Estimated breakpoint:", round(breakpoint_estimate_ph, 2), "km\n")
  cat("  Standard error:", round(breakpoint_se_ph, 2), "km\n")
  cat("  95% CI: [", round(breakpoint_estimate_ph - 1.96*breakpoint_se_ph, 2), ",",
      round(breakpoint_estimate_ph + 1.96*breakpoint_se_ph, 2), "]\n\n")
  
  # Fit lme with estimated breakpoint
  plant_height_final <- plant_height_final |>
    mutate(distance_tundra_estimated = pmin(distance_to_boundary_km - breakpoint_estimate_ph, 0),
           distance_boreal_estimated = pmax(distance_to_boundary_km - breakpoint_estimate_ph, 0))
  
  model_ph_estimated <- lme(log(CleanedTraitValue) ~ distance_tundra_estimated + distance_boreal_estimated,
                            random = list(StandardSpeciesName = ~ 1, site_name = ~ 1),
                            data = plant_height_final,
                            method = "REML")
  
  cat("  Model with estimated breakpoint:\n")
  cat("    AIC:", round(AIC(model_ph_estimated), 2), "\n")
  coefs_estimated_ph <- fixef(model_ph_estimated)
  cat("    Tundra slope:", round(coefs_estimated_ph["distance_tundra_estimated"], 5), "\n")
  cat("    Boreal slope:", round(coefs_estimated_ph["distance_boreal_estimated"], 5), "\n\n")
  
  estimated_converged_ph <- TRUE
  
}, error = function(e) {
  cat("  ✗ Segmented model failed to converge\n")
  cat("  Error:", e$message, "\n\n")
  estimated_converged_ph <<- FALSE
  breakpoint_estimate_ph <<- NA
}) #breakpoint estimated at 6.6km

# Compare models
AICtab(model_plant_height, model_ph_estimated)
# dAIC df
# model_plant_height  0.0 6 
# model_ph_estimated 12.1 6 

# 4. SLA -----------------------------------------------------------------------

## 4.1. Prepare SLA data -------------------------------------------------------

# Filter to SLA records
sla_final <- cleaned_traits_final |>
  filter(TraitNameNew == "SLA") |>
  filter(!is.na(CleanedTraitValue))

# Check number of records, species and sited
nrow(sla_final) # 18023
sum(sla_final$biome == "boreal") # 2769
sum(sla_final$biome == "tundra") # 15254
length(unique(sla_final$StandardSpeciesName)) # 401
length(unique(sla_final$site_name)) # 188

## 4.2. Fit SLA breakpoint model -----------------------------------------------

# Run model
model_sla <- lme(log(CleanedTraitValue) ~ distance_tundra_side + distance_boreal_side,
                 random = list(StandardSpeciesName = ~ 1, site_name = ~ 1),
                 data = sla_final,
                 method = "REML")

# Check diagnostics
round(AIC(model_sla), 2)
round(BIC(model_sla), 2)
round(logLik(model_sla), 2)

# Check fixed effects
coefs_sla <- fixef(model_sla)
coefs_sla

# Check random effects
print(VarCorr(model_sla))

## 4.3. SLA diagnostics --------------------------------------------------------

# Extract diagnostics
residuals_sla <- residuals(model_sla, type = "normalized")

# Save diagnostic plots
png(here("figures", "RQ3_SLA_breakpoint_diagnostics.png"),
    width = 10, height = 8, units = "in", res = 300)

par(mfrow = c(2, 2))

plot(fitted(model_sla), residuals_sla,
     xlab = "Fitted values", ylab = "Normalized residuals",
     main = "SLA - Residuals vs Fitted")
abline(h = 0, col = "red", lty = 2)

qqnorm(residuals_sla, main = "SLA - Q-Q Plot")
qqline(residuals_sla, col = "red")

plot(fitted(model_sla), sqrt(abs(residuals_sla)),
     xlab = "Fitted values", ylab = "√|Normalized residuals|",
     main = "SLA - Scale-Location")

plot(sla_final$distance_to_boundary_km, residuals_sla,
     xlab = "Distance to boundary (km)", ylab = "Normalized residuals",
     main = "SLA - Residuals vs Distance")
abline(h = 0, col = "red", lty = 2)
abline(v = 0, col = "blue", lty = 2)

par(mfrow = c(1, 1))
dev.off()

## 4.4. SLA Visualisation ------------------------------------------------------

# Create prediction data
pred_data_sla <- data.frame(distance_to_boundary_km = seq(min(sla_final$distance_to_boundary_km),
                                                          max(sla_final$distance_to_boundary_km),
                                                          length.out = 200)) |>
  mutate(distance_tundra_side = pmin(distance_to_boundary_km, 0),
         distance_boreal_side = pmax(distance_to_boundary_km, 0))

# Get predictions
pred_data_sla$predicted <- predict(model_sla, newdata = pred_data_sla, level = 0)
pred_data_sla$predicted_original <- exp(pred_data_sla$predicted)

# Plot
plot_sla <- ggplot() +
  geom_point(data = sla_final,
             aes(x = distance_to_boundary_km, y = CleanedTraitValue, color = biome),
             alpha = 0.2, size = 1) +
  geom_line(data = pred_data_sla,
            aes(x = distance_to_boundary_km, y = predicted_original),
            color = "#dcd0ff", linewidth = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_color_manual(values = c("boreal" = "darkgreen", "tundra" = "darkblue"),
                     name = "Biome") +
  scale_y_log10() +
  labs(x = "Distance to Biome Boundary (km)",
       y = "SLA (mm²/mg, log scale)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"))

# Check plot
print(plot_sla)

# Save
ggsave(here("figures", "Figure4b_SLA_breakpoint_fitted.png"),
       plot = plot_sla, width = 10, height = 6, dpi = 600)

## 4.5. Alternative model with estimated breakpoint ----------------------------

# Fit simple lm to find breakpoint
simple_lm_sla <- lm(log(CleanedTraitValue) ~ distance_to_boundary_km,
                    data = sla_final)

# Fit segmented model
tryCatch({
  seg_model_sla <- segmented(simple_lm_sla, 
                             seg.Z = ~ distance_to_boundary_km,
                             psi = 0)
  
  breakpoint_estimate_sla <- seg_model_sla$psi[, "Est."]
  breakpoint_se_sla <- seg_model_sla$psi[, "St.Err"]
  
  cat("  Estimated breakpoint:", round(breakpoint_estimate_sla, 2), "km\n")
  cat("  Standard error:", round(breakpoint_se_sla, 2), "km\n")
  cat("  95% CI: [", round(breakpoint_estimate_sla - 1.96*breakpoint_se_sla, 2), ",",
      round(breakpoint_estimate_sla + 1.96*breakpoint_se_sla, 2), "]\n\n")
  
  sla_final <- sla_final %>%
    mutate(distance_tundra_estimated = pmin(distance_to_boundary_km - breakpoint_estimate_sla, 0),
           distance_boreal_estimated = pmax(distance_to_boundary_km - breakpoint_estimate_sla, 0))
  
  model_sla_estimated <- lme(log(CleanedTraitValue) ~ distance_tundra_estimated + distance_boreal_estimated,
                             random = list(StandardSpeciesName = ~ 1, site_name = ~ 1),
                             data = sla_final,
                             method = "REML")
  
  cat("  Model with estimated breakpoint:\n")
  cat("    AIC:", round(AIC(model_sla_estimated), 2), "\n")
  coefs_estimated_sla <- fixef(model_sla_estimated)
  cat("    Tundra slope:", round(coefs_estimated_sla["distance_tundra_estimated"], 5), "\n")
  cat("    Boreal slope:", round(coefs_estimated_sla["distance_boreal_estimated"], 5), "\n\n")
  
  estimated_converged_sla <- TRUE
  
}, error = function(e) {
  cat("  ✗ Segmented model failed to converge\n\n")
  estimated_converged_sla <<- FALSE
  breakpoint_estimate_sla <<- NA
}) # breakpoint estimated at 184.11km

# Compare models
AICtab(model_sla, model_sla_estimated)
# dAIC df
# model_sla            0.0 6 
# model_sla_estimated 35.9 6 

# 5. LEAF N --------------------------------------------------------------------

## 5.1. Prepare Leaf N data ----------------------------------------------------

# Filter to LeafN records
leafn_final <- cleaned_traits_final |>
  filter(TraitNameNew == "LeafN") |>
  filter(!is.na(CleanedTraitValue))

# Check number of records, species and sites
nrow(leafn_final) # 7239
sum(leafn_final$biome == "boreal") # 1766
sum(leafn_final$biome == "tundra") # 5473
length(unique(leafn_final$StandardSpeciesName)) # 423
length(unique(leafn_final$site_name)) # 210

## 5.2. Fit breakpoint model ---------------------------------------------------

# Define model
model_leafn <- lme(log(CleanedTraitValue) ~ distance_tundra_side + distance_boreal_side,
                   random = list(StandardSpeciesName = ~ 1, site_name = ~ 1),
                   data = leafn_final,
                   method = "REML")

# Check model diagnostics
round(AIC(model_leafn), 2)
round(BIC(model_leafn), 2)
round(logLik(model_leafn), 2)

# Check fixed effects
coefs_leafn <- fixef(model_leafn)
coefs_leafn

# Check random effects
print(VarCorr(model_leafn))

## 5.3. Leaf N diagnostics -----------------------------------------------------

# Extract residuals
residuals_leafn <- residuals(model_leafn, type = "normalized")

# Diagnostic plots
png(here("figures", "RQ3_LeafN_breakpoint_diagnostics.png"),
    width = 10, height = 8, units = "in", res = 300)

par(mfrow = c(2, 2))

plot(fitted(model_leafn), residuals_leafn,
     xlab = "Fitted values", ylab = "Normalized residuals",
     main = "Leaf N - Residuals vs Fitted")
abline(h = 0, col = "red", lty = 2)

qqnorm(residuals_leafn, main = "Leaf N - Q-Q Plot")
qqline(residuals_leafn, col = "red")

plot(fitted(model_leafn), sqrt(abs(residuals_leafn)),
     xlab = "Fitted values", ylab = "√|Normalized residuals|",
     main = "Leaf N - Scale-Location")

plot(leafn_final$distance_to_boundary_km, residuals_leafn,
     xlab = "Distance to boundary (km)", ylab = "Normalized residuals",
     main = "Leaf N - Residuals vs Distance")
abline(h = 0, col = "red", lty = 2)
abline(v = 0, col = "blue", lty = 2)

par(mfrow = c(1, 1))
dev.off()

## 5.4. Leaf N visualisation ---------------------------------------------------

# Create prediction data
pred_data_leafn <- data.frame(distance_to_boundary_km = seq(min(leafn_final$distance_to_boundary_km),
                                                            max(leafn_final$distance_to_boundary_km),
                                                            length.out = 200)) |>
  mutate(distance_tundra_side = pmin(distance_to_boundary_km, 0),
         distance_boreal_side = pmax(distance_to_boundary_km, 0))

# Get predictions
pred_data_leafn$predicted <- predict(model_leafn, newdata = pred_data_leafn, level = 0)
pred_data_leafn$predicted_original <- exp(pred_data_leafn$predicted)

# Plot
plot_leafn <- ggplot() +
  geom_point(data = leafn_final,
             aes(x = distance_to_boundary_km, y = CleanedTraitValue, color = biome),
             alpha = 0.2, size = 1) +
  geom_line(data = pred_data_leafn,
            aes(x = distance_to_boundary_km, y = predicted_original),
            color = "#dcd0ff", linewidth = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_color_manual(values = c("boreal" = "darkgreen", "tundra" = "darkblue"),
                     name = "Biome") +
  scale_y_log10() +
  labs(x = "Distance to Biome Boundary (km)",
       y = "Leaf N (mg/g, log scale)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"))

# Check model
print(plot_leafn)

# Save figure
ggsave(here("figures", "Figure4c_LeafN_breakpoint_fitted.png"),
       plot = plot_leafn, width = 10, height = 6, dpi = 600)

## 5.5. Alternative model with estimated breakpoint ----------------------------

# Run simple lm to find breakpoint
simple_lm_leafn <- lm(log(CleanedTraitValue) ~ distance_to_boundary_km,
                      data = leafn_final)

# Fit segmented model
tryCatch({
  seg_model_leafn <- segmented(simple_lm_leafn, 
                               seg.Z = ~ distance_to_boundary_km,
                               psi = 0)
  
  breakpoint_estimate_leafn <- seg_model_leafn$psi[, "Est."]
  breakpoint_se_leafn <- seg_model_leafn$psi[, "St.Err"]
  
  cat("  Estimated breakpoint:", round(breakpoint_estimate_leafn, 2), "km\n")
  cat("  Standard error:", round(breakpoint_se_leafn, 2), "km\n")
  cat("  95% CI: [", round(breakpoint_estimate_leafn - 1.96*breakpoint_se_leafn, 2), ",",
      round(breakpoint_estimate_leafn + 1.96*breakpoint_se_leafn, 2), "]\n\n")
  
  leafn_final <- leafn_final %>%
    mutate(distance_tundra_estimated = pmin(distance_to_boundary_km - breakpoint_estimate_leafn, 0),
           distance_boreal_estimated = pmax(distance_to_boundary_km - breakpoint_estimate_leafn, 0))
  
  model_leafn_estimated <- lme(log(CleanedTraitValue) ~ distance_tundra_estimated + distance_boreal_estimated,
                               random = list(StandardSpeciesName = ~ 1, site_name = ~ 1),
                               data = leafn_final,
                               method = "REML")
  
  cat("  Model with estimated breakpoint:\n")
  cat("    AIC:", round(AIC(model_leafn_estimated), 2), "\n")
  coefs_estimated_leafn <- fixef(model_leafn_estimated)
  cat("    Tundra slope:", round(coefs_estimated_leafn["distance_tundra_estimated"], 5), "\n")
  cat("    Boreal slope:", round(coefs_estimated_leafn["distance_boreal_estimated"], 5), "\n\n")
  
  estimated_converged_leafn <- TRUE
  
}, error = function(e) {
  cat("  ✗ Segmented model failed to converge\n\n")
  estimated_converged_leafn <<- FALSE
  breakpoint_estimate_leafn <<- NA
}) # breakpoint estimated at -107.11km

# Compare models
AICtab(model_leafn, model_leafn_estimated)
# dAIC df
# model_leafn_estimated  0.0 6 
# model_leafn           12.6 6

# 6. SEED MASS -----------------------------------------------------------------

## 6.1. Prepare Seed Mass data -------------------------------------------------

# Filter to SeedMass records
seedmass_final <- cleaned_traits_final |>
  filter(TraitNameNew == "SeedMass") |>
  filter(!is.na(CleanedTraitValue))

# Check the number of sites, species and records
nrow(seedmass_final) # 814
sum(seedmass_final$biome == "boreal") # 23
sum(seedmass_final$biome == "tundra") # 791
length(unique(seedmass_final$StandardSpeciesName)) # 123
length(unique(seedmass_final$site_name)) # 36

## 6.2. Fit Seed Mass model ----------------------------------------------------

# Define model
model_seedmass <- lme(log(CleanedTraitValue) ~ distance_tundra_side + distance_boreal_side,
                      random = list(StandardSpeciesName = ~ 1, site_name = ~ 1),
                      data = seedmass_final,
                      method = "REML")

# Check model diagnostics
round(AIC(model_seedmass), 2)
round(BIC(model_seedmass), 2)
round(logLik(model_seedmass), 2)

# Check fixed effects
coefs_seedmass <- fixef(model_seedmass)
coefs_seedmass

# Check random effects
print(VarCorr(model_seedmass))

## 6.3. Seed Mass diagnostics --------------------------------------------------

# Extract residuals
residuals_seedmass <- residuals(model_seedmass, type = "normalized")

# Plot model diagnostics
png(here("figures", "RQ3_SeedMass_diagnostics.png"),
    width = 10, height = 8, units = "in", res = 300)

par(mfrow = c(2, 2))

plot(fitted(model_seedmass), residuals_seedmass,
     xlab = "Fitted values", ylab = "Normalized residuals",
     main = "Seed Mass - Residuals vs Fitted")
abline(h = 0, col = "red", lty = 2)

qqnorm(residuals_seedmass, main = "Seed Mass - Q-Q Plot")
qqline(residuals_seedmass, col = "red")

plot(fitted(model_seedmass), sqrt(abs(residuals_seedmass)),
     xlab = "Fitted values", ylab = "√|Normalized residuals|",
     main = "Seed Mass - Scale-Location")

plot(seedmass_final$distance_to_boundary_km, residuals_seedmass,
     xlab = "Distance to boundary (km)", ylab = "Normalized residuals",
     main = "Seed Mass - Residuals vs Distance")
abline(h = 0, col = "red", lty = 2)
abline(v = 0, col = "blue", lty = 2)

par(mfrow = c(1, 1))
dev.off()

## 6.4. Seed Mass visualisation ------------------------------------------------

# Create prediction data
pred_data_seedmass <- data.frame(distance_to_boundary_km = seq(min(seedmass_final$distance_to_boundary_km),
                                                               max(seedmass_final$distance_to_boundary_km),
                                                               length.out = 200)) |>
  mutate(distance_tundra_side = pmin(distance_to_boundary_km, 0),
         distance_boreal_side = pmax(distance_to_boundary_km, 0))

# Get predictions
pred_data_seedmass$predicted <- predict(model_seedmass, newdata = pred_data_seedmass, level = 0)
pred_data_seedmass$predicted_original <- exp(pred_data_seedmass$predicted)

# Plot
plot_seedmass <- ggplot() +
  geom_point(data = seedmass_final,
             aes(x = distance_to_boundary_km, y = CleanedTraitValue, color = biome),
             alpha = 0.2, size = 1) +
  geom_line(data = pred_data_seedmass,
            aes(x = distance_to_boundary_km, y = predicted_original),
            color = "#dcd0ff", linewidth = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_color_manual(values = c("boreal" = "darkgreen", "tundra" = "darkblue"),
                     name = "Biome") +
  scale_y_log10() +
  labs(x = "Distance to Biome Boundary (km)",
       y = "Seed Mass (mg, log scale)") +
  theme_classic() +
  theme(legend.position = "right",
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

# Check plot
print(plot_seedmass)

# Save plot
ggsave(here("figures", "Figure3d_SeedMass_breakpoint_fitted.png"),
       plot = plot_seedmass, width = 10, height = 6, dpi = 600)

# Combine all model plots into a single one
trait_breakpoints <- plot_grid(plot_ph, plot_sla, plot_leafn, plot_seedmass,
                               labels = c("a)", "b)", "c)", "d"),
                               nrow = 2)

# Save combined plot
ggsave(here("figures", "Figure4_traits_breakpoint_fitted.png"),
       plot = trait_breakpoints, width = 15, height = 10, dpi = 600)
ggsave(here("figures", "Figure4_traits_breakpoint_fitted.pdf"),
       plot = trait_breakpoints, width = 15, height = 10, dpi = 600)

## 6.5. Alternative model with estimated breakpoint ----------------------------

# Fit simple lm to find breakpoint
simple_lm_seedmass <- lm(log(CleanedTraitValue) ~ distance_to_boundary_km,
                         data = seedmass_final)

# Fit segmented model
tryCatch({
  seg_model_seedmass <- segmented(simple_lm_seedmass, 
                                  seg.Z = ~ distance_to_boundary_km,
                                  psi = 0)
  
  breakpoint_estimate_seedmass <- seg_model_seedmass$psi[, "Est."]
  breakpoint_se_seedmass <- seg_model_seedmass$psi[, "St.Err"]
  
  cat("  Estimated breakpoint:", round(breakpoint_estimate_seedmass, 2), "km\n")
  cat("  Standard error:", round(breakpoint_se_seedmass, 2), "km\n")
  cat("  95% CI: [", round(breakpoint_estimate_seedmass - 1.96*breakpoint_se_seedmass, 2), ",",
      round(breakpoint_estimate_seedmass + 1.96*breakpoint_se_seedmass, 2), "]\n\n")
  
  seedmass_final <- seedmass_final %>%
    mutate(distance_tundra_estimated = pmin(distance_to_boundary_km - breakpoint_estimate_seedmass, 0),
           distance_boreal_estimated = pmax(distance_to_boundary_km - breakpoint_estimate_seedmass, 0))
  
  model_seedmass_estimated <- lme(log(CleanedTraitValue) ~ distance_tundra_estimated + distance_boreal_estimated,
                                  random = list(StandardSpeciesName = ~ 1, site_name = ~ 1),
                                  data = seedmass_final,
                                  method = "REML")
  
  cat("  Model with estimated breakpoint:\n")
  cat("    AIC:", round(AIC(model_seedmass_estimated), 2), "\n")
  coefs_estimated_seedmass <- fixef(model_seedmass_estimated)
  cat("    Tundra slope:", round(coefs_estimated_seedmass["distance_tundra_estimated"], 5), "\n")
  cat("    Boreal slope:", round(coefs_estimated_seedmass["distance_boreal_estimated"], 5), "\n\n")
  
  estimated_converged_seedmass <- TRUE
  
}, error = function(e) {
  cat("  ✗ Segmented model failed to converge\n\n")
  estimated_converged_seedmass <<- FALSE
  breakpoint_estimate_seedmass <<- NA
}) #estimated breakpoint at -0.33

# Compare models
AICtab(model_seedmass, model_seedmass_estimated) # equally good?

# 7. GAMs ----------------------------------------------------------------------

## 7.1. Plant Height GAMs ------------------------------------------------------

# Fit GAM 
gam_full <- mgcv::gamm(log(CleanedTraitValue) ~ 
                         s(distance_to_boundary_km, k=200) +  
                         s(LON_site, LAT_site, bs="tp", k=100),  
                       random = list(StandardSpeciesName = ~1),
                       data = plant_height_final,
                       method = "REML")

# Check the model
mgcv::gam.check(gam_full$gam)

# The edf for distance to biome boundary differs based on the value of k given to'
# the smoother that deals with the latitude and longitude 

# END OF SCRIPT ----------------------------------------------------------------