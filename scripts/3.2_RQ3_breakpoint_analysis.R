##----------------------------------------------------------------------------##
# PAPER 3: BOREAL AND ARCTIC PLANT SPECIES TRAIT SPACES 
# 3.2_RQ3_breakpoint_analysis
# This script contains code for breakpoint regression analysis testing whether
# trait-distance relationships differ on either side of the biome boundary
##----------------------------------------------------------------------------##

# 1. SETUP ---------------------------------------------------------------------

## 1.1. Settings ---------------------------------------------------------------

# Grid size (degrees) for grouping nearby plots into locations
location_grid_deg <- 0.1

# Starting value of breakpoint search (km)
seg_start_psi <- 0

# Number of restarts from each starting value
seg_restarts  <- 20

# Spatial autocorrelation: distance classes of the correlogram (km)
correlogram_step_km <- 100
correlogram_max_km  <- 2000

# Figure settings (same as script 3.1)
biome_colours <- c("boreal" = "darkgreen", "tundra" = "darkblue")
biome_labels <- c("boreal" = "Boreal", "tundra" = "Tundra")
fit_line_colour <- "#dcd0ff"
segmented_label <- "Segmented (no random effects)"
fig_axis_text <- 16
fig_axis_title <- 18
fig_legend_text <- 16
fig_legend_title <- 18
fig_annot_text <- 14
fig_panel_label <- 22

## 1.2. Load packages and data -------------------------------------------------

# Load packages
library(here)
source(here("scripts", "0_setup.R"))
library(segmented)   # breakpoint regression
library(ape)         # Moran's I (spatial autocorrelation)

# Load the cleaned traits
load(here("data", "derived_data", "TRY_traits_cleaned_July2025.RData"))
cleaned_traits <- cleaned_traits_July2025

# Get biomes
global_biomes <- st_read(here("data", "raw_data", "biomes", "wwf_terr_ecos.shp"))

# CAFF quality check
caff_check <- read.xlsx(here("data", "derived_data", "caff_quality_check.xlsx"),
                        sheet = 1, skipEmptyRows = TRUE)

## 1.3. Prepare biomes and the treeline ----------------------------------------

# Load Boreal Forest (BIOME = 6)
boreal_forest <- st_union(global_biomes[global_biomes$BIOME == 6, ])

# Load Tundra (BIOME = 11)
tundra <- st_union(global_biomes[global_biomes$BIOME == 11 &
                                   (global_biomes$REALM == "PA" | global_biomes$REALM == "NA"), ])

# Re-project biomes to North Pole Lambert Azimuthal Equal Area (EPSG: 3574)
boreal_sf <- st_transform(st_make_valid(boreal_forest), "EPSG:3574")
tundra_sf <- st_transform(st_make_valid(tundra), "EPSG:3574")

# Get biome boundaries for distance calculations
boreal_boundary <- st_boundary(boreal_sf)
tundra_boundary <- st_boundary(tundra_sf)

# Extract the "treeline", i.e. the parts where the boreal and tundra biomes touch
shared_edge <- st_intersection(boreal_boundary, st_buffer(tundra_sf, 1000))
cat("  Length of treeline:", round(as.numeric(sum(st_length(shared_edge))) / 1000), "km\n") #93846

# 2. PREPARE THE DATA ----------------------------------------------------------

## 2.1. Biome and distance to the treeline for each record ----------------------

# Keep only records with valid lat and long
traits_with_coords <- cleaned_traits |>
  filter(!is.na(LON_site), !is.na(LAT_site),
         LON_site >= -180, LON_site <= 180,
         LAT_site >= -90, LAT_site <= 90)

# Convert traits to spatial object
traits_sf <- st_as_sf(traits_with_coords, coords = c("LON_site", "LAT_site"),
                      crs = 4326, remove = FALSE) |>
  st_transform(crs = "EPSG:3574")

# Check which biome each record is in
in_boreal <- st_intersects(traits_sf, boreal_sf, sparse = FALSE)[, 1]
in_tundra <- st_intersects(traits_sf, tundra_sf, sparse = FALSE)[, 1]

# Add the biome
traits_with_coords$biome <- NA_character_
traits_with_coords$biome[in_boreal] <- "boreal"
traits_with_coords$biome[in_tundra] <- "tundra"

# Distance to the treeline (NEGATIVE = boreal, POSITIVE = tundra)
traits_with_coords$distance_to_boundary_km <- NA_real_
traits_with_coords$distance_to_boundary_km[in_boreal] <-
  -as.numeric(st_distance(traits_sf[in_boreal, ], shared_edge)) / 1000
traits_with_coords$distance_to_boundary_km[in_tundra] <-
  as.numeric(st_distance(traits_sf[in_tundra, ], shared_edge)) / 1000

# Projected coordinates in km (EPSG:3574), used to calculate distances between
# sites for the spatial autocorrelation check
coords_km <- st_coordinates(traits_sf) / 1000
traits_with_coords$x_km <- coords_km[, "X"]
traits_with_coords$y_km <- coords_km[, "Y"]

# Site identifier (records at the same coordinates)
traits_with_coords <- traits_with_coords |>
  mutate(site_name = paste0(round(LON_site, 3), "_", round(LAT_site, 3)))

## 2.2. Compare shared edge to biome edge --------------------------------------

# Distance to each biome's OWN edge includes coastlines and the southern boreal
# edge, which puts e.g. coastal Arctic sites "near the boundary" although they
# are far from any treeline.
traits_with_coords$distance_biome_edge_km <- NA_real_
traits_with_coords$distance_biome_edge_km[in_boreal] <-
  -as.numeric(st_distance(traits_sf[in_boreal, ], boreal_boundary)) / 1000
traits_with_coords$distance_biome_edge_km[in_tundra] <-
  as.numeric(st_distance(traits_sf[in_tundra, ], tundra_boundary)) / 1000

# Use one row per site with both distances
site_distances <- traits_with_coords |>
  filter(!is.na(biome)) |>
  distinct(site_name, biome, distance_biome_edge_km, distance_to_boundary_km)

# Check the number of sites that the alternative distance would wrongly place near the boundary
cat("Sites < 50 km from their own biome edge but > 200 km from the treeline:",
    sum(abs(site_distances$distance_biome_edge_km) < 50 &
          abs(site_distances$distance_to_boundary_km) > 200),
    "of", nrow(site_distances), "\n")

# Define a theme used by all figures in this script
theme_rq3 <- theme_classic() +
  theme(axis.text = element_text(size = fig_axis_text),
        axis.title = element_text(size = fig_axis_title, face = "bold"),
        legend.text = element_text(size = fig_legend_text),
        legend.title = element_text(size = fig_legend_title, face = "bold"),
        plot.title = element_text(size = fig_annot_text, face = "bold"))

# Create a supplementary figure (one point per site)
plot_distance_comparison <- ggplot(site_distances,
                                   aes(x = distance_biome_edge_km, y = distance_to_boundary_km,
                                       colour = biome)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(alpha = 0.6, size = 2) +
  scale_colour_manual(values = biome_colours, labels = biome_labels, name = "Biome") +
  labs(x = "Distance to own biome edge (km)", y = "Distance to treeline (km)") +
  theme_rq3

# Save figure to file
ggsave(here("figures", "FigureS9_RQ3_distance_definitions_comparison.png"),
       plot = plot_distance_comparison, width = 8, height = 7, dpi = 300)

## 2.3. CAFF quality control ---------------------------------------------------

# Keep only species with the CAFF classification that is not "remove",
# and only records that lie in the boreal or tundra biome
caff_check <- caff_check |>
  rename(StandardSpeciesName = SPECIES_CLEAN) |>
  dplyr::select(StandardSpeciesName, final.category)

# Remove records with NA for biome, final category, and distance to boundary
cleaned_traits_final <- traits_with_coords |>
  left_join(caff_check, by = "StandardSpeciesName") |>
  filter(!is.na(final.category), final.category != "remove") |>
  filter(!is.na(biome), !is.na(distance_to_boundary_km))

# Check how many species are left after filtering
cat("Species after CAFF filtering:", n_distinct(cleaned_traits_final$StandardSpeciesName),
    "| records:", nrow(cleaned_traits_final), "\n")
# Species after CAFF filtering: 643 | records: 57964

## 2.4. Use standardised trait values (StdValue) ------------------------------

# N.B! All analyses use StdValue.

# List the traits we are focusing on
key_traits <- c("PlantHeight", "SLA", "LeafN", "SeedMass")

# Convert StdValue to numeric
cleaned_traits_final <- cleaned_traits_final |>
  mutate(StdValue = as.numeric(StdValue))

# Check: one standardised unit per trait, and no missing standardised values
print(cleaned_traits_final |>
        filter(TraitNameNew %in% key_traits) |>
        group_by(TraitNameNew, UnitName) |>
        summarise(n_records = n(),
                  n_original_units = n_distinct(OrigUnitStr),
                  n_missing_StdValue = sum(is.na(StdValue)),
                  .groups = "drop"))
# TraitNameNew UnitName n_records n_original_units n_missing_StdValue
# <chr>        <chr>        <int>            <int>              <int>
# 1 LeafN        mg/g          7239                5                  0
# 2 PlantHeight  m            29601                4                  0
# 3 SLA          mm2 mg-1     18023               12                  0
# 4 SeedMass     mg             814                1                  0

## 2.5. Group nearby plots into locations --------------------------------------

# Group plots into locations to account for the fact that their proximity means they likely have similar
# climate (and other abiotic conditions)
cleaned_traits_final <- cleaned_traits_final |>
  mutate(location = paste0(round(LON_site / location_grid_deg), "_",
                           round(LAT_site / location_grid_deg)))

# 3. MODEL FORMULA -------------------------------------------------------------

# A straight line (lm, no random effects) to which segmented() adds one
# breakpoint in distance to the treeline
lm_formula <- log_trait ~ distance_to_boundary_km

# Distance classes of the correlogram (step .3 of each trait)
correlogram_breaks <- seq(0, correlogram_max_km, by = correlogram_step_km)

# 4. PLANT HEIGHT ------------------------------------------------------------

## 4.1. Prepare the data -------------------------------------------------------

# Keep only records that can be logged (> 0)
plant_height_data <- cleaned_traits_final |>
  filter(TraitNameNew == "PlantHeight", !is.na(StdValue), StdValue > 0)

# Outliers per biome: values more than 5 SD from the biome mean on the LOG
# scale (the scale used in the models)
plant_height_data <- plant_height_data |>
  group_by(biome) |>
  mutate(lower_bound = exp(mean(log(StdValue)) - 5 * sd(log(StdValue))),
         upper_bound = exp(mean(log(StdValue)) + 5 * sd(log(StdValue)))) |>
  ungroup()

# Check how many outliers are removed
cat("Plant height outliers removed:", sum(plant_height_data$StdValue < plant_height_data$lower_bound |
                                            plant_height_data$StdValue > plant_height_data$upper_bound), "\n")

# Remove the outliers, log-transform the trait, and make the grouping
# variables factors (needed for the random effects in the GAM script)
plant_height_data <- plant_height_data |>
  filter(StdValue >= lower_bound, StdValue <= upper_bound) |>
  dplyr::select(-lower_bound, -upper_bound) |>
  mutate(log_trait = log(StdValue),
         species_f = factor(StandardSpeciesName),
         location_f = factor(location),
         site_f = factor(site_name))

# Check how much data there is on each side of the treeline
print(plant_height_data |>
        group_by(biome) |>
        summarise(n_records = n(),
                  n_species = n_distinct(StandardSpeciesName),
                  n_sites = n_distinct(site_name),
                  n_locations = n_distinct(location)))

# biome  n_records n_species n_sites n_locations
# <chr>      <int>     <int>   <int>       <int>
# 1 boreal     14072       302      85          50
# 2 tundra     15529       243     215          99

## 4.2. Segmented regression without random effects ----------------------------
# A straight line (lm) is fitted first; segmented() then adds one breakpoint
# and estimates where it lies, starting the search at seg_start_psi (the treeline) and
# restarting seg_restarts times from other values. The model treats every
# record as independent (no random effects), so its CIs and p-values are too
# narrow if the residuals are spatially autocorrelated (checked in step .3).

# N.B. We are using an lm because of convergence issues when using species as a random effect (see email sent)

# Straight line
lm_ph <- lm(lm_formula, data = plant_height_data)

# Breakpoint model
seg_ph <- segmented(lm_ph, seg.Z = ~ distance_to_boundary_km, psi = seg_start_psi,
                    control = seg.control(n.boot = seg_restarts))
print(summary(seg_ph))
# ***Regression Model with Segmented Relationship(s)***
#   
#   Call: 
#   segmented.lm(obj = lm_ph, seg.Z = ~distance_to_boundary_km, psi = seg_start_psi, 
#                control = seg.control(n.boot = seg_restarts))
# 
# Estimated Break-Point(s):
#   Est. St.Err
# psi1.distance_to_boundary_km -44.892  5.049
# 
# Coefficients of the linear terms:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                -2.527e+00  1.578e-02 -160.11   <2e-16 ***
#   distance_to_boundary_km    -3.747e-03  8.484e-05  -44.17   <2e-16 ***
#   U1.distance_to_boundary_km  3.324e-03  8.970e-05   37.05       NA    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.153 on 29597 degrees of freedom
# Multiple R-Squared: 0.0891,  Adjusted R-squared: 0.08901 
# 
# Boot restarting based on 6 samples. Last fit:
#   Convergence attained in 2 iterations (rel. change 3.0545e-12)


# Note on the summary: the coefficient 'distance_to_boundary_km' is the slope
# BEFORE the breakpoint; 'U1.distance_to_boundary_km' is the CHANGE in slope at
# the breakpoint (slope after = the two added up; slope() below gives both
# slopes directly). Need to use the Davies test to check the p value

# Davies test: is there a change in slope anywhere along the gradient?
# (tests the straight line against a breakpoint model, allowing for the fact
# that the breakpoint location is estimated)
davies_ph <- davies.test(lm_ph, seg.Z = ~ distance_to_boundary_km)
print(davies_ph)
# Davies' test for a change in the slope
# 
# data:  formula = log_trait ~ distance_to_boundary_km ,   method = lm 
# model = gaussian , link = identity  
# segmented variable = distance_to_boundary_km
# 'best' at = -149.22, n.points = 8, p-value < 2.2e-16
# alternative hypothesis: two.sided

# Breakpoint and its 95% CI (first column = estimate, then lower and upper limit)
breakpoint_ci_ph <- confint(seg_ph)
breakpoint_ph <- breakpoint_ci_ph[1, 1]
breakpoint_lower_ph <- breakpoint_ci_ph[1, 2]
breakpoint_upper_ph <- breakpoint_ci_ph[1, 3]

# Slopes before (row 1) and after (row 2) the breakpoint, with 95% CIs
# (columns: estimate, SE, t value, lower and upper CI limit)
slopes_ph <- slope(seg_ph)$distance_to_boundary_km
print(slopes_ph)
# Est.    St.Err. t value   CI(95%).l  CI(95%).u
# slope1 -0.00374720 8.4837e-05 -44.170 -0.00391350 -0.0035809
# slope2 -0.00042361 2.9137e-05 -14.539 -0.00048073 -0.0003665

# Delta AIC: positive = the breakpoint model fits better than the straight line
daic_seg_ph <- AIC(lm_ph) - AIC(seg_ph)

cat("\nPlant height - breakpoint:", round(breakpoint_ph), "km (95% CI",
    round(breakpoint_lower_ph), "to", round(breakpoint_upper_ph), "km)",
    ifelse(breakpoint_lower_ph <= 0 & breakpoint_upper_ph >= 0,
           "- CI includes the treeline", "- CI does not include the treeline"), "\n",
    "  Davies test p =", signif(davies_ph$p.value, 3),
    "| delta AIC vs straight line =", round(daic_seg_ph, 1), "\n")
# Plant height - breakpoint: -45 km (95% CI -55 to -35 km) - CI does not include the treeline 
# Davies test p = 1.01e-292 | delta AIC vs straight line = 1446.9 

## 4.3. Spatial autocorrelation of the residuals --------------------------------
# Are the residuals of nearby sites more similar than those of distant sites?
# If so, the records are not independent and the model's CIs and p-values
# are too narrow. Residuals are averaged per site (sites with many records
# would otherwise dominate), then:
#   - Moran's I over all sites, with weights = 1 / distance between sites
#     (I > 0 = nearby sites are similar; p < 0.05 = significant autocorrelation)
#   - a correlogram: Moran's I among pairs of sites within each distance class
#     (distance BETWEEN sites, not distance to the treeline), showing up to
#     what distance sites remain similar

# Extract residuals from the segmented regression
plant_height_data$resid_seg <- residuals(seg_ph)

# Mean residual per site, with the site's projected coordinates
sites_ph <- plant_height_data |>
  group_by(site_name) |>
  summarise(x_km = first(x_km),
            y_km = first(y_km),
            resid_seg = mean(resid_seg),
            n_records = n(),
            .groups = "drop")

# Distances between all pairs of sites (km)
site_dist_ph <- as.matrix(dist(cbind(sites_ph$x_km, sites_ph$y_km)))

# Moran's I over all sites (inverse-distance weights)
weights_ph <- 1 / site_dist_ph
weights_ph[!is.finite(weights_ph)] <- 0   # a site's distance to itself
moran_seg_ph <- Moran.I(sites_ph$resid_seg, weights_ph)

cat("\nPlant height - Moran's I of the segmented residuals (", nrow(sites_ph), " sites): I =",
    round(moran_seg_ph$observed, 3), ", p =", signif(moran_seg_ph$p.value, 3), "\n")
#Plant height - Moran's I of the segmented residuals ( 300  sites): I = 0.572 , p = 0

# Correlogram: Moran's I among pairs of sites within each distance class
correlogram_ph <- data.frame()
for (i in seq_len(length(correlogram_breaks) - 1)) {
  
  # Weight 1 for pairs of sites within this distance class, 0 otherwise
  class_weights <- (site_dist_ph > correlogram_breaks[i] &
                      site_dist_ph <= correlogram_breaks[i + 1]) * 1
  
  # Skip classes without any pairs of sites
  if (sum(class_weights) == 0) next
  
  moran_class <- Moran.I(sites_ph$resid_seg, class_weights)
  
  correlogram_ph <- bind_rows(correlogram_ph, data.frame(
    Trait = "PlantHeight",
    distance_km = (correlogram_breaks[i] + correlogram_breaks[i + 1]) / 2,
    n_pairs = sum(class_weights) / 2,
    Model = segmented_label,
    morans_I = moran_class$observed,
    p = moran_class$p.value))
}

## 4.4. Model validation ---------------------------------------------------------

# Residuals vs fitted values (no pattern and an even spread around 0) and a
# Q-Q plot (points close to the line = roughly normal residuals)
# N.B. some misfit always shows up in segmented regressions because of the nature of the model itself 
# that being said, the QQ plot doesnt look as bad as expected

png(here("figures", "RQ3_PlantHeight_validation_segmented.png"),
    width = 12, height = 5, units = "in", res = 200)
par(mfrow = c(1, 2))

plot(fitted(seg_ph), plant_height_data$resid_seg, pch = 16, cex = 0.4, col = rgb(0, 0, 0, 0.2),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Plant height - Segmented: residuals vs fitted")
abline(h = 0, col = "red", lty = 2)

qqnorm(plant_height_data$resid_seg, pch = 16, cex = 0.4, main = "Plant height - Segmented: Q-Q plot")
qqline(plant_height_data$resid_seg, col = "red")

par(mfrow = c(1, 1))
dev.off()

## 4.5. Figure -------------------------------------------------------------------

# All records with the fitted segmented line, the estimated breakpoint (solid
# line) and its 95% CI (grey band), and the treeline (dashed line)

# Combine into sigle df
seg_line_ph <- data.frame(distance_to_boundary_km = seq(min(plant_height_data$distance_to_boundary_km),
                                                        max(plant_height_data$distance_to_boundary_km),
                                                        length.out = 300))
# Predict values
seg_line_ph$predicted <- predict(seg_ph, newdata = seg_line_ph)

# Slope labels as % change per km
label_before_ph <- paste0("Before breakpoint: ",
                          sprintf("%+.3f", 100 * (exp(slopes_ph[1, 1]) - 1)), "% per km")
label_after_ph <- paste0("After breakpoint: ",
                         sprintf("%+.3f", 100 * (exp(slopes_ph[2, 1]) - 1)), "% per km")
title_seg_ph <- paste0("Breakpoint: ", round(breakpoint_ph), " km (95% CI ",
                       round(breakpoint_lower_ph), " to ", round(breakpoint_upper_ph),
                       "); Davies p ", ifelse(davies_ph$p.value < 0.001, "< 0.001",
                                              paste0("= ", sprintf("%.3f", davies_ph$p.value))),
                       ", \u0394AIC = ", sprintf("%.1f", daic_seg_ph))

# Plot figure
(plot_seg_ph <- ggplot() +
  annotate("rect", xmin = breakpoint_lower_ph, xmax = breakpoint_upper_ph,
           ymin = -Inf, ymax = Inf, fill = "grey70", alpha = 0.3) +
  geom_point(data = plant_height_data, aes(x = distance_to_boundary_km, y = log_trait, colour = biome),
             alpha = 0.2, size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray30", linewidth = 0.7) +
  geom_vline(xintercept = breakpoint_ph, colour = "black", linewidth = 0.7) +
  geom_line(data = seg_line_ph, aes(x = distance_to_boundary_km, y = predicted),
            colour = fit_line_colour, linewidth = 1.5) +
  annotate("text", x = -Inf, y = Inf, label = label_before_ph,
           hjust = -0.03, vjust = 1.2, size = fig_annot_text / .pt, fontface = "bold") +
  annotate("text", x = Inf, y = Inf, label = label_after_ph,
           hjust = 1.03, vjust = 1.2, size = fig_annot_text / .pt, fontface = "bold") +
  scale_colour_manual(values = biome_colours, labels = biome_labels, name = "Biome") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +   # room for the labels
  labs(x = "Distance to Biome Boundary (km)", y = "Plant Height (m, log scale)", title = title_seg_ph) +
  theme_rq3 +
  theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))))


# Save figure to file 
ggsave(here("figures", "Figure4a_PlantHeight_segmented.png"),
       plot = plot_seg_ph, width = 10, height = 6, dpi = 600)

# Collect the results of the segmented regression into single df to save later
results_seg_ph <- data.frame(Trait = "PlantHeight",
                             n_records = nrow(plant_height_data),
                             n_species = n_distinct(plant_height_data$StandardSpeciesName),
                             n_sites = n_distinct(plant_height_data$site_name),
                             n_locations = n_distinct(plant_height_data$location),
                             breakpoint_km = breakpoint_ph,
                             breakpoint_lower = breakpoint_lower_ph,
                             breakpoint_upper = breakpoint_upper_ph,
                             ci_includes_treeline = breakpoint_lower_ph <= 0 & breakpoint_upper_ph >= 0,
                             davies_p = davies_ph$p.value,
                             dAIC_vs_straight_line = daic_seg_ph,
                             slope_before = slopes_ph[1, 1],
                             slope_before_lower = slopes_ph[1, 4],
                             slope_before_upper = slopes_ph[1, 5],
                             slope_after = slopes_ph[2, 1],
                             slope_after_lower = slopes_ph[2, 4],
                             slope_after_upper = slopes_ph[2, 5],
                             n_sites_autocorrelation = nrow(sites_ph),
                             morans_I = moran_seg_ph$observed,
                             morans_I_p = moran_seg_ph$p.value)

# 5. SLA ------------------------------------------------------------

## 5.1. Prepare the data -------------------------------------------------------

# Keep only records that can be logged (> 0)
sla_data <- cleaned_traits_final |>
  filter(TraitNameNew == "SLA", !is.na(StdValue), StdValue > 0)

# Outliers per biome: values more than 5 SD from the biome mean on the LOG
# scale (the scale used in the models).
# SLA: the lower bound is at least 0.1 mm2/mg, because values below this are
# implausible for vascular plants
sla_data <- sla_data |>
  group_by(biome) |>
  mutate(lower_bound = pmax(0.1, exp(mean(log(StdValue)) - 5 * sd(log(StdValue)))),
         upper_bound = exp(mean(log(StdValue)) + 5 * sd(log(StdValue)))) |>
  ungroup()

# Check how many outliers are removed
cat("SLA outliers removed:", sum(sla_data$StdValue < sla_data$lower_bound |
                                   sla_data$StdValue > sla_data$upper_bound), "\n") #6

# Remove the outliers, log-transform the trait, and make the grouping
# variables factors (needed for the random effects in the GAM script)
sla_data <- sla_data |>
  filter(StdValue >= lower_bound, StdValue <= upper_bound) |>
  dplyr::select(-lower_bound, -upper_bound) |>
  mutate(log_trait = log(StdValue),
         species_f = factor(StandardSpeciesName),
         location_f = factor(location),
         site_f = factor(site_name))

# Check how much data there is on each side of the treeline
print(sla_data |>
        group_by(biome) |>
        summarise(n_records = n(),
                  n_species = n_distinct(StandardSpeciesName),
                  n_sites = n_distinct(site_name),
                  n_locations = n_distinct(location)))

# biome  n_records n_species n_sites n_locations
# <chr>      <int>     <int>   <int>       <int>
# 1 boreal      2769       146      55          38
# 2 tundra     15248       295     133          67

## 5.2. Segmented regression without random effects ----------------------------

# Fit the lm
lm_sla <- lm(lm_formula, data = sla_data)

# Fit the breakpoint model
seg_sla <- segmented(lm_sla, seg.Z = ~ distance_to_boundary_km, psi = seg_start_psi,
                     control = seg.control(n.boot = seg_restarts))

# Quick check of the segmented regression summary
print(summary(seg_sla))
# Call: 
#   segmented.lm(obj = lm_sla, seg.Z = ~distance_to_boundary_km, 
#                psi = seg_start_psi, control = seg.control(n.boot = seg_restarts))
# 
# Estimated Break-Point(s):
#   Est. St.Err
# psi1.distance_to_boundary_km 67.696 16.976
# 
# Coefficients of the linear terms:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                 2.525e+00  5.717e-03  441.71   <2e-16 ***
#   distance_to_boundary_km     4.987e-04  2.500e-05   19.94   <2e-16 ***
#   U1.distance_to_boundary_km -7.643e-04  3.084e-05  -24.79       NA    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.6014 on 18013 degrees of freedom
# Multiple R-Squared: 0.03406,  Adjusted R-squared: 0.0339 
# 
# Boot restarting based on 11 samples. Last fit:
#   Convergence attained in 1 iterations (rel. change 3.5167e-06)

# Davies test: is there a change in slope anywhere along the gradient?
# (tests the straight line against a breakpoint model, allowing for the fact
# that the breakpoint location is estimated)
davies_sla <- davies.test(lm_sla, seg.Z = ~ distance_to_boundary_km)
print(davies_sla)
# data:  formula = log_trait ~ distance_to_boundary_km ,   method = lm 
# model = gaussian , link = identity  
# segmented variable = distance_to_boundary_km
# 'best' at = 225.19, n.points = 8, p-value < 2.2e-16
# alternative hypothesis: two.sided

# Breakpoint and its 95% CI (first column = estimate, then lower and upper limit)
breakpoint_ci_sla <- confint(seg_sla)
breakpoint_sla <- breakpoint_ci_sla[1, 1]
breakpoint_lower_sla <- breakpoint_ci_sla[1, 2]
breakpoint_upper_sla <- breakpoint_ci_sla[1, 3]

# Slopes before (row 1) and after (row 2) the breakpoint, with 95% CIs
# (columns: estimate, SE, t value, lower and upper CI limit)
slopes_sla <- slope(seg_sla)$distance_to_boundary_km
print(slopes_sla)
# Est.    St.Err. t value   CI(95%).l   CI(95%).u
# slope1  0.00049867 2.5005e-05  19.943  0.00044965  0.00054768
# slope2 -0.00026562 1.8046e-05 -14.719 -0.00030099 -0.00023025

# Delta AIC: positive = the breakpoint model fits better than the straight line
daic_seg_sla <- AIC(lm_sla) - AIC(seg_sla)

cat("\nSLA - breakpoint:", round(breakpoint_sla), "km (95% CI",
    round(breakpoint_lower_sla), "to", round(breakpoint_upper_sla), "km)",
    ifelse(breakpoint_lower_sla <= 0 & breakpoint_upper_sla >= 0,
           "- CI includes the treeline", "- CI does not include the treeline"), "\n",
    "  Davies test p =", signif(davies_sla$p.value, 3),
    "| delta AIC vs straight line =", round(daic_seg_sla, 1), "\n")
# SLA - breakpoint: 68 km (95% CI 34 to 101 km) - CI does not include the treeline 
# Davies test p = 3.71e-130 | delta AIC vs straight line = 610.4

## 5.3. Spatial autocorrelation of the residuals --------------------------------

# Extract residuals
sla_data$resid_seg <- residuals(seg_sla)

# Get mean residual per site, with the site's projected coordinates
sites_sla <- sla_data |>
  group_by(site_name) |>
  summarise(x_km = first(x_km),
            y_km = first(y_km),
            resid_seg = mean(resid_seg),
            n_records = n(),
            .groups = "drop")

# Distances between all pairs of sites (km)
site_dist_sla <- as.matrix(dist(cbind(sites_sla$x_km, sites_sla$y_km)))

# Moran's I over all sites (inverse-distance weights)
weights_sla <- 1 / site_dist_sla
weights_sla[!is.finite(weights_sla)] <- 0   # a site's distance to itself
moran_seg_sla <- Moran.I(sites_sla$resid_seg, weights_sla)

cat("\nSLA - Moran's I of the segmented residuals (", nrow(sites_sla), " sites): I =",
    round(moran_seg_sla$observed, 3), ", p =", signif(moran_seg_sla$p.value, 3), "\n") 
# SLA - Moran's I of the segmented residuals ( 188  sites): I = 0.276 , p = 3.32e-11

# Correlogram: Moran's I among pairs of sites within each distance class
correlogram_sla <- data.frame()
for (i in seq_len(length(correlogram_breaks) - 1)) {
  
  # Weight 1 for pairs of sites within this distance class, 0 otherwise
  class_weights <- (site_dist_sla > correlogram_breaks[i] &
                      site_dist_sla <= correlogram_breaks[i + 1]) * 1
  
  # Skip classes without any pairs of sites
  if (sum(class_weights) == 0) next
  
  moran_class <- Moran.I(sites_sla$resid_seg, class_weights)
  
  correlogram_sla <- bind_rows(correlogram_sla, data.frame(
    Trait = "SLA",
    distance_km = (correlogram_breaks[i] + correlogram_breaks[i + 1]) / 2,
    n_pairs = sum(class_weights) / 2,
    Model = segmented_label,
    morans_I = moran_class$observed,
    p = moran_class$p.value))
}

## 5.4. Model validation ---------------------------------------------------------

# Residuals vs fitted values (no pattern and an even spread around 0) and a
# Q-Q plot (points close to the line = roughly normal residuals)
png(here("figures", "RQ3_SLA_validation_segmented.png"),
    width = 12, height = 5, units = "in", res = 200)
par(mfrow = c(1, 2))

plot(fitted(seg_sla), sla_data$resid_seg, pch = 16, cex = 0.4, col = rgb(0, 0, 0, 0.2),
     xlab = "Fitted values", ylab = "Residuals",
     main = "SLA - Segmented: residuals vs fitted")
abline(h = 0, col = "red", lty = 2)

qqnorm(sla_data$resid_seg, pch = 16, cex = 0.4, main = "SLA - Segmented: Q-Q plot")
qqline(sla_data$resid_seg, col = "red")

par(mfrow = c(1, 1))
dev.off()

## 5.5. Figure -------------------------------------------------------------------

# Combine everything into a single df
seg_line_sla <- data.frame(distance_to_boundary_km = seq(min(sla_data$distance_to_boundary_km),
                                                         max(sla_data$distance_to_boundary_km),
                                                         length.out = 300))
seg_line_sla$predicted <- predict(seg_sla, newdata = seg_line_sla)

# Slope labels as % change per km
label_before_sla <- paste0("Before breakpoint: ",
                           sprintf("%+.3f", 100 * (exp(slopes_sla[1, 1]) - 1)), "% per km")
label_after_sla <- paste0("After breakpoint: ",
                          sprintf("%+.3f", 100 * (exp(slopes_sla[2, 1]) - 1)), "% per km")
title_seg_sla <- paste0("Breakpoint: ", round(breakpoint_sla), " km (95% CI ",
                        round(breakpoint_lower_sla), " to ", round(breakpoint_upper_sla),
                        "); Davies p ", ifelse(davies_sla$p.value < 0.001, "< 0.001",
                                               paste0("= ", sprintf("%.3f", davies_sla$p.value))),
                        ", \u0394AIC = ", sprintf("%.1f", daic_seg_sla))

# Plot figure
(plot_seg_sla <- ggplot() +
  annotate("rect", xmin = breakpoint_lower_sla, xmax = breakpoint_upper_sla,
           ymin = -Inf, ymax = Inf, fill = "grey70", alpha = 0.3) +
  geom_point(data = sla_data, aes(x = distance_to_boundary_km, y = log_trait, colour = biome),
             alpha = 0.2, size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray30", linewidth = 0.7) +
  geom_vline(xintercept = breakpoint_sla, colour = "black", linewidth = 0.7) +
  geom_line(data = seg_line_sla, aes(x = distance_to_boundary_km, y = predicted),
            colour = fit_line_colour, linewidth = 1.5) +
  annotate("text", x = -Inf, y = Inf, label = label_before_sla,
           hjust = -0.03, vjust = 1.2, size = fig_annot_text / .pt, fontface = "bold") +
  annotate("text", x = Inf, y = Inf, label = label_after_sla,
           hjust = 1.03, vjust = 1.2, size = fig_annot_text / .pt, fontface = "bold") +
  scale_colour_manual(values = biome_colours, labels = biome_labels, name = "Biome") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +   # room for the labels
  labs(x = "Distance to Biome Boundary (km)", y = expression(bold(paste("SLA (mm"^2, " mg"^-1, ", log scale)"))), title = title_seg_sla) +
  theme_rq3 +
  theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))))

# Save figure to file
ggsave(here("figures", "Figure4b_SLA_segmented.png"),
       plot = plot_seg_sla, width = 10, height = 6, dpi = 600)

# Collect results into signle df
results_seg_sla <- data.frame(Trait = "SLA",
                              n_records = nrow(sla_data),
                              n_species = n_distinct(sla_data$StandardSpeciesName),
                              n_sites = n_distinct(sla_data$site_name),
                              n_locations = n_distinct(sla_data$location),
                              breakpoint_km = breakpoint_sla,
                              breakpoint_lower = breakpoint_lower_sla,
                              breakpoint_upper = breakpoint_upper_sla,
                              ci_includes_treeline = breakpoint_lower_sla <= 0 & breakpoint_upper_sla >= 0,
                              davies_p = davies_sla$p.value,
                              dAIC_vs_straight_line = daic_seg_sla,
                              slope_before = slopes_sla[1, 1],
                              slope_before_lower = slopes_sla[1, 4],
                              slope_before_upper = slopes_sla[1, 5],
                              slope_after = slopes_sla[2, 1],
                              slope_after_lower = slopes_sla[2, 4],
                              slope_after_upper = slopes_sla[2, 5],
                              n_sites_autocorrelation = nrow(sites_sla),
                              morans_I = moran_seg_sla$observed,
                              morans_I_p = moran_seg_sla$p.value)

# 6. LEAF N ------------------------------------------------------------

## 6.1. Prepare the data -------------------------------------------------------

# Keep only records that can be logged (> 0)
leafn_data <- cleaned_traits_final |>
  filter(TraitNameNew == "LeafN", !is.na(StdValue), StdValue > 0)

# Outliers per biome: values more than 5 SD from the biome mean on the LOG
# scale (the scale used in the models)
leafn_data <- leafn_data |>
  group_by(biome) |>
  mutate(lower_bound = exp(mean(log(StdValue)) - 5 * sd(log(StdValue))),
         upper_bound = exp(mean(log(StdValue)) + 5 * sd(log(StdValue)))) |>
  ungroup()

# Check how many outliers are removed
cat("Leaf N outliers removed:", sum(leafn_data$StdValue < leafn_data$lower_bound |
                                      leafn_data$StdValue > leafn_data$upper_bound), "\n") #12

# Remove the outliers
leafn_data <- leafn_data |>
  filter(StdValue >= lower_bound, StdValue <= upper_bound) |>
  dplyr::select(-lower_bound, -upper_bound) |>
  mutate(log_trait = log(StdValue),
         species_f = factor(StandardSpeciesName),
         location_f = factor(location),
         site_f = factor(site_name))

# Check how much data there is on each side of the treeline
print(leafn_data |>
        group_by(biome) |>
        summarise(n_records = n(),
                  n_species = n_distinct(StandardSpeciesName),
                  n_sites = n_distinct(site_name),
                  n_locations = n_distinct(location)))
# biome  n_records n_species n_sites n_locations
# <chr>      <int>     <int>   <int>       <int>
# 1 boreal      1766       252     113          99
# 2 tundra      5441       240      97          60

## 6.2. Segmented regression without random effects ----------------------------

# Straight line
lm_leafn <- lm(lm_formula, data = leafn_data)

# Breakpoint model
seg_leafn <- segmented(lm_leafn, seg.Z = ~ distance_to_boundary_km, psi = seg_start_psi,
                       control = seg.control(n.boot = seg_restarts))

# Check output of the segmented regression
print(summary(seg_leafn))
# Estimated Break-Point(s):
#   Est. St.Err
# psi1.distance_to_boundary_km 752.935 37.818
# 
# Coefficients of the linear terms:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                 2.928e+00  6.509e-03  449.76   <2e-16 ***
#   distance_to_boundary_km     4.013e-04  2.151e-05   18.66   <2e-16 ***
#   U1.distance_to_boundary_km -1.048e-03  7.106e-05  -14.75       NA    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.5162 on 7203 degrees of freedom
# Multiple R-Squared: 0.05749,  Adjusted R-squared: 0.0571 

# Davies test: is there a change in slope anywhere along the gradient?
davies_leafn <- davies.test(lm_leafn, seg.Z = ~ distance_to_boundary_km)
print(davies_leafn)
# data:  formula = log_trait ~ distance_to_boundary_km ,   method = lm 
# model = gaussian , link = identity  
# segmented variable = distance_to_boundary_km
# 'best' at = 601.61, n.points = 8, p-value < 2.2e-16
# alternative hypothesis: two.sided

# Breakpoint and its 95% CI (first column = estimate, then lower and upper limit)
breakpoint_ci_leafn <- confint(seg_leafn)
breakpoint_leafn <- breakpoint_ci_leafn[1, 1]
breakpoint_lower_leafn <- breakpoint_ci_leafn[1, 2]
breakpoint_upper_leafn <- breakpoint_ci_leafn[1, 3]

# Slopes before (row 1) and after (row 2) the breakpoint
slopes_leafn <- slope(seg_leafn)$distance_to_boundary_km
print(slopes_leafn)
# Est.    St.Err. t value   CI(95%).l   CI(95%).u
# slope1  0.00040131 2.1510e-05 18.6560  0.00035914  0.00044347
# slope2 -0.00064659 6.7726e-05 -9.5472 -0.00077935 -0.00051383

# Delta AIC: positive = the breakpoint model fits better than the straight line
daic_seg_leafn <- AIC(lm_leafn) - AIC(seg_leafn)

cat("\nLeaf N - breakpoint:", round(breakpoint_leafn), "km (95% CI",
    round(breakpoint_lower_leafn), "to", round(breakpoint_upper_leafn), "km)",
    ifelse(breakpoint_lower_leafn <= 0 & breakpoint_upper_leafn >= 0,
           "- CI includes the treeline", "- CI does not include the treeline"), "\n",
    "  Davies test p =", signif(davies_leafn$p.value, 3),
    "| delta AIC vs straight line =", round(daic_seg_leafn, 1), "\n")
# Leaf N - breakpoint: 753 km (95% CI 679 to 827 km) - CI does not include the treeline 
# Davies test p = 1.05e-75 | delta AIC vs straight line = 350.1 

## 6.3. Spatial autocorrelation of the residuals --------------------------------

# GEt residuals
leafn_data$resid_seg <- residuals(seg_leafn)

# Mean residual per site, with the site's projected coordinates
sites_leafn <- leafn_data |>
  group_by(site_name) |>
  summarise(x_km = first(x_km),
            y_km = first(y_km),
            resid_seg = mean(resid_seg),
            n_records = n(),
            .groups = "drop")

# Distances between all pairs of sites (km)
site_dist_leafn <- as.matrix(dist(cbind(sites_leafn$x_km, sites_leafn$y_km)))

# Moran's I over all sites (inverse-distance weights)
weights_leafn <- 1 / site_dist_leafn
weights_leafn[!is.finite(weights_leafn)] <- 0   # a site's distance to itself
moran_seg_leafn <- Moran.I(sites_leafn$resid_seg, weights_leafn)

cat("\nLeaf N - Moran's I of the segmented residuals (", nrow(sites_leafn), " sites): I =",
    round(moran_seg_leafn$observed, 3), ", p =", signif(moran_seg_leafn$p.value, 3), "\n")
# Leaf N - Moran's I of the segmented residuals ( 210  sites): I = 0.28 , p = 8.73e-14

# Correlogram: Moran's I among pairs of sites within each distance class
correlogram_leafn <- data.frame()
for (i in seq_len(length(correlogram_breaks) - 1)) {
  
  # Weight 1 for pairs of sites within this distance class, 0 otherwise
  class_weights <- (site_dist_leafn > correlogram_breaks[i] &
                      site_dist_leafn <= correlogram_breaks[i + 1]) * 1
  
  # Skip classes without any pairs of sites
  if (sum(class_weights) == 0) next
  
  moran_class <- Moran.I(sites_leafn$resid_seg, class_weights)
  
  correlogram_leafn <- bind_rows(correlogram_leafn, data.frame(
    Trait = "LeafN",
    distance_km = (correlogram_breaks[i] + correlogram_breaks[i + 1]) / 2,
    n_pairs = sum(class_weights) / 2,
    Model = segmented_label,
    morans_I = moran_class$observed,
    p = moran_class$p.value))
}

## 6.4. Model validation ---------------------------------------------------------

# Residuals vs fitted values and Q-Q plot 
png(here("figures", "RQ3_LeafN_validation_segmented.png"),
    width = 12, height = 5, units = "in", res = 200)
par(mfrow = c(1, 2))

plot(fitted(seg_leafn), leafn_data$resid_seg, pch = 16, cex = 0.4, col = rgb(0, 0, 0, 0.2),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Leaf N - Segmented: residuals vs fitted")
abline(h = 0, col = "red", lty = 2)

qqnorm(leafn_data$resid_seg, pch = 16, cex = 0.4, main = "Leaf N - Segmented: Q-Q plot")
qqline(leafn_data$resid_seg, col = "red")

par(mfrow = c(1, 1))
dev.off()

## 6.5. Figure -------------------------------------------------------------------

# Put results into single df
seg_line_leafn <- data.frame(distance_to_boundary_km = seq(min(leafn_data$distance_to_boundary_km),
                                                           max(leafn_data$distance_to_boundary_km),
                                                           length.out = 300))
# Predict values
seg_line_leafn$predicted <- predict(seg_leafn, newdata = seg_line_leafn)

# Slope labels as % change per km
label_before_leafn <- paste0("Before breakpoint: ",
                             sprintf("%+.3f", 100 * (exp(slopes_leafn[1, 1]) - 1)), "% per km")
label_after_leafn <- paste0("After breakpoint: ",
                            sprintf("%+.3f", 100 * (exp(slopes_leafn[2, 1]) - 1)), "% per km")
title_seg_leafn <- paste0("Breakpoint: ", round(breakpoint_leafn), " km (95% CI ",
                          round(breakpoint_lower_leafn), " to ", round(breakpoint_upper_leafn),
                          "); Davies p ", ifelse(davies_leafn$p.value < 0.001, "< 0.001",
                                                 paste0("= ", sprintf("%.3f", davies_leafn$p.value))),
                          ", \u0394AIC = ", sprintf("%.1f", daic_seg_leafn))

# Plot figure
(plot_seg_leafn <- ggplot() +
  annotate("rect", xmin = breakpoint_lower_leafn, xmax = breakpoint_upper_leafn,
           ymin = -Inf, ymax = Inf, fill = "grey70", alpha = 0.3) +
  geom_point(data = leafn_data, aes(x = distance_to_boundary_km, y = log_trait, colour = biome),
             alpha = 0.2, size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray30", linewidth = 0.7) +
  geom_vline(xintercept = breakpoint_leafn, colour = "black", linewidth = 0.7) +
  geom_line(data = seg_line_leafn, aes(x = distance_to_boundary_km, y = predicted),
            colour = fit_line_colour, linewidth = 1.5) +
  annotate("text", x = -Inf, y = Inf, label = label_before_leafn,
           hjust = -0.03, vjust = 1.2, size = fig_annot_text / .pt, fontface = "bold") +
  annotate("text", x = Inf, y = Inf, label = label_after_leafn,
           hjust = 1.03, vjust = 1.2, size = fig_annot_text / .pt, fontface = "bold") +
  scale_colour_manual(values = biome_colours, labels = biome_labels, name = "Biome") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +   # room for the labels
  labs(x = "Distance to Biome Boundary (km)", y = "Leaf N (mg/g, log scale)", title = title_seg_leafn) +
  theme_rq3 +
  theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))))

# Save figure
ggsave(here("figures", "Figure4c_LeafN_segmented.png"),
       plot = plot_seg_leafn, width = 10, height = 6, dpi = 600)

# Collect all results into df
results_seg_leafn <- data.frame(Trait = "LeafN",
                                n_records = nrow(leafn_data),
                                n_species = n_distinct(leafn_data$StandardSpeciesName),
                                n_sites = n_distinct(leafn_data$site_name),
                                n_locations = n_distinct(leafn_data$location),
                                breakpoint_km = breakpoint_leafn,
                                breakpoint_lower = breakpoint_lower_leafn,
                                breakpoint_upper = breakpoint_upper_leafn,
                                ci_includes_treeline = breakpoint_lower_leafn <= 0 & breakpoint_upper_leafn >= 0,
                                davies_p = davies_leafn$p.value,
                                dAIC_vs_straight_line = daic_seg_leafn,
                                slope_before = slopes_leafn[1, 1],
                                slope_before_lower = slopes_leafn[1, 4],
                                slope_before_upper = slopes_leafn[1, 5],
                                slope_after = slopes_leafn[2, 1],
                                slope_after_lower = slopes_leafn[2, 4],
                                slope_after_upper = slopes_leafn[2, 5],
                                n_sites_autocorrelation = nrow(sites_leafn),
                                morans_I = moran_seg_leafn$observed,
                                morans_I_p = moran_seg_leafn$p.value)

# 7. SEED MASS ------------------------------------------------------------
# N.B.: very few boreal seed mass records (3 locations in total), so results
# for the boreal side are very uncertain.

## 7.1. Prepare the data -------------------------------------------------------

# Keep only records that can be logged (> 0)
seedmass_data <- cleaned_traits_final |>
  filter(TraitNameNew == "SeedMass", !is.na(StdValue), StdValue > 0)

# Outliers per biome: values more than 5 SD from the biome mean on the LOG
# scale (the scale used in the models)
seedmass_data <- seedmass_data |>
  group_by(biome) |>
  mutate(lower_bound = exp(mean(log(StdValue)) - 5 * sd(log(StdValue))),
         upper_bound = exp(mean(log(StdValue)) + 5 * sd(log(StdValue)))) |>
  ungroup()

# Check how many outliers are removed
cat("Seed mass outliers removed:", sum(seedmass_data$StdValue < seedmass_data$lower_bound |
                                         seedmass_data$StdValue > seedmass_data$upper_bound), "\n") #0

# Remove the outliers, log-transform the trait, and make the grouping
# variables factors (needed for the random effects in the GAM script)
seedmass_data <- seedmass_data |>
  filter(StdValue >= lower_bound, StdValue <= upper_bound) |>
  dplyr::select(-lower_bound, -upper_bound) |>
  mutate(log_trait = log(StdValue),
         species_f = factor(StandardSpeciesName),
         location_f = factor(location),
         site_f = factor(site_name))

# Check how much data there is on each side of the treeline
print(seedmass_data |>
        group_by(biome) |>
        summarise(n_records = n(),
                  n_species = n_distinct(StandardSpeciesName),
                  n_sites = n_distinct(site_name),
                  n_locations = n_distinct(location)))
# biome  n_records n_species n_sites n_locations
# <chr>      <int>     <int>   <int>       <int>
# 1 boreal        23         9       3           3
# 2 tundra       791       117      33          17

## 7.2. Segmented regression without random effects ----------------------------

# Straight line
lm_seedmass <- lm(lm_formula, data = seedmass_data)

# Breakpoint model
seg_seedmass <- segmented(lm_seedmass, seg.Z = ~ distance_to_boundary_km, psi = median(seedmass_data$distance_to_boundary_km),
                          control = seg.control(n.boot = seg_restarts))

# Check the summary of the segmented model
print(summary(seg_seedmass))
# Estimated Break-Point(s):
#   Est. St.Err
# psi1.distance_to_boundary_km 97.597 34.507
# 
# Coefficients of the linear terms:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                -1.0930806  0.0535458 -20.414  < 2e-16 ***
#   distance_to_boundary_km     0.0062413  0.0008657   7.210 1.29e-12 ***
#   U1.distance_to_boundary_km -0.0063084  0.0009015  -6.998       NA    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.262 on 810 degrees of freedom
# Multiple R-Squared: 0.06659,  Adjusted R-squared: 0.06314 

# Davies test: is there a change in slope anywhere along the gradient?
davies_seedmass <- davies.test(lm_seedmass, seg.Z = ~ distance_to_boundary_km)
print(davies_seedmass)
# data:  formula = log_trait ~ distance_to_boundary_km ,   method = lm 
# model = gaussian , link = identity  
# segmented variable = distance_to_boundary_km
# 'best' at = 186.16, n.points = 8, p-value = 1.073e-08
# alternative hypothesis: two.sided

# Breakpoint and its 95% CI (first column = estimate, then lower and upper limit)
breakpoint_ci_seedmass <- confint(seg_seedmass)
breakpoint_seedmass <- breakpoint_ci_seedmass[1, 1]
breakpoint_lower_seedmass <- breakpoint_ci_seedmass[1, 2]
breakpoint_upper_seedmass <- breakpoint_ci_seedmass[1, 3]

# Slopes before (row 1) and after (row 2) the breakpoint, with 95% CIs
# (columns: estimate, SE, t value, lower and upper CI limit)
slopes_seedmass <- slope(seg_seedmass)$distance_to_boundary_km
print(slopes_seedmass)
# Est.    St.Err.  t value   CI(95%).l  CI(95%).u
# slope1  6.2413e-03 0.00086567  7.20980  0.00454210 0.00794060
# slope2 -6.7043e-05 0.00025146 -0.26661 -0.00056063 0.00042655

# Delta AIC: positive = the breakpoint model fits better than the straight line
daic_seg_seedmass <- AIC(lm_seedmass) - AIC(seg_seedmass)

cat("\nSeed mass - breakpoint:", round(breakpoint_seedmass), "km (95% CI",
    round(breakpoint_lower_seedmass), "to", round(breakpoint_upper_seedmass), "km)",
    ifelse(breakpoint_lower_seedmass <= 0 & breakpoint_upper_seedmass >= 0,
           "- CI includes the treeline", "- CI does not include the treeline"), "\n",
    "  Davies test p =", signif(davies_seedmass$p.value, 3),
    "| delta AIC vs straight line =", round(daic_seg_seedmass, 1), "\n")
# Seed mass - breakpoint: 98 km (95% CI 30 to 165 km) - CI does not include the treeline 
# Davies test p = 1.07e-08 | delta AIC vs straight line = 44.1 

## 7.3. Spatial autocorrelation of the residuals --------------------------------

# Get residuals
seedmass_data$resid_seg <- residuals(seg_seedmass)

# Mean residual per site, with the site's projected coordinates
sites_seedmass <- seedmass_data |>
  group_by(site_name) |>
  summarise(x_km = first(x_km),
            y_km = first(y_km),
            resid_seg = mean(resid_seg),
            n_records = n(),
            .groups = "drop")

# Distances between all pairs of sites (km)
site_dist_seedmass <- as.matrix(dist(cbind(sites_seedmass$x_km, sites_seedmass$y_km)))

# Moran's I over all sites (inverse-distance weights)
weights_seedmass <- 1 / site_dist_seedmass
weights_seedmass[!is.finite(weights_seedmass)] <- 0   # a site's distance to itself
moran_seg_seedmass <- Moran.I(sites_seedmass$resid_seg, weights_seedmass)

cat("\nSeed mass - Moran's I of the segmented residuals (", nrow(sites_seedmass), " sites): I =",
    round(moran_seg_seedmass$observed, 3), ", p =", signif(moran_seg_seedmass$p.value, 3), "\n")
#Seed mass - Moran's I of the segmented residuals ( 36  sites): I = 0.009 , p = 0.698

# Correlogram: Moran's I among pairs of sites within each distance class
correlogram_seedmass <- data.frame()
for (i in seq_len(length(correlogram_breaks) - 1)) {
  
  # Weight 1 for pairs of sites within this distance class, 0 otherwise
  class_weights <- (site_dist_seedmass > correlogram_breaks[i] &
                      site_dist_seedmass <= correlogram_breaks[i + 1]) * 1
  
  # Skip classes without any pairs of sites
  if (sum(class_weights) == 0) next
  
  moran_class <- Moran.I(sites_seedmass$resid_seg, class_weights)
  
  correlogram_seedmass <- bind_rows(correlogram_seedmass, data.frame(
    Trait = "SeedMass",
    distance_km = (correlogram_breaks[i] + correlogram_breaks[i + 1]) / 2,
    n_pairs = sum(class_weights) / 2,
    Model = segmented_label,
    morans_I = moran_class$observed,
    p = moran_class$p.value))
}
 
## 7.4. Model validation ---------------------------------------------------------

# Residuals vs fitted values and Q-Q plot
png(here("figures", "RQ3_SeedMass_validation_segmented.png"),
    width = 12, height = 5, units = "in", res = 200)
par(mfrow = c(1, 2))

plot(fitted(seg_seedmass), seedmass_data$resid_seg, pch = 16, cex = 0.4, col = rgb(0, 0, 0, 0.2),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Seed mass - Segmented: residuals vs fitted")
abline(h = 0, col = "red", lty = 2)

qqnorm(seedmass_data$resid_seg, pch = 16, cex = 0.4, main = "Seed mass - Segmented: Q-Q plot")
qqline(seedmass_data$resid_seg, col = "red")

par(mfrow = c(1, 1))
dev.off()

## 7.5. Figure -------------------------------------------------------------------

# All results in one df
seg_line_seedmass <- data.frame(distance_to_boundary_km = seq(min(seedmass_data$distance_to_boundary_km),
                                                              max(seedmass_data$distance_to_boundary_km),
                                                              length.out = 300))
# Predict values
seg_line_seedmass$predicted <- predict(seg_seedmass, newdata = seg_line_seedmass)

# Slope labels as % change per km
label_before_seedmass <- paste0("Before breakpoint: ",
                                sprintf("%+.3f", 100 * (exp(slopes_seedmass[1, 1]) - 1)), "% per km")
label_after_seedmass <- paste0("After breakpoint: ",
                               sprintf("%+.3f", 100 * (exp(slopes_seedmass[2, 1]) - 1)), "% per km")
title_seg_seedmass <- paste0("Breakpoint: ", round(breakpoint_seedmass), " km (95% CI ",
                             round(breakpoint_lower_seedmass), " to ", round(breakpoint_upper_seedmass),
                             "); Davies p ", ifelse(davies_seedmass$p.value < 0.001, "< 0.001",
                                                    paste0("= ", sprintf("%.3f", davies_seedmass$p.value))),
                             ", \u0394AIC = ", sprintf("%.1f", daic_seg_seedmass))

# Plot figure
(plot_seg_seedmass <- ggplot() +
  annotate("rect", xmin = breakpoint_lower_seedmass, xmax = breakpoint_upper_seedmass,
           ymin = -Inf, ymax = Inf, fill = "grey70", alpha = 0.3) +
  geom_point(data = seedmass_data, aes(x = distance_to_boundary_km, y = log_trait, colour = biome),
             alpha = 0.2, size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray30", linewidth = 0.7) +
  geom_vline(xintercept = breakpoint_seedmass, colour = "black", linewidth = 0.7) +
  geom_line(data = seg_line_seedmass, aes(x = distance_to_boundary_km, y = predicted),
            colour = fit_line_colour, linewidth = 1.5) +
  annotate("text", x = -Inf, y = Inf, label = label_before_seedmass,
           hjust = -0.03, vjust = 1.2, size = fig_annot_text / .pt, fontface = "bold") +
  annotate("text", x = Inf, y = Inf, label = label_after_seedmass,
           hjust = 1.03, vjust = 1.2, size = fig_annot_text / .pt, fontface = "bold") +
  scale_colour_manual(values = biome_colours, labels = biome_labels, name = "Biome") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +   # room for the labels
  labs(x = "Distance to Biome Boundary (km)", y = "Seed Mass (mg, log scale)", title = title_seg_seedmass) +
  theme_rq3 +
  theme(legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))))

# Save to file
ggsave(here("figures", "Figure4d_SeedMass_segmented.png"),
       plot = plot_seg_seedmass, width = 10, height = 6, dpi = 600)

# Collect all records into single df
results_seg_seedmass <- data.frame(Trait = "SeedMass",
                                   n_records = nrow(seedmass_data),
                                   n_species = n_distinct(seedmass_data$StandardSpeciesName),
                                   n_sites = n_distinct(seedmass_data$site_name),
                                   n_locations = n_distinct(seedmass_data$location),
                                   breakpoint_km = breakpoint_seedmass,
                                   breakpoint_lower = breakpoint_lower_seedmass,
                                   breakpoint_upper = breakpoint_upper_seedmass,
                                   ci_includes_treeline = breakpoint_lower_seedmass <= 0 & breakpoint_upper_seedmass >= 0,
                                   davies_p = davies_seedmass$p.value,
                                   dAIC_vs_straight_line = daic_seg_seedmass,
                                   slope_before = slopes_seedmass[1, 1],
                                   slope_before_lower = slopes_seedmass[1, 4],
                                   slope_before_upper = slopes_seedmass[1, 5],
                                   slope_after = slopes_seedmass[2, 1],
                                   slope_after_lower = slopes_seedmass[2, 4],
                                   slope_after_upper = slopes_seedmass[2, 5],
                                   n_sites_autocorrelation = nrow(sites_seedmass),
                                   morans_I = moran_seg_seedmass$observed,
                                   morans_I_p = moran_seg_seedmass$p.value)


# 8. RESULTS TABLES, FIGURES AND SAVED DATA ------------------------------------

## 8.1. Results table ----------------------------------------------------------

# Combine all the results dfs
segmented_results <- bind_rows(results_seg_ph, results_seg_sla,
                               results_seg_leafn, results_seg_seedmass) |>
  mutate(pct_before = 100 * (exp(slope_before) - 1),
         pct_before_lower = 100 * (exp(slope_before_lower) - 1),
         pct_before_upper = 100 * (exp(slope_before_upper) - 1),
         pct_after = 100 * (exp(slope_after) - 1),
         pct_after_lower = 100 * (exp(slope_after_lower) - 1),
         pct_after_upper = 100 * (exp(slope_after_upper) - 1))

# Save table
write.csv(segmented_results, here("data", "derived_data", "RQ3_segmented_results_full.csv"),
          row.names = FALSE)

# Tidy up the table to make it pretty for the Supplementary
segmented_table <- segmented_results |>
  transmute(Trait,
            `Records / species / sites / locations` =
              paste0(n_records, " / ", n_species, " / ", n_sites, " / ", n_locations),
            `Breakpoint, km (95% CI)` = paste0(round(breakpoint_km), " (", round(breakpoint_lower),
                                               " to ", round(breakpoint_upper), ")"),
            `CI includes treeline` = ifelse(ci_includes_treeline, "yes", "no"),
            `Davies test p` = format.pval(davies_p, digits = 3, eps = 1e-4),
            `ΔAIC vs straight line` = sprintf("%.1f", dAIC_vs_straight_line),
            `Before breakpoint: % per km (95% CI)` =
              paste0(sprintf("%.3f", pct_before), " (", sprintf("%.3f", pct_before_lower),
                     " to ", sprintf("%.3f", pct_before_upper), ")"),
            `After breakpoint: % per km (95% CI)` =
              paste0(sprintf("%.3f", pct_after), " (", sprintf("%.3f", pct_after_lower),
                     " to ", sprintf("%.3f", pct_after_upper), ")"),
            `Residuals: Moran's I (p)` =
              paste0(sprintf("%.3f", morans_I), " (",
                     format.pval(morans_I_p, digits = 3, eps = 1e-4), ")"))

# Quick check of the table
print(tibble::as_tibble(segmented_table), width = Inf)

# Save the nice table
write.csv(segmented_table, here("data", "derived_data", "TableS_RQ3_segmented.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# Caption notes:
# - Segmented regression (segmented::segmented) of log trait values on distance
#   to the treeline, without random effects. Breakpoint CI from the standard
#   error; Davies test of a change in slope anywhere along the gradient;
#   ΔAIC = AIC(straight line) - AIC(breakpoint model), positive = breakpoint
#   model better. Slopes as % change per km towards the tundra.
# - Moran's I of the mean residual per site (inverse-distance weights; I > 0 =
#   nearby sites have similar residuals). Significant spatial autocorrelation
#   means the records are not independent, so the CIs and p-values above are
#   too narrow.

## 8.2. Combined figures --------------------------------------------------------

# Combine all panels of the segmented regression into single figure
all_traits_segmented <- plot_grid(plot_seg_ph, plot_seg_sla, plot_seg_leafn, plot_seg_seedmass,
                                  labels = c("a)", "b)", "c)", "d)"),
                                  label_size = fig_panel_label, nrow = 2)

# Save as both .png and .pdf
ggsave(here("figures", "Figure4_traits_segmented.png"),
       plot = all_traits_segmented, width = 20, height = 12, dpi = 600)
ggsave(here("figures", "Figure4_traits_segmented.pdf"),
       plot = all_traits_segmented, width = 20, height = 12, dpi = 600)

## 8.3. Save the prepared data and the correlograms ------------------------------

# Used by 3.3_RQ3_GAM.R: the prepared trait data (without the segmented
# residuals) and the correlograms of the segmented residuals, so that the GAM
# script can compare them with its own residuals.
plant_height_data <- plant_height_data |> dplyr::select(-resid_seg)
sla_data <- sla_data |> dplyr::select(-resid_seg)
leafn_data <- leafn_data |> dplyr::select(-resid_seg)
seedmass_data <- seedmass_data |> dplyr::select(-resid_seg)

# Save into one file you can then quickly load
save(plant_height_data, sla_data, leafn_data, seedmass_data,
     file = here("data", "derived_data", "RQ3_trait_data.RData"))

# Save the correlograms 
write.csv(bind_rows(correlogram_ph, correlogram_sla, correlogram_leafn, correlogram_seedmass),
          here("data", "derived_data", "RQ3_correlograms_segmented.csv"), row.names = FALSE)

# END OF SCRIPT ----------------------------------------------------------------