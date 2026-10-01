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

# Minimum number of locations needed on each side of the treeline to fit the models
# (trying to keep as many as possible)
min_locations_per_side <- 2

# Breakpoint search: candidates between these quantiles of distance, first
# every coarse_step_km, then every fine_step_km around the best coarse value
profile_trim <- 0.05
coarse_step_km <- 20
fine_step_km <- 1
fine_window_km <- 30

# Figure settings (same as script 3.1)
biome_colours <- c("boreal" = "darkgreen", "tundra" = "darkblue")
biome_labels <- c("boreal" = "Boreal", "tundra" = "Tundra")
fit_line_colour <- "#dcd0ff"
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
library(lme4)
library(lmerTest) 

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

# Re-project biomes to North Pole Lamvert Azimuthal Equal Area (EPSG: 3574)
boreal_sf <- st_transform(st_make_valid(boreal_forest), "EPSG:3574")
tundra_sf <- st_transform(st_make_valid(tundra), "EPSG:3574")

# Get biome boundaries for distance calculations
boreal_boundary <- st_boundary(boreal_sf)
tundra_boundary <- st_boundary(tundra_sf)

# Extract "treeline", i.e. the parts where the boreal and tundra biomes are touch
shared_edge <- st_intersection(boreal_boundary, st_buffer(tundra_sf, 1000))
cat("  Length of treeline:", round(as.numeric(sum(st_length(shared_edge))) / 1000), "km\n")

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

# Initialize distance column
traits_with_coords$biome <- NA_character_
traits_with_coords$biome[in_boreal] <- "boreal"
traits_with_coords$biome[in_tundra] <- "tundra"

# Distance to the treeline (NEGATIVE = boreal, POSITIVE = tundra)
traits_with_coords$distance_to_boundary_km <- NA_real_
traits_with_coords$distance_to_boundary_km[in_boreal] <-
  -as.numeric(st_distance(traits_sf[in_boreal, ], shared_edge)) / 1000
traits_with_coords$distance_to_boundary_km[in_tundra] <-
  as.numeric(st_distance(traits_sf[in_tundra, ], shared_edge)) / 1000

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

## 2.4. Use standardised trait values (StdValue) ------------------------------

# All analyses use StdValue.

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

## 2.5. Group nearby plots into locations --------------------------------------

# Group plots into locations to account for the fact that their proximity means they likely have similar
# climate (and other abiotic conditions)
cleaned_traits_final <- cleaned_traits_final |>
  mutate(location = paste0(round(LON_site / location_grid_deg), "_",
                           round(LAT_site / location_grid_deg)))

# 3. MODEL FORMULAS ------------------------------------------------------------

# Random effects: species, location, and site (plot) within location.
# Distance is in km (distance_to_boundary_km): slopes are the change in the
# log trait value per km towards the tundra.
# In the breakpoint model, for a breakpoint at psi (km):
#   below = pmin(distance_to_boundary_km - psi, 0) -> slope BELOW the breakpoint (boreal side)
#   above = pmax(distance_to_boundary_km - psi, 0) -> slope ABOVE the breakpoint (tundra side)

# No breakpoint: one straight line
linear_formula <- log_trait ~ distance_to_boundary_km +
  (1 | StandardSpeciesName) + (1 | location) + (1 | site_name)

# Breakpoint: separate slopes on each side
breakpoint_formula <- log_trait ~ below + above +
  (1 | StandardSpeciesName) + (1 | location) + (1 | site_name)

# Robustness check: one side of the treeline on its own
side_formula <- log_trait ~ distance_to_boundary_km +
  (1 | StandardSpeciesName) + (1 | location) + (1 | site_name)

# 4. PLANT HEIGHT ------------------------------------------------------------

## 4.1. Prepare the data -------------------------------------------------------

# Keep only records that can be logged (>0)
plant_height_data <- cleaned_traits_final |>
  filter(TraitNameNew == "PlantHeight", !is.na(StdValue), StdValue > 0)

# Calculate outliers per biome: values more than 5 SD from the biome mean on the
# LOG scale (the scale used in the models). On the original scale, tall trees
# and shrubs look extreme next to many small plants and would be removed.
plant_height_data <- plant_height_data |>
  group_by(biome) |>
  mutate(lower_bound = exp(mean(log(StdValue)) - 5 * sd(log(StdValue))),
         upper_bound = exp(mean(log(StdValue)) + 5 * sd(log(StdValue)))) |>
  ungroup()

# Check how many outliers were removed
cat("Plant height outliers removed:", sum(plant_height_data$StdValue < plant_height_data$lower_bound |
                                            plant_height_data$StdValue > plant_height_data$upper_bound), "\n")

# Remove the outliers and log-transform the trait
plant_height_data <- plant_height_data |>
  filter(StdValue >= lower_bound, StdValue <= upper_bound) |>
  dplyr::select(-lower_bound, -upper_bound) |>
  # log-transform the trait
  mutate(log_trait = log(StdValue))

# Check how much data there is on each side of the treeline
print(plant_height_data |>
        group_by(biome) |>
        summarise(n_records = n(),
                  n_species = n_distinct(StandardSpeciesName),
                  n_sites = n_distinct(site_name),
                  n_locations = n_distinct(location)))

# Number of independent locations on each side
n_locations_boreal_ph <- n_distinct(plant_height_data$location[plant_height_data$distance_to_boundary_km < 0])
n_locations_tundra_ph <- n_distinct(plant_height_data$location[plant_height_data$distance_to_boundary_km >= 0])

# Check if there is enough locations
enough_data_ph <- n_locations_boreal_ph >= min_locations_per_side &
  n_locations_tundra_ph >= min_locations_per_side

if (!enough_data_ph) {
  cat("NOTE: fewer than", min_locations_per_side, "locations on one side of the treeline.",
      "The breakpoint models cannot be fitted for Plant height.\n")
}

## 4.2. Question 1: is there a breakpoint at the treeline (0 km)? --------------

# Compare two models:
#   - one straight line along the whole gradient (no breakpoint)
#   - a line that can change slope at the treeline (breakpoint at 0 km)
# Both are fitted with maximum likelihood (REML = FALSE), which is needed to
# compare models with different fixed effects. Because the breakpoint location
# (0 km) is fixed in advance, the likelihood ratio test with 1 df is valid.

# Add the breakpoint variables for a breakpoint at 0 km
data_treeline_ph <- plant_height_data |>
  mutate(below = pmin(distance_to_boundary_km, 0),    # boreal side
         above = pmax(distance_to_boundary_km, 0))    # tundra side

# Run the linear model
model_linear_ph <- lmerTest::lmer(linear_formula, data = data_treeline_ph, REML = FALSE)

# Run the model with breakpoint at "treeline" (0km)
model_treeline_ml_ph <- lmerTest::lmer(breakpoint_formula, data = data_treeline_ph,
                                       REML = FALSE)

# Use a likelihood ratio test to check if the breakpoint model fits significantly better
lr_treeline_ph <- 2 * (as.numeric(logLik(model_treeline_ml_ph)) -
                         as.numeric(logLik(model_linear_ph)))
p_treeline_ph <- pchisq(lr_treeline_ph, df = 1, lower.tail = FALSE)

# Extract delta AIC
daic_treeline_ph <- AIC(model_linear_ph) - AIC(model_treeline_ml_ph)

# Quick output of Q1 result
cat("\nQUESTION 1 - breakpoint at the treeline: LR =", round(lr_treeline_ph, 2),
    ", p =", signif(p_treeline_ph, 3),
    ", delta AIC =", round(daic_treeline_ph, 2),
    "(positive = breakpoint model fits better)\n")

# Fit slopes on each side of the "treeline"the same breakpoint model refitted with
# REML, which gives better estimates of the random effects for reporting.
# 'below' = boreal-side slope, 'above' = tundra-side slope (per km)
model_treeline_ph <- lmerTest::lmer(breakpoint_formula, data = data_treeline_ph,
                                    REML = TRUE)

# Print summary 
print(summary(model_treeline_ph))

# Get the coefficients
coefs_treeline_ph <- summary(model_treeline_ph)$coefficients

# Model diagnostics: residuals vs fitted values (should show no pattern) and
# a Q-Q plot (points should lie close to the line)
residuals_ph <- residuals(model_treeline_ph, type = "pearson")
png(here("figures", "RQ3_PlantHeight_model_diagnostics.png"),
    width = 10, height = 5, units = "in", res = 300)
par(mfrow = c(1, 2))
plot(fitted(model_treeline_ph), residuals_ph,
     xlab = "Fitted values", ylab = "Pearson residuals", main = "Plant height")
abline(h = 0, col = "red", lty = 2)
qqnorm(residuals_ph, main = "Plant height - Q-Q Plot")
qqline(residuals_ph, col = "red")
par(mfrow = c(1, 1))
dev.off()

## 4.3. Question 2: is there a breakpoint somewhere else? ----------------------
# The breakpoint model is fitted with the breakpoint at many candidate
# locations. The best-fitting one is compared with (a) no breakpoint and
# (b) a breakpoint at the treeline. Candidates whose fit is not significantly
# worse than the best (deviance difference <= 3.84) form the 95% CI.
# Only trust the location if the profile has one clear dip and the CI does not
# reach the edge of the search range.

# Range of candidates: the central 90% of the data
search_range_ph <- quantile(plant_height_data$distance_to_boundary_km, c(profile_trim, 1 - profile_trim))

# Step 1: coarse grid of candidates; fit the breakpoint model at each one and
# store how well it fits (log-likelihood)
coarse_psi_ph <- seq(ceiling(search_range_ph[1]), floor(search_range_ph[2]),
                     by = coarse_step_km)
coarse_loglik_ph <- numeric(length(coarse_psi_ph))
cat("Fitting", length(coarse_psi_ph), "coarse candidate breakpoints...\n")

for (i in seq_along(coarse_psi_ph)) {
  candidate_data <- plant_height_data |>
    mutate(below = pmin(distance_to_boundary_km - coarse_psi_ph[i], 0),
           above = pmax(distance_to_boundary_km - coarse_psi_ph[i], 0))
  candidate_model <- suppressMessages(suppressWarnings(
    lmerTest::lmer(breakpoint_formula, data = candidate_data, REML = FALSE,
                   control = lmerControl(calc.derivs = FALSE))))
  coarse_loglik_ph[i] <- as.numeric(logLik(candidate_model))
}

# Step 2: finer grid of candidates around the best coarse candidate
best_coarse_ph <- coarse_psi_ph[which.max(coarse_loglik_ph)]
fine_psi_ph <- seq(max(ceiling(search_range_ph[1]), best_coarse_ph - fine_window_km),
                   min(floor(search_range_ph[2]), best_coarse_ph + fine_window_km),
                   by = fine_step_km)
fine_psi_ph <- setdiff(fine_psi_ph, coarse_psi_ph)   # skip candidates already fitted
fine_loglik_ph <- numeric(length(fine_psi_ph))
cat("Fitting", length(fine_psi_ph), "fine candidate breakpoints...\n")

for (i in seq_along(fine_psi_ph)) {
  candidate_data <- plant_height_data |>
    mutate(below = pmin(distance_to_boundary_km - fine_psi_ph[i], 0),
           above = pmax(distance_to_boundary_km - fine_psi_ph[i], 0))
  candidate_model <- suppressMessages(suppressWarnings(
    lmerTest::lmer(breakpoint_formula, data = candidate_data, REML = FALSE,
                   control = lmerControl(calc.derivs = FALSE))))
  fine_loglik_ph[i] <- as.numeric(logLik(candidate_model))
}

# All candidates, sorted by distance. delta_deviance = how much worse each
# candidate fits than the best one (0 = the best candidate)
profile_ph <- data.frame(psi_km = c(coarse_psi_ph, fine_psi_ph),
                         logLik = c(coarse_loglik_ph, fine_loglik_ph)) |>
  arrange(psi_km) |>
  mutate(delta_deviance = 2 * (max(logLik) - logLik))

# Best breakpoint
best_index_ph <- which.max(profile_ph$logLik)
breakpoint_ph <- profile_ph$psi_km[best_index_ph]

# 95% CI: starting at the best candidate, move outwards in both directions for
# as long as the candidates fit almost as well (delta_deviance <= 3.84)
ci_threshold <- qchisq(0.95, df = 1)

lower_index_ph <- best_index_ph
while (lower_index_ph > 1 &&
       profile_ph$delta_deviance[lower_index_ph - 1] <= ci_threshold) {
  lower_index_ph <- lower_index_ph - 1
}

upper_index_ph <- best_index_ph
while (upper_index_ph < nrow(profile_ph) &&
       profile_ph$delta_deviance[upper_index_ph + 1] <= ci_threshold) {
  upper_index_ph <- upper_index_ph + 1
}

breakpoint_lower_ph <- profile_ph$psi_km[lower_index_ph]
breakpoint_upper_ph <- profile_ph$psi_km[upper_index_ph]

# Warning signs that the location cannot be trusted:
#   - the CI reaches the first or last candidate (the edge of the search range)
#   - a second, separate group of candidates also fits almost as well
ci_at_edge_ph <- lower_index_ph == 1 | upper_index_ph == nrow(profile_ph)
second_dip_ph <- any(profile_ph$delta_deviance <= ci_threshold &
                       (seq_len(nrow(profile_ph)) < lower_index_ph |
                          seq_len(nrow(profile_ph)) > upper_index_ph))

# (a) Best breakpoint vs no breakpoint: delta AIC, counting the breakpoint
#     location as an extra estimated parameter (positive = breakpoint better)
daic_best_vs_linear_ph <- 2 * (max(profile_ph$logLik) -
                                 as.numeric(logLik(model_linear_ph))) - 2 * 2

# (b) Best breakpoint vs breakpoint at the treeline: deviance difference
#     (below 3.84 = the treeline fits about as well as the best location)
deviance_treeline_vs_best_ph <- 2 * (max(profile_ph$logLik) -
                                       as.numeric(logLik(model_treeline_ml_ph)))

cat("\nQUESTION 2 - best breakpoint:", breakpoint_ph, "km (95% CI",
    breakpoint_lower_ph, "to", breakpoint_upper_ph, "km)\n",
    "  delta AIC vs no breakpoint:", round(daic_best_vs_linear_ph, 2), "\n",
    "  treeline vs best location, deviance difference:",
    round(deviance_treeline_vs_best_ph, 2),
    ifelse(deviance_treeline_vs_best_ph <= ci_threshold,
           "(treeline fits about as well)", "(treeline fits significantly worse)"), "\n",
    ifelse(ci_at_edge_ph, "  WARNING: CI reaches the edge of the search range\n", ""),
    ifelse(second_dip_ph, "  WARNING: second dip in the profile\n", ""))

# Supplementary figure: how well each candidate fits. Candidates below the red
# line lie within the 95% CI; the solid line marks the best candidate and the
# dashed grey line the treeline.
plot_profile_ph <- ggplot(profile_ph, aes(x = psi_km, y = delta_deviance)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1) +
  geom_hline(yintercept = ci_threshold, linetype = "dashed", colour = "red") +
  geom_vline(xintercept = breakpoint_ph, colour = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  labs(x = "Candidate breakpoint (km from the treeline)",
       y = "Deviance relative to best breakpoint",
       title = "Plant height: points below the red line lie within the 95% CI") +
  theme_rq3

# Save figure 
ggsave(here("figures", "FigureS10_RQ3_PlantHeight_breakpoint_profile.png"),
       plot = plot_profile_ph, width = 9, height = 6, dpi = 300)


## 4.4. Figure -------------------------------------------------------------------

# All records (log scale) against distance to the treeline, with the fitted
# line from the breakpoint model at the treeline (step .2). The title gives the
# result of Question 1; the labels give the slope on each side.

# Slope labels, converted from log-scale slopes to % change per km
label_boreal_ph <- paste0("Boreal side: ",
                          sprintf("%+.2f", 100 * (exp(coefs_treeline_ph["below", "Estimate"]) - 1)),
                          "% per km\np = ",
                          sprintf("%.3f", coefs_treeline_ph["below", "Pr(>|t|)"]))
label_tundra_ph <- paste0("Tundra side: ",
                          sprintf("%+.2f", 100 * (exp(coefs_treeline_ph["above", "Estimate"]) - 1)),
                          "% per km\np = ",
                          sprintf("%.3f", coefs_treeline_ph["above", "Pr(>|t|)"]))

# Title: result of the test for a breakpoint at the treeline
title_ph <- paste0("Breakpoint at the treeline: LR = ", sprintf("%.2f", lr_treeline_ph),
                   ", p = ", sprintf("%.3f", p_treeline_ph),
                   ", \u0394AIC = ", sprintf("%.1f", daic_treeline_ph))

# Fitted line for an average species, location and site (random effects left
# out with re.form = NA)
line_ph <- data.frame(distance_to_boundary_km = seq(min(plant_height_data$distance_to_boundary_km),
                                                    max(plant_height_data$distance_to_boundary_km),
                                                    length.out = 200)) |>
  mutate(below = pmin(distance_to_boundary_km, 0),
         above = pmax(distance_to_boundary_km, 0))
line_ph$predicted <- predict(model_treeline_ph, newdata = line_ph, re.form = NA)

# Plot figure
(plot_ph <- ggplot() +
  geom_point(data = plant_height_data, aes(x = distance_to_boundary_km, y = log_trait, colour = biome),
             alpha = 0.2, size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray30", linewidth = 0.7) +
  geom_line(data = line_ph, aes(x = distance_to_boundary_km, y = predicted),
            colour = fit_line_colour, linewidth = 1.5) +
  annotate("text", x = -Inf, y = Inf, label = label_boreal_ph,
           hjust = -0.03, vjust = 1.2, size = fig_annot_text / .pt, fontface = "bold") +
  annotate("text", x = Inf, y = Inf, label = label_tundra_ph,
           hjust = 1.03, vjust = 1.2, size = fig_annot_text / .pt, fontface = "bold") +
  scale_colour_manual(values = biome_colours, labels = biome_labels, name = "Biome") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +   # room for the labels
  labs(x = "Distance to treeline (km)", y = "Plant height (m, log scale)", title = title_ph) +
  theme_rq3 +
  theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))))

# Save figure to file
ggsave(here("figures", "Figure4a_PlantHeight_treeline.png"), 
       plot = plot_ph, width = 10, height = 6, dpi = 600)

## 4.5. Robustness check: separate models for each side --------------------------
# The slopes in step .2 come from ONE model in which the boreal and tundra lines
# meet at the treeline and share the random effects. Here each side is
# analysed on its own, without those assumptions. If the slopes are similar,
# the reported slopes do not depend on the lines meeting at the treeline.
# Note: these models do NOT test whether the two slopes differ; that is the
# test in step .2.

# Extract the boreal and tundra sides
boreal_side_ph <- plant_height_data |> filter(distance_to_boundary_km < 0)
tundra_side_ph <- plant_height_data |> filter(distance_to_boundary_km >= 0)

# Test the boreal side & get coefficients
model_boreal_side_ph <- lmerTest::lmer(side_formula, data = boreal_side_ph, REML = TRUE)
coef_boreal_side_ph <- summary(model_boreal_side_ph)$coefficients["distance_to_boundary_km", ]

# Test the tundra side & get coefficients
model_tundra_side_ph <- lmerTest::lmer(side_formula, data = tundra_side_ph, REML = TRUE)
coef_tundra_side_ph <- summary(model_tundra_side_ph)$coefficients["distance_to_boundary_km", ]

# Slopes from the one model (step .2) next to the separate models, as % per km
print(data.frame(Side = c("Boreal side", "Tundra side"),
                 one_model_pct_per_km = 100 * (exp(coefs_treeline_ph[c("below", "above"), "Estimate"]) - 1),
                 one_model_p = coefs_treeline_ph[c("below", "above"), "Pr(>|t|)"],
                 separate_pct_per_km = 100 * (exp(c(coef_boreal_side_ph[["Estimate"]],
                                     coef_tundra_side_ph[["Estimate"]])) - 1),
                 separate_p = c(coef_boreal_side_ph[["Pr(>|t|)"]], coef_tundra_side_ph[["Pr(>|t|)"]]),
                 row.names = NULL
), digits = 3)

## 4.6. Collect results ------------------------------------------------------------

# One row with everything for this trait; combined across traits in section 8

results_ph <- data.frame(Trait = "PlantHeight",
                         n_records = nrow(plant_height_data),
                         n_species = n_distinct(plant_height_data$StandardSpeciesName),
                         n_sites = n_distinct(plant_height_data$site_name),
                         n_locations_boreal = n_locations_boreal_ph,
                         n_locations_tundra = n_locations_tundra_ph,
                         # Question 1: breakpoint at the treeline, and the slopes on each side
                         LR_treeline = lr_treeline_ph,
                         p_treeline = p_treeline_ph,
                         dAIC_treeline = daic_treeline_ph,
                         slope_boreal = coefs_treeline_ph["below", "Estimate"],
                         SE_boreal = coefs_treeline_ph["below", "Std. Error"],
                         df_boreal = coefs_treeline_ph["below", "df"],
                         p_boreal = coefs_treeline_ph["below", "Pr(>|t|)"],
                         slope_tundra = coefs_treeline_ph["above", "Estimate"],
                         SE_tundra = coefs_treeline_ph["above", "Std. Error"],
                         df_tundra = coefs_treeline_ph["above", "df"],
                         p_tundra = coefs_treeline_ph["above", "Pr(>|t|)"],
                         # Robustness check: separate models per side
                         sep_n_locations_boreal = n_distinct(boreal_side_ph$location),
                         sep_slope_boreal = coef_boreal_side_ph[["Estimate"]],
                         sep_SE_boreal = coef_boreal_side_ph[["Std. Error"]],
                         sep_df_boreal = coef_boreal_side_ph[["df"]],
                         sep_p_boreal = coef_boreal_side_ph[["Pr(>|t|)"]],
                         sep_n_locations_tundra = n_distinct(tundra_side_ph$location),
                         sep_slope_tundra = coef_tundra_side_ph[["Estimate"]],
                         sep_SE_tundra = coef_tundra_side_ph[["Std. Error"]],
                         sep_df_tundra = coef_tundra_side_ph[["df"]],
                         sep_p_tundra = coef_tundra_side_ph[["Pr(>|t|)"]],
                         # Question 2: best breakpoint anywhere along the gradient
                         best_breakpoint_km = breakpoint_ph,
                         best_breakpoint_lower95_km = breakpoint_lower_ph,
                         best_breakpoint_upper95_km = breakpoint_upper_ph,
                         ci_at_search_edge = ci_at_edge_ph,
                         second_dip = second_dip_ph,
                         dAIC_best_vs_linear = daic_best_vs_linear_ph,
                         deviance_treeline_vs_best = deviance_treeline_vs_best_ph)

# 5. SLA ------------------------------------------------------------

## 5.1. Prepare the data -------------------------------------------------------

# Records of this trait with a positive standardised value (log() needs > 0)
sla_data <- cleaned_traits_final |>
  filter(TraitNameNew == "SLA", !is.na(StdValue), StdValue > 0)

# Outliers: values more than 5 SD from their biome's mean, calculated on the
# LOG scale (the scale used in the models). On the original scale, e.g. tall
# trees look extreme next to many small plants and would be wrongly removed.
# SLA: the lower bound is at least 0.1 mm2/mg, because values below this are
# implausible for vascular plants.
sla_data <- sla_data |>
  group_by(biome) |>
  mutate(lower_bound = pmax(0.1, exp(mean(log(StdValue)) - 5 * sd(log(StdValue)))),
         upper_bound = exp(mean(log(StdValue)) + 5 * sd(log(StdValue)))) |>
  ungroup()

cat("SLA outliers removed:", sum(sla_data$StdValue < sla_data$lower_bound |
                                   sla_data$StdValue > sla_data$upper_bound), "\n")

# Remove the outliers and log-transform the trait
sla_data <- sla_data |>
  filter(StdValue >= lower_bound, StdValue <= upper_bound) |>
  dplyr::select(-lower_bound, -upper_bound) |>
  mutate(log_trait = log(StdValue))

# How much data there is on each side of the treeline
print(sla_data |>
        group_by(biome) |>
        summarise(n_records = n(),
                  n_species = n_distinct(StandardSpeciesName),
                  n_sites = n_distinct(site_name),
                  n_locations = n_distinct(location)))

# Number of independent locations on each side (reported in the results table)
n_locations_boreal_sla <- n_distinct(sla_data$location[sla_data$distance_to_boundary_km < 0])
n_locations_tundra_sla <- n_distinct(sla_data$location[sla_data$distance_to_boundary_km >= 0])

## 5.2. Question 1: is there a breakpoint at the treeline (0 km)? --------------
# Two models are compared:
#   - one straight line along the whole gradient (no breakpoint)
#   - a line that can change slope at the treeline (breakpoint at 0 km)
# Both are fitted with maximum likelihood (REML = FALSE), which is needed to
# compare models with different fixed effects. Because the breakpoint location
# (0 km) is fixed in advance, the likelihood ratio test with 1 df is valid.

# Add the breakpoint variables for a breakpoint at 0 km
data_treeline_sla <- sla_data |>
  mutate(below = pmin(distance_to_boundary_km, 0),    # boreal side
         above = pmax(distance_to_boundary_km, 0))    # tundra side

model_linear_sla <- lmerTest::lmer(linear_formula, data = data_treeline_sla, REML = FALSE)
model_treeline_ml_sla <- lmerTest::lmer(breakpoint_formula, data = data_treeline_sla,
                                        REML = FALSE)

# Likelihood ratio test: does the breakpoint model fit significantly better?
lr_treeline_sla <- 2 * (as.numeric(logLik(model_treeline_ml_sla)) -
                          as.numeric(logLik(model_linear_sla)))
p_treeline_sla <- pchisq(lr_treeline_sla, df = 1, lower.tail = FALSE)

# Delta AIC: positive = the breakpoint model fits better (~2 = weak, >4 = clearer support)
daic_treeline_sla <- AIC(model_linear_sla) - AIC(model_treeline_ml_sla)

cat("\nQUESTION 1 - breakpoint at the treeline: LR =", round(lr_treeline_sla, 2),
    ", p =", signif(p_treeline_sla, 3),
    ", delta AIC =", round(daic_treeline_sla, 2),
    "(positive = breakpoint model fits better)\n")

# Slopes on each side of the treeline: the same breakpoint model refitted with
# REML, which gives better estimates of the random effects for reporting.
# 'below' = boreal-side slope, 'above' = tundra-side slope (per km)
model_treeline_sla <- lmerTest::lmer(breakpoint_formula, data = data_treeline_sla,
                                     REML = TRUE)
print(summary(model_treeline_sla))
coefs_treeline_sla <- summary(model_treeline_sla)$coefficients

# Model diagnostics: residuals vs fitted values (should show no pattern) and
# a Q-Q plot (points should lie close to the line)
residuals_sla <- residuals(model_treeline_sla, type = "pearson")
png(here("figures", "RQ3_SLA_model_diagnostics.png"),
    width = 10, height = 5, units = "in", res = 300)
par(mfrow = c(1, 2))
plot(fitted(model_treeline_sla), residuals_sla,
     xlab = "Fitted values", ylab = "Pearson residuals", main = "SLA")
abline(h = 0, col = "red", lty = 2)
qqnorm(residuals_sla, main = "SLA - Q-Q Plot")
qqline(residuals_sla, col = "red")
par(mfrow = c(1, 1))
dev.off()

## 5.3. Question 2: is there a breakpoint somewhere else? ----------------------
# The breakpoint model is fitted with the breakpoint placed at many candidate
# locations along the gradient ("profile likelihood"). The best-fitting
# candidate is then compared with
#   (a) no breakpoint (delta AIC), and
#   (b) a breakpoint at the treeline (deviance difference).
# Candidates that fit almost as well as the best one (deviance difference
# <= 3.84) form the 95% confidence interval (CI) of the breakpoint location.
# Only trust the location if the profile has one clear dip and the CI does not
# reach the edge of the search range (warnings are printed if not).

# Range of candidates: the central 90% of the data
search_range_sla <- quantile(sla_data$distance_to_boundary_km, c(profile_trim, 1 - profile_trim))

# Step 1: coarse grid of candidates; fit the breakpoint model at each one and
# store how well it fits (log-likelihood)
coarse_psi_sla <- seq(ceiling(search_range_sla[1]), floor(search_range_sla[2]),
                      by = coarse_step_km)
coarse_loglik_sla <- numeric(length(coarse_psi_sla))
cat("Fitting", length(coarse_psi_sla), "coarse candidate breakpoints...\n")

for (i in seq_along(coarse_psi_sla)) {
  candidate_data <- sla_data |>
    mutate(below = pmin(distance_to_boundary_km - coarse_psi_sla[i], 0),
           above = pmax(distance_to_boundary_km - coarse_psi_sla[i], 0))
  candidate_model <- suppressMessages(suppressWarnings(
    lmerTest::lmer(breakpoint_formula, data = candidate_data, REML = FALSE,
                   control = lmerControl(calc.derivs = FALSE))))
  coarse_loglik_sla[i] <- as.numeric(logLik(candidate_model))
}

# Step 2: finer grid of candidates around the best coarse candidate
best_coarse_sla <- coarse_psi_sla[which.max(coarse_loglik_sla)]
fine_psi_sla <- seq(max(ceiling(search_range_sla[1]), best_coarse_sla - fine_window_km),
                    min(floor(search_range_sla[2]), best_coarse_sla + fine_window_km),
                    by = fine_step_km)
fine_psi_sla <- setdiff(fine_psi_sla, coarse_psi_sla)   # skip candidates already fitted
fine_loglik_sla <- numeric(length(fine_psi_sla))
cat("Fitting", length(fine_psi_sla), "fine candidate breakpoints...\n")

for (i in seq_along(fine_psi_sla)) {
  candidate_data <- sla_data |>
    mutate(below = pmin(distance_to_boundary_km - fine_psi_sla[i], 0),
           above = pmax(distance_to_boundary_km - fine_psi_sla[i], 0))
  candidate_model <- suppressMessages(suppressWarnings(
    lmerTest::lmer(breakpoint_formula, data = candidate_data, REML = FALSE,
                   control = lmerControl(calc.derivs = FALSE))))
  fine_loglik_sla[i] <- as.numeric(logLik(candidate_model))
}

# All candidates, sorted by distance. delta_deviance = how much worse each
# candidate fits than the best one (0 = the best candidate)
profile_sla <- data.frame(psi_km = c(coarse_psi_sla, fine_psi_sla),
                          logLik = c(coarse_loglik_sla, fine_loglik_sla)) |>
  arrange(psi_km) |>
  mutate(delta_deviance = 2 * (max(logLik) - logLik))

# Best breakpoint
best_index_sla <- which.max(profile_sla$logLik)
breakpoint_sla <- profile_sla$psi_km[best_index_sla]

# 95% CI: starting at the best candidate, move outwards in both directions for
# as long as the candidates fit almost as well (delta_deviance <= 3.84)
ci_threshold <- qchisq(0.95, df = 1)

lower_index_sla <- best_index_sla
while (lower_index_sla > 1 &&
       profile_sla$delta_deviance[lower_index_sla - 1] <= ci_threshold) {
  lower_index_sla <- lower_index_sla - 1
}

upper_index_sla <- best_index_sla
while (upper_index_sla < nrow(profile_sla) &&
       profile_sla$delta_deviance[upper_index_sla + 1] <= ci_threshold) {
  upper_index_sla <- upper_index_sla + 1
}

breakpoint_lower_sla <- profile_sla$psi_km[lower_index_sla]
breakpoint_upper_sla <- profile_sla$psi_km[upper_index_sla]

# Warning signs that the location cannot be trusted:
#   - the CI reaches the first or last candidate (the edge of the search range)
#   - a second, separate group of candidates also fits almost as well
ci_at_edge_sla <- lower_index_sla == 1 | upper_index_sla == nrow(profile_sla)
second_dip_sla <- any(profile_sla$delta_deviance <= ci_threshold &
                        (seq_len(nrow(profile_sla)) < lower_index_sla |
                           seq_len(nrow(profile_sla)) > upper_index_sla))

# (a) Best breakpoint vs no breakpoint: delta AIC, counting the breakpoint
#     location as an extra estimated parameter (positive = breakpoint better)
daic_best_vs_linear_sla <- 2 * (max(profile_sla$logLik) -
                                  as.numeric(logLik(model_linear_sla))) - 2 * 2

# (b) Best breakpoint vs breakpoint at the treeline: deviance difference
#     (below 3.84 = the treeline fits about as well as the best location)
deviance_treeline_vs_best_sla <- 2 * (max(profile_sla$logLik) -
                                        as.numeric(logLik(model_treeline_ml_sla)))

cat("\nQUESTION 2 - best breakpoint:", breakpoint_sla, "km (95% CI",
    breakpoint_lower_sla, "to", breakpoint_upper_sla, "km)\n",
    "  delta AIC vs no breakpoint:", round(daic_best_vs_linear_sla, 2), "\n",
    "  treeline vs best location, deviance difference:",
    round(deviance_treeline_vs_best_sla, 2),
    ifelse(deviance_treeline_vs_best_sla <= ci_threshold,
           "(treeline fits about as well)", "(treeline fits significantly worse)"), "\n",
    ifelse(ci_at_edge_sla, "  WARNING: CI reaches the edge of the search range\n", ""),
    ifelse(second_dip_sla, "  WARNING: second dip in the profile\n", ""))

# Supplementary figure: how well each candidate fits. Candidates below the red
# line lie within the 95% CI; the solid line marks the best candidate and the
# dashed grey line the treeline.
plot_profile_sla <- ggplot(profile_sla, aes(x = psi_km, y = delta_deviance)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1) +
  geom_hline(yintercept = ci_threshold, linetype = "dashed", colour = "red") +
  geom_vline(xintercept = breakpoint_sla, colour = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  labs(x = "Candidate breakpoint (km from the treeline)",
       y = "Deviance relative to best breakpoint",
       title = "SLA: points below the red line lie within the 95% CI") +
  theme_rq3

ggsave(here("figures", "FigureS_RQ3_SLA_breakpoint_profile.png"),
       plot = plot_profile_sla, width = 9, height = 6, dpi = 300)

## 5.4. Figure -------------------------------------------------------------------
# All records (log scale) against distance to the treeline, with the fitted
# line from the breakpoint model at the treeline (step .2). The title gives the
# result of Question 1; the labels give the slope on each side.

# Slope labels, converted from log-scale slopes to % change per km
label_boreal_sla <- paste0("Boreal side: ",
                           sprintf("%+.2f", 100 * (exp(coefs_treeline_sla["below", "Estimate"]) - 1)),
                           "% per km\np = ",
                           sprintf("%.3f", coefs_treeline_sla["below", "Pr(>|t|)"]))
label_tundra_sla <- paste0("Tundra side: ",
                           sprintf("%+.2f", 100 * (exp(coefs_treeline_sla["above", "Estimate"]) - 1)),
                           "% per km\np = ",
                           sprintf("%.3f", coefs_treeline_sla["above", "Pr(>|t|)"]))

# Title: result of the test for a breakpoint at the treeline
title_sla <- paste0("Breakpoint at the treeline: LR = ", sprintf("%.2f", lr_treeline_sla),
                    ", p = ", sprintf("%.3f", p_treeline_sla),
                    ", \u0394AIC = ", sprintf("%.1f", daic_treeline_sla))

# Fitted line for an average species, location and site (random effects left
# out with re.form = NA)
line_sla <- data.frame(distance_to_boundary_km = seq(min(sla_data$distance_to_boundary_km),
                                                     max(sla_data$distance_to_boundary_km),
                                                     length.out = 200)) |>
  mutate(below = pmin(distance_to_boundary_km, 0),
         above = pmax(distance_to_boundary_km, 0))
line_sla$predicted <- predict(model_treeline_sla, newdata = line_sla, re.form = NA)

plot_sla <- ggplot() +
  geom_point(data = sla_data, aes(x = distance_to_boundary_km, y = log_trait, colour = biome),
             alpha = 0.2, size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray30", linewidth = 0.7) +
  geom_line(data = line_sla, aes(x = distance_to_boundary_km, y = predicted),
            colour = fit_line_colour, linewidth = 1.5) +
  annotate("text", x = -Inf, y = Inf, label = label_boreal_sla,
           hjust = -0.03, vjust = 1.2, size = fig_annot_text / .pt, fontface = "bold") +
  annotate("text", x = Inf, y = Inf, label = label_tundra_sla,
           hjust = 1.03, vjust = 1.2, size = fig_annot_text / .pt, fontface = "bold") +
  scale_colour_manual(values = biome_colours, labels = biome_labels, name = "Biome") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +   # room for the labels
  labs(x = "Distance to treeline (km)", y = expression(bold(paste("SLA (mm"^2, " mg"^-1, ", log scale)"))), title = title_sla) +
  theme_rq3 +
  theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))

print(plot_sla)
ggsave(here("figures", "Figure4b_SLA_treeline.png"), plot = plot_sla, width = 10, height = 6, dpi = 600)

## 5.5. Robustness check: separate models for each side --------------------------
# The slopes in step .2 come from ONE model in which the boreal and tundra lines
# meet at the treeline and share the random effects. Here each side is
# analysed on its own, without those assumptions. If the slopes are similar,
# the reported slopes do not depend on the lines meeting at the treeline.
# Note: these models do NOT test whether the two slopes differ; that is the
# test in step .2.

boreal_side_sla <- sla_data |> filter(distance_to_boundary_km < 0)
tundra_side_sla <- sla_data |> filter(distance_to_boundary_km >= 0)

model_boreal_side_sla <- lmerTest::lmer(side_formula, data = boreal_side_sla, REML = TRUE)
coef_boreal_side_sla <- summary(model_boreal_side_sla)$coefficients["distance_to_boundary_km", ]

model_tundra_side_sla <- lmerTest::lmer(side_formula, data = tundra_side_sla, REML = TRUE)
coef_tundra_side_sla <- summary(model_tundra_side_sla)$coefficients["distance_to_boundary_km", ]

# Slopes from the one model (step .2) next to the separate models, as % per km
print(data.frame(
  Side = c("Boreal side", "Tundra side"),
  one_model_pct_per_km = 100 * (exp(coefs_treeline_sla[c("below", "above"), "Estimate"]) - 1),
  one_model_p = coefs_treeline_sla[c("below", "above"), "Pr(>|t|)"],
  separate_pct_per_km = 100 * (exp(c(coef_boreal_side_sla[["Estimate"]],
                                     coef_tundra_side_sla[["Estimate"]])) - 1),
  separate_p = c(coef_boreal_side_sla[["Pr(>|t|)"]], coef_tundra_side_sla[["Pr(>|t|)"]]),
  row.names = NULL
), digits = 3)

## 5.6. Collect the results ---------------------------------------------------------
# One row with everything for this trait; combined across traits in section 8

results_sla <- data.frame(
  Trait = "SLA",
  n_records = nrow(sla_data),
  n_species = n_distinct(sla_data$StandardSpeciesName),
  n_sites = n_distinct(sla_data$site_name),
  n_locations_boreal = n_locations_boreal_sla,
  n_locations_tundra = n_locations_tundra_sla,
  # Question 1: breakpoint at the treeline, and the slopes on each side
  LR_treeline = lr_treeline_sla,
  p_treeline = p_treeline_sla,
  dAIC_treeline = daic_treeline_sla,
  slope_boreal = coefs_treeline_sla["below", "Estimate"],
  SE_boreal = coefs_treeline_sla["below", "Std. Error"],
  df_boreal = coefs_treeline_sla["below", "df"],
  p_boreal = coefs_treeline_sla["below", "Pr(>|t|)"],
  slope_tundra = coefs_treeline_sla["above", "Estimate"],
  SE_tundra = coefs_treeline_sla["above", "Std. Error"],
  df_tundra = coefs_treeline_sla["above", "df"],
  p_tundra = coefs_treeline_sla["above", "Pr(>|t|)"],
  # Robustness check: separate models per side
  sep_n_locations_boreal = n_distinct(boreal_side_sla$location),
  sep_slope_boreal = coef_boreal_side_sla[["Estimate"]],
  sep_SE_boreal = coef_boreal_side_sla[["Std. Error"]],
  sep_df_boreal = coef_boreal_side_sla[["df"]],
  sep_p_boreal = coef_boreal_side_sla[["Pr(>|t|)"]],
  sep_n_locations_tundra = n_distinct(tundra_side_sla$location),
  sep_slope_tundra = coef_tundra_side_sla[["Estimate"]],
  sep_SE_tundra = coef_tundra_side_sla[["Std. Error"]],
  sep_df_tundra = coef_tundra_side_sla[["df"]],
  sep_p_tundra = coef_tundra_side_sla[["Pr(>|t|)"]],
  # Question 2: best breakpoint anywhere along the gradient
  best_breakpoint_km = breakpoint_sla,
  best_breakpoint_lower95_km = breakpoint_lower_sla,
  best_breakpoint_upper95_km = breakpoint_upper_sla,
  ci_at_search_edge = ci_at_edge_sla,
  second_dip = second_dip_sla,
  dAIC_best_vs_linear = daic_best_vs_linear_sla,
  deviance_treeline_vs_best = deviance_treeline_vs_best_sla
)

# 6. LEAF N ------------------------------------------------------------

## 6.1. Prepare the data -------------------------------------------------------

# Records of this trait with a positive standardised value (log() needs > 0)
leafn_data <- cleaned_traits_final |>
  filter(TraitNameNew == "LeafN", !is.na(StdValue), StdValue > 0)

# Outliers: values more than 5 SD from their biome's mean, calculated on the
# LOG scale (the scale used in the models). On the original scale, e.g. tall
# trees look extreme next to many small plants and would be wrongly removed.
leafn_data <- leafn_data |>
  group_by(biome) |>
  mutate(lower_bound = exp(mean(log(StdValue)) - 5 * sd(log(StdValue))),
         upper_bound = exp(mean(log(StdValue)) + 5 * sd(log(StdValue)))) |>
  ungroup()

cat("Leaf N outliers removed:", sum(leafn_data$StdValue < leafn_data$lower_bound |
                                      leafn_data$StdValue > leafn_data$upper_bound), "\n")

# Remove the outliers and log-transform the trait
leafn_data <- leafn_data |>
  filter(StdValue >= lower_bound, StdValue <= upper_bound) |>
  dplyr::select(-lower_bound, -upper_bound) |>
  mutate(log_trait = log(StdValue))

# How much data there is on each side of the treeline
print(leafn_data |>
        group_by(biome) |>
        summarise(n_records = n(),
                  n_species = n_distinct(StandardSpeciesName),
                  n_sites = n_distinct(site_name),
                  n_locations = n_distinct(location)))

# Number of independent locations on each side (reported in the results table)
n_locations_boreal_leafn <- n_distinct(leafn_data$location[leafn_data$distance_to_boundary_km < 0])
n_locations_tundra_leafn <- n_distinct(leafn_data$location[leafn_data$distance_to_boundary_km >= 0])

## 6.2. Question 1: is there a breakpoint at the treeline (0 km)? --------------
# Two models are compared:
#   - one straight line along the whole gradient (no breakpoint)
#   - a line that can change slope at the treeline (breakpoint at 0 km)
# Both are fitted with maximum likelihood (REML = FALSE), which is needed to
# compare models with different fixed effects. Because the breakpoint location
# (0 km) is fixed in advance, the likelihood ratio test with 1 df is valid.

# Add the breakpoint variables for a breakpoint at 0 km
data_treeline_leafn <- leafn_data |>
  mutate(below = pmin(distance_to_boundary_km, 0),    # boreal side
         above = pmax(distance_to_boundary_km, 0))    # tundra side

model_linear_leafn <- lmerTest::lmer(linear_formula, data = data_treeline_leafn, REML = FALSE)
model_treeline_ml_leafn <- lmerTest::lmer(breakpoint_formula, data = data_treeline_leafn,
                                          REML = FALSE)

# Likelihood ratio test: does the breakpoint model fit significantly better?
lr_treeline_leafn <- 2 * (as.numeric(logLik(model_treeline_ml_leafn)) -
                            as.numeric(logLik(model_linear_leafn)))
p_treeline_leafn <- pchisq(lr_treeline_leafn, df = 1, lower.tail = FALSE)

# Delta AIC: positive = the breakpoint model fits better (~2 = weak, >4 = clearer support)
daic_treeline_leafn <- AIC(model_linear_leafn) - AIC(model_treeline_ml_leafn)

cat("\nQUESTION 1 - breakpoint at the treeline: LR =", round(lr_treeline_leafn, 2),
    ", p =", signif(p_treeline_leafn, 3),
    ", delta AIC =", round(daic_treeline_leafn, 2),
    "(positive = breakpoint model fits better)\n")

# Slopes on each side of the treeline: the same breakpoint model refitted with
# REML, which gives better estimates of the random effects for reporting.
# 'below' = boreal-side slope, 'above' = tundra-side slope (per km)
model_treeline_leafn <- lmerTest::lmer(breakpoint_formula, data = data_treeline_leafn,
                                       REML = TRUE)
print(summary(model_treeline_leafn))
coefs_treeline_leafn <- summary(model_treeline_leafn)$coefficients

# Model diagnostics: residuals vs fitted values (should show no pattern) and
# a Q-Q plot (points should lie close to the line)
residuals_leafn <- residuals(model_treeline_leafn, type = "pearson")
png(here("figures", "RQ3_LeafN_model_diagnostics.png"),
    width = 10, height = 5, units = "in", res = 300)
par(mfrow = c(1, 2))
plot(fitted(model_treeline_leafn), residuals_leafn,
     xlab = "Fitted values", ylab = "Pearson residuals", main = "Leaf N")
abline(h = 0, col = "red", lty = 2)
qqnorm(residuals_leafn, main = "Leaf N - Q-Q Plot")
qqline(residuals_leafn, col = "red")
par(mfrow = c(1, 1))
dev.off()

## 6.3. Question 2: is there a breakpoint somewhere else? ----------------------
# The breakpoint model is fitted with the breakpoint placed at many candidate
# locations along the gradient ("profile likelihood"). The best-fitting
# candidate is then compared with
#   (a) no breakpoint (delta AIC), and
#   (b) a breakpoint at the treeline (deviance difference).
# Candidates that fit almost as well as the best one (deviance difference
# <= 3.84) form the 95% confidence interval (CI) of the breakpoint location.
# Only trust the location if the profile has one clear dip and the CI does not
# reach the edge of the search range (warnings are printed if not).

# Range of candidates: the central 90% of the data
search_range_leafn <- quantile(leafn_data$distance_to_boundary_km, c(profile_trim, 1 - profile_trim))

# Step 1: coarse grid of candidates; fit the breakpoint model at each one and
# store how well it fits (log-likelihood)
coarse_psi_leafn <- seq(ceiling(search_range_leafn[1]), floor(search_range_leafn[2]),
                        by = coarse_step_km)
coarse_loglik_leafn <- numeric(length(coarse_psi_leafn))
cat("Fitting", length(coarse_psi_leafn), "coarse candidate breakpoints...\n")

for (i in seq_along(coarse_psi_leafn)) {
  candidate_data <- leafn_data |>
    mutate(below = pmin(distance_to_boundary_km - coarse_psi_leafn[i], 0),
           above = pmax(distance_to_boundary_km - coarse_psi_leafn[i], 0))
  candidate_model <- suppressMessages(suppressWarnings(
    lmerTest::lmer(breakpoint_formula, data = candidate_data, REML = FALSE,
                   control = lmerControl(calc.derivs = FALSE))))
  coarse_loglik_leafn[i] <- as.numeric(logLik(candidate_model))
}

# Step 2: finer grid of candidates around the best coarse candidate
best_coarse_leafn <- coarse_psi_leafn[which.max(coarse_loglik_leafn)]
fine_psi_leafn <- seq(max(ceiling(search_range_leafn[1]), best_coarse_leafn - fine_window_km),
                      min(floor(search_range_leafn[2]), best_coarse_leafn + fine_window_km),
                      by = fine_step_km)
fine_psi_leafn <- setdiff(fine_psi_leafn, coarse_psi_leafn)   # skip candidates already fitted
fine_loglik_leafn <- numeric(length(fine_psi_leafn))
cat("Fitting", length(fine_psi_leafn), "fine candidate breakpoints...\n")

for (i in seq_along(fine_psi_leafn)) {
  candidate_data <- leafn_data |>
    mutate(below = pmin(distance_to_boundary_km - fine_psi_leafn[i], 0),
           above = pmax(distance_to_boundary_km - fine_psi_leafn[i], 0))
  candidate_model <- suppressMessages(suppressWarnings(
    lmerTest::lmer(breakpoint_formula, data = candidate_data, REML = FALSE,
                   control = lmerControl(calc.derivs = FALSE))))
  fine_loglik_leafn[i] <- as.numeric(logLik(candidate_model))
}

# All candidates, sorted by distance. delta_deviance = how much worse each
# candidate fits than the best one (0 = the best candidate)
profile_leafn <- data.frame(psi_km = c(coarse_psi_leafn, fine_psi_leafn),
                            logLik = c(coarse_loglik_leafn, fine_loglik_leafn)) |>
  arrange(psi_km) |>
  mutate(delta_deviance = 2 * (max(logLik) - logLik))

# Best breakpoint
best_index_leafn <- which.max(profile_leafn$logLik)
breakpoint_leafn <- profile_leafn$psi_km[best_index_leafn]

# 95% CI: starting at the best candidate, move outwards in both directions for
# as long as the candidates fit almost as well (delta_deviance <= 3.84)
ci_threshold <- qchisq(0.95, df = 1)

lower_index_leafn <- best_index_leafn
while (lower_index_leafn > 1 &&
       profile_leafn$delta_deviance[lower_index_leafn - 1] <= ci_threshold) {
  lower_index_leafn <- lower_index_leafn - 1
}

upper_index_leafn <- best_index_leafn
while (upper_index_leafn < nrow(profile_leafn) &&
       profile_leafn$delta_deviance[upper_index_leafn + 1] <= ci_threshold) {
  upper_index_leafn <- upper_index_leafn + 1
}

breakpoint_lower_leafn <- profile_leafn$psi_km[lower_index_leafn]
breakpoint_upper_leafn <- profile_leafn$psi_km[upper_index_leafn]

# Warning signs that the location cannot be trusted:
#   - the CI reaches the first or last candidate (the edge of the search range)
#   - a second, separate group of candidates also fits almost as well
ci_at_edge_leafn <- lower_index_leafn == 1 | upper_index_leafn == nrow(profile_leafn)
second_dip_leafn <- any(profile_leafn$delta_deviance <= ci_threshold &
                          (seq_len(nrow(profile_leafn)) < lower_index_leafn |
                             seq_len(nrow(profile_leafn)) > upper_index_leafn))

# (a) Best breakpoint vs no breakpoint: delta AIC, counting the breakpoint
#     location as an extra estimated parameter (positive = breakpoint better)
daic_best_vs_linear_leafn <- 2 * (max(profile_leafn$logLik) -
                                    as.numeric(logLik(model_linear_leafn))) - 2 * 2

# (b) Best breakpoint vs breakpoint at the treeline: deviance difference
#     (below 3.84 = the treeline fits about as well as the best location)
deviance_treeline_vs_best_leafn <- 2 * (max(profile_leafn$logLik) -
                                          as.numeric(logLik(model_treeline_ml_leafn)))

cat("\nQUESTION 2 - best breakpoint:", breakpoint_leafn, "km (95% CI",
    breakpoint_lower_leafn, "to", breakpoint_upper_leafn, "km)\n",
    "  delta AIC vs no breakpoint:", round(daic_best_vs_linear_leafn, 2), "\n",
    "  treeline vs best location, deviance difference:",
    round(deviance_treeline_vs_best_leafn, 2),
    ifelse(deviance_treeline_vs_best_leafn <= ci_threshold,
           "(treeline fits about as well)", "(treeline fits significantly worse)"), "\n",
    ifelse(ci_at_edge_leafn, "  WARNING: CI reaches the edge of the search range\n", ""),
    ifelse(second_dip_leafn, "  WARNING: second dip in the profile\n", ""))

# Supplementary figure: how well each candidate fits. Candidates below the red
# line lie within the 95% CI; the solid line marks the best candidate and the
# dashed grey line the treeline.
plot_profile_leafn <- ggplot(profile_leafn, aes(x = psi_km, y = delta_deviance)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1) +
  geom_hline(yintercept = ci_threshold, linetype = "dashed", colour = "red") +
  geom_vline(xintercept = breakpoint_leafn, colour = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  labs(x = "Candidate breakpoint (km from the treeline)",
       y = "Deviance relative to best breakpoint",
       title = "Leaf N: points below the red line lie within the 95% CI") +
  theme_rq3

ggsave(here("figures", "FigureS_RQ3_LeafN_breakpoint_profile.png"),
       plot = plot_profile_leafn, width = 9, height = 6, dpi = 300)

## 6.4. Figure -------------------------------------------------------------------
# All records (log scale) against distance to the treeline, with the fitted
# line from the breakpoint model at the treeline (step .2). The title gives the
# result of Question 1; the labels give the slope on each side.

# Slope labels, converted from log-scale slopes to % change per km
label_boreal_leafn <- paste0("Boreal side: ",
                             sprintf("%+.2f", 100 * (exp(coefs_treeline_leafn["below", "Estimate"]) - 1)),
                             "% per km\np = ",
                             sprintf("%.3f", coefs_treeline_leafn["below", "Pr(>|t|)"]))
label_tundra_leafn <- paste0("Tundra side: ",
                             sprintf("%+.2f", 100 * (exp(coefs_treeline_leafn["above", "Estimate"]) - 1)),
                             "% per km\np = ",
                             sprintf("%.3f", coefs_treeline_leafn["above", "Pr(>|t|)"]))

# Title: result of the test for a breakpoint at the treeline
title_leafn <- paste0("Breakpoint at the treeline: LR = ", sprintf("%.2f", lr_treeline_leafn),
                      ", p = ", sprintf("%.3f", p_treeline_leafn),
                      ", \u0394AIC = ", sprintf("%.1f", daic_treeline_leafn))

# Fitted line for an average species, location and site (random effects left
# out with re.form = NA)
line_leafn <- data.frame(distance_to_boundary_km = seq(min(leafn_data$distance_to_boundary_km),
                                                       max(leafn_data$distance_to_boundary_km),
                                                       length.out = 200)) |>
  mutate(below = pmin(distance_to_boundary_km, 0),
         above = pmax(distance_to_boundary_km, 0))
line_leafn$predicted <- predict(model_treeline_leafn, newdata = line_leafn, re.form = NA)

plot_leafn <- ggplot() +
  geom_point(data = leafn_data, aes(x = distance_to_boundary_km, y = log_trait, colour = biome),
             alpha = 0.2, size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray30", linewidth = 0.7) +
  geom_line(data = line_leafn, aes(x = distance_to_boundary_km, y = predicted),
            colour = fit_line_colour, linewidth = 1.5) +
  annotate("text", x = -Inf, y = Inf, label = label_boreal_leafn,
           hjust = -0.03, vjust = 1.2, size = fig_annot_text / .pt, fontface = "bold") +
  annotate("text", x = Inf, y = Inf, label = label_tundra_leafn,
           hjust = 1.03, vjust = 1.2, size = fig_annot_text / .pt, fontface = "bold") +
  scale_colour_manual(values = biome_colours, labels = biome_labels, name = "Biome") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +   # room for the labels
  labs(x = "Distance to treeline (km)", y = "Leaf N (mg/g, log scale)", title = title_leafn) +
  theme_rq3 +
  theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))

print(plot_leafn)
ggsave(here("figures", "Figure4c_LeafN_treeline.png"), plot = plot_leafn, width = 10, height = 6, dpi = 600)

## 6.5. Robustness check: separate models for each side --------------------------
# The slopes in step .2 come from ONE model in which the boreal and tundra lines
# meet at the treeline and share the random effects. Here each side is
# analysed on its own, without those assumptions. If the slopes are similar,
# the reported slopes do not depend on the lines meeting at the treeline.
# Note: these models do NOT test whether the two slopes differ; that is the
# test in step .2.

boreal_side_leafn <- leafn_data |> filter(distance_to_boundary_km < 0)
tundra_side_leafn <- leafn_data |> filter(distance_to_boundary_km >= 0)

model_boreal_side_leafn <- lmerTest::lmer(side_formula, data = boreal_side_leafn, REML = TRUE)
coef_boreal_side_leafn <- summary(model_boreal_side_leafn)$coefficients["distance_to_boundary_km", ]

model_tundra_side_leafn <- lmerTest::lmer(side_formula, data = tundra_side_leafn, REML = TRUE)
coef_tundra_side_leafn <- summary(model_tundra_side_leafn)$coefficients["distance_to_boundary_km", ]

# Slopes from the one model (step .2) next to the separate models, as % per km
print(data.frame(
  Side = c("Boreal side", "Tundra side"),
  one_model_pct_per_km = 100 * (exp(coefs_treeline_leafn[c("below", "above"), "Estimate"]) - 1),
  one_model_p = coefs_treeline_leafn[c("below", "above"), "Pr(>|t|)"],
  separate_pct_per_km = 100 * (exp(c(coef_boreal_side_leafn[["Estimate"]],
                                     coef_tundra_side_leafn[["Estimate"]])) - 1),
  separate_p = c(coef_boreal_side_leafn[["Pr(>|t|)"]], coef_tundra_side_leafn[["Pr(>|t|)"]]),
  row.names = NULL
), digits = 3)

## 6.6. Collect the results ---------------------------------------------------------
# One row with everything for this trait; combined across traits in section 8

results_leafn <- data.frame(
  Trait = "LeafN",
  n_records = nrow(leafn_data),
  n_species = n_distinct(leafn_data$StandardSpeciesName),
  n_sites = n_distinct(leafn_data$site_name),
  n_locations_boreal = n_locations_boreal_leafn,
  n_locations_tundra = n_locations_tundra_leafn,
  # Question 1: breakpoint at the treeline, and the slopes on each side
  LR_treeline = lr_treeline_leafn,
  p_treeline = p_treeline_leafn,
  dAIC_treeline = daic_treeline_leafn,
  slope_boreal = coefs_treeline_leafn["below", "Estimate"],
  SE_boreal = coefs_treeline_leafn["below", "Std. Error"],
  df_boreal = coefs_treeline_leafn["below", "df"],
  p_boreal = coefs_treeline_leafn["below", "Pr(>|t|)"],
  slope_tundra = coefs_treeline_leafn["above", "Estimate"],
  SE_tundra = coefs_treeline_leafn["above", "Std. Error"],
  df_tundra = coefs_treeline_leafn["above", "df"],
  p_tundra = coefs_treeline_leafn["above", "Pr(>|t|)"],
  # Robustness check: separate models per side
  sep_n_locations_boreal = n_distinct(boreal_side_leafn$location),
  sep_slope_boreal = coef_boreal_side_leafn[["Estimate"]],
  sep_SE_boreal = coef_boreal_side_leafn[["Std. Error"]],
  sep_df_boreal = coef_boreal_side_leafn[["df"]],
  sep_p_boreal = coef_boreal_side_leafn[["Pr(>|t|)"]],
  sep_n_locations_tundra = n_distinct(tundra_side_leafn$location),
  sep_slope_tundra = coef_tundra_side_leafn[["Estimate"]],
  sep_SE_tundra = coef_tundra_side_leafn[["Std. Error"]],
  sep_df_tundra = coef_tundra_side_leafn[["df"]],
  sep_p_tundra = coef_tundra_side_leafn[["Pr(>|t|)"]],
  # Question 2: best breakpoint anywhere along the gradient
  best_breakpoint_km = breakpoint_leafn,
  best_breakpoint_lower95_km = breakpoint_lower_leafn,
  best_breakpoint_upper95_km = breakpoint_upper_leafn,
  ci_at_search_edge = ci_at_edge_leafn,
  second_dip = second_dip_leafn,
  dAIC_best_vs_linear = daic_best_vs_linear_leafn,
  deviance_treeline_vs_best = deviance_treeline_vs_best_leafn
)

# 7. SEED MASS ------------------------------------------------------------

# N.B: very few boreal seed mass records (3 locations in total), so the
# breakpoint models will most likely be marked 'not fitted'.

## 7.1. Prepare the data -------------------------------------------------------

# Records of this trait with a positive standardised value (log() needs > 0)
seedmass_data <- cleaned_traits_final |>
  filter(TraitNameNew == "SeedMass", !is.na(StdValue), StdValue > 0)

# Outliers: values more than 5 SD from their biome's mean, calculated on the
# LOG scale (the scale used in the models). On the original scale, e.g. tall
# trees look extreme next to many small plants and would be wrongly removed.
seedmass_data <- seedmass_data |>
  group_by(biome) |>
  mutate(lower_bound = exp(mean(log(StdValue)) - 5 * sd(log(StdValue))),
         upper_bound = exp(mean(log(StdValue)) + 5 * sd(log(StdValue)))) |>
  ungroup()

cat("Seed mass outliers removed:", sum(seedmass_data$StdValue < seedmass_data$lower_bound |
                                         seedmass_data$StdValue > seedmass_data$upper_bound), "\n")

# Remove the outliers and log-transform the trait
seedmass_data <- seedmass_data |>
  filter(StdValue >= lower_bound, StdValue <= upper_bound) |>
  dplyr::select(-lower_bound, -upper_bound) |>
  mutate(log_trait = log(StdValue))

# How much data there is on each side of the treeline
print(seedmass_data |>
        group_by(biome) |>
        summarise(n_records = n(),
                  n_species = n_distinct(StandardSpeciesName),
                  n_sites = n_distinct(site_name),
                  n_locations = n_distinct(location)))

# Number of independent locations on each side (reported in the results table)
n_locations_boreal_seedmass <- n_distinct(seedmass_data$location[seedmass_data$distance_to_boundary_km < 0])
n_locations_tundra_seedmass <- n_distinct(seedmass_data$location[seedmass_data$distance_to_boundary_km >= 0])

## 7.2. Question 1: is there a breakpoint at the treeline (0 km)? --------------
# Two models are compared:
#   - one straight line along the whole gradient (no breakpoint)
#   - a line that can change slope at the treeline (breakpoint at 0 km)
# Both are fitted with maximum likelihood (REML = FALSE), which is needed to
# compare models with different fixed effects. Because the breakpoint location
# (0 km) is fixed in advance, the likelihood ratio test with 1 df is valid.

# Add the breakpoint variables for a breakpoint at 0 km
data_treeline_seedmass <- seedmass_data |>
  mutate(below = pmin(distance_to_boundary_km, 0),    # boreal side
         above = pmax(distance_to_boundary_km, 0))    # tundra side

model_linear_seedmass <- lmerTest::lmer(linear_formula, data = data_treeline_seedmass, REML = FALSE)
model_treeline_ml_seedmass <- lmerTest::lmer(breakpoint_formula, data = data_treeline_seedmass,
                                             REML = FALSE)

# Likelihood ratio test: does the breakpoint model fit significantly better?
lr_treeline_seedmass <- 2 * (as.numeric(logLik(model_treeline_ml_seedmass)) -
                               as.numeric(logLik(model_linear_seedmass)))
p_treeline_seedmass <- pchisq(lr_treeline_seedmass, df = 1, lower.tail = FALSE)

# Delta AIC: positive = the breakpoint model fits better (~2 = weak, >4 = clearer support)
daic_treeline_seedmass <- AIC(model_linear_seedmass) - AIC(model_treeline_ml_seedmass)

cat("\nQUESTION 1 - breakpoint at the treeline: LR =", round(lr_treeline_seedmass, 2),
    ", p =", signif(p_treeline_seedmass, 3),
    ", delta AIC =", round(daic_treeline_seedmass, 2),
    "(positive = breakpoint model fits better)\n")

# Slopes on each side of the treeline: the same breakpoint model refitted with
# REML, which gives better estimates of the random effects for reporting.
# 'below' = boreal-side slope, 'above' = tundra-side slope (per km)
model_treeline_seedmass <- lmerTest::lmer(breakpoint_formula, data = data_treeline_seedmass,
                                          REML = TRUE)
print(summary(model_treeline_seedmass))
coefs_treeline_seedmass <- summary(model_treeline_seedmass)$coefficients

# Model diagnostics: residuals vs fitted values (should show no pattern) and
# a Q-Q plot (points should lie close to the line)
residuals_seedmass <- residuals(model_treeline_seedmass, type = "pearson")
png(here("figures", "RQ3_SeedMass_model_diagnostics.png"),
    width = 10, height = 5, units = "in", res = 300)
par(mfrow = c(1, 2))
plot(fitted(model_treeline_seedmass), residuals_seedmass,
     xlab = "Fitted values", ylab = "Pearson residuals", main = "Seed mass")
abline(h = 0, col = "red", lty = 2)
qqnorm(residuals_seedmass, main = "Seed mass - Q-Q Plot")
qqline(residuals_seedmass, col = "red")
par(mfrow = c(1, 1))
dev.off()

## 7.3. Question 2: is there a breakpoint somewhere else? ----------------------
# The breakpoint model is fitted with the breakpoint placed at many candidate
# locations along the gradient ("profile likelihood"). The best-fitting
# candidate is then compared with
#   (a) no breakpoint (delta AIC), and
#   (b) a breakpoint at the treeline (deviance difference).
# Candidates that fit almost as well as the best one (deviance difference
# <= 3.84) form the 95% confidence interval (CI) of the breakpoint location.
# Only trust the location if the profile has one clear dip and the CI does not
# reach the edge of the search range (warnings are printed if not).

# Range of candidates: the central 90% of the data
search_range_seedmass <- quantile(seedmass_data$distance_to_boundary_km, c(profile_trim, 1 - profile_trim))

# Step 1: coarse grid of candidates; fit the breakpoint model at each one and
# store how well it fits (log-likelihood)
coarse_psi_seedmass <- seq(ceiling(search_range_seedmass[1]), floor(search_range_seedmass[2]),
                           by = coarse_step_km)
coarse_loglik_seedmass <- numeric(length(coarse_psi_seedmass))
cat("Fitting", length(coarse_psi_seedmass), "coarse candidate breakpoints...\n")

for (i in seq_along(coarse_psi_seedmass)) {
  candidate_data <- seedmass_data |>
    mutate(below = pmin(distance_to_boundary_km - coarse_psi_seedmass[i], 0),
           above = pmax(distance_to_boundary_km - coarse_psi_seedmass[i], 0))
  candidate_model <- suppressMessages(suppressWarnings(
    lmerTest::lmer(breakpoint_formula, data = candidate_data, REML = FALSE,
                   control = lmerControl(calc.derivs = FALSE))))
  coarse_loglik_seedmass[i] <- as.numeric(logLik(candidate_model))
}

# Step 2: finer grid of candidates around the best coarse candidate
best_coarse_seedmass <- coarse_psi_seedmass[which.max(coarse_loglik_seedmass)]
fine_psi_seedmass <- seq(max(ceiling(search_range_seedmass[1]), best_coarse_seedmass - fine_window_km),
                         min(floor(search_range_seedmass[2]), best_coarse_seedmass + fine_window_km),
                         by = fine_step_km)
fine_psi_seedmass <- setdiff(fine_psi_seedmass, coarse_psi_seedmass)   # skip candidates already fitted
fine_loglik_seedmass <- numeric(length(fine_psi_seedmass))
cat("Fitting", length(fine_psi_seedmass), "fine candidate breakpoints...\n")

for (i in seq_along(fine_psi_seedmass)) {
  candidate_data <- seedmass_data |>
    mutate(below = pmin(distance_to_boundary_km - fine_psi_seedmass[i], 0),
           above = pmax(distance_to_boundary_km - fine_psi_seedmass[i], 0))
  candidate_model <- suppressMessages(suppressWarnings(
    lmerTest::lmer(breakpoint_formula, data = candidate_data, REML = FALSE,
                   control = lmerControl(calc.derivs = FALSE))))
  fine_loglik_seedmass[i] <- as.numeric(logLik(candidate_model))
}

# All candidates, sorted by distance. delta_deviance = how much worse each
# candidate fits than the best one (0 = the best candidate)
profile_seedmass <- data.frame(psi_km = c(coarse_psi_seedmass, fine_psi_seedmass),
                               logLik = c(coarse_loglik_seedmass, fine_loglik_seedmass)) |>
  arrange(psi_km) |>
  mutate(delta_deviance = 2 * (max(logLik) - logLik))

# Best breakpoint
best_index_seedmass <- which.max(profile_seedmass$logLik)
breakpoint_seedmass <- profile_seedmass$psi_km[best_index_seedmass]

# 95% CI: starting at the best candidate, move outwards in both directions for
# as long as the candidates fit almost as well (delta_deviance <= 3.84)
ci_threshold <- qchisq(0.95, df = 1)

lower_index_seedmass <- best_index_seedmass
while (lower_index_seedmass > 1 &&
       profile_seedmass$delta_deviance[lower_index_seedmass - 1] <= ci_threshold) {
  lower_index_seedmass <- lower_index_seedmass - 1
}

upper_index_seedmass <- best_index_seedmass
while (upper_index_seedmass < nrow(profile_seedmass) &&
       profile_seedmass$delta_deviance[upper_index_seedmass + 1] <= ci_threshold) {
  upper_index_seedmass <- upper_index_seedmass + 1
}

breakpoint_lower_seedmass <- profile_seedmass$psi_km[lower_index_seedmass]
breakpoint_upper_seedmass <- profile_seedmass$psi_km[upper_index_seedmass]

# Warning signs that the location cannot be trusted:
#   - the CI reaches the first or last candidate (the edge of the search range)
#   - a second, separate group of candidates also fits almost as well
ci_at_edge_seedmass <- lower_index_seedmass == 1 | upper_index_seedmass == nrow(profile_seedmass)
second_dip_seedmass <- any(profile_seedmass$delta_deviance <= ci_threshold &
                             (seq_len(nrow(profile_seedmass)) < lower_index_seedmass |
                                seq_len(nrow(profile_seedmass)) > upper_index_seedmass))

# (a) Best breakpoint vs no breakpoint: delta AIC, counting the breakpoint
#     location as an extra estimated parameter (positive = breakpoint better)
daic_best_vs_linear_seedmass <- 2 * (max(profile_seedmass$logLik) -
                                       as.numeric(logLik(model_linear_seedmass))) - 2 * 2

# (b) Best breakpoint vs breakpoint at the treeline: deviance difference
#     (below 3.84 = the treeline fits about as well as the best location)
deviance_treeline_vs_best_seedmass <- 2 * (max(profile_seedmass$logLik) -
                                             as.numeric(logLik(model_treeline_ml_seedmass)))

cat("\nQUESTION 2 - best breakpoint:", breakpoint_seedmass, "km (95% CI",
    breakpoint_lower_seedmass, "to", breakpoint_upper_seedmass, "km)\n",
    "  delta AIC vs no breakpoint:", round(daic_best_vs_linear_seedmass, 2), "\n",
    "  treeline vs best location, deviance difference:",
    round(deviance_treeline_vs_best_seedmass, 2),
    ifelse(deviance_treeline_vs_best_seedmass <= ci_threshold,
           "(treeline fits about as well)", "(treeline fits significantly worse)"), "\n",
    ifelse(ci_at_edge_seedmass, "  WARNING: CI reaches the edge of the search range\n", ""),
    ifelse(second_dip_seedmass, "  WARNING: second dip in the profile\n", ""))

# Supplementary figure: how well each candidate fits. Candidates below the red
# line lie within the 95% CI; the solid line marks the best candidate and the
# dashed grey line the treeline.
plot_profile_seedmass <- ggplot(profile_seedmass, aes(x = psi_km, y = delta_deviance)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1) +
  geom_hline(yintercept = ci_threshold, linetype = "dashed", colour = "red") +
  geom_vline(xintercept = breakpoint_seedmass, colour = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  labs(x = "Candidate breakpoint (km from the treeline)",
       y = "Deviance relative to best breakpoint",
       title = "Seed mass: points below the red line lie within the 95% CI") +
  theme_rq3

ggsave(here("figures", "FigureS_RQ3_SeedMass_breakpoint_profile.png"),
       plot = plot_profile_seedmass, width = 9, height = 6, dpi = 300)

## 7.4. Figure -------------------------------------------------------------------
# All records (log scale) against distance to the treeline, with the fitted
# line from the breakpoint model at the treeline (step .2). The title gives the
# result of Question 1; the labels give the slope on each side.

# Slope labels, converted from log-scale slopes to % change per km
label_boreal_seedmass <- paste0("Boreal side: ",
                                sprintf("%+.2f", 100 * (exp(coefs_treeline_seedmass["below", "Estimate"]) - 1)),
                                "% per km\np = ",
                                sprintf("%.3f", coefs_treeline_seedmass["below", "Pr(>|t|)"]))
label_tundra_seedmass <- paste0("Tundra side: ",
                                sprintf("%+.2f", 100 * (exp(coefs_treeline_seedmass["above", "Estimate"]) - 1)),
                                "% per km\np = ",
                                sprintf("%.3f", coefs_treeline_seedmass["above", "Pr(>|t|)"]))

# Title: result of the test for a breakpoint at the treeline
title_seedmass <- paste0("Breakpoint at the treeline: LR = ", sprintf("%.2f", lr_treeline_seedmass),
                         ", p = ", sprintf("%.3f", p_treeline_seedmass),
                         ", \u0394AIC = ", sprintf("%.1f", daic_treeline_seedmass))

# Fitted line for an average species, location and site (random effects left
# out with re.form = NA)
line_seedmass <- data.frame(distance_to_boundary_km = seq(min(seedmass_data$distance_to_boundary_km),
                                                          max(seedmass_data$distance_to_boundary_km),
                                                          length.out = 200)) |>
  mutate(below = pmin(distance_to_boundary_km, 0),
         above = pmax(distance_to_boundary_km, 0))
line_seedmass$predicted <- predict(model_treeline_seedmass, newdata = line_seedmass, re.form = NA)

plot_seedmass <- ggplot() +
  geom_point(data = seedmass_data, aes(x = distance_to_boundary_km, y = log_trait, colour = biome),
             alpha = 0.2, size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray30", linewidth = 0.7) +
  geom_line(data = line_seedmass, aes(x = distance_to_boundary_km, y = predicted),
            colour = fit_line_colour, linewidth = 1.5) +
  annotate("text", x = -Inf, y = Inf, label = label_boreal_seedmass,
           hjust = -0.03, vjust = 1.2, size = fig_annot_text / .pt, fontface = "bold") +
  annotate("text", x = Inf, y = Inf, label = label_tundra_seedmass,
           hjust = 1.03, vjust = 1.2, size = fig_annot_text / .pt, fontface = "bold") +
  scale_colour_manual(values = biome_colours, labels = biome_labels, name = "Biome") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +   # room for the labels
  labs(x = "Distance to treeline (km)", y = "Seed mass (mg, log scale)", title = title_seedmass) +
  theme_rq3 +
  theme(legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))

print(plot_seedmass)
ggsave(here("figures", "Figure4d_SeedMass_treeline.png"), plot = plot_seedmass, width = 10, height = 6, dpi = 600)

## 7.5. Robustness check: separate models for each side --------------------------
# The slopes in step .2 come from ONE model in which the boreal and tundra lines
# meet at the treeline and share the random effects. Here each side is
# analysed on its own, without those assumptions. If the slopes are similar,
# the reported slopes do not depend on the lines meeting at the treeline.
# Note: these models do NOT test whether the two slopes differ; that is the
# test in step .2.

boreal_side_seedmass <- seedmass_data |> filter(distance_to_boundary_km < 0)
tundra_side_seedmass <- seedmass_data |> filter(distance_to_boundary_km >= 0)

model_boreal_side_seedmass <- lmerTest::lmer(side_formula, data = boreal_side_seedmass, REML = TRUE)
coef_boreal_side_seedmass <- summary(model_boreal_side_seedmass)$coefficients["distance_to_boundary_km", ]

model_tundra_side_seedmass <- lmerTest::lmer(side_formula, data = tundra_side_seedmass, REML = TRUE)
coef_tundra_side_seedmass <- summary(model_tundra_side_seedmass)$coefficients["distance_to_boundary_km", ]

# Slopes from the one model (step .2) next to the separate models, as % per km
print(data.frame(
  Side = c("Boreal side", "Tundra side"),
  one_model_pct_per_km = 100 * (exp(coefs_treeline_seedmass[c("below", "above"), "Estimate"]) - 1),
  one_model_p = coefs_treeline_seedmass[c("below", "above"), "Pr(>|t|)"],
  separate_pct_per_km = 100 * (exp(c(coef_boreal_side_seedmass[["Estimate"]],
                                     coef_tundra_side_seedmass[["Estimate"]])) - 1),
  separate_p = c(coef_boreal_side_seedmass[["Pr(>|t|)"]], coef_tundra_side_seedmass[["Pr(>|t|)"]]),
  row.names = NULL
), digits = 3)

## 7.6. Collect the results ---------------------------------------------------------
# One row with everything for this trait; combined across traits in section 8

results_seedmass <- data.frame(
  Trait = "SeedMass",
  n_records = nrow(seedmass_data),
  n_species = n_distinct(seedmass_data$StandardSpeciesName),
  n_sites = n_distinct(seedmass_data$site_name),
  n_locations_boreal = n_locations_boreal_seedmass,
  n_locations_tundra = n_locations_tundra_seedmass,
  # Question 1: breakpoint at the treeline, and the slopes on each side
  LR_treeline = lr_treeline_seedmass,
  p_treeline = p_treeline_seedmass,
  dAIC_treeline = daic_treeline_seedmass,
  slope_boreal = coefs_treeline_seedmass["below", "Estimate"],
  SE_boreal = coefs_treeline_seedmass["below", "Std. Error"],
  df_boreal = coefs_treeline_seedmass["below", "df"],
  p_boreal = coefs_treeline_seedmass["below", "Pr(>|t|)"],
  slope_tundra = coefs_treeline_seedmass["above", "Estimate"],
  SE_tundra = coefs_treeline_seedmass["above", "Std. Error"],
  df_tundra = coefs_treeline_seedmass["above", "df"],
  p_tundra = coefs_treeline_seedmass["above", "Pr(>|t|)"],
  # Robustness check: separate models per side
  sep_n_locations_boreal = n_distinct(boreal_side_seedmass$location),
  sep_slope_boreal = coef_boreal_side_seedmass[["Estimate"]],
  sep_SE_boreal = coef_boreal_side_seedmass[["Std. Error"]],
  sep_df_boreal = coef_boreal_side_seedmass[["df"]],
  sep_p_boreal = coef_boreal_side_seedmass[["Pr(>|t|)"]],
  sep_n_locations_tundra = n_distinct(tundra_side_seedmass$location),
  sep_slope_tundra = coef_tundra_side_seedmass[["Estimate"]],
  sep_SE_tundra = coef_tundra_side_seedmass[["Std. Error"]],
  sep_df_tundra = coef_tundra_side_seedmass[["df"]],
  sep_p_tundra = coef_tundra_side_seedmass[["Pr(>|t|)"]],
  # Question 2: best breakpoint anywhere along the gradient
  best_breakpoint_km = breakpoint_seedmass,
  best_breakpoint_lower95_km = breakpoint_lower_seedmass,
  best_breakpoint_upper95_km = breakpoint_upper_seedmass,
  ci_at_search_edge = ci_at_edge_seedmass,
  second_dip = second_dip_seedmass,
  dAIC_best_vs_linear = daic_best_vs_linear_seedmass,
  deviance_treeline_vs_best = deviance_treeline_vs_best_seedmass
)

# 8. RESULTS TABLE AND COMBINED FIGURE -----------------------------------------

## 8.1. Main results table -----------------------------------------------------

# Combine the four traits; add Holm-corrected p-values (across the four
# traits, for reference) and convert slopes and their 95% CIs to % per km
results <- bind_rows(results_ph, results_sla, results_leafn, results_seedmass) |>
  mutate(p_treeline_holm = p.adjust(p_treeline, method = "holm"),
         pct_boreal = 100 * (exp(slope_boreal) - 1),
         pct_boreal_lower = 100 * (exp(slope_boreal - qt(0.975, df_boreal) * SE_boreal) - 1),
         pct_boreal_upper = 100 * (exp(slope_boreal + qt(0.975, df_boreal) * SE_boreal) - 1),
         pct_tundra = 100 * (exp(slope_tundra) - 1),
         pct_tundra_lower = 100 * (exp(slope_tundra - qt(0.975, df_tundra) * SE_tundra) - 1),
         pct_tundra_upper = 100 * (exp(slope_tundra + qt(0.975, df_tundra) * SE_tundra) - 1))

# All numbers, unrounded
write.csv(results, here("data", "derived_data", "RQ3_results_full.csv"), row.names = FALSE)

# Formatted table for the supplement
results_table <- results |>
  transmute(
    Trait,
    `Records / species / sites` = paste0(n_records, " / ", n_species, " / ", n_sites),
    `Locations boreal / tundra` = paste0(n_locations_boreal, " / ", n_locations_tundra),
    # Question 1
    `Breakpoint at treeline: LR` = sprintf("%.2f", LR_treeline),
    p = format.pval(p_treeline, digits = 3, eps = 1e-4),
    `p (Holm)` = format.pval(p_treeline_holm, digits = 3, eps = 1e-4),
    `ΔAIC` = sprintf("%.1f", dAIC_treeline),
    `Boreal side: % per km (95% CI)` =
      paste0(sprintf("%.2f", pct_boreal), " (", sprintf("%.2f", pct_boreal_lower),
             " to ", sprintf("%.2f", pct_boreal_upper), ")"),
    `Tundra side: % per km (95% CI)` =
      paste0(sprintf("%.2f", pct_tundra), " (", sprintf("%.2f", pct_tundra_lower),
             " to ", sprintf("%.2f", pct_tundra_upper), ")"),
    # Question 2
    `Best breakpoint, km (95% CI)` =
      paste0(best_breakpoint_km, " (", best_breakpoint_lower95_km, " to ",
             best_breakpoint_upper95_km, ")",
             ifelse(ci_at_search_edge, " †", ""),
             ifelse(second_dip, " ‡", "")),
    `Best breakpoint: ΔAIC vs no breakpoint` = sprintf("%.1f", dAIC_best_vs_linear),
    `Treeline vs best: deviance difference` = sprintf("%.1f", deviance_treeline_vs_best))

print(tibble::as_tibble(results_table), width = Inf)

write.csv(results_table, here("data", "derived_data", "TableS_RQ3_breakpoints.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# Caption notes:
# - All records along the gradient; distance = distance to the treeline
#   (shared boreal-tundra boundary). Linear mixed models of log trait values
#   with random effects for species, location and site.
# - Breakpoint at treeline: likelihood ratio test of a model with a change in
#   slope at 0 km against a single straight line (1 df). ΔAIC = AIC(straight
#   line) - AIC(breakpoint model); positive = breakpoint model fits better.
#   Holm-corrected p across the four traits given for reference.
# - Slopes: % change per km towards the tundra on each side of the treeline
#   (from the log-scale slope b: 100 * (exp(b) - 1)).
# - Best breakpoint: location found by fitting the breakpoint at many candidate
#   locations (profile likelihood). ΔAIC vs no breakpoint counts the location
#   as an extra parameter. Treeline vs best: deviance difference between a
#   breakpoint at 0 km and at the best location; < 3.84 = the treeline fits
#   about as well. † = CI reaches the edge of the search range; ‡ = second dip
#   in the profile; in either case the location should not be interpreted.

## 8.2. Robustness check table: one model vs separate models per side ----------

# Convert the separate-model slopes and their 95% CIs to % per km
side_check <- results |>
  mutate(sep_pct_boreal = 100 * (exp(sep_slope_boreal) - 1),
         sep_pct_boreal_lower = 100 * (exp(sep_slope_boreal - qt(0.975, sep_df_boreal) * sep_SE_boreal) - 1),
         sep_pct_boreal_upper = 100 * (exp(sep_slope_boreal + qt(0.975, sep_df_boreal) * sep_SE_boreal) - 1),
         sep_pct_tundra = 100 * (exp(sep_slope_tundra) - 1),
         sep_pct_tundra_lower = 100 * (exp(sep_slope_tundra - qt(0.975, sep_df_tundra) * sep_SE_tundra) - 1),
         sep_pct_tundra_upper = 100 * (exp(sep_slope_tundra + qt(0.975, sep_df_tundra) * sep_SE_tundra) - 1))

# One row per trait and side
side_check_table <- bind_rows(
  side_check |>
    transmute(Trait, Side = "Boreal side",
              Locations = sep_n_locations_boreal,
              one_pct = pct_boreal, one_lower = pct_boreal_lower, one_upper = pct_boreal_upper,
              one_p = p_boreal,
              sep_pct = sep_pct_boreal, sep_lower = sep_pct_boreal_lower, sep_upper = sep_pct_boreal_upper,
              sep_p = sep_p_boreal),
  side_check |>
    transmute(Trait, Side = "Tundra side",
              Locations = sep_n_locations_tundra,
              one_pct = pct_tundra, one_lower = pct_tundra_lower, one_upper = pct_tundra_upper,
              one_p = p_tundra,
              sep_pct = sep_pct_tundra, sep_lower = sep_pct_tundra_lower, sep_upper = sep_pct_tundra_upper,
              sep_p = sep_p_tundra)
) |>
  arrange(Trait, Side) |>
  transmute(Trait, Side, Locations,
            `One model: % per km (95% CI)` =
              paste0(sprintf("%.2f", one_pct), " (", sprintf("%.2f", one_lower),
                     " to ", sprintf("%.2f", one_upper), ")"),
            `One model: p` = format.pval(one_p, digits = 3, eps = 1e-4),
            `Separate model: % per km (95% CI)` =
              paste0(sprintf("%.2f", sep_pct), " (", sprintf("%.2f", sep_lower),
                     " to ", sprintf("%.2f", sep_upper), ")"),
            `Separate model: p` = format.pval(sep_p, digits = 3, eps = 1e-4))

print(tibble::as_tibble(side_check_table), width = Inf)

write.csv(side_check_table, here("data", "derived_data", "TableS_RQ3_side_slopes_check.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# Caption notes:
# - One model: slopes from the breakpoint model at the treeline (the lines on
#   both sides meet at 0 km; random effects estimated from all data).
# - Separate model: each side analysed on its own (no assumption that the lines
#   meet; random effects estimated per side). Similar values show that the
#   reported slopes do not depend on the lines meeting at the treeline.
# - These p-values test whether each slope differs from zero, not whether the
#   two sides differ; that is the breakpoint test in the main table.

## 8.3. Combined figure ---------------------------------------------------------

all_traits_plot <- plot_grid(plot_ph, plot_sla, plot_leafn, plot_seedmass,
                             labels = c("a)", "b)", "c)", "d)"),
                             label_size = fig_panel_label, nrow = 2)

ggsave(here("figures", "Figure4_traits_treeline.png"),
       plot = all_traits_plot, width = 20, height = 12, dpi = 600)
ggsave(here("figures", "Figure4_traits_treeline.pdf"),
       plot = all_traits_plot, width = 20, height = 12, dpi = 600)

# END OF SCRIPT ----------------------------------------------------------------