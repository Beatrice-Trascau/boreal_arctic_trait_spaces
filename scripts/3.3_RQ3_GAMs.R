##----------------------------------------------------------------------------##
# PAPER 3: BOREAL AND ARCTIC PLANT SPECIES TRAIT SPACES 
# 3.3_RQ3_GAMs
# This script contains code to fit GAMs on the relationships between traits
# and distance to biome boundary
##----------------------------------------------------------------------------##

# 1. SETUP ---------------------------------------------------------------------

## 1.1. Settings ---------------------------------------------------------------

# Maximum flexibility (basis dimension) of the smooths (the distance smooth is
# only used for SLA; see section 2)
gam_k_distance  <- 10
gam_k_longitude <- 10

# Spatial autocorrelation: distance classes of the correlogram (km); use the
# same values as in 3.2a so the correlograms can be compared
correlogram_step_km <- 100
correlogram_max_km  <- 2000

# Figure settings (same as script 3.1)
biome_colours <- c("boreal" = "darkgreen", "tundra" = "darkblue")
biome_labels <- c("boreal" = "Boreal", "tundra" = "Tundra")
gam_label <- "GAM (distance + longitude, random effects)"
residual_colours <- c("Segmented (no random effects)" = "#E69F00",
                      "GAM (distance + longitude, random effects)" = "black")
fig_axis_text <- 16
fig_axis_title <- 18
fig_legend_text <- 16
fig_legend_title <- 18
fig_annot_text <- 14
fig_panel_label <- 22

## 1.2. Load packages and data -------------------------------------------------

library(here)
source(here("scripts", "0_setup.R"))
library(mgcv)   # GAMs with random effects
library(ape)    # Moran's I (spatial autocorrelation)

# Prepared trait data from 3.2a (plant_height_data, sla_data, leafn_data,
# seedmass_data): outliers removed, log-transformed, grouping factors added
load(here("data", "derived_data", "RQ3_trait_data.RData"))

# Correlograms of the segmented residuals from 3.2a (for comparison)
correlograms_segmented <- read.csv(here("data", "derived_data", "RQ3_correlograms_segmented.csv"))

# Figure theme used throughout
theme_rq3 <- theme_classic() +
  theme(axis.text = element_text(size = fig_axis_text),
        axis.title = element_text(size = fig_axis_title, face = "bold"),
        legend.text = element_text(size = fig_legend_text),
        legend.title = element_text(size = fig_legend_title, face = "bold"),
        plot.title = element_text(size = fig_annot_text, face = "bold"))

# 2. MODEL FORMULAS ------------------------------------------------------------

# GAM: distance to the treeline, a smooth of longitude, and random intercepts
# for species, location and site 

# Longitude wraps around (-180 and +180 degrees are the same line, e.g. across
# the Bering Strait), so its smooth is cyclic (bs = "cc"): the curve joins up
# at +/-180. The ends are set with gam_knots_longitude in the bam() calls.

# Distance enters in one of two ways:
#   - as a SMOOTH for SLA: in a GAM with a distance smooth, the SLA
#     relationship was clearly curved (edf = 2.60)
#   - as a LINEAR term for plant height, leaf N and seed mass: there the
#     distance smooth reduced to (almost) a straight line (edf = 1.00, 1.41
#     and 1.00), so one slope describes it and is simpler to report
gam_formula_smooth <- log_trait ~ s(distance_to_boundary_km, k = gam_k_distance) +
  s(LON_site, bs = "cc", k = gam_k_longitude) +
  s(species_f, bs = "re") + s(location_f, bs = "re") + s(site_f, bs = "re")

gam_formula_linear <- log_trait ~ distance_to_boundary_km +
  s(LON_site, bs = "cc", k = gam_k_longitude) +
  s(species_f, bs = "re") + s(location_f, bs = "re") + s(site_f, bs = "re")

# Ends of the cyclic longitude smooth
gam_knots_longitude <- list(LON_site = c(-180, 180))

# Distance classes of the correlogram (step .3 of each trait)
correlogram_breaks <- seq(0, correlogram_max_km, by = correlogram_step_km)

# 3. PLANT HEIGHT --------------------------------------------------------------

## 3.1. GAM: distance + longitude, with random effects -------------------------

# Distance to the treeline as a LINEAR term (one slope along the gradient; see section 2).
# Plus a cyclic smooth of longitude and random intercepts for species,
# location and site. For a smooth, edf (effective degrees of freedom)
# describes its shape: 1 = a straight line, higher = more curved; its p-value
# tests whether it differs from a flat line. The longitude smooth captures
# east-west (regional) differences.
gam_ph <- bam(gam_formula_linear, data = plant_height_data, knots = gam_knots_longitude,
              method = "fREML", discrete = TRUE)

# Check model summary
print(summary(gam_ph))
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             -1.6168847  0.0805494 -20.073   <2e-16 ***
#   distance_to_boundary_km -0.0001712  0.0001565  -1.094    0.274    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df       F p-value    
# s(LON_site)     2.806      8 1006047  0.0259 *  
#   s(species_f)  388.510    440   20327  <2e-16 ***
#   s(location_f) 103.761    148   17052  0.0630 .  
# s(site_f)     160.328    298    9930  0.0109 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =   0.81   Deviance explained = 81.4%
# fREML =  24236  Scale est. = 0.27723   n = 29601

# Distance effect: slope per km (estimate, SE, t value, p-value)
gam_distance_ph <- summary(gam_ph)$p.table["distance_to_boundary_km", ]

# Longitude smooth (edf, test statistic, p-value)
gam_longitude_ph <- summary(gam_ph)$s.table["s(LON_site)", ]

# Concurvity between distance and longitude: how much the two overlap
# (0 = not at all, 1 = completely; above ~0.8 their effects cannot be cleanly
# separated). Pairwise (full = FALSE), so the random effects are not
# included: distance and longitude only vary between locations, so the
# location random effect would always "explain" them completely.
# With distance as a linear term, it is part of the parametric terms
# ("para", together with the intercept) in the concurvity output.
concurvity_ph <- tryCatch(concurvity(gam_ph, full = FALSE)$worst, error = function(e) NULL)
print(concurvity_ph)
#                    para s(LON_site) s(species_f) s(location_f) s(site_f)
# para          1.0000000   0.5060442     1.000000             1         1
# s(LON_site)   0.5060442   1.0000000     0.888102             1         1
# s(species_f)  1.0000000   0.8881020     1.000000             1         1
# s(location_f) 1.0000000   1.0000000     1.000000             1         1
# s(site_f)     1.0000000   1.0000000     1.000000             1         1

# Predicted distance effect for an average species, location, site and
# longitude (all other terms left out of the prediction)
gam_curve_ph <- data.frame(distance_to_boundary_km = seq(min(plant_height_data$distance_to_boundary_km),
                                max(plant_height_data$distance_to_boundary_km), length.out = 300),
                           LON_site = median(plant_height_data$LON_site),
                           species_f = factor(levels(plant_height_data$species_f)[1], levels = levels(plant_height_data$species_f)),
                           location_f = factor(levels(plant_height_data$location_f)[1], levels = levels(plant_height_data$location_f)),
                           site_f = factor(levels(plant_height_data$site_f)[1], levels = levels(plant_height_data$site_f)))
gam_fit_ph <- predict(gam_ph, newdata = gam_curve_ph, se.fit = TRUE,
                      exclude = c("s(LON_site)", "s(species_f)", "s(location_f)", "s(site_f)"))
gam_curve_ph <- gam_curve_ph |>
  mutate(fit = as.numeric(gam_fit_ph$fit),
         lower = fit - 1.96 * as.numeric(gam_fit_ph$se.fit),
         upper = fit + 1.96 * as.numeric(gam_fit_ph$se.fit))

# Predicted longitude effect (distance at its median, random effects left
# out), over the range of longitudes in the data
gam_lon_curve_ph <- data.frame(LON_site = seq(min(plant_height_data$LON_site), max(plant_height_data$LON_site), length.out = 300),
                               distance_to_boundary_km = median(plant_height_data$distance_to_boundary_km),
                               species_f = factor(levels(plant_height_data$species_f)[1], levels = levels(plant_height_data$species_f)),
                               location_f = factor(levels(plant_height_data$location_f)[1], levels = levels(plant_height_data$location_f)),
                               site_f = factor(levels(plant_height_data$site_f)[1], levels = levels(plant_height_data$site_f)))
gam_lon_fit_ph <- predict(gam_ph, newdata = gam_lon_curve_ph, se.fit = TRUE,
                          exclude = c("s(species_f)", "s(location_f)", "s(site_f)"))
gam_lon_curve_ph <- gam_lon_curve_ph |>
  mutate(fit = as.numeric(gam_lon_fit_ph$fit),
         lower = fit - 1.96 * as.numeric(gam_lon_fit_ph$se.fit),
         upper = fit + 1.96 * as.numeric(gam_lon_fit_ph$se.fit))

## 3.2. Data coverage ----------------------------------------------------------

# One point per location, by longitude and distance to the treeline (dashed
# line), coloured by biome and sized by the number of records: shows where
# along the treeline, and how far from it, the data come from.
coverage_ph <- plant_height_data |>
  group_by(location, biome) |>
  summarise(longitude = median(LON_site),
            distance_km = median(distance_to_boundary_km),
            n_records = n(),
            .groups = "drop")

# Plot figure
(plot_coverage_ph <- ggplot(coverage_ph, aes(x = longitude, y = distance_km, colour = biome)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "gray30") +
  geom_point(aes(size = n_records), alpha = 0.6) +
  scale_colour_manual(values = biome_colours, labels = biome_labels, name = "Biome") +
  scale_x_continuous(limits = c(-180, 180), breaks = seq(-180, 180, by = 60)) +
  labs(x = "Longitude (\u00b0E)", y = "Distance to Biome Boundary (km)", size = "Records",
       title = "Plant height") +
  theme_rq3)

# Save figure
ggsave(here("figures", "FigureS_RQ3_PlantHeight_data_coverage.png"),
       plot = plot_coverage_ph, width = 10, height = 6, dpi = 300)

## 3.3. Spatial autocorrelation of the residuals -------------------------------

# Same check as for the segmented regression (3.2a): are the residuals of
# nearby sites more similar than those of distant sites? Residuals are
# averaged per site, then Moran's I over all sites (weights = 1 / distance
# between sites) and a correlogram (Moran's I within each class of distance
# BETWEEN sites).

# Get residuals
plant_height_data$resid_gam <- residuals(gam_ph)

# Mean residual per site, with the site's projected coordinates
sites_ph <- plant_height_data |>
  group_by(site_name) |>
  summarise(x_km = first(x_km),
            y_km = first(y_km),
            resid_gam = mean(resid_gam),
            n_records = n(),
            .groups = "drop")

# Distances between all pairs of sites (km)
site_dist_ph <- as.matrix(dist(cbind(sites_ph$x_km, sites_ph$y_km)))

# Moran's I over all sites (inverse-distance weights)
weights_ph <- 1 / site_dist_ph
weights_ph[!is.finite(weights_ph)] <- 0   # a site's distance to itself
moran_gam_ph <- Moran.I(sites_ph$resid_gam, weights_ph)

cat("\nPlant height - Moran's I of the GAM residuals (", nrow(sites_ph), " sites): I =",
    round(moran_gam_ph$observed, 3), ", p =", signif(moran_gam_ph$p.value, 3), "\n")
# Plant height - Moran's I of the GAM residuals ( 300  sites): I = 0.008 , p = 0.747

# Correlogram: Moran's I among pairs of sites within each distance class
correlogram_ph <- data.frame()
for (i in seq_len(length(correlogram_breaks) - 1)) {
  
  # Weight 1 for pairs of sites within this distance class, 0 otherwise
  class_weights <- (site_dist_ph > correlogram_breaks[i] &
                      site_dist_ph <= correlogram_breaks[i + 1]) * 1
  
  # Skip classes without any pairs of sites
  if (sum(class_weights) == 0) next
  
  moran_class <- Moran.I(sites_ph$resid_gam, class_weights)
  
  correlogram_ph <- bind_rows(correlogram_ph, data.frame(
    Trait = "PlantHeight",
    distance_km = (correlogram_breaks[i] + correlogram_breaks[i + 1]) / 2,
    n_pairs = sum(class_weights) / 2,
    Model = gam_label,
    morans_I = moran_class$observed,
    p = moran_class$p.value))
}

## 3.4. Model validation -------------------------------------------------------

# Residuals vs fitted values and a Q-Q plot
png(here("figures", "RQ3_PlantHeight_validation_GAM.png"),
    width = 12, height = 5, units = "in", res = 200)
par(mfrow = c(1, 2))

plot(fitted(gam_ph), plant_height_data$resid_gam, pch = 16, cex = 0.4, col = rgb(0, 0, 0, 0.2),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Plant height - GAM: residuals vs fitted")
abline(h = 0, col = "red", lty = 2)

qqnorm(plant_height_data$resid_gam, pch = 16, cex = 0.4, main = "Plant height - GAM: Q-Q plot")
qqline(plant_height_data$resid_gam, col = "red")

par(mfrow = c(1, 1))
dev.off()

## 3.5. Figures ----------------------------------------------------------------

# Plot GAM effects: distance effect (left) and longitude effect (right) averaged 
# across species, location and site with 95% CI
# LEft side
(plot_gam_distance_ph <- ggplot(gam_curve_ph, aes(x = distance_to_boundary_km)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70", alpha = 0.5) +
  geom_line(aes(y = fit), linewidth = 1.2) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray30") +
  labs(x = "Distance to Biome Boundary (km)", y = "Plant Height (m, log scale)",
       title = paste0("Distance (linear): ",
                      sprintf("%+.3f", 100 * (exp(gam_distance_ph[["Estimate"]]) - 1)),
                      "% per km, p ", ifelse(gam_distance_ph[["Pr(>|t|)"]] < 0.001, "< 0.001",
                                             paste0("= ", sprintf("%.3f", gam_distance_ph[["Pr(>|t|)"]]))))) +
  theme_rq3)

# Right side
(plot_gam_longitude_ph <- ggplot(gam_lon_curve_ph, aes(x = LON_site)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70", alpha = 0.5) +
  geom_line(aes(y = fit), linewidth = 1.2) +
  labs(x = "Longitude (\u00b0E)", y = "Plant Height (m, log scale)",
       title = paste0("Longitude: edf = ", sprintf("%.1f", gam_longitude_ph[["edf"]]),
                      ", p ", ifelse(gam_longitude_ph[["p-value"]] < 0.001, "< 0.001",
                                     paste0("= ", sprintf("%.3f", gam_longitude_ph[["p-value"]]))))) +
  theme_rq3)

# Combine into single plot 
(plot_gam_ph <- plot_grid(plot_gam_distance_ph, plot_gam_longitude_ph, nrow = 1))

# Save figure
ggsave(here("figures", "FigureS_RQ3_PlantHeight_GAM.png"),
       plot = plot_gam_ph, width = 16, height = 6, dpi = 300)

# Plot correlogram for the segmented and GAM residuals
# Moran's I of the site residuals per class of distance between sites, for the
# segmented regression (from 3.2a) and the GAM. Filled points = significant
# (p < 0.05). Values near 0 = no autocorrelation.

# Combine into a single df
correlogram_both_ph <- bind_rows(correlograms_segmented |> filter(Trait == "PlantHeight"),
                                 correlogram_ph)

# Plot the correlogram
(plot_correlogram_ph <- ggplot(correlogram_both_ph,
                              aes(x = distance_km, y = morans_I, colour = Model)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_line(linewidth = 0.8) +
  geom_point(aes(shape = p < 0.05), size = 3) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), name = "p < 0.05") +
  scale_colour_manual(values = residual_colours, name = NULL) +
  labs(x = "Distance between sites (km)", y = "Moran's I of residuals", title = "Plant height") +
  theme_rq3 +
  theme(legend.position = "bottom", legend.direction = "vertical"))

# Save to file
ggsave(here("figures", "FigureS_RQ3_PlantHeight_correlogram.png"),
       plot = plot_correlogram_ph, width = 10, height = 7, dpi = 300)

# Collect all results into single df
results_gam_ph <- data.frame(Trait = "PlantHeight",
                             distance_term = "linear",
                             distance_edf = NA,
                             distance_slope = gam_distance_ph[["Estimate"]],
                             distance_slope_lower = gam_distance_ph[["Estimate"]] - 1.96 * gam_distance_ph[["Std. Error"]],
                             distance_slope_upper = gam_distance_ph[["Estimate"]] + 1.96 * gam_distance_ph[["Std. Error"]],
                             distance_statistic = gam_distance_ph[["t value"]],   # t value; full CSV only
                             distance_p = gam_distance_ph[["Pr(>|t|)"]],
                             longitude_edf = gam_longitude_ph[["edf"]],
                             longitude_F = gam_longitude_ph[["F"]],   
                             longitude_p = gam_longitude_ph[["p-value"]],
                             concurvity_distance_longitude = ifelse(is.null(concurvity_ph), NA,
                                         concurvity_ph["para", "s(LON_site)"]),
                             deviance_explained_pct = 100 * summary(gam_ph)$dev.expl,
                             n_sites_autocorrelation = nrow(sites_ph),
                             morans_I = moran_gam_ph$observed,
                             morans_I_p = moran_gam_ph$p.value)

# 4. SLA -----------------------------------------------------------------------

## 4.1. GAM: distance + longitude, with random effects -------------------------

# Run GAM (using smooher on distance because previous testing revealed that edf > 1 for distance)
gam_sla <- bam(gam_formula_smooth, data = sla_data, knots = gam_knots_longitude,
               method = "fREML", discrete = TRUE)

# Get summary
print(summary(gam_sla))
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  2.82576    0.03698   76.42   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf  Ref.df         F p-value    
# s(distance_to_boundary_km) 2.598e+00   2.643     4.606 0.00353 ** 
#   s(LON_site)                1.264e-03   8.000     0.000 0.90680    
# s(species_f)               3.761e+02 400.000 10303.780 < 2e-16 ***
#   s(location_f)              6.578e+01 103.000  2458.292 0.06319 .  
# s(site_f)                  8.064e+01 186.000  1247.891 0.00380 ** 

# Distance effect: smooth (edf, test statistic, p-value)
gam_distance_sla <- summary(gam_sla)$s.table["s(distance_to_boundary_km)", ]

# Longitude smooth (edf, test statistic, p-value)
gam_longitude_sla <- summary(gam_sla)$s.table["s(LON_site)", ]

# Concurvity between distance and longitude: how much they overlap
concurvity_sla <- tryCatch(concurvity(gam_sla, full = FALSE)$worst, error = function(e) NULL)
print(concurvity_sla)
#                            para s(distance_to_boundary_km) s(LON_site) s(species_f) s(location_f) s(site_f)
# para                       1.0000000                  0.4247996   0.6102720    1.0000000             1         1
# s(distance_to_boundary_km) 0.4247996                  1.0000000   0.9086047    0.9319017             1         1
# s(LON_site)                0.6102720                  0.9086047   1.0000000    0.8863979             1         1
# s(species_f)               1.0000000                  0.9319017   0.8863979    1.0000000             1         1
# s(location_f)              1.0000000                  1.0000000   1.0000000    1.0000000             1         1
# s(site_f)                  1.0000000                  1.0000000   1.0000000    1.0000000             1         1

# Predicted distance effect for an average species, location, site and
# longitude (all other terms left out of the prediction)
gam_curve_sla <- data.frame(distance_to_boundary_km = seq(min(sla_data$distance_to_boundary_km),
                                max(sla_data$distance_to_boundary_km), length.out = 300),
                            LON_site = median(sla_data$LON_site),
                            species_f = factor(levels(sla_data$species_f)[1], levels = levels(sla_data$species_f)),
                            location_f = factor(levels(sla_data$location_f)[1], levels = levels(sla_data$location_f)),
                            site_f = factor(levels(sla_data$site_f)[1], levels = levels(sla_data$site_f)))
gam_fit_sla <- predict(gam_sla, newdata = gam_curve_sla, se.fit = TRUE,
                       exclude = c("s(LON_site)", "s(species_f)", "s(location_f)", "s(site_f)"))
gam_curve_sla <- gam_curve_sla |>
  mutate(fit = as.numeric(gam_fit_sla$fit),
         lower = fit - 1.96 * as.numeric(gam_fit_sla$se.fit),
         upper = fit + 1.96 * as.numeric(gam_fit_sla$se.fit))

# Predicted longitude effect (distance at its median, random effects left
# out), over the range of longitudes in the data
gam_lon_curve_sla <- data.frame(
  LON_site = seq(min(sla_data$LON_site), max(sla_data$LON_site), length.out = 300),
  distance_to_boundary_km = median(sla_data$distance_to_boundary_km),
  species_f = factor(levels(sla_data$species_f)[1], levels = levels(sla_data$species_f)),
  location_f = factor(levels(sla_data$location_f)[1], levels = levels(sla_data$location_f)),
  site_f = factor(levels(sla_data$site_f)[1], levels = levels(sla_data$site_f)))
gam_lon_fit_sla <- predict(gam_sla, newdata = gam_lon_curve_sla, se.fit = TRUE,
                           exclude = c("s(distance_to_boundary_km)", "s(species_f)", "s(location_f)", "s(site_f)"))
gam_lon_curve_sla <- gam_lon_curve_sla |>
  mutate(fit = as.numeric(gam_lon_fit_sla$fit),
         lower = fit - 1.96 * as.numeric(gam_lon_fit_sla$se.fit),
         upper = fit + 1.96 * as.numeric(gam_lon_fit_sla$se.fit))

## 4.2. Data coverage ----------------------------------------------------------

# Check data coverage
coverage_sla <- sla_data |>
  group_by(location, biome) |>
  summarise(longitude = median(LON_site),
            distance_km = median(distance_to_boundary_km),
            n_records = n(),
            .groups = "drop")

# Plot data coverage
(plot_coverage_sla <- ggplot(coverage_sla, aes(x = longitude, y = distance_km, colour = biome)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "gray30") +
  geom_point(aes(size = n_records), alpha = 0.6) +
  scale_colour_manual(values = biome_colours, labels = biome_labels, name = "Biome") +
  scale_x_continuous(limits = c(-180, 180), breaks = seq(-180, 180, by = 60)) +
  labs(x = "Longitude (\u00b0E)", y = "Distance to Biome Boundary (km)", size = "Records",
       title = "SLA") +
  theme_rq3)

# Saev to file
ggsave(here("figures", "FigureS_RQ3_SLA_data_coverage.png"),
       plot = plot_coverage_sla, width = 10, height = 6, dpi = 300)

## 4.3. Spatial autocorrelation of the residuals -------------------------------

# Get residuals
sla_data$resid_gam <- residuals(gam_sla)

# Mean residual per site, with the site's projected coordinates
sites_sla <- sla_data |>
  group_by(site_name) |>
  summarise(x_km = first(x_km),
            y_km = first(y_km),
            resid_gam = mean(resid_gam),
            n_records = n(),
            .groups = "drop")

# Distances between all pairs of sites (km)
site_dist_sla <- as.matrix(dist(cbind(sites_sla$x_km, sites_sla$y_km)))

# Moran's I over all sites (inverse-distance weights)
weights_sla <- 1 / site_dist_sla
weights_sla[!is.finite(weights_sla)] <- 0   # a site's distance to itself
moran_gam_sla <- Moran.I(sites_sla$resid_gam, weights_sla)

# Check Moran's I
cat("\nSLA - Moran's I of the GAM residuals (", nrow(sites_sla), " sites): I =",
    round(moran_gam_sla$observed, 3), ", p =", signif(moran_gam_sla$p.value, 3), "\n")
# SLA - Moran's I of the GAM residuals ( 188  sites): I = -0.117 , p = 0.0071 

# Correlogram: Moran's I among pairs of sites within each distance class
correlogram_sla <- data.frame()
for (i in seq_len(length(correlogram_breaks) - 1)) {
  
  # Weight 1 for pairs of sites within this distance class, 0 otherwise
  class_weights <- (site_dist_sla > correlogram_breaks[i] &
                      site_dist_sla <= correlogram_breaks[i + 1]) * 1
  
  # Skip classes without any pairs of sites
  if (sum(class_weights) == 0) next
  
  moran_class <- Moran.I(sites_sla$resid_gam, class_weights)
  
  correlogram_sla <- bind_rows(correlogram_sla, data.frame(
    Trait = "SLA",
    distance_km = (correlogram_breaks[i] + correlogram_breaks[i + 1]) / 2,
    n_pairs = sum(class_weights) / 2,
    Model = gam_label,
    morans_I = moran_class$observed,
    p = moran_class$p.value))
}

## 4.4. Model validation -------------------------------------------------------

# Plot residuals vs fitted values and qq plot
png(here("figures", "RQ3_SLA_validation_GAM.png"),
    width = 12, height = 5, units = "in", res = 200)
par(mfrow = c(1, 2))

plot(fitted(gam_sla), sla_data$resid_gam, pch = 16, cex = 0.4, col = rgb(0, 0, 0, 0.2),
     xlab = "Fitted values", ylab = "Residuals",
     main = "SLA - GAM: residuals vs fitted")
abline(h = 0, col = "red", lty = 2)

qqnorm(sla_data$resid_gam, pch = 16, cex = 0.4, main = "SLA - GAM: Q-Q plot")
qqline(sla_data$resid_gam, col = "red")

par(mfrow = c(1, 1))
dev.off()

## 4.5. Figures ----------------------------------------------------------------

# Plot GAM output for distance effect
(plot_gam_distance_sla <- ggplot(gam_curve_sla, aes(x = distance_to_boundary_km)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70", alpha = 0.5) +
  geom_line(aes(y = fit), linewidth = 1.2) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray30") +
  labs(x = "Distance to Biome Boundary (km)", y = expression(bold(paste("SLA (mm"^2, " mg"^-1, ", log scale)"))),
       title = paste0("Distance (smooth): edf = ", sprintf("%.1f", gam_distance_sla[["edf"]]),
                      ", p ", ifelse(gam_distance_sla[["p-value"]] < 0.001, "< 0.001",
                                     paste0("= ", sprintf("%.3f", gam_distance_sla[["p-value"]]))))) +
  theme_rq3)

# Plot GAM output for longitude effect
plot_gam_longitude_sla <- ggplot(gam_lon_curve_sla, aes(x = LON_site)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70", alpha = 0.5) +
  geom_line(aes(y = fit), linewidth = 1.2) +
  labs(x = "Longitude (\u00b0E)", y = expression(bold(paste("SLA (mm"^2, " mg"^-1, ", log scale)"))),
       title = paste0("Longitude: edf = ", sprintf("%.1f", gam_longitude_sla[["edf"]]),
                      ", p ", ifelse(gam_longitude_sla[["p-value"]] < 0.001, "< 0.001",
                                     paste0("= ", sprintf("%.3f", gam_longitude_sla[["p-value"]]))))) +
  theme_rq3

# Combine panels
(plot_gam_sla <- plot_grid(plot_gam_distance_sla, plot_gam_longitude_sla, nrow = 1))

# Save to file
ggsave(here("figures", "FigureS_RQ3_SLA_GAM.png"),
       plot = plot_gam_sla, width = 16, height = 6, dpi = 300)


# Combine correlogram results for segmented regression and GAM
correlogram_both_sla <- bind_rows(correlograms_segmented |> filter(Trait == "SLA"),
                                  correlogram_sla)

# Plot correlogram for segmented regression and GAM
(plot_correlogram_sla <- ggplot(correlogram_both_sla,
                               aes(x = distance_km, y = morans_I, colour = Model)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_line(linewidth = 0.8) +
  geom_point(aes(shape = p < 0.05), size = 3) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), name = "p < 0.05") +
  scale_colour_manual(values = residual_colours, name = NULL) +
  labs(x = "Distance between sites (km)", y = "Moran's I of residuals", title = "SLA") +
  theme_rq3 +
  theme(legend.position = "bottom", legend.direction = "vertical"))

# Save to file
ggsave(here("figures", "FigureS_RQ3_SLA_correlogram.png"),
       plot = plot_correlogram_sla, width = 10, height = 7, dpi = 300)

# Collect all results into a single df
results_gam_sla <- data.frame(Trait = "SLA",
                              distance_term = "smooth",
                              distance_edf = gam_distance_sla[["edf"]],
                              distance_slope = NA,
                              distance_slope_lower = NA,
                              distance_slope_upper = NA,
                              distance_statistic = gam_distance_sla[["F"]],
                              distance_p = gam_distance_sla[["p-value"]],
                              longitude_edf = gam_longitude_sla[["edf"]],
                              longitude_F = gam_longitude_sla[["F"]],
                              longitude_p = gam_longitude_sla[["p-value"]],
                              concurvity_distance_longitude = ifelse(is.null(concurvity_sla), NA,
                                         concurvity_sla["s(distance_to_boundary_km)", "s(LON_site)"]),
                              deviance_explained_pct = 100 * summary(gam_sla)$dev.expl,
                              n_sites_autocorrelation = nrow(sites_sla),
                              morans_I = moran_gam_sla$observed,
                              morans_I_p = moran_gam_sla$p.value)

# 5. LEAF N --------------------------------------------------------------------

## 5.1. GAM: distance + longitude, with random effects -------------------------

# Run GAM
gam_leafn <- bam(gam_formula_linear, data = leafn_data, knots = gam_knots_longitude,
                 method = "fREML", discrete = TRUE)

# Get GAM summary
print(summary(gam_leafn))
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             3.034e+00  2.939e-02 103.238   <2e-16 ***
#   distance_to_boundary_km 3.899e-05  4.367e-05   0.893    0.372    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df       F  p-value    
# s(LON_site)     3.893      8 44183.4  0.00251 ** 
#   s(species_f)  344.733    422   259.9  < 2e-16 ***
#   s(location_f)  20.904    158   135.0  0.01406 *  
#   s(site_f)     126.393    208   634.8 5.62e-05 ***

# Distance effect: slope per km (estimate, SE, t value, p-value)
gam_distance_leafn <- summary(gam_leafn)$p.table["distance_to_boundary_km", ]

# Longitude smooth (edf, test statistic, p-value)
gam_longitude_leafn <- summary(gam_leafn)$s.table["s(LON_site)", ]

# Concurvity between distance and longitude: how much the two overlap
concurvity_leafn <- tryCatch(concurvity(gam_leafn, full = FALSE)$worst, error = function(e) NULL)
print(concurvity_leafn)
# para s(LON_site) s(species_f) s(location_f) s(site_f)
# para          1.0000000   0.6538495    1.0000000             1         1
# s(LON_site)   0.6538495   1.0000000    0.9098405             1         1
# s(species_f)  1.0000000   0.9098405    1.0000000             1         1
# s(location_f) 1.0000000   1.0000000    1.0000000             1         1
# s(site_f)     1.0000000   1.0000000    1.0000000             1         1

# Predicted distance effect for an average species, location, site and
# longitude (all other terms left out of the prediction)
gam_curve_leafn <- data.frame(
  distance_to_boundary_km = seq(min(leafn_data$distance_to_boundary_km),
                                max(leafn_data$distance_to_boundary_km), length.out = 300),
  LON_site = median(leafn_data$LON_site),
  species_f = factor(levels(leafn_data$species_f)[1], levels = levels(leafn_data$species_f)),
  location_f = factor(levels(leafn_data$location_f)[1], levels = levels(leafn_data$location_f)),
  site_f = factor(levels(leafn_data$site_f)[1], levels = levels(leafn_data$site_f)))
gam_fit_leafn <- predict(gam_leafn, newdata = gam_curve_leafn, se.fit = TRUE,
                         exclude = c("s(LON_site)", "s(species_f)", "s(location_f)", "s(site_f)"))
gam_curve_leafn <- gam_curve_leafn |>
  mutate(fit = as.numeric(gam_fit_leafn$fit),
         lower = fit - 1.96 * as.numeric(gam_fit_leafn$se.fit),
         upper = fit + 1.96 * as.numeric(gam_fit_leafn$se.fit))

# Predicted longitude effect (distance at its median, random effects left
# out), over the range of longitudes in the data
gam_lon_curve_leafn <- data.frame(
  LON_site = seq(min(leafn_data$LON_site), max(leafn_data$LON_site), length.out = 300),
  distance_to_boundary_km = median(leafn_data$distance_to_boundary_km),
  species_f = factor(levels(leafn_data$species_f)[1], levels = levels(leafn_data$species_f)),
  location_f = factor(levels(leafn_data$location_f)[1], levels = levels(leafn_data$location_f)),
  site_f = factor(levels(leafn_data$site_f)[1], levels = levels(leafn_data$site_f)))
gam_lon_fit_leafn <- predict(gam_leafn, newdata = gam_lon_curve_leafn, se.fit = TRUE,
                             exclude = c("s(species_f)", "s(location_f)", "s(site_f)"))
gam_lon_curve_leafn <- gam_lon_curve_leafn |>
  mutate(fit = as.numeric(gam_lon_fit_leafn$fit),
         lower = fit - 1.96 * as.numeric(gam_lon_fit_leafn$se.fit),
         upper = fit + 1.96 * as.numeric(gam_lon_fit_leafn$se.fit))

## 5.2. Data coverage ----------------------------------------------------------

# Check coverage
coverage_leafn <- leafn_data |>
  group_by(location, biome) |>
  summarise(longitude = median(LON_site),
            distance_km = median(distance_to_boundary_km),
            n_records = n(),
            .groups = "drop")

# Plot coverage
(plot_coverage_leafn <- ggplot(coverage_leafn, aes(x = longitude, y = distance_km, colour = biome)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "gray30") +
  geom_point(aes(size = n_records), alpha = 0.6) +
  scale_colour_manual(values = biome_colours, labels = biome_labels, name = "Biome") +
  scale_x_continuous(limits = c(-180, 180), breaks = seq(-180, 180, by = 60)) +
  labs(x = "Longitude (\u00b0E)", y = "Distance to treeline (km)", size = "Records",
       title = "Leaf N") +
  theme_rq3)

# Save to file
ggsave(here("figures", "FigureS_RQ3_LeafN_data_coverage.png"),
       plot = plot_coverage_leafn, width = 10, height = 6, dpi = 300)

## 5.3. Spatial autocorrelation of the residuals -------------------------------

# Get residuals from the GAM
leafn_data$resid_gam <- residuals(gam_leafn)

# Mean residual per site, with the site's projected coordinates
sites_leafn <- leafn_data |>
  group_by(site_name) |>
  summarise(x_km = first(x_km),
            y_km = first(y_km),
            resid_gam = mean(resid_gam),
            n_records = n(),
            .groups = "drop")

# Distances between all pairs of sites (km)
site_dist_leafn <- as.matrix(dist(cbind(sites_leafn$x_km, sites_leafn$y_km)))

# Moran's I over all sites (inverse-distance weights)
weights_leafn <- 1 / site_dist_leafn
weights_leafn[!is.finite(weights_leafn)] <- 0   # a site's distance to itself
moran_gam_leafn <- Moran.I(sites_leafn$resid_gam, weights_leafn)

# Check Moran's I output
cat("\nLeaf N - Moran's I of the GAM residuals (", nrow(sites_leafn), " sites): I =",
    round(moran_gam_leafn$observed, 3), ", p =", signif(moran_gam_leafn$p.value, 3), "\n")
# Leaf N - Moran's I of the GAM residuals ( 210  sites): I = 0.019 , p = 0.525

# Correlogram: Moran's I among pairs of sites within each distance class
correlogram_leafn <- data.frame()
for (i in seq_len(length(correlogram_breaks) - 1)) {
  
  # Weight 1 for pairs of sites within this distance class, 0 otherwise
  class_weights <- (site_dist_leafn > correlogram_breaks[i] &
                      site_dist_leafn <= correlogram_breaks[i + 1]) * 1
  
  # Skip classes without any pairs of sites
  if (sum(class_weights) == 0) next
  
  moran_class <- Moran.I(sites_leafn$resid_gam, class_weights)
  
  correlogram_leafn <- bind_rows(correlogram_leafn, data.frame(
    Trait = "LeafN",
    distance_km = (correlogram_breaks[i] + correlogram_breaks[i + 1]) / 2,
    n_pairs = sum(class_weights) / 2,
    Model = gam_label,
    morans_I = moran_class$observed,
    p = moran_class$p.value))
}

## 5.4. Model validation -------------------------------------------------------

# Plot residuals vs fitted and QQ plot
png(here("figures", "RQ3_LeafN_validation_GAM.png"),
    width = 12, height = 5, units = "in", res = 200)
par(mfrow = c(1, 2))

plot(fitted(gam_leafn), leafn_data$resid_gam, pch = 16, cex = 0.4, col = rgb(0, 0, 0, 0.2),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Leaf N - GAM: residuals vs fitted")
abline(h = 0, col = "red", lty = 2)

qqnorm(leafn_data$resid_gam, pch = 16, cex = 0.4, main = "Leaf N - GAM: Q-Q plot")
qqline(leafn_data$resid_gam, col = "red")

par(mfrow = c(1, 1))
dev.off()

## 5.5. Figures ----------------------------------------------------------------

# Plot GAM effects of distance to biome boundary
(plot_gam_distance_leafn <- ggplot(gam_curve_leafn, aes(x = distance_to_boundary_km)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70", alpha = 0.5) +
  geom_line(aes(y = fit), linewidth = 1.2) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray30") +
  labs(x = "Distance to Biome Boundary (km)", y = "Leaf N (mg/g, log scale)",
       title = paste0("Distance (linear): ",
                      sprintf("%+.3f", 100 * (exp(gam_distance_leafn[["Estimate"]]) - 1)),
                      "% per km, p ", ifelse(gam_distance_leafn[["Pr(>|t|)"]] < 0.001, "< 0.001",
                                             paste0("= ", sprintf("%.3f", gam_distance_leafn[["Pr(>|t|)"]]))))) +
  theme_rq3)

# Plot GAM effects of longitude
(plot_gam_longitude_leafn <- ggplot(gam_lon_curve_leafn, aes(x = LON_site)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70", alpha = 0.5) +
  geom_line(aes(y = fit), linewidth = 1.2) +
  labs(x = "Longitude (\u00b0E)", y = "Leaf N (mg/g, log scale)",
       title = paste0("Longitude: edf = ", sprintf("%.1f", gam_longitude_leafn[["edf"]]),
                      ", p ", ifelse(gam_longitude_leafn[["p-value"]] < 0.001, "< 0.001",
                                     paste0("= ", sprintf("%.3f", gam_longitude_leafn[["p-value"]]))))) +
  theme_rq3)

# Combine the two panles into a signle figure
plot_gam_leafn <- plot_grid(plot_gam_distance_leafn, plot_gam_longitude_leafn, nrow = 1)

# Save to file
ggsave(here("figures", "FigureS_RQ3_LeafN_GAM.png"),
       plot = plot_gam_leafn, width = 16, height = 6, dpi = 300)

# Combine correlograms for segmented regression and GAM into a single df
correlogram_both_leafn <- bind_rows(correlograms_segmented |> filter(Trait == "LeafN"),
                                    correlogram_leafn)

# Plot correlogram
(plot_correlogram_leafn <- ggplot(correlogram_both_leafn,
                                 aes(x = distance_km, y = morans_I, colour = Model)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_line(linewidth = 0.8) +
  geom_point(aes(shape = p < 0.05), size = 3) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), name = "p < 0.05") +
  scale_colour_manual(values = residual_colours, name = NULL) +
  labs(x = "Distance between sites (km)", y = "Moran's I of residuals", title = "Leaf N") +
  theme_rq3 +
  theme(legend.position = "bottom", legend.direction = "vertical"))

# Save to file
ggsave(here("figures", "FigureS_RQ3_LeafN_correlogram.png"),
       plot = plot_correlogram_leafn, width = 10, height = 7, dpi = 300)

# Combine all results into a single df
results_gam_leafn <- data.frame(Trait = "LeafN",
                                distance_term = "linear",
                                distance_edf = NA,
                                distance_slope = gam_distance_leafn[["Estimate"]],
                                distance_slope_lower = gam_distance_leafn[["Estimate"]] - 1.96 * gam_distance_leafn[["Std. Error"]],
                                distance_slope_upper = gam_distance_leafn[["Estimate"]] + 1.96 * gam_distance_leafn[["Std. Error"]],
                                distance_statistic = gam_distance_leafn[["t value"]],
                                distance_p = gam_distance_leafn[["Pr(>|t|)"]],
                                longitude_edf = gam_longitude_leafn[["edf"]],
                                longitude_F = gam_longitude_leafn[["F"]],
                                longitude_p = gam_longitude_leafn[["p-value"]],
                                concurvity_distance_longitude = ifelse(is.null(concurvity_leafn), NA,
                                         concurvity_leafn["para", "s(LON_site)"]),
                                deviance_explained_pct = 100 * summary(gam_leafn)$dev.expl,
                                n_sites_autocorrelation = nrow(sites_leafn),
                                morans_I = moran_gam_leafn$observed,
                                morans_I_p = moran_gam_leafn$p.value)

# 6. SEED MASS -----------------------------------------------------------------

# N.B.: very few boreal seed mass records (3 locations in total), so results
# for the boreal side are very uncertain.

## 6.1. GAM: distance + longitude, with random effects -------------------------

# Run GAM without a smoother on the distance to biome boundary
gam_seedmass <- bam(gam_formula_linear, data = seedmass_data, knots = gam_knots_longitude,
                    method = "fREML", discrete = TRUE)

# Check summary
print(summary(gam_seedmass))
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             -0.7015475  0.1877819  -3.736 0.000203 ***
#   distance_to_boundary_km  0.0002994  0.0003761   0.796 0.426343    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df     F  p-value    
# s(LON_site)     0.2058      7  42.7 0.259202    
# s(species_f)  112.4555    122 126.3 1.31e-06 ***
#   s(location_f)   4.1419     19 111.7 0.124914    
# s(site_f)      18.7993     34 444.9 0.000154 ***

# Distance effect: slope per km (estimate, SE, t value, p-value)
gam_distance_seedmass <- summary(gam_seedmass)$p.table["distance_to_boundary_km", ]

# Longitude smooth (edf, test statistic, p-value)
gam_longitude_seedmass <- summary(gam_seedmass)$s.table["s(LON_site)", ]

# Concurvity between distance and longitude: how much the two overlap
concurvity_seedmass <- tryCatch(concurvity(gam_seedmass, full = FALSE)$worst, error = function(e) NULL)
print(concurvity_seedmass)
#                    para s(LON_site) s(species_f) s(location_f) s(site_f)
# para          1.0000000   0.5690017    1.0000000             1         1
# s(LON_site)   0.5690017   1.0000000    0.9681766             1         1
# s(species_f)  1.0000000   0.9681766    1.0000000             1         1
# s(location_f) 1.0000000   1.0000000    1.0000000             1         1
# s(site_f)     1.0000000   1.0000000    1.0000000             1         1

# Predicted distance effect for an average species, location, site and
# longitude (all other terms left out of the prediction)
gam_curve_seedmass <- data.frame(
  distance_to_boundary_km = seq(min(seedmass_data$distance_to_boundary_km),
                                max(seedmass_data$distance_to_boundary_km), length.out = 300),
  LON_site = median(seedmass_data$LON_site),
  species_f = factor(levels(seedmass_data$species_f)[1], levels = levels(seedmass_data$species_f)),
  location_f = factor(levels(seedmass_data$location_f)[1], levels = levels(seedmass_data$location_f)),
  site_f = factor(levels(seedmass_data$site_f)[1], levels = levels(seedmass_data$site_f)))
gam_fit_seedmass <- predict(gam_seedmass, newdata = gam_curve_seedmass, se.fit = TRUE,
                            exclude = c("s(LON_site)", "s(species_f)", "s(location_f)", "s(site_f)"))
gam_curve_seedmass <- gam_curve_seedmass |>
  mutate(fit = as.numeric(gam_fit_seedmass$fit),
         lower = fit - 1.96 * as.numeric(gam_fit_seedmass$se.fit),
         upper = fit + 1.96 * as.numeric(gam_fit_seedmass$se.fit))

# Predicted longitude effect (distance at its median, random effects left
# out), over the range of longitudes in the data
gam_lon_curve_seedmass <- data.frame(
  LON_site = seq(min(seedmass_data$LON_site), max(seedmass_data$LON_site), length.out = 300),
  distance_to_boundary_km = median(seedmass_data$distance_to_boundary_km),
  species_f = factor(levels(seedmass_data$species_f)[1], levels = levels(seedmass_data$species_f)),
  location_f = factor(levels(seedmass_data$location_f)[1], levels = levels(seedmass_data$location_f)),
  site_f = factor(levels(seedmass_data$site_f)[1], levels = levels(seedmass_data$site_f)))
gam_lon_fit_seedmass <- predict(gam_seedmass, newdata = gam_lon_curve_seedmass, se.fit = TRUE,
                                exclude = c("s(species_f)", "s(location_f)", "s(site_f)"))
gam_lon_curve_seedmass <- gam_lon_curve_seedmass |>
  mutate(fit = as.numeric(gam_lon_fit_seedmass$fit),
         lower = fit - 1.96 * as.numeric(gam_lon_fit_seedmass$se.fit),
         upper = fit + 1.96 * as.numeric(gam_lon_fit_seedmass$se.fit))

## 6.2. Data coverage ----------------------------------------------------------

# Check data coverage
coverage_seedmass <- seedmass_data |>
  group_by(location, biome) |>
  summarise(longitude = median(LON_site),
            distance_km = median(distance_to_boundary_km),
            n_records = n(),
            .groups = "drop")

# Plot data coverage figutre
(plot_coverage_seedmass <- ggplot(coverage_seedmass, aes(x = longitude, y = distance_km, colour = biome)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "gray30") +
  geom_point(aes(size = n_records), alpha = 0.6) +
  scale_colour_manual(values = biome_colours, labels = biome_labels, name = "Biome") +
  scale_x_continuous(limits = c(-180, 180), breaks = seq(-180, 180, by = 60)) +
  labs(x = "Longitude (\u00b0E)", y = "Distance to treeline (km)", size = "Records",
       title = "Seed mass") +
  theme_rq3)

# Save figure
ggsave(here("figures", "FigureS_RQ3_SeedMass_data_coverage.png"),
       plot = plot_coverage_seedmass, width = 10, height = 6, dpi = 300)

## 6.3. Spatial autocorrelation of the residuals -------------------------------

# Get residuals
seedmass_data$resid_gam <- residuals(gam_seedmass)

# Mean residual per site, with the site's projected coordinates
sites_seedmass <- seedmass_data |>
  group_by(site_name) |>
  summarise(x_km = first(x_km),
            y_km = first(y_km),
            resid_gam = mean(resid_gam),
            n_records = n(),
            .groups = "drop")

# Distances between all pairs of sites (km)
site_dist_seedmass <- as.matrix(dist(cbind(sites_seedmass$x_km, sites_seedmass$y_km)))

# Moran's I over all sites (inverse-distance weights)
weights_seedmass <- 1 / site_dist_seedmass
weights_seedmass[!is.finite(weights_seedmass)] <- 0   # a site's distance to itself
moran_gam_seedmass <- Moran.I(sites_seedmass$resid_gam, weights_seedmass)

# Check Moran's I
cat("\nSeed mass - Moran's I of the GAM residuals (", nrow(sites_seedmass), " sites): I =",
    round(moran_gam_seedmass$observed, 3), ", p =", signif(moran_gam_seedmass$p.value, 3), "\n")
# Seed mass - Moran's I of the GAM residuals ( 36  sites): I = -0.05 , p = 0.81 

# Correlogram: Moran's I among pairs of sites within each distance class
correlogram_seedmass <- data.frame()
for (i in seq_len(length(correlogram_breaks) - 1)) {
  
  # Weight 1 for pairs of sites within this distance class, 0 otherwise
  class_weights <- (site_dist_seedmass > correlogram_breaks[i] &
                      site_dist_seedmass <= correlogram_breaks[i + 1]) * 1
  
  # Skip classes without any pairs of sites
  if (sum(class_weights) == 0) next
  
  moran_class <- Moran.I(sites_seedmass$resid_gam, class_weights)
  
  correlogram_seedmass <- bind_rows(correlogram_seedmass, data.frame(
    Trait = "SeedMass",
    distance_km = (correlogram_breaks[i] + correlogram_breaks[i + 1]) / 2,
    n_pairs = sum(class_weights) / 2,
    Model = gam_label,
    morans_I = moran_class$observed,
    p = moran_class$p.value))
}

## 6.4. Model validation -------------------------------------------------------

# Plot Residuals vs fitted values and Q-Q plot
png(here("figures", "RQ3_SeedMass_validation_GAM.png"),
    width = 12, height = 5, units = "in", res = 200)
par(mfrow = c(1, 2))

plot(fitted(gam_seedmass), seedmass_data$resid_gam, pch = 16, cex = 0.4, col = rgb(0, 0, 0, 0.2),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Seed mass - GAM: residuals vs fitted")
abline(h = 0, col = "red", lty = 2)

qqnorm(seedmass_data$resid_gam, pch = 16, cex = 0.4, main = "Seed mass - GAM: Q-Q plot")
qqline(seedmass_data$resid_gam, col = "red")

par(mfrow = c(1, 1))
dev.off()

## 6.5. Figures ----------------------------------------------------------------

# Plot GAM effect of distance to biome boundary
(plot_gam_distance_seedmass <- ggplot(gam_curve_seedmass, aes(x = distance_to_boundary_km)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70", alpha = 0.5) +
  geom_line(aes(y = fit), linewidth = 1.2) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray30") +
  labs(x = "Distance to Biome Boundary (km)", y = "Seed mass (mg, log scale)",
       title = paste0("Distance (linear): ",
                      sprintf("%+.3f", 100 * (exp(gam_distance_seedmass[["Estimate"]]) - 1)),
                      "% per km, p ", ifelse(gam_distance_seedmass[["Pr(>|t|)"]] < 0.001, "< 0.001",
                                             paste0("= ", sprintf("%.3f", gam_distance_seedmass[["Pr(>|t|)"]]))))) +
  theme_rq3)

# Plot GAM effect of longitude
(plot_gam_longitude_seedmass <- ggplot(gam_lon_curve_seedmass, aes(x = LON_site)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70", alpha = 0.5) +
  geom_line(aes(y = fit), linewidth = 1.2) +
  labs(x = "Longitude (\u00b0E)", y = "Seed mass (mg, log scale)",
       title = paste0("Longitude: edf = ", sprintf("%.1f", gam_longitude_seedmass[["edf"]]),
                      ", p ", ifelse(gam_longitude_seedmass[["p-value"]] < 0.001, "< 0.001",
                                     paste0("= ", sprintf("%.3f", gam_longitude_seedmass[["p-value"]]))))) +
  theme_rq3)

# Combine into single figure
(plot_gam_seedmass <- plot_grid(plot_gam_distance_seedmass, 
                                plot_gam_longitude_seedmass, nrow = 1))

# Save to file
ggsave(here("figures", "FigureS_RQ3_SeedMass_GAM.png"),
       plot = plot_gam_seedmass, width = 16, height = 6, dpi = 300)

# Combine correlogram for segmented regression and GAM
correlogram_both_seedmass <- bind_rows(correlograms_segmented |> 
                                         filter(Trait == "SeedMass"),
                                       correlogram_seedmass)

# Plot correlogram
(plot_correlogram_seedmass <- ggplot(correlogram_both_seedmass,
                                    aes(x = distance_km, y = morans_I, colour = Model)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_line(linewidth = 0.8) +
  geom_point(aes(shape = p < 0.05), size = 3) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), name = "p < 0.05") +
  scale_colour_manual(values = residual_colours, name = NULL) +
  labs(x = "Distance between sites (km)", y = "Moran's I of residuals", title = "Seed mass") +
  theme_rq3 +
  theme(legend.position = "bottom", legend.direction = "vertical"))

# Save supplementary figure
ggsave(here("figures", "FigureS_RQ3_SeedMass_correlogram.png"),
       plot = plot_correlogram_seedmass, width = 10, height = 7, dpi = 300)

# Collect all results into single df
results_gam_seedmass <- data.frame(Trait = "SeedMass",
                                   distance_term = "linear",
                                   distance_edf = NA,
                                   distance_slope = gam_distance_seedmass[["Estimate"]],
                                   distance_slope_lower = gam_distance_seedmass[["Estimate"]] - 1.96 * gam_distance_seedmass[["Std. Error"]],
                                   distance_slope_upper = gam_distance_seedmass[["Estimate"]] + 1.96 * gam_distance_seedmass[["Std. Error"]],
                                   distance_statistic = gam_distance_seedmass[["t value"]],
                                   distance_p = gam_distance_seedmass[["Pr(>|t|)"]],
                                   longitude_edf = gam_longitude_seedmass[["edf"]],
                                   longitude_F = gam_longitude_seedmass[["F"]],
                                   longitude_p = gam_longitude_seedmass[["p-value"]],
                                   concurvity_distance_longitude = ifelse(is.null(concurvity_seedmass), NA,
                                         concurvity_seedmass["para", "s(LON_site)"]),
                                   deviance_explained_pct = 100 * summary(gam_seedmass)$dev.expl,
                                   n_sites_autocorrelation = nrow(sites_seedmass),
                                   morans_I = moran_gam_seedmass$observed,
                                   morans_I_p = moran_gam_seedmass$p.value)

# 7. RESULTS TABLES AND COMBINED FIGURES ---------------------------------------

## 7.1. GAM results table ------------------------------------------------------

# Combine all GAM results
gam_results <- bind_rows(results_gam_ph, results_gam_sla,
                         results_gam_leafn, results_gam_seedmass)

# Save to file
write.csv(gam_results, here("data", "derived_data", "RQ3_GAM_results_full.csv"), 
          row.names = FALSE)

# Tidy uup the results table to use in supplementary
gam_table <- gam_results |>
  transmute(Trait,
            `Distance term` = distance_term,
            `Distance: slope, % per km (95% CI), or edf` =
              ifelse(distance_term == "linear",
                     paste0(sprintf("%.3f", 100 * (exp(distance_slope) - 1)), " (",
                            sprintf("%.3f", 100 * (exp(distance_slope_lower) - 1)), " to ",
                            sprintf("%.3f", 100 * (exp(distance_slope_upper) - 1)), ")"),
                     paste0("edf = ", sprintf("%.2f", distance_edf))),
            `Distance: p` = format.pval(distance_p, digits = 3, eps = 1e-4),
            `Longitude: edf` = sprintf("%.2f", longitude_edf),
            `Longitude: p` = format.pval(longitude_p, digits = 3, eps = 1e-4),
            `Concurvity, distance vs longitude` = sprintf("%.2f", concurvity_distance_longitude),
            `Deviance explained (%)` = sprintf("%.1f", deviance_explained_pct),
            `Residuals: Moran's I (p)` =
              paste0(sprintf("%.3f", morans_I), " (",
                     format.pval(morans_I_p, digits = 3, eps = 1e-4), ")"))

print(tibble::as_tibble(gam_table), width = Inf)
# Trait       `Distance term` `Distance: slope, % per km (95% CI), or edf` `Distance: p` `Longitude: edf`
# <chr>       <chr>           <chr>                                        <chr>         <chr>           
#   1 PlantHeight linear          -0.017 (-0.048 to 0.014)                     0.27403       2.81            
# 2 SLA         smooth          edf = 2.60                                   0.00353       0.00            
# 3 LeafN       linear          0.004 (-0.005 to 0.012)                      0.37204       3.89            
# 4 SeedMass    linear          0.030 (-0.044 to 0.104)                      0.42634       0.21            
# `Longitude: p` `Concurvity, distance vs longitude` `Deviance explained (%)` `Residuals: Moran's I (p)`
# <chr>          <chr>                               <chr>                    <chr>                     
#   1 0.02591        0.51                                81.4                     0.008 (0.7465)            
# 2 0.90680        0.91                                84.6                     -0.117 (0.0071)           
# 3 0.00251        0.65                                79.6                     0.019 (0.5255)            
# 4 0.25920        0.57                                85.4                     -0.050 (0.8098)         

# Save to file
write.csv(gam_table, here("data", "derived_data", "TableS_RQ3_GAM.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# Caption notes:
# - GAM (mgcv::bam) of log trait values with distance to the treeline, a
#   cyclic smooth of longitude, and random intercepts for species, location
#   and site. Distance was a linear term (slope as % change per km towards the
#   tundra) for plant height, leaf N and seed mass, where a distance smooth
#   reduced to (almost) a straight line, and a smooth for SLA, where it was
#   curved. edf = effective degrees of freedom (1 = straight line).
# - Test statistics of the smooths are not shown: mgcv's statistic for
#   penalised terms cannot be read as a classical F when terms overlap (here
#   longitude and the location random effect); edf and the (approximate)
#   p-values are reported. The statistics are in RQ3_GAM_results_full.csv.
# - Concurvity ("worst" measure, pairwise): overlap between distance and
#   longitude; values above ~0.8 mean their effects cannot be cleanly
#   separated.
# - Deviance explained includes the random effects.
# - Moran's I of the mean residual per site (inverse-distance weights).

## 7.2. Combined figures -------------------------------------------------------

# Combine all the data coverage figures
all_traits_coverage <- plot_grid(plot_coverage_ph, plot_coverage_sla,
                                 plot_coverage_leafn, plot_coverage_seedmass,
                                 labels = c("a)", "b)", "c)", "d)"),
                                 label_size = fig_panel_label, nrow = 2)

# Save combined figure
ggsave(here("figures", "FigureS_RQ3_data_coverage.png"),
       plot = all_traits_coverage, width = 20, height = 12, dpi = 300)

# Combine all correlograms
all_traits_correlograms <- plot_grid(plot_correlogram_ph, plot_correlogram_sla,
                                     plot_correlogram_leafn, plot_correlogram_seedmass,
                                     labels = c("a)", "b)", "c)", "d)"),
                                     label_size = fig_panel_label, nrow = 2)

# Save the combined correlograms
ggsave(here("figures", "FigureS_RQ3_correlograms.png"),
       plot = all_traits_correlograms, width = 20, height = 14, dpi = 300)

# END OF SCRIPT ----------------------------------------------------------------