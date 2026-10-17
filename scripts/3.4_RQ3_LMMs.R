##----------------------------------------------------------------------------##
# PAPER 3: BOREAL AND ARCTIC PLANT SPECIES TRAIT SPACES 
# 3.4_RQ3_LMMs
# This script contains code to fit linear mixed models on the relationships 
# between traits and distance to biome boundary
##----------------------------------------------------------------------------##

# 1. SETUP ---------------------------------------------------------------------

## 1.1. Settings ---------------------------------------------------------------

# Spatial autocorrelation: distance classes of the correlogram (km); same
# values as in 3.2a and 3.2b
correlogram_step_km <- 100
correlogram_max_km  <- 2000

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

# 1.2. Load packages and data --------------------------------------------------

library(here)
source(here("scripts", "0_setup.R"))
library(lme4)       # mixed models
library(lmerTest)   # p-values for mixed models (Satterthwaite degrees of freedom)
library(ape)        # Moran's I (spatial autocorrelation)

# Prepared trait data from 3.2a (plant_height_data, sla_data, leafn_data,
# seedmass_data): outliers removed, log-transformed, grouping factors added
load(here("data", "derived_data", "RQ3_trait_data.RData"))

# Biome as a factor with boreal as the reference level, so that the
# coefficients are tundra relative to boreal
plant_height_data$biome <- factor(plant_height_data$biome, levels = c("boreal", "tundra"))
sla_data$biome <- factor(sla_data$biome, levels = c("boreal", "tundra"))
leafn_data$biome <- factor(leafn_data$biome, levels = c("boreal", "tundra"))
seedmass_data$biome <- factor(seedmass_data$biome, levels = c("boreal", "tundra"))

# Figure theme used throughout
theme_rq3 <- theme_classic() +
  theme(axis.text = element_text(size = fig_axis_text),
        axis.title = element_text(size = fig_axis_title, face = "bold"),
        legend.text = element_text(size = fig_legend_text),
        legend.title = element_text(size = fig_legend_title, face = "bold"),
        plot.title = element_text(size = fig_annot_text, face = "bold"))

# 2. MODEL FORMULAS ------------------------------------------------------------

# Main model: distance, biome and their interaction, with random intercepts
# for species and site
lmm_formula <- log_trait ~ distance_to_boundary_km * biome +
  (1 | species_f) + (1 | location_f) + (1 | site_f)

# The same model written so that each biome gets its own slope coefficient
# (biome / distance = a slope within each biome, gives the same fit but used to report the slope on each side with its CI)
lmm_formula_slopes <- log_trait ~ biome / distance_to_boundary_km +
  (1 | species_f) + (1 | location_f) + (1 | site_f)

# Distance classes of the correlogram (step .2 of each trait)
correlogram_breaks <- seq(0, correlogram_max_km, by = correlogram_step_km)

# 3. PLANT HEIGHT --------------------------------------------------------------

## 3.1. Linear mixed model -----------------------------------------------------

# Fitted model with REML
lmm_ph <- lmer(lmm_formula, data = plant_height_data, REML = TRUE)
print(summary(lmm_ph))
# REML criterion at convergence: 48494.6
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -6.9595 -0.5296  0.0234  0.5770  5.4529 
# 
# Random effects:
#   Groups     Name        Variance Std.Dev.
# species_f  (Intercept) 0.8234   0.9074  
# site_f     (Intercept) 0.1464   0.3826  
# location_f (Intercept) 0.4714   0.6866  
# Residual               0.2773   0.5266  
# Number of obs: 29601, groups:  species_f, 441; site_f, 300; location_f, 149
# 
# Fixed effects:
#   Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                         -1.487e+00  1.531e-01  1.621e+02  -9.713   <2e-16 ***
#   distance_to_boundary_km             -2.549e-04  3.337e-04  1.548e+02  -0.764    0.446    
# biometundra                         -7.764e-02  1.740e-01  1.341e+02  -0.446    0.656    
# distance_to_boundary_km:biometundra -7.652e-05  3.956e-04  1.407e+02  -0.193    0.847    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) dst___ bmtndr
# dstnc_t_bn_  0.558              
# biometundra -0.804 -0.483       
# dstnc_t_b_: -0.471 -0.843  0.265

# Tests of the fixed effects (F-tests with Satterthwaite degrees of freedom).
# The key rows: 'biome' (jump at the treeline) and
# 'distance_to_boundary_km:biome' (difference in slope between the two sides)
anova_ph <- anova(lmm_ph)
print(anova_ph)
# Type III Analysis of Variance Table with Satterthwaite's method
#                                Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# distance_to_boundary_km       0.60737 0.60737     1 141.51  2.1903 0.1411
# biome                         0.05523 0.05523     1 134.14  0.1992 0.6561
# distance_to_boundary_km:biome 0.01038 0.01038     1 140.70  0.0374 0.8469 

# Coefficients of the main model
coefs_ph <- summary(lmm_ph)$coefficients

# Slope on each side of the treeline, with 95% CIs (same model, written with
# one slope per biome)
lmm_slopes_ph <- lmer(lmm_formula_slopes, data = plant_height_data, REML = TRUE)
slopes_ph <- summary(lmm_slopes_ph)$coefficients
slope_boreal_ph <- slopes_ph["biomeboreal:distance_to_boundary_km", ]
slope_tundra_ph <- slopes_ph["biometundra:distance_to_boundary_km", ]
print(slopes_ph)
# Estimate   Std. Error       df    t value     Pr(>|t|)
# (Intercept)                         -1.4870321450 0.1530897916 162.0924 -9.7134638 7.282500e-18
# biometundra                         -0.0776365127 0.1739669756 134.1414 -0.4462716 6.561207e-01
# biomeboreal:distance_to_boundary_km -0.0002548562 0.0003336889 154.7862 -0.7637539 4.461765e-01
# biometundra:distance_to_boundary_km -0.0003313766 0.0002129542 113.8764 -1.5560935 1.224608e-01

# Check random effects (i.e. how much varaition is due to species and sites)
print(VarCorr(lmm_ph))
# Groups     Name        Std.Dev.
# species_f  (Intercept) 0.90744 
# site_f     (Intercept) 0.38259 
# location_f (Intercept) 0.68655 
# Residual               0.52659 

# Check results
cat("\nPlant height:",
    "\n  boreal slope:", sprintf("%+.3f", 100 * (exp(slope_boreal_ph[["Estimate"]]) - 1)),
    "% per km, p =", signif(slope_boreal_ph[["Pr(>|t|)"]], 3),
    "\n  tundra slope:", sprintf("%+.3f", 100 * (exp(slope_tundra_ph[["Estimate"]]) - 1)),
    "% per km, p =", signif(slope_tundra_ph[["Pr(>|t|)"]], 3),
    "\n  slope difference: p =", signif(anova_ph["distance_to_boundary_km:biome", "Pr(>F)"], 3),
    "\n  jump at the treeline:", sprintf("%+.1f", 100 * (exp(coefs_ph["biometundra", "Estimate"]) - 1)),
    "%, p =", signif(coefs_ph["biometundra", "Pr(>|t|)"], 3), "\n")
# Plant height: 
#   boreal slope: -0.025 % per km, p = 0.446 
# tundra slope: -0.033 % per km, p = 0.122 
# slope difference: p = 0.847 
# jump at the treeline: -7.5 %, p = 0.656 

## 3.2. Spatial autocorrelation of the residuals -------------------------------

# Get residuals
plant_height_data$resid_lmm <- residuals(lmm_ph)

# Mean residual per site, with the site's projected coordinates
sites_ph <- plant_height_data |>
  group_by(site_name) |>
  summarise(x_km = first(x_km),
            y_km = first(y_km),
            resid_lmm = mean(resid_lmm),
            n_records = n(),
            .groups = "drop")

# Distances between all pairs of sites (km)
site_dist_ph <- as.matrix(dist(cbind(sites_ph$x_km, sites_ph$y_km)))

# Moran's I over all sites (inverse-distance weights)
weights_ph <- 1 / site_dist_ph
weights_ph[!is.finite(weights_ph)] <- 0   # a site's distance to itself
moran_ph <- Moran.I(sites_ph$resid_lmm, weights_ph)

# Check Moran's I
cat("\nPlant height - Moran's I of the residuals (", nrow(sites_ph), " sites): I =",
    round(moran_ph$observed, 3), ", p =", signif(moran_ph$p.value, 3), "\n")
# Plant height - Moran's I of the residuals ( 300  sites): I = 0.003 , p = 0.861 

# Correlogram: Moran's I among pairs of sites within each distance class
correlogram_ph <- data.frame()
for (i in seq_len(length(correlogram_breaks) - 1)) {
  
  # Weight 1 for pairs of sites within this distance class, 0 otherwise
  class_weights <- (site_dist_ph > correlogram_breaks[i] &
                      site_dist_ph <= correlogram_breaks[i + 1]) * 1
  
  # Skip classes without any pairs of sites
  if (sum(class_weights) == 0) next
  
  moran_class <- Moran.I(sites_ph$resid_lmm, class_weights)
  
  correlogram_ph <- bind_rows(correlogram_ph, data.frame(
    Trait = "PlantHeight",
    distance_km = (correlogram_breaks[i] + correlogram_breaks[i + 1]) / 2,
    n_pairs = sum(class_weights) / 2,
    morans_I = moran_class$observed,
    p = moran_class$p.value))
}

## 3.3. Model validation -------------------------------------------------------

# Plot Residuals vs fitted values and Q-Q plot 
png(here("figures", "RQ3_PlantHeight_validation_LMM.png"),
    width = 12, height = 5, units = "in", res = 200)
par(mfrow = c(1, 2))

plot(fitted(lmm_ph), plant_height_data$resid_lmm, pch = 16, cex = 0.4, col = rgb(0, 0, 0, 0.2),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Plant height - LMM: residuals vs fitted")
abline(h = 0, col = "red", lty = 2)

qqnorm(plant_height_data$resid_lmm, pch = 16, cex = 0.4, main = "Plant height - LMM: Q-Q plot")
qqline(plant_height_data$resid_lmm, col = "red")

par(mfrow = c(1, 1))
dev.off()

## 3.4. Figures ----------------------------------------------------------------

# Fitted lines (fixed effects only: re.form = NA)
line_ph <- bind_rows(data.frame(distance_to_boundary_km = seq(min(plant_height_data$distance_to_boundary_km), 0, length.out = 100),
                                biome = factor("boreal", levels = c("boreal", "tundra"))),
                     data.frame(distance_to_boundary_km = seq(0, max(plant_height_data$distance_to_boundary_km), length.out = 100),
                                biome = factor("tundra", levels = c("boreal", "tundra"))))

# Predict values
line_ph$predicted <- predict(lmm_ph, newdata = line_ph, re.form = NA)

# Labels: slope on each side as % change per km
label_boreal_ph <- paste0("Boreal: ",
                          sprintf("%+.3f", 100 * (exp(slope_boreal_ph[["Estimate"]]) - 1)),
                          "% per km\np ", ifelse(slope_boreal_ph[["Pr(>|t|)"]] < 0.001, "< 0.001",
                                                 paste0("= ", sprintf("%.3f", slope_boreal_ph[["Pr(>|t|)"]]))))
label_tundra_ph <- paste0("Tundra: ",
                          sprintf("%+.3f", 100 * (exp(slope_tundra_ph[["Estimate"]]) - 1)),
                          "% per km\np ", ifelse(slope_tundra_ph[["Pr(>|t|)"]] < 0.001, "< 0.001",
                                                 paste0("= ", sprintf("%.3f", slope_tundra_ph[["Pr(>|t|)"]]))))

# Title: slope difference and jump at the treeline
p_slope_diff_ph <- anova_ph["distance_to_boundary_km:biome", "Pr(>F)"]
p_jump_ph <- coefs_ph["biometundra", "Pr(>|t|)"]
title_lmm_ph <- paste0("Slope difference: p ", ifelse(p_slope_diff_ph < 0.001, "< 0.001",
                                                      paste0("= ", sprintf("%.3f", p_slope_diff_ph))),
                       "; jump at treeline: p ", ifelse(p_jump_ph < 0.001, "< 0.001",
                                                        paste0("= ", sprintf("%.3f", p_jump_ph))))

# Plot fitted lines
(plot_lmm_ph <- ggplot() +
  geom_point(data = plant_height_data, aes(x = distance_to_boundary_km, y = log_trait, colour = biome),
             alpha = 0.2, size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray30", linewidth = 0.7) +
  geom_line(data = line_ph, aes(x = distance_to_boundary_km, y = predicted, group = biome),
            colour = fit_line_colour, linewidth = 1.5) +
  annotate("text", x = -Inf, y = Inf, label = label_boreal_ph,
           hjust = -0.03, vjust = 1.2, size = fig_annot_text / .pt, fontface = "bold") +
  annotate("text", x = Inf, y = Inf, label = label_tundra_ph,
           hjust = 1.03, vjust = 1.2, size = fig_annot_text / .pt, fontface = "bold") +
  scale_colour_manual(values = biome_colours, labels = biome_labels, name = "Biome") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +   # room for the labels
  labs(x = "Distance to Biome Boundary (km)", y = "Plant Height (m, log scale)", title = title_lmm_ph) +
  theme_rq3 +
  theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))))

# Save to file
ggsave(here("figures", "FigureS_RQ3_PlantHeight_LMM.png"),
       plot = plot_lmm_ph, width = 10, height = 6, dpi = 600)

# Plot correlogram
(plot_correlogram_ph <- ggplot(correlogram_ph, aes(x = distance_km, y = morans_I)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_line(linewidth = 0.8) +
  geom_point(aes(shape = p < 0.05), size = 3) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), name = "p < 0.05") +
  labs(x = "Distance between sites (km)", y = "Moran's I of residuals", title = "Plant height") +
  theme_rq3)

# Save correlogram for supplementary information
ggsave(here("figures", "FigureS_RQ3_PlantHeight_correlogram_LMM.png"),
       plot = plot_correlogram_ph, width = 10, height = 6, dpi = 300)

# Collect results in single df
results_lmm_ph <- data.frame(Trait = "PlantHeight",
                             n_records = nrow(plant_height_data),
                             n_species = n_distinct(plant_height_data$StandardSpeciesName),
                             n_sites_boreal = n_distinct(plant_height_data$site_name[plant_height_data$biome == "boreal"]),
                             n_sites_tundra = n_distinct(plant_height_data$site_name[plant_height_data$biome == "tundra"]),
                             # Slope on each side (log scale, per km)
                             slope_boreal = slope_boreal_ph[["Estimate"]],
                             se_boreal = slope_boreal_ph[["Std. Error"]],
                             df_boreal = slope_boreal_ph[["df"]],
                             p_boreal = slope_boreal_ph[["Pr(>|t|)"]],
                             slope_tundra = slope_tundra_ph[["Estimate"]],
                             se_tundra = slope_tundra_ph[["Std. Error"]],
                             df_tundra = slope_tundra_ph[["df"]],
                             p_tundra = slope_tundra_ph[["Pr(>|t|)"]],
                             # Difference in slope (tundra minus boreal) and its F-test
                             slope_difference = coefs_ph["distance_to_boundary_km:biometundra", "Estimate"],
                             F_slope_difference = anova_ph["distance_to_boundary_km:biome", "F value"],
                             p_slope_difference = anova_ph["distance_to_boundary_km:biome", "Pr(>F)"],
                             # Jump at the treeline (tundra minus boreal at 0 km, log scale)
                             jump = coefs_ph["biometundra", "Estimate"],
                             se_jump = coefs_ph["biometundra", "Std. Error"],
                             df_jump = coefs_ph["biometundra", "df"],
                             p_jump = coefs_ph["biometundra", "Pr(>|t|)"],
                             # Spatial autocorrelation of the residuals
                             morans_I = moran_ph$observed,
                             morans_I_p = moran_ph$p.value)

# 4. SLA -----------------------------------------------------------------------

## 4.1. Linear mixed model -----------------------------------------------------

# Run model
lmm_sla <- lmer(lmm_formula, data = sla_data, REML = TRUE)

# Get summary
print(summary(lmm_sla))
# REML criterion at convergence: 2316.6
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -8.4712 -0.5582  0.0062  0.5624  5.9427 
# 
# Random effects:
#   Groups     Name        Variance Std.Dev.
# species_f  (Intercept) 0.25974  0.5097  
# site_f     (Intercept) 0.01830  0.1353  
# location_f (Intercept) 0.03541  0.1882  
# Residual               0.05945  0.2438  
# Number of obs: 18017, groups:  species_f, 401; site_f, 188; location_f, 104
# 
# Fixed effects:
#   Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                          3.020e+00  7.354e-02  1.557e+02  41.063  < 2e-16 ***
#   distance_to_boundary_km             -3.155e-05  1.521e-04  1.026e+02  -0.207 0.836099    
# biometundra                         -2.888e-01  7.748e-02  1.064e+02  -3.727 0.000312 ***
#   distance_to_boundary_km:biometundra  4.857e-05  1.645e-04  9.706e+01   0.295 0.768504    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) dst___ bmtndr
# dstnc_t_bn_  0.695              
# biometundra -0.835 -0.664       
# dstnc_t_b_: -0.646 -0.924  0.521

# Tests of the fixed effects
anova_sla <- anova(lmm_sla)
print(anova_sla)

# Coefficients of the main model
coefs_sla <- summary(lmm_sla)$coefficients
# Type III Analysis of Variance Table with Satterthwaite's method
#                                Sum Sq Mean Sq NumDF   DenDF F value    Pr(>F)    
# distance_to_boundary_km       0.00046 0.00046     1  97.842  0.0078 0.9298546    
# biome                         0.82595 0.82595     1 106.358 13.8923 0.0003121 ***
# distance_to_boundary_km:biome 0.00518 0.00518     1  97.059  0.0871 0.7685042    

# Slope on each side of the treeline, with 95% CIs 
lmm_slopes_sla <- lmer(lmm_formula_slopes, data = sla_data, REML = TRUE)
slopes_sla <- summary(lmm_slopes_sla)$coefficients
slope_boreal_sla <- slopes_sla["biomeboreal:distance_to_boundary_km", ]
slope_tundra_sla <- slopes_sla["biometundra:distance_to_boundary_km", ]
print(slopes_sla)
# Estimate   Std. Error        df    t value     Pr(>|t|)
# (Intercept)                          3.019863e+00 7.354244e-02 155.74712 41.0628632 1.894591e-85
# biometundra                         -2.888028e-01 7.748453e-02 106.35762 -3.7272315 3.120988e-04
# biomeboreal:distance_to_boundary_km -3.154871e-05 1.521065e-04 102.56526 -0.2074120 8.360994e-01
# biometundra:distance_to_boundary_km  1.701933e-05 6.287643e-05  73.57781  0.2706791 7.873958e-01

# Random effects: how much variation is due to species, locations and sites
print(VarCorr(lmm_sla))
# Groups     Name        Std.Dev.
# species_f  (Intercept) 0.50965 
# site_f     (Intercept) 0.13529 
# location_f (Intercept) 0.18817 
# Residual               0.24383

# Get summary
cat("\nSLA:",
    "\n  boreal slope:", sprintf("%+.3f", 100 * (exp(slope_boreal_sla[["Estimate"]]) - 1)),
    "% per km, p =", signif(slope_boreal_sla[["Pr(>|t|)"]], 3),
    "\n  tundra slope:", sprintf("%+.3f", 100 * (exp(slope_tundra_sla[["Estimate"]]) - 1)),
    "% per km, p =", signif(slope_tundra_sla[["Pr(>|t|)"]], 3),
    "\n  slope difference: p =", signif(anova_sla["distance_to_boundary_km:biome", "Pr(>F)"], 3),
    "\n  jump at the treeline:", sprintf("%+.1f", 100 * (exp(coefs_sla["biometundra", "Estimate"]) - 1)),
    "%, p =", signif(coefs_sla["biometundra", "Pr(>|t|)"], 3), "\n")
# SLA: 
#   boreal slope: -0.003 % per km, p = 0.836 
# tundra slope: +0.002 % per km, p = 0.787 
# slope difference: p = 0.769 
# jump at the treeline: -25.1 %, p = 0.000312 

## 4.2. Spatial autocorrelation of the residuals -------------------------------

# Get residuals
sla_data$resid_lmm <- residuals(lmm_sla)

# Mean residual per site, with the site's projected coordinates
sites_sla <- sla_data |>
  group_by(site_name) |>
  summarise(x_km = first(x_km),
            y_km = first(y_km),
            resid_lmm = mean(resid_lmm),
            n_records = n(),
            .groups = "drop")

# Distances between all pairs of sites (km)
site_dist_sla <- as.matrix(dist(cbind(sites_sla$x_km, sites_sla$y_km)))

# Moran's I over all sites (inverse-distance weights)
weights_sla <- 1 / site_dist_sla
weights_sla[!is.finite(weights_sla)] <- 0   # a site's distance to itself
moran_sla <- Moran.I(sites_sla$resid_lmm, weights_sla)

# Get the spatial autocorrelation check
cat("\nSLA - Moran's I of the residuals (", nrow(sites_sla), " sites): I =",
    round(moran_sla$observed, 3), ", p =", signif(moran_sla$p.value, 3), "\n")
# SLA - Moran's I of the residuals ( 188  sites): I = -0.104 , p = 0.0181

# Correlogram: Moran's I among pairs of sites within each distance class
correlogram_sla <- data.frame()
for (i in seq_len(length(correlogram_breaks) - 1)) {
  
  # Weight 1 for pairs of sites within this distance class, 0 otherwise
  class_weights <- (site_dist_sla > correlogram_breaks[i] &
                      site_dist_sla <= correlogram_breaks[i + 1]) * 1
  
  # Skip classes without any pairs of sites
  if (sum(class_weights) == 0) next
  
  moran_class <- Moran.I(sites_sla$resid_lmm, class_weights)
  
  correlogram_sla <- bind_rows(correlogram_sla, data.frame(
    Trait = "SLA",
    distance_km = (correlogram_breaks[i] + correlogram_breaks[i + 1]) / 2,
    n_pairs = sum(class_weights) / 2,
    morans_I = moran_class$observed,
    p = moran_class$p.value))
}

## 4.3. Model validation -------------------------------------------------------

# Residuals vs fitted values Q-Q plot
png(here("figures", "RQ3_SLA_validation_LMM.png"),
    width = 12, height = 5, units = "in", res = 200)
par(mfrow = c(1, 2))

plot(fitted(lmm_sla), sla_data$resid_lmm, pch = 16, cex = 0.4, col = rgb(0, 0, 0, 0.2),
     xlab = "Fitted values", ylab = "Residuals",
     main = "SLA - LMM: residuals vs fitted")
abline(h = 0, col = "red", lty = 2)

qqnorm(sla_data$resid_lmm, pch = 16, cex = 0.4, main = "SLA - LMM: Q-Q plot")
qqline(sla_data$resid_lmm, col = "red")

par(mfrow = c(1, 1))
dev.off()

## 4.4. Figures ----------------------------------------------------------------

# Get fitted lines (fixed effects only: re.form = NA)
line_sla <- bind_rows(data.frame(distance_to_boundary_km = seq(min(sla_data$distance_to_boundary_km), 0, length.out = 100),
             biome = factor("boreal", levels = c("boreal", "tundra"))),
             data.frame(distance_to_boundary_km = seq(0, max(sla_data$distance_to_boundary_km), length.out = 100),
             biome = factor("tundra", levels = c("boreal", "tundra"))))
line_sla$predicted <- predict(lmm_sla, newdata = line_sla, re.form = NA)

# Labels: slope on each side as % change per km
label_boreal_sla <- paste0("Boreal: ",
                           sprintf("%+.3f", 100 * (exp(slope_boreal_sla[["Estimate"]]) - 1)),
                           "% per km\np ", ifelse(slope_boreal_sla[["Pr(>|t|)"]] < 0.001, "< 0.001",
                                                  paste0("= ", sprintf("%.3f", slope_boreal_sla[["Pr(>|t|)"]]))))
label_tundra_sla <- paste0("Tundra: ",
                           sprintf("%+.3f", 100 * (exp(slope_tundra_sla[["Estimate"]]) - 1)),
                           "% per km\np ", ifelse(slope_tundra_sla[["Pr(>|t|)"]] < 0.001, "< 0.001",
                                                  paste0("= ", sprintf("%.3f", slope_tundra_sla[["Pr(>|t|)"]]))))

# Title: slope difference and jump at the treeline
p_slope_diff_sla <- anova_sla["distance_to_boundary_km:biome", "Pr(>F)"]
p_jump_sla <- coefs_sla["biometundra", "Pr(>|t|)"]
title_lmm_sla <- paste0("Slope difference: p ", ifelse(p_slope_diff_sla < 0.001, "< 0.001",
                                                       paste0("= ", sprintf("%.3f", p_slope_diff_sla))),
                        "; jump at treeline: p ", ifelse(p_jump_sla < 0.001, "< 0.001",
                                                         paste0("= ", sprintf("%.3f", p_jump_sla))))

# Plot figure
(plot_lmm_sla <- ggplot() +
  geom_point(data = sla_data, aes(x = distance_to_boundary_km, y = log_trait, colour = biome),
             alpha = 0.2, size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray30", linewidth = 0.7) +
  geom_line(data = line_sla, aes(x = distance_to_boundary_km, y = predicted, group = biome),
            colour = fit_line_colour, linewidth = 1.5) +
  annotate("text", x = -Inf, y = Inf, label = label_boreal_sla,
           hjust = -0.03, vjust = 1.2, size = fig_annot_text / .pt, fontface = "bold") +
  annotate("text", x = Inf, y = Inf, label = label_tundra_sla,
           hjust = 1.03, vjust = 1.2, size = fig_annot_text / .pt, fontface = "bold") +
  scale_colour_manual(values = biome_colours, labels = biome_labels, name = "Biome") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +   # room for the labels
  labs(x = "Distance to Biome Boundary (km)", y = expression(bold(paste("SLA (mm"^2, " mg"^-1, ", log scale)"))), title = title_lmm_sla) +
  theme_rq3 +
  theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))))

# Save figure
ggsave(here("figures", "FigureS_RQ3_SLA_LMM.png"),
       plot = plot_lmm_sla, width = 10, height = 6, dpi = 600)

# Plot correlogram
(plot_correlogram_sla <- ggplot(correlogram_sla, aes(x = distance_km, y = morans_I)) +
    geom_hline(yintercept = 0, colour = "grey50") +
    geom_line(linewidth = 0.8) +
    geom_point(aes(shape = p < 0.05), size = 3) +
    scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), name = "p < 0.05") +
    labs(x = "Distance between sites (km)", y = "Moran's I of residuals", title = "SLA") +
    theme_rq3)

# Save figure
ggsave(here("figures", "FigureS_RQ3_SLA_correlogram_LMM.png"),
       plot = plot_correlogram_sla, width = 10, height = 6, dpi = 300)

# Collect all results into single df
results_lmm_sla <- data.frame(Trait = "SLA",
                              n_records = nrow(sla_data),
                              n_species = n_distinct(sla_data$StandardSpeciesName),
                              n_sites_boreal = n_distinct(sla_data$site_name[sla_data$biome == "boreal"]),
                              n_sites_tundra = n_distinct(sla_data$site_name[sla_data$biome == "tundra"]),
                              # Slope on each side (log scale, per km)
                              slope_boreal = slope_boreal_sla[["Estimate"]],
                              se_boreal = slope_boreal_sla[["Std. Error"]],
                              df_boreal = slope_boreal_sla[["df"]],
                              p_boreal = slope_boreal_sla[["Pr(>|t|)"]],
                              slope_tundra = slope_tundra_sla[["Estimate"]],
                              se_tundra = slope_tundra_sla[["Std. Error"]],
                              df_tundra = slope_tundra_sla[["df"]],
                              p_tundra = slope_tundra_sla[["Pr(>|t|)"]],
                              # Difference in slope (tundra minus boreal) and its F-test
                              slope_difference = coefs_sla["distance_to_boundary_km:biometundra", "Estimate"],
                              F_slope_difference = anova_sla["distance_to_boundary_km:biome", "F value"],
                              p_slope_difference = anova_sla["distance_to_boundary_km:biome", "Pr(>F)"],
                              # Jump at the treeline (tundra minus boreal at 0 km, log scale)
                              jump = coefs_sla["biometundra", "Estimate"],
                              se_jump = coefs_sla["biometundra", "Std. Error"],
                              df_jump = coefs_sla["biometundra", "df"],
                              p_jump = coefs_sla["biometundra", "Pr(>|t|)"],
                              # Spatial autocorrelation of the residuals
                              morans_I = moran_sla$observed,
                              morans_I_p = moran_sla$p.value)

# 5. LEAF N --------------------------------------------------------------------

## 5.1. Linear mixed model -----------------------------------------------------

# Fit model
lmm_leafn <- lmer(lmm_formula, data = leafn_data, REML = TRUE)

# Get model summary
print(summary(lmm_leafn))
# REML criterion at convergence: 1802.2
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -6.1491 -0.4310  0.0473  0.4834  6.9429 
# 
# Random effects:
#   Groups     Name        Variance Std.Dev.
# species_f  (Intercept) 0.117782 0.34319 
# site_f     (Intercept) 0.038048 0.19506 
# location_f (Intercept) 0.009662 0.09829 
# Residual               0.061887 0.24877 
# Number of obs: 7207, groups:  species_f, 423; site_f, 210; location_f, 159
# 
# Fixed effects:
#   Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                          2.985e+00  5.040e-02  3.977e+02  59.220   <2e-16 ***
#   distance_to_boundary_km             -1.946e-04  9.464e-05  2.013e+02  -2.056   0.0411 *  
#   biometundra                          2.243e-02  5.577e-02  2.009e+02   0.402   0.6879    
# distance_to_boundary_km:biometundra  2.528e-04  1.098e-04  1.855e+02   2.304   0.0223 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) dst___ bmtndr
# dstnc_t_bn_  0.632              
# biometundra -0.715 -0.539       
# dstnc_t_b_: -0.577 -0.869  0.304

# Tests of the fixed effects 
anova_leafn <- anova(lmm_leafn)
print(anova_leafn)
# Type III Analysis of Variance Table with Satterthwaite's method
#                                Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# distance_to_boundary_km       0.09775 0.09775     1 184.98  1.5795 0.21042  
# biome                         0.01001 0.01001     1 200.94  0.1618 0.68794  
# distance_to_boundary_km:biome 0.32845 0.32845     1 185.49  5.3074 0.02234 *

# Coefficients of the main model
coefs_leafn <- summary(lmm_leafn)$coefficients

# Slope on each side of the treeline, with 95% CIs
lmm_slopes_leafn <- lmer(lmm_formula_slopes, data = leafn_data, REML = TRUE)
slopes_leafn <- summary(lmm_slopes_leafn)$coefficients
slope_boreal_leafn <- slopes_leafn["biomeboreal:distance_to_boundary_km", ]
slope_tundra_leafn <- slopes_leafn["biometundra:distance_to_boundary_km", ]
print(slopes_leafn)
# Estimate   Std. Error       df    t value      Pr(>|t|)
# (Intercept)                          2.984782e+00 5.040176e-02 397.6723 59.2198026 2.334468e-199
# biometundra                          2.243146e-02 5.576759e-02 200.9369  0.4022311  6.879415e-01
# biomeboreal:distance_to_boundary_km -1.945590e-04 9.463998e-05 201.2983 -2.0557802  4.109432e-02
# biometundra:distance_to_boundary_km  5.829354e-05 5.426239e-05 141.3640  1.0742898  2.845231e-01

# Random effects: how much variation is due to species, locations and sites
print(VarCorr(lmm_leafn))
# Groups     Name        Std.Dev.
# species_f  (Intercept) 0.343194
# site_f     (Intercept) 0.195060
# location_f (Intercept) 0.098294
# Residual               0.248770

# Print summary
cat("\nLeaf N:",
    "\n  boreal slope:", sprintf("%+.3f", 100 * (exp(slope_boreal_leafn[["Estimate"]]) - 1)),
    "% per km, p =", signif(slope_boreal_leafn[["Pr(>|t|)"]], 3),
    "\n  tundra slope:", sprintf("%+.3f", 100 * (exp(slope_tundra_leafn[["Estimate"]]) - 1)),
    "% per km, p =", signif(slope_tundra_leafn[["Pr(>|t|)"]], 3),
    "\n  slope difference: p =", signif(anova_leafn["distance_to_boundary_km:biome", "Pr(>F)"], 3),
    "\n  jump at the treeline:", sprintf("%+.1f", 100 * (exp(coefs_leafn["biometundra", "Estimate"]) - 1)),
    "%, p =", signif(coefs_leafn["biometundra", "Pr(>|t|)"], 3), "\n")
# Leaf N: 
# boreal slope: -0.019 % per km, p = 0.0411 
# tundra slope: +0.006 % per km, p = 0.285 
# slope difference: p = 0.0223 
# jump at the treeline: +2.3 %, p = 0.688 

## 5.2. Spatial autocorrelation of the residuals -------------------------------

# Get residuals
leafn_data$resid_lmm <- residuals(lmm_leafn)

# Mean residual per site, with the site's projected coordinates
sites_leafn <- leafn_data |>
  group_by(site_name) |>
  summarise(x_km = first(x_km),
            y_km = first(y_km),
            resid_lmm = mean(resid_lmm),
            n_records = n(),
            .groups = "drop")

# Distances between all pairs of sites (km)
site_dist_leafn <- as.matrix(dist(cbind(sites_leafn$x_km, sites_leafn$y_km)))

# Moran's I over all sites (inverse-distance weights)
weights_leafn <- 1 / site_dist_leafn
weights_leafn[!is.finite(weights_leafn)] <- 0   # a site's distance to itself
moran_leafn <- Moran.I(sites_leafn$resid_lmm, weights_leafn)

# Print Moran's I
cat("\nLeaf N - Moran's I of the residuals (", nrow(sites_leafn), " sites): I =",
    round(moran_leafn$observed, 3), ", p =", signif(moran_leafn$p.value, 3), "\n")
# Leaf N - Moran's I of the residuals ( 210  sites): I = 0.022 , p = 0.475 

# Correlogram: Moran's I among pairs of sites within each distance class
correlogram_leafn <- data.frame()
for (i in seq_len(length(correlogram_breaks) - 1)) {
  
  # Weight 1 for pairs of sites within this distance class, 0 otherwise
  class_weights <- (site_dist_leafn > correlogram_breaks[i] &
                      site_dist_leafn <= correlogram_breaks[i + 1]) * 1
  
  # Skip classes without any pairs of sites
  if (sum(class_weights) == 0) next
  
  moran_class <- Moran.I(sites_leafn$resid_lmm, class_weights)
  
  correlogram_leafn <- bind_rows(correlogram_leafn, data.frame(
    Trait = "LeafN",
    distance_km = (correlogram_breaks[i] + correlogram_breaks[i + 1]) / 2,
    n_pairs = sum(class_weights) / 2,
    morans_I = moran_class$observed,
    p = moran_class$p.value))
}

## 5.3. Model validation ---------------------------------------------------------

# Plot residuals vs fitted values and Q-Q plot
png(here("figures", "RQ3_LeafN_validation_LMM.png"),
    width = 12, height = 5, units = "in", res = 200)
par(mfrow = c(1, 2))

plot(fitted(lmm_leafn), leafn_data$resid_lmm, pch = 16, cex = 0.4, col = rgb(0, 0, 0, 0.2),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Leaf N - LMM: residuals vs fitted")
abline(h = 0, col = "red", lty = 2)

qqnorm(leafn_data$resid_lmm, pch = 16, cex = 0.4, main = "Leaf N - LMM: Q-Q plot")
qqline(leafn_data$resid_lmm, col = "red")

par(mfrow = c(1, 1))
dev.off()

## 5.4. Figures ----------------------------------------------------------------

# Fitted lines (fixed effects only: re.form = NA)
line_leafn <- bind_rows(data.frame(distance_to_boundary_km = seq(min(leafn_data$distance_to_boundary_km), 0, length.out = 100),
                                   biome = factor("boreal", levels = c("boreal", "tundra"))),
                        data.frame(distance_to_boundary_km = seq(0, max(leafn_data$distance_to_boundary_km), length.out = 100),
                                   biome = factor("tundra", levels = c("boreal", "tundra"))))
line_leafn$predicted <- predict(lmm_leafn, newdata = line_leafn, re.form = NA)

# Labels: slope on each side as % change per km
label_boreal_leafn <- paste0("Boreal: ",
                             sprintf("%+.3f", 100 * (exp(slope_boreal_leafn[["Estimate"]]) - 1)),
                             "% per km\np ", ifelse(slope_boreal_leafn[["Pr(>|t|)"]] < 0.001, "< 0.001",
                                                    paste0("= ", sprintf("%.3f", slope_boreal_leafn[["Pr(>|t|)"]]))))
label_tundra_leafn <- paste0("Tundra: ",
                             sprintf("%+.3f", 100 * (exp(slope_tundra_leafn[["Estimate"]]) - 1)),
                             "% per km\np ", ifelse(slope_tundra_leafn[["Pr(>|t|)"]] < 0.001, "< 0.001",
                                                    paste0("= ", sprintf("%.3f", slope_tundra_leafn[["Pr(>|t|)"]]))))

# Title: slope difference and jump at the treeline
p_slope_diff_leafn <- anova_leafn["distance_to_boundary_km:biome", "Pr(>F)"]
p_jump_leafn <- coefs_leafn["biometundra", "Pr(>|t|)"]
title_lmm_leafn <- paste0("Slope difference: p ", ifelse(p_slope_diff_leafn < 0.001, "< 0.001",
                                                         paste0("= ", sprintf("%.3f", p_slope_diff_leafn))),
                          "; jump at treeline: p ", ifelse(p_jump_leafn < 0.001, "< 0.001",
                                                           paste0("= ", sprintf("%.3f", p_jump_leafn))))

# Plot model fit figure
(plot_lmm_leafn <- ggplot() +
  geom_point(data = leafn_data, aes(x = distance_to_boundary_km, y = log_trait, colour = biome),
             alpha = 0.2, size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray30", linewidth = 0.7) +
  geom_line(data = line_leafn, aes(x = distance_to_boundary_km, y = predicted, group = biome),
            colour = fit_line_colour, linewidth = 1.5) +
  annotate("text", x = -Inf, y = Inf, label = label_boreal_leafn,
           hjust = -0.03, vjust = 1.2, size = fig_annot_text / .pt, fontface = "bold") +
  annotate("text", x = Inf, y = Inf, label = label_tundra_leafn,
           hjust = 1.03, vjust = 1.2, size = fig_annot_text / .pt, fontface = "bold") +
  scale_colour_manual(values = biome_colours, labels = biome_labels, name = "Biome") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +   # room for the labels
  labs(x = "Distance to Biome Boundary (km)", y = "Leaf N (mg/g, log scale)", title = title_lmm_leafn) +
  theme_rq3 +
  theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))))


# Save figure
ggsave(here("figures", "FigureS_RQ3_LeafN_LMM.png"),
       plot = plot_lmm_leafn, width = 10, height = 6, dpi = 600)

# Plot correlogram
(plot_correlogram_leafn <- ggplot(correlogram_leafn, aes(x = distance_km, y = morans_I)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_line(linewidth = 0.8) +
  geom_point(aes(shape = p < 0.05), size = 3) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), name = "p < 0.05") +
  labs(x = "Distance between sites (km)", y = "Moran's I of residuals", title = "Leaf N") +
  theme_rq3)

# Save correlogram
ggsave(here("figures", "FigureS_RQ3_LeafN_correlogram_LMM.png"),
       plot = plot_correlogram_leafn, width = 10, height = 6, dpi = 300)

# Collect all the results into a single df
results_lmm_leafn <- data.frame(Trait = "LeafN",
                                n_records = nrow(leafn_data),
                                n_species = n_distinct(leafn_data$StandardSpeciesName),
                                n_sites_boreal = n_distinct(leafn_data$site_name[leafn_data$biome == "boreal"]),
                                n_sites_tundra = n_distinct(leafn_data$site_name[leafn_data$biome == "tundra"]),
                                # Slope on each side (log scale, per km)
                                slope_boreal = slope_boreal_leafn[["Estimate"]],
                                se_boreal = slope_boreal_leafn[["Std. Error"]],
                                df_boreal = slope_boreal_leafn[["df"]],
                                p_boreal = slope_boreal_leafn[["Pr(>|t|)"]],
                                slope_tundra = slope_tundra_leafn[["Estimate"]],
                                se_tundra = slope_tundra_leafn[["Std. Error"]],
                                df_tundra = slope_tundra_leafn[["df"]],
                                p_tundra = slope_tundra_leafn[["Pr(>|t|)"]],
                                # Difference in slope (tundra minus boreal) and its F-test
                                slope_difference = coefs_leafn["distance_to_boundary_km:biometundra", "Estimate"],
                                F_slope_difference = anova_leafn["distance_to_boundary_km:biome", "F value"],
                                p_slope_difference = anova_leafn["distance_to_boundary_km:biome", "Pr(>F)"],
                                # Jump at the treeline (tundra minus boreal at 0 km, log scale)
                                jump = coefs_leafn["biometundra", "Estimate"],
                                se_jump = coefs_leafn["biometundra", "Std. Error"],
                                df_jump = coefs_leafn["biometundra", "df"],
                                p_jump = coefs_leafn["biometundra", "Pr(>|t|)"],
                                # Spatial autocorrelation of the residuals
                                morans_I = moran_leafn$observed,
                                morans_I_p = moran_leafn$p.value)

# 6. SEED MASS -----------------------------------------------------------------
# N.B.: very few boreal seed mass records (3 locations in total), so the
# boreal slope and the jump at the treeline are very uncertain.

## 6.1. Linear mixed model -----------------------------------------------------

# Fit linear model
lmm_seedmass <- lmer(lmm_formula, data = seedmass_data, REML = TRUE)

# Get model summary
print(summary(lmm_seedmass))
# REML criterion at convergence: 1810.7
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -4.4472 -0.3335  0.0386  0.3878  5.7802 
# 
# Random effects:
#   Groups     Name        Variance Std.Dev.
# species_f  (Intercept) 1.90326  1.3796  
# site_f     (Intercept) 0.30737  0.5544  
# location_f (Intercept) 0.09141  0.3023  
# Residual               0.29830  0.5462  
# Number of obs: 814, groups:  species_f, 123; site_f, 36; location_f, 20
# 
# Fixed effects:
#   Estimate Std. Error        df t value Pr(>|t|)
# (Intercept)                         -0.953926   0.635853 42.059868  -1.500    0.141
# distance_to_boundary_km              0.001290   0.004298 33.445777   0.300    0.766
# biometundra                          0.306951   0.646379 37.355293   0.475    0.638
# distance_to_boundary_km:biometundra -0.001089   0.004313 33.433683  -0.252    0.802
# 
# Correlation of Fixed Effects:
#   (Intr) dst___ bmtndr
# dstnc_t_bn_  0.671              
# biometundra -0.952 -0.662       
# dstnc_t_b_: -0.669 -0.996  0.652

# Tests of the fixed effects
anova_seedmass <- anova(lmm_seedmass)
print(anova_seedmass)
# Type III Analysis of Variance Table with Satterthwaite's method
#                                 Sum Sq  Mean Sq NumDF  DenDF F value Pr(>F)
# distance_to_boundary_km       0.035611 0.035611     1 33.392  0.1194 0.7319
# biome                         0.067269 0.067269     1 37.355  0.2255 0.6376
# distance_to_boundary_km:biome 0.019015 0.019015     1 33.434  0.0637 0.8022

# Coefficients of the main model
coefs_seedmass <- summary(lmm_seedmass)$coefficients

# Slope on each side of the treeline, with 95% CIs
lmm_slopes_seedmass <- lmer(lmm_formula_slopes, data = seedmass_data, REML = TRUE)
slopes_seedmass <- summary(lmm_slopes_seedmass)$coefficients
slope_boreal_seedmass <- slopes_seedmass["biomeboreal:distance_to_boundary_km", ]
slope_tundra_seedmass <- slopes_seedmass["biometundra:distance_to_boundary_km", ]
print(slopes_seedmass)

# Random effects: how much variation is due to species, locations and sites
print(VarCorr(lmm_seedmass))
# Groups     Name        Std.Dev.
# species_f  (Intercept) 1.37959 
# site_f     (Intercept) 0.55441 
# location_f (Intercept) 0.30235 
# Residual               0.54617 

# Print quick summary
cat("\nSeed mass:",
    "\n  boreal slope:", sprintf("%+.3f", 100 * (exp(slope_boreal_seedmass[["Estimate"]]) - 1)),
    "% per km, p =", signif(slope_boreal_seedmass[["Pr(>|t|)"]], 3),
    "\n  tundra slope:", sprintf("%+.3f", 100 * (exp(slope_tundra_seedmass[["Estimate"]]) - 1)),
    "% per km, p =", signif(slope_tundra_seedmass[["Pr(>|t|)"]], 3),
    "\n  slope difference: p =", signif(anova_seedmass["distance_to_boundary_km:biome", "Pr(>F)"], 3),
    "\n  jump at the treeline:", sprintf("%+.1f", 100 * (exp(coefs_seedmass["biometundra", "Estimate"]) - 1)),
    "%, p =", signif(coefs_seedmass["biometundra", "Pr(>|t|)"], 3), "\n")
# Seed mass: 
#   boreal slope: +0.129 % per km, p = 0.766 
# tundra slope: +0.020 % per km, p = 0.605 
# slope difference: p = 0.802 
# jump at the treeline: +35.9 %, p = 0.638 

## 6.2. Spatial autocorrelation of the residuals -------------------------------

# Get residuals
seedmass_data$resid_lmm <- residuals(lmm_seedmass)

# Mean residual per site, with the site's projected coordinates
sites_seedmass <- seedmass_data |>
  group_by(site_name) |>
  summarise(x_km = first(x_km),
            y_km = first(y_km),
            resid_lmm = mean(resid_lmm),
            n_records = n(),
            .groups = "drop")

# Distances between all pairs of sites (km)
site_dist_seedmass <- as.matrix(dist(cbind(sites_seedmass$x_km, sites_seedmass$y_km)))

# Moran's I over all sites (inverse-distance weights)
weights_seedmass <- 1 / site_dist_seedmass
weights_seedmass[!is.finite(weights_seedmass)] <- 0   # a site's distance to itself
moran_seedmass <- Moran.I(sites_seedmass$resid_lmm, weights_seedmass)

# Print out Moran's I summary
cat("\nSeed mass - Moran's I of the residuals (", nrow(sites_seedmass), " sites): I =",
    round(moran_seedmass$observed, 3), ", p =", signif(moran_seedmass$p.value, 3), "\n")
# Seed mass - Moran's I of the residuals ( 36  sites): I = -0.046 , p = 0.839 

# Correlogram: Moran's I among pairs of sites within each distance class
correlogram_seedmass <- data.frame()
for (i in seq_len(length(correlogram_breaks) - 1)) {
  
  # Weight 1 for pairs of sites within this distance class, 0 otherwise
  class_weights <- (site_dist_seedmass > correlogram_breaks[i] &
                      site_dist_seedmass <= correlogram_breaks[i + 1]) * 1
  
  # Skip classes without any pairs of sites
  if (sum(class_weights) == 0) next
  
  moran_class <- Moran.I(sites_seedmass$resid_lmm, class_weights)
  
  correlogram_seedmass <- bind_rows(correlogram_seedmass, data.frame(
    Trait = "SeedMass",
    distance_km = (correlogram_breaks[i] + correlogram_breaks[i + 1]) / 2,
    n_pairs = sum(class_weights) / 2,
    morans_I = moran_class$observed,
    p = moran_class$p.value))
}

## 6.3. Model validation ---------------------------------------------------------

# Plot Residuals vs fitted values and Q-Q plot 
png(here("figures", "RQ3_SeedMass_validation_LMM.png"),
    width = 12, height = 5, units = "in", res = 200)
par(mfrow = c(1, 2))

plot(fitted(lmm_seedmass), seedmass_data$resid_lmm, pch = 16, cex = 0.4, col = rgb(0, 0, 0, 0.2),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Seed mass - LMM: residuals vs fitted")
abline(h = 0, col = "red", lty = 2)

qqnorm(seedmass_data$resid_lmm, pch = 16, cex = 0.4, main = "Seed mass - LMM: Q-Q plot")
qqline(seedmass_data$resid_lmm, col = "red")

par(mfrow = c(1, 1))
dev.off()


## 6.4. Figures ----------------------------------------------------------------

# Fitted lines (fixed effects only: re.form = NA)
line_seedmass <- bind_rows(
  data.frame(distance_to_boundary_km = seq(min(seedmass_data$distance_to_boundary_km), 0, length.out = 100),
             biome = factor("boreal", levels = c("boreal", "tundra"))),
  data.frame(distance_to_boundary_km = seq(0, max(seedmass_data$distance_to_boundary_km), length.out = 100),
             biome = factor("tundra", levels = c("boreal", "tundra"))))
line_seedmass$predicted <- predict(lmm_seedmass, newdata = line_seedmass, re.form = NA)

# Labels: slope on each side as % change per km
label_boreal_seedmass <- paste0("Boreal: ",
                                sprintf("%+.3f", 100 * (exp(slope_boreal_seedmass[["Estimate"]]) - 1)),
                                "% per km\np ", ifelse(slope_boreal_seedmass[["Pr(>|t|)"]] < 0.001, "< 0.001",
                                                       paste0("= ", sprintf("%.3f", slope_boreal_seedmass[["Pr(>|t|)"]]))))
label_tundra_seedmass <- paste0("Tundra: ",
                                sprintf("%+.3f", 100 * (exp(slope_tundra_seedmass[["Estimate"]]) - 1)),
                                "% per km\np ", ifelse(slope_tundra_seedmass[["Pr(>|t|)"]] < 0.001, "< 0.001",
                                                       paste0("= ", sprintf("%.3f", slope_tundra_seedmass[["Pr(>|t|)"]]))))

# Title: slope difference and jump at the treeline
p_slope_diff_seedmass <- anova_seedmass["distance_to_boundary_km:biome", "Pr(>F)"]
p_jump_seedmass <- coefs_seedmass["biometundra", "Pr(>|t|)"]
title_lmm_seedmass <- paste0("Slope difference: p ", ifelse(p_slope_diff_seedmass < 0.001, "< 0.001",
                                                            paste0("= ", sprintf("%.3f", p_slope_diff_seedmass))),
                             "; jump at treeline: p ", ifelse(p_jump_seedmass < 0.001, "< 0.001",
                                                              paste0("= ", sprintf("%.3f", p_jump_seedmass))))

# Plot figure
(plot_lmm_seedmass <- ggplot() +
  geom_point(data = seedmass_data, aes(x = distance_to_boundary_km, y = log_trait, colour = biome),
             alpha = 0.2, size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray30", linewidth = 0.7) +
  geom_line(data = line_seedmass, aes(x = distance_to_boundary_km, y = predicted, group = biome),
            colour = fit_line_colour, linewidth = 1.5) +
  annotate("text", x = -Inf, y = Inf, label = label_boreal_seedmass,
           hjust = -0.03, vjust = 1.2, size = fig_annot_text / .pt, fontface = "bold") +
  annotate("text", x = Inf, y = Inf, label = label_tundra_seedmass,
           hjust = 1.03, vjust = 1.2, size = fig_annot_text / .pt, fontface = "bold") +
  scale_colour_manual(values = biome_colours, labels = biome_labels, name = "Biome") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +   # room for the labels
  labs(x = "Distance to Biome Boundary (km)", y = "Seed Mass (mg, log scale)", title = title_lmm_seedmass) +
  theme_rq3 +
  theme(legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))))

# Save figure
ggsave(here("figures", "FigureS_RQ3_SeedMass_LMM.png"),
       plot = plot_lmm_seedmass, width = 10, height = 6, dpi = 600)

# Plot correlogram
(plot_correlogram_seedmass <- ggplot(correlogram_seedmass, aes(x = distance_km, y = morans_I)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_line(linewidth = 0.8) +
  geom_point(aes(shape = p < 0.05), size = 3) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), name = "p < 0.05") +
  labs(x = "Distance between sites (km)", y = "Moran's I of residuals", title = "Seed mass") +
  theme_rq3)

# Save figure
ggsave(here("figures", "FigureS_RQ3_SeedMass_correlogram_LMM.png"),
       plot = plot_correlogram_seedmass, width = 10, height = 6, dpi = 300)

# Combine all model results into a single df
results_lmm_seedmass <- data.frame(Trait = "SeedMass",
                                   n_records = nrow(seedmass_data),
                                   n_species = n_distinct(seedmass_data$StandardSpeciesName),
                                   n_sites_boreal = n_distinct(seedmass_data$site_name[seedmass_data$biome == "boreal"]),
                                   n_sites_tundra = n_distinct(seedmass_data$site_name[seedmass_data$biome == "tundra"]),
                                   # Slope on each side (log scale, per km)
                                   slope_boreal = slope_boreal_seedmass[["Estimate"]],
                                   se_boreal = slope_boreal_seedmass[["Std. Error"]],
                                   df_boreal = slope_boreal_seedmass[["df"]],
                                   p_boreal = slope_boreal_seedmass[["Pr(>|t|)"]],
                                   slope_tundra = slope_tundra_seedmass[["Estimate"]],
                                   se_tundra = slope_tundra_seedmass[["Std. Error"]],
                                   df_tundra = slope_tundra_seedmass[["df"]],
                                   p_tundra = slope_tundra_seedmass[["Pr(>|t|)"]],
                                   # Difference in slope (tundra minus boreal) and its F-test
                                   slope_difference = coefs_seedmass["distance_to_boundary_km:biometundra", "Estimate"],
                                   F_slope_difference = anova_seedmass["distance_to_boundary_km:biome", "F value"],
                                   p_slope_difference = anova_seedmass["distance_to_boundary_km:biome", "Pr(>F)"],
                                   # Jump at the treeline (tundra minus boreal at 0 km, log scale)
                                   jump = coefs_seedmass["biometundra", "Estimate"],
                                   se_jump = coefs_seedmass["biometundra", "Std. Error"],
                                   df_jump = coefs_seedmass["biometundra", "df"],
                                   p_jump = coefs_seedmass["biometundra", "Pr(>|t|)"],
                                   # Spatial autocorrelation of the residuals
                                   morans_I = moran_seedmass$observed,
                                   morans_I_p = moran_seedmass$p.value)

# 7. RESULTS TABLE AND COMBINED FIGURES ----------------------------------------

## 7.1. Results table ----------------------------------------------------------

# Combine everything into single df
lmm_results <- bind_rows(results_lmm_ph, results_lmm_sla,
                         results_lmm_leafn, results_lmm_seedmass) |>
  mutate(pct_boreal = 100 * (exp(slope_boreal) - 1),
         pct_boreal_lower = 100 * (exp(slope_boreal - qt(0.975, df_boreal) * se_boreal) - 1),
         pct_boreal_upper = 100 * (exp(slope_boreal + qt(0.975, df_boreal) * se_boreal) - 1),
         pct_tundra = 100 * (exp(slope_tundra) - 1),
         pct_tundra_lower = 100 * (exp(slope_tundra - qt(0.975, df_tundra) * se_tundra) - 1),
         pct_tundra_upper = 100 * (exp(slope_tundra + qt(0.975, df_tundra) * se_tundra) - 1),
         pct_jump = 100 * (exp(jump) - 1),
         pct_jump_lower = 100 * (exp(jump - qt(0.975, df_jump) * se_jump) - 1),
         pct_jump_upper = 100 * (exp(jump + qt(0.975, df_jump) * se_jump) - 1))

# Save df to use in the supplementary
write.csv(lmm_results, here("data", "derived_data", "RQ3_LMM_results_full.csv"), row.names = FALSE)

# Tidy up table
lmm_table <- lmm_results |>
  transmute(Trait,
            `Records / species` = paste0(n_records, " / ", n_species),
            `Sites boreal / tundra` = paste0(n_sites_boreal, " / ", n_sites_tundra),
            `Boreal slope: % per km (95% CI)` =
              paste0(sprintf("%.3f", pct_boreal), " (", sprintf("%.3f", pct_boreal_lower),
                     " to ", sprintf("%.3f", pct_boreal_upper), ")"),
            `Boreal slope: p` = format.pval(p_boreal, digits = 3, eps = 1e-4),
            `Tundra slope: % per km (95% CI)` =
              paste0(sprintf("%.3f", pct_tundra), " (", sprintf("%.3f", pct_tundra_lower),
                     " to ", sprintf("%.3f", pct_tundra_upper), ")"),
            `Tundra slope: p` = format.pval(p_tundra, digits = 3, eps = 1e-4),
            `Slope difference: F (p)` = paste0(sprintf("%.2f", F_slope_difference), " (",
                                               format.pval(p_slope_difference, digits = 3, eps = 1e-4), ")"),
            `Jump at treeline: % (95% CI)` =
              paste0(sprintf("%.1f", pct_jump), " (", sprintf("%.1f", pct_jump_lower),
                     " to ", sprintf("%.1f", pct_jump_upper), ")"),
            `Jump at treeline: p` = format.pval(p_jump, digits = 3, eps = 1e-4),
            `Residuals: Moran's I (p)` = paste0(sprintf("%.3f", morans_I), " (",
                                                format.pval(morans_I_p, digits = 3, eps = 1e-4), ")"))

# Quick look at the table
print(tibble::as_tibble(lmm_table), width = Inf)
# A tibble: 4 × 11
# Trait       `Records / species` `Sites boreal / tundra` `Boreal slope: % per km (95% CI)` `Boreal slope: p`
# <chr>       <chr>               <chr>                   <chr>                             <chr>            
#   1 PlantHeight 29601 / 441         85 / 215                -0.025 (-0.091 to 0.040)          0.4462           
# 2 SLA         18017 / 401         55 / 133                -0.003 (-0.033 to 0.027)          0.8361           
# 3 LeafN       7207 / 423          113 / 97                -0.019 (-0.038 to -0.001)         0.0411           
# 4 SeedMass    814 / 123           3 / 33                  0.129 (-0.742 to 1.008)           0.7659           
# `Tundra slope: % per km (95% CI)` `Tundra slope: p` `Slope difference: F (p)` `Jump at treeline: % (95% CI)`
# <chr>                             <chr>             <chr>                     <chr>                         
#   1 -0.033 (-0.075 to 0.009)          0.122             0.04 (0.8469)             -7.5 (-34.4 to 30.5)          
# 2 0.002 (-0.011 to 0.014)           0.787             0.09 (0.7685)             -25.1 (-35.8 to -12.6)        
# 3 0.006 (-0.005 to 0.017)           0.285             5.31 (0.0223)             2.3 (-8.4 to 14.2)            
# 4 0.020 (-0.059 to 0.099)           0.605             0.06 (0.8022)             35.9 (-63.3 to 403.4)         
# `Jump at treeline: p` `Residuals: Moran's I (p)`
# <chr>                 <chr>                     
#   1 0.656121              0.003 (0.8613)            
# 2 0.000312              -0.104 (0.0181)           
# 3 0.687942              0.022 (0.4749)            
# 4 0.637639              -0.046 (0.8385)

# Save pretty tabke
write.csv(lmm_table, here("data", "derived_data", "TableS_RQ3_LMM.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# Caption notes:
# - Linear mixed models (lme4/lmerTest) of log trait values:
#   log(trait) ~ distance * biome + (1 | species) + (1 | location) + (1 | site),
#   with sites (plots) nested within locations. One straight line on each side
#   of the treeline; the lines do not have to meet.
# - Slopes as % change per km towards the tundra, with 95% CIs (Satterthwaite
#   degrees of freedom). Slope difference: F-test of the distance x biome
#   interaction. Jump at treeline: difference between the tundra and boreal
#   lines at 0 km, as % (tundra relative to boreal).

## 7.2. Combined figures --------------------------------------------------------

# Combine all traits figures into a single one
all_traits_lmm <- plot_grid(plot_lmm_ph, plot_lmm_sla, plot_lmm_leafn, plot_lmm_seedmass,
                            labels = c("a)", "b)", "c)", "d)"),
                            label_size = fig_panel_label, nrow = 2)

# Save to file
ggsave(here("figures", "FigureS_RQ3_traits_LMM.png"),
       plot = all_traits_lmm, width = 20, height = 12, dpi = 600)

# Combine all correlograms together
all_traits_correlograms_lmm <- plot_grid(plot_correlogram_ph, plot_correlogram_sla,
                                         plot_correlogram_leafn, plot_correlogram_seedmass,
                                         labels = c("a)", "b)", "c)", "d)"),
                                         label_size = fig_panel_label, nrow = 2)

# Save the combined correlograms
ggsave(here("figures", "FigureS_RQ3_correlograms_LMM.png"),
       plot = all_traits_correlograms_lmm, width = 20, height = 12, dpi = 300)

# END OF SCRIPT ----------------------------------------------------------------
# - Moran's I of the mean residual per site (inverse-distance weights).