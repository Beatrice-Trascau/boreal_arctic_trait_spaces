# =============================================================================
# 2_3_recalculate_biome_boundary_distance.R
#Script to recalculate the biome assignment and distance to biome boundary for 
# every grid cell used in script 2.1, WITHOUT re-downloading GBIF occurrences.
#
# This is a fix for an issue in script 2.1 (lines 323-332) where cells 
# overlapping both biomes were always treated as boreal
# I fixed this in this script by assigning the biome which overlapped the
# majority of the cell's area
#
# This script is reusing the saved cell-level output of 2.1 (cell_id, species, 
# count, in_boreal, in_tundra). No redownload of occurrences was needed.
# =============================================================================

# 0. SETUP ---------------------------------------------------------------------

# Load libraries
library(here)
library(dplyr)
library(sf)
library(terra)
library(ggplot2)

# Load data
input_cells <- here("data", "derived_data", "dist_to_biome_boundary_June25.rds")
input_summaries <- here("data", "derived_data", "species_summaries_dist_to_biome_boundary_June25.rds")
out_dir <- here("data", "derived_data")
fig_dir <- here("figures")
suffix <- "corrected_Sept2026"

# Create a 1m tolerance when comparing the new distances to the old ones
tol_km <- 0.001

# Set a maximum distance from a 50 km cell's centroid to any point in it (~35.36 km)
half_diag_km <- 25 * sqrt(2)  

# Create directory if needed
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# 1. LOAD SAVED DATA -----------------------------------------------------------

## 1.1. Cell-level data from script 2.1 ----------------------------------------
cell_data <- as.data.frame(readRDS(input_cells))
summaries_original <- as.data.frame(readRDS(input_summaries))

required_cols <- c("cell_id", "species", "count", "in_boreal", "in_tundra",
                   "distance_to_boundary_km", "biome")
missing_cols <- setdiff(required_cols, names(cell_data))
if (length(missing_cols) > 0) {
  stop("Saved cell data is missing columns: ", paste(missing_cols, collapse = ", "))
}

cat("Loaded", nrow(cell_data), "cell-species records for",
    length(unique(cell_data$species)), "species\n")

## 1.2. One row per cell with the ORIGINAL values ------------------------------
original_cell_lookup <- cell_data |>
  distinct(cell_id, in_boreal, in_tundra, distance_to_boundary_km, biome)

if (any(duplicated(original_cell_lookup$cell_id))) {
  stop("Some cells have more than one distance/biome value in the saved data. ",
       "Inspect before continuing.")
}

cat("Unique cells in saved data:", nrow(original_cell_lookup), "\n")
cat("  of which overlap both biomes:",
    sum(original_cell_lookup$in_boreal & original_cell_lookup$in_tundra), "\n\n")


# 2. REBUILD BIOMES AND GRID (identical to script 2.1, lines 38-97) -----------

cat("Rebuilding biomes and grid...\n")

global_biomes <- st_read(here("data", "raw_data", "biomes", "wwf_terr_ecos.shp"), quiet = TRUE)

boreal_forest <- st_union(global_biomes[global_biomes$BIOME == 6, ])
tundra <- st_union(global_biomes[global_biomes$BIOME == 11 &
                                   (global_biomes$REALM == "PA" | global_biomes$REALM == "NA"), ])

boreal_forest <- st_make_valid(boreal_forest)
tundra <- st_make_valid(tundra)

boreal_sf <- st_transform(boreal_forest, "EPSG:3574")
tundra_sf <- st_transform(tundra, "EPSG:3574")

combined_biomes <- st_union(boreal_sf, tundra_sf)
combined_extent <- st_bbox(combined_biomes)

grid <- rast(extent = combined_extent,
             resolution = c(50000, 50000),
             crs = "EPSG:3574")

polygrid <- as.polygons(grid)
polygrid$cell_id <- 1:nrow(polygrid)
polygrid_sf <- st_as_sf(polygrid)

cells_in_boreal <- st_intersects(polygrid_sf, boreal_sf, sparse = FALSE)[, 1]
cells_in_tundra <- st_intersects(polygrid_sf, tundra_sf, sparse = FALSE)[, 1]
cells_in_biomes <- cells_in_boreal | cells_in_tundra

polygrid_filtered <- polygrid_sf[cells_in_biomes, ]
polygrid_filtered$in_boreal <- cells_in_boreal[cells_in_biomes]
polygrid_filtered$in_tundra <- cells_in_tundra[cells_in_biomes]

cat("Grid rebuilt:", nrow(polygrid_filtered), "cells within biomes",
    "(compare with 'Grid created:' in the original 2.1 output)\n\n")


# 3. CHECK 1: REBUILT GRID MATCHES THE ORIGINAL --------------------------------

cat("--- CHECK 1: grid identity ---\n")

grid_flags <- st_drop_geometry(polygrid_filtered)[, c("cell_id", "in_boreal", "in_tundra")]

check1 <- original_cell_lookup |>
  select(cell_id, in_boreal, in_tundra) |>
  left_join(grid_flags, by = "cell_id", suffix = c("_saved", "_new"))

n_missing  <- sum(is.na(check1$in_boreal_new))
n_mismatch <- sum(check1$in_boreal_saved != check1$in_boreal_new |
                    check1$in_tundra_saved != check1$in_tundra_new, na.rm = TRUE)

cat("  Saved cells not found in rebuilt grid:", n_missing, "\n")
cat("  Cells with mismatching biome flags:  ", n_mismatch, "\n")

if (n_missing > 0 || n_mismatch > 0) {
  write.csv(check1, file.path(out_dir, paste0("CHECK1_failed_", suffix, ".csv")), row.names = FALSE)
  stop("CHECK 1 FAILED: the rebuilt grid does not match the original. cell_ids cannot be trusted.")
}
cat("  PASSED\n\n")

# Keep only the cells that actually carry occurrence data
cells_used <- polygrid_filtered[polygrid_filtered$cell_id %in% original_cell_lookup$cell_id, ]
st_agr(cells_used) <- "constant"


# 4. MAJORITY OVERLAP FOR CELLS IN BOTH BIOMES ---------------------------------

overlap_polys <- cells_used[cells_used$in_boreal & cells_used$in_tundra, "cell_id"]
st_agr(overlap_polys) <- "constant"

cat("Calculating boreal/tundra area for", nrow(overlap_polys), "overlap cells",
    "(may take a few minutes)...\n")

intersect_area <- function(polys, biome_polygon) {
  x <- st_intersection(polys, biome_polygon)
  data.frame(cell_id = x$cell_id, area_m2 = as.numeric(st_area(x))) |>
    group_by(cell_id) |>
    summarise(area_m2 = sum(area_m2), .groups = "drop")
}

boreal_areas <- intersect_area(overlap_polys, boreal_sf)
tundra_areas <- intersect_area(overlap_polys, tundra_sf)

overlap_lookup <- st_drop_geometry(overlap_polys) |>
  left_join(rename(boreal_areas, boreal_area_m2 = area_m2), by = "cell_id") |>
  left_join(rename(tundra_areas, tundra_area_m2 = area_m2), by = "cell_id") |>
  mutate(boreal_area_m2 = coalesce(boreal_area_m2, 0),
         tundra_area_m2 = coalesce(tundra_area_m2, 0),
         cell_boreal_area_pct = 100 * boreal_area_m2 / (boreal_area_m2 + tundra_area_m2),
         cell_tundra_area_pct = 100 - cell_boreal_area_pct,
         # ties (very unlikely) go to boreal, as in the original script
         majority_biome = ifelse(boreal_area_m2 >= tundra_area_m2, "boreal", "tundra"))

cat("  Majority boreal:", sum(overlap_lookup$majority_biome == "boreal"), "\n")
cat("  Majority tundra:", sum(overlap_lookup$majority_biome == "tundra"), "\n")
cat("  Exact ties:     ", sum(overlap_lookup$boreal_area_m2 == overlap_lookup$tundra_area_m2), "\n\n")


# 5. RECALCULATE DISTANCE FOR EVERY CELL ---------------------------------------

cat("Recalculating distance to boundary for all", nrow(cells_used), "cells...\n")

cell_lookup_new <- st_drop_geometry(cells_used)[, c("cell_id", "in_boreal", "in_tundra")] |>
  left_join(select(overlap_lookup, cell_id, cell_boreal_area_pct,
                   cell_tundra_area_pct, majority_biome), by = "cell_id") |>
  mutate(biome_new = case_when(in_boreal & in_tundra ~ majority_biome,
                               in_boreal             ~ "boreal",
                               in_tundra             ~ "tundra"))

boreal_boundary <- st_boundary(boreal_sf)
tundra_boundary <- st_boundary(tundra_sf)

centroids <- st_centroid(cells_used)
centroids <- centroids[match(cell_lookup_new$cell_id, centroids$cell_id), ]

is_boreal <- cell_lookup_new$biome_new == "boreal"
dist_m <- numeric(nrow(cell_lookup_new))
dist_m[is_boreal]  <-  as.numeric(st_distance(centroids[is_boreal, ],  boreal_boundary))
dist_m[!is_boreal] <- -as.numeric(st_distance(centroids[!is_boreal, ], tundra_boundary))

cell_lookup_new$distance_to_boundary_m  <- dist_m
cell_lookup_new$distance_to_boundary_km <- dist_m / 1000

cat("  Done\n\n")


# 6. CHECK 2: ONLY THE INTENDED CELLS CHANGE ----------------------------------

cat("--- CHECK 2: reproduction of original values ---\n")

check2 <- original_cell_lookup |>
  select(cell_id, biome_original = biome, distance_original_km = distance_to_boundary_km) |>
  left_join(cell_lookup_new, by = "cell_id") |>
  mutate(overlap         = in_boreal & in_tundra,
         expected_change = overlap & biome_new == "tundra",
         changed         = abs(distance_to_boundary_km - distance_original_km) > tol_km |
           biome_new != biome_original)

check2_summary <- check2 |>
  mutate(cell_type = case_when(!overlap                 ~ "single biome",
                               biome_new == "boreal"    ~ "overlap, majority boreal",
                               TRUE                     ~ "overlap, majority tundra")) |>
  count(cell_type, changed)
print(check2_summary)

unexpected <- filter(check2, changed & !expected_change)
cat("\n  Cells that changed but should not have:", nrow(unexpected), "\n")
cat("  Tundra-majority overlap cells that changed:",
    sum(check2$changed & check2$expected_change), "of", sum(check2$expected_change), "\n")

if (nrow(unexpected) > 0) {
  write.csv(unexpected, file.path(out_dir, paste0("CHECK2_failed_", suffix, ".csv")), row.names = FALSE)
  stop("CHECK 2 FAILED: cells outside the tundra-majority overlap group changed. ",
       "The recalculation does not reproduce the original method.")
}
cat("  PASSED\n\n")


# 7. SPECIES SUMMARY FUNCTION (identical to script 2.1, section 2.3) -----------

calculate_species_summaries <- function(results_df) {
  results_df |>
    group_by(species) |>
    summarise(total_cells = n(),
              total_occurrences = sum(count),
              cells_in_boreal = sum(biome == "boreal"),
              cells_in_tundra = sum(biome == "tundra"),
              pct_cells_boreal = round(100 * sum(biome == "boreal") / n(), 1),
              pct_cells_tundra = round(100 * sum(biome == "tundra") / n(), 1),
              mean_distance_km = round(mean(distance_to_boundary_km, na.rm = TRUE), 2),
              median_distance_km = round(median(distance_to_boundary_km, na.rm = TRUE), 2),
              sd_distance_km = round(sd(distance_to_boundary_km, na.rm = TRUE), 2),
              min_distance_km = round(min(distance_to_boundary_km, na.rm = TRUE), 2),
              max_distance_km = round(max(distance_to_boundary_km, na.rm = TRUE), 2),
              mean_distance_boreal_km = round(mean(distance_to_boundary_km[biome == "boreal"],
                                                   na.rm = TRUE), 2),
              median_distance_boreal_km = round(median(distance_to_boundary_km[biome == "boreal"],
                                                       na.rm = TRUE), 2),
              mean_distance_tundra_km = round(mean(distance_to_boundary_km[biome == "tundra"],
                                                   na.rm = TRUE), 2),
              median_distance_tundra_km = round(median(distance_to_boundary_km[biome == "tundra"],
                                                       na.rm = TRUE), 2),
              weighted_mean_distance_km = round(weighted.mean(distance_to_boundary_km, count,
                                                              na.rm = TRUE), 2),
              latitudinal_range_km = round(max_distance_km - min_distance_km, 2),
              .groups = "drop") |>
    mutate(across(where(is.numeric), ~ ifelse(is.nan(.x), NA, .x))) |>
    as.data.frame()
}

## CHECK 2b: saved summaries can be rebuilt from the saved cell data -----------
cat("--- CHECK 2b: saved summaries match saved cell data ---\n")

summaries_rebuilt <- calculate_species_summaries(cell_data)
shared_cols <- intersect(names(summaries_original), names(summaries_rebuilt))

a <- summaries_original |> arrange(species) |> select(all_of(shared_cols))
b <- summaries_rebuilt  |> arrange(species) |> select(all_of(shared_cols))

if (isTRUE(all.equal(a, b, check.attributes = FALSE, tolerance = 1e-6))) {
  cat("  PASSED\n\n")
} else {
  warning("CHECK 2b: the saved species summaries could not be exactly rebuilt from the saved ",
          "cell data. The before/after comparison below uses the saved summaries as baseline; ",
          "treat it with care.")
  cat("  WARNING (see warnings())\n\n")
}


# 8. CHECK 3: VALUE RANGES -----------------------------------------------------

cat("--- CHECK 3: value ranges ---\n")

range_table <- cell_lookup_new |>
  mutate(cell_type = ifelse(in_boreal & in_tundra, "overlap", "single biome")) |>
  group_by(cell_type, biome_new) |>
  summarise(n = n(),
            min_km = round(min(distance_to_boundary_km), 2),
            max_km = round(max(distance_to_boundary_km), 2),
            .groups = "drop")
print(range_table)

sign_problems <- filter(cell_lookup_new,
                        (biome_new == "boreal" & distance_to_boundary_km < 0) |
                          (biome_new == "tundra" & distance_to_boundary_km > 0))

bound_problems <- filter(cell_lookup_new,
                         in_boreal & in_tundra,
                         abs(distance_to_boundary_km) > half_diag_km + tol_km)

cat("\n  Cells with wrong sign:", nrow(sign_problems), "\n")
cat("  Overlap cells further than", round(half_diag_km, 2), "km from boundary:",
    nrow(bound_problems), "\n")

if (nrow(sign_problems) > 0) {
  write.csv(sign_problems, file.path(out_dir, paste0("CHECK3_sign_failed_", suffix, ".csv")), row.names = FALSE)
  stop("CHECK 3 FAILED: some cells have a distance sign that does not match their biome.")
}
if (nrow(bound_problems) > 0) {
  write.csv(bound_problems, file.path(out_dir, paste0("CHECK3_range_warning_", suffix, ".csv")), row.names = FALSE)
  warning("CHECK 3: some overlap cells lie further from the boundary than geometrically expected. ",
          "This can happen if the biome polygons overlap each other. Inspect the saved CSV.")
} else {
  cat("  PASSED\n\n")
}


# 9. CHECK 4: MAP OF NEW DISTANCES ---------------------------------------------

cat("--- CHECK 4: saving map for visual inspection ---\n")

map_points <- centroids
map_points$distance_km <- cell_lookup_new$distance_to_boundary_km
map_points$changed <- map_points$cell_id %in% check2$cell_id[check2$changed]

p_map <- ggplot() +
  geom_sf(data = boreal_sf, fill = "#556B2F", alpha = 0.25, colour = NA) +
  geom_sf(data = tundra_sf, fill = "#B8860B", alpha = 0.25, colour = NA) +
  geom_sf(data = map_points, aes(colour = distance_km), size = 0.5) +
  geom_sf(data = map_points[map_points$changed, ], shape = 1, colour = "red", size = 1.2) +
  scale_colour_gradient2(low = "#8B4513", mid = "grey90", high = "#1B4D1B", midpoint = 0,
                         name = "Distance (km)\n+ boreal / - tundra") +
  labs(title = "Recalculated distance to biome boundary per grid cell",
       subtitle = "Red circles = cells whose value changed (tundra-majority overlap cells)") +
  theme_minimal()

ggsave(file.path(fig_dir, paste0("cell_distances_", suffix, ".png")), p_map,
       width = 10, height = 10, dpi = 300)
cat("  Saved figures/cell_distances_", suffix, ".png\n\n", sep = "")


# 10. BUILD CORRECTED DATA AND SUMMARIES ---------------------------------------

cell_data_corrected <- cell_data |>
  rename(biome_original = biome,
         distance_original_km = distance_to_boundary_km) |>
  select(-any_of("distance_to_boundary_m")) |>
  left_join(cell_lookup_new |>
              select(cell_id,
                     biome = biome_new,
                     distance_to_boundary_m,
                     distance_to_boundary_km,
                     cell_boreal_area_pct,
                     cell_tundra_area_pct),
            by = "cell_id")

summaries_corrected <- calculate_species_summaries(cell_data_corrected)
cat("Corrected summaries calculated for", nrow(summaries_corrected), "species\n\n")


# 11. COMPARE ORIGINAL VS CORRECTED SPECIES VALUES -----------------------------

cat("=== ORIGINAL VS CORRECTED ===\n\n")

compare_cols <- c("species", "total_cells", "pct_cells_boreal", "mean_distance_km",
                  "median_distance_km", "weighted_mean_distance_km")

comparison <- summaries_original |>
  select(all_of(compare_cols)) |>
  full_join(select(summaries_corrected, all_of(compare_cols)),
            by = "species", suffix = c("_original", "_corrected")) |>
  mutate(change_median_km        = median_distance_km_corrected - median_distance_km_original,
         change_mean_km          = mean_distance_km_corrected - mean_distance_km_original,
         change_weighted_mean_km = weighted_mean_distance_km_corrected - weighted_mean_distance_km_original,
         change_pct_boreal       = pct_cells_boreal_corrected - pct_cells_boreal_original,
         median_changed          = abs(change_median_km) > 0,
         any_change              = median_changed | abs(change_mean_km) > 0 |
           abs(change_weighted_mean_km) > 0 | change_pct_boreal != 0,
         median_sign_flip        = sign(median_distance_km_original) != sign(median_distance_km_corrected),
         mean_sign_flip          = sign(mean_distance_km_original) != sign(mean_distance_km_corrected)) |>
  arrange(desc(abs(change_median_km)))

if (any(comparison$total_cells_original != comparison$total_cells_corrected, na.rm = TRUE)) {
  warning("Number of cells per species differs between original and corrected summaries.")
}

cat("Species compared:                    ", nrow(comparison), "\n")
cat("Species with any change:             ", sum(comparison$any_change, na.rm = TRUE), "\n")
cat("Species whose MEDIAN changed:        ", sum(comparison$median_changed, na.rm = TRUE), "\n")
cat("Species whose median changed sign:   ", sum(comparison$median_sign_flip, na.rm = TRUE), "\n")
cat("Species whose mean changed sign:     ", sum(comparison$mean_sign_flip, na.rm = TRUE), "\n\n")

cat("Absolute change in median distance (species with a change):\n")
print(summary(abs(comparison$change_median_km[comparison$median_changed])))
cat("\nAbsolute change in mean distance (species with a change):\n")
print(summary(abs(comparison$change_mean_km[comparison$any_change])))

cat("\nTop 20 species by change in median distance:\n")
print(comparison |>
        filter(median_changed) |>
        select(species, median_distance_km_original, median_distance_km_corrected,
               change_median_km, mean_distance_km_original, mean_distance_km_corrected,
               pct_cells_boreal_original, pct_cells_boreal_corrected) |>
        head(20))

sign_flippers <- filter(comparison, median_sign_flip | mean_sign_flip)
if (nrow(sign_flippers) > 0) {
  cat("\nSpecies whose median or mean moved to the other side of the boundary:\n")
  print(select(sign_flippers, species, median_distance_km_original, median_distance_km_corrected,
               mean_distance_km_original, mean_distance_km_corrected))
}

p_compare <- ggplot(comparison, aes(median_distance_km_original, median_distance_km_corrected)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
  geom_point(aes(colour = median_changed), alpha = 0.6, size = 1.8) +
  scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "darkred"),
                      labels = c(`FALSE` = "Unchanged", `TRUE` = "Changed"), name = NULL) +
  labs(title = "Species median distance to biome boundary",
       x = "Original median (km)", y = "Corrected median (km)") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, paste0("species_median_original_vs_", suffix, ".png")), p_compare,
       width = 8, height = 8, dpi = 300)


# 12. SAVE OUTPUTS -------------------------------------------------------------

cat("\nSaving outputs...\n")

saveRDS(cell_data_corrected, file.path(out_dir, paste0("dist_to_biome_boundary_", suffix, ".rds")))
write.csv(cell_data_corrected, file.path(out_dir, paste0("dist_to_biome_boundary_", suffix, ".csv")),
          row.names = FALSE)

saveRDS(summaries_corrected,
        file.path(out_dir, paste0("species_summaries_dist_to_biome_boundary_", suffix, ".rds")))
write.csv(summaries_corrected,
          file.path(out_dir, paste0("species_summaries_dist_to_biome_boundary_", suffix, ".csv")),
          row.names = FALSE)

# Per-cell audit trail: original and new biome/distance for every cell
write.csv(check2, file.path(out_dir, paste0("cell_level_audit_", suffix, ".csv")), row.names = FALSE)

write.csv(comparison, file.path(out_dir, paste0("species_comparison_", suffix, ".csv")), row.names = FALSE)
write.csv(filter(comparison, any_change),
          file.path(out_dir, paste0("species_changed_", suffix, ".csv")), row.names = FALSE)

# Save the grid so it never has to be rebuilt again
saveRDS(polygrid_filtered, file.path(out_dir, "analysis_grid_50km.rds"))

cat("\n=== DONE ===\n")
cat("Use species_summaries_dist_to_biome_boundary_", suffix, ".rds in downstream scripts ",
    "(e.g. update the file name in script 1.4).\n", sep = "")

# END OF SCRIPT  ---------------------------------------------------------------