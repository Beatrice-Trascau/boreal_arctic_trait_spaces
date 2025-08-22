##----------------------------------------------------------------------------##
# PAPER 3: BOREAL AND ARCTIC PLANT SPECIES TRAIT SPACES 
# 3.3_traits_biome_NMDS
# This script contains code which run multivariate analysis on the species
# trait values
##----------------------------------------------------------------------------##

# The original NMDS analysis is in script 3.2 with plans to move it here when I get more time

# Here is the code used to label the outliers in the NMDS plot
# Function to identify points outside ellipses (no spatial packages needed!)
identify_outliers <- function(nmds_data, confidence_level = 0.95) {
  
  outliers <- data.frame()
  
  for(biome in unique(nmds_data$biome_category)) {
    # Get data for this biome
    biome_data <- nmds_data[nmds_data$biome_category == biome, ]
    
    if(nrow(biome_data) < 3) next  # Need at least 3 points for ellipse
    
    # Calculate ellipse parameters
    center <- c(mean(biome_data$MDS1), mean(biome_data$MDS2))
    cov_matrix <- cov(biome_data[, c("MDS1", "MDS2")])
    
    # Calculate Mahalanobis distances
    mahal_dist <- mahalanobis(biome_data[, c("MDS1", "MDS2")], 
                              center = center, 
                              cov = cov_matrix)
    
    # Critical value for confidence ellipse (chi-square distribution with 2 df)
    critical_value <- qchisq(confidence_level, df = 2)
    
    # Identify outliers (points outside the ellipse)
    is_outlier <- mahal_dist > critical_value
    
    if(any(is_outlier)) {
      outlier_species <- biome_data[is_outlier, ]
      outliers <- rbind(outliers, outlier_species)
    }
  }
  
  return(outliers)
}

# Identify outlier species
outlier_species <- identify_outliers(nmds_plot_data)

# Create the plot with outlier labels using ggrepel
library(ggrepel)

nmds_plot_with_outlier_labels <- nmds_plot_with_hulls +
  geom_segment(data = trait_vectors, 
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "black", 
               size = 1,
               inherit.aes = FALSE) +
  geom_text(data = trait_vectors,
            aes(x = NMDS1 * 1.1, y = NMDS2 * 1.1, label = trait),
            color = "black", 
            size = 3, 
            fontface = "bold",
            inherit.aes = FALSE) +
  # Add repelling labels for outlier species
  geom_text_repel(data = outlier_species,
                  aes(x = MDS1, y = MDS2, label = StandardSpeciesName),
                  size = 2.5,
                  color = "red",
                  fontface = "italic",
                  box.padding = 0.5,
                  point.padding = 0.3,
                  max.overlaps = Inf,
                  inherit.aes = FALSE)

# Alternative: Identify extreme points (most distant from centroid)
calculate_distance_from_center <- function(data) {
  center_x <- mean(data$MDS1)
  center_y <- mean(data$MDS2)
  data$distance_from_center <- sqrt((data$MDS1 - center_x)^2 + (data$MDS2 - center_y)^2)
  return(data)
}

# Add distance calculations
nmds_plot_data_with_dist <- calculate_distance_from_center(nmds_plot_data)

# Get top 10% most extreme points
extreme_threshold <- quantile(nmds_plot_data_with_dist$distance_from_center, 0.9)
extreme_points <- nmds_plot_data_with_dist[
  nmds_plot_data_with_dist$distance_from_center > extreme_threshold, ]

# Alternative plot with extreme points
nmds_plot_extreme_points <- nmds_plot_with_hulls +
  geom_segment(data = trait_vectors, 
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "black", 
               size = 1,
               inherit.aes = FALSE) +
  geom_text(data = trait_vectors,
            aes(x = NMDS1 * 1.1, y = NMDS2 * 1.1, label = trait),
            color = "black", 
            size = 3, 
            fontface = "bold",
            inherit.aes = FALSE) +
  geom_text_repel(data = extreme_points,
                  aes(x = MDS1, y = MDS2, label = StandardSpeciesName),
                  size = 2.5,
                  color = "darkred",
                  fontface = "italic",
                  box.padding = 0.5,
                  point.padding = 0.3,
                  max.overlaps = Inf,
                  inherit.aes = FALSE)

# Print information about outliers
cat("Statistical outliers (outside 95% confidence ellipses):\n")
if(nrow(outlier_species) > 0) {
  for(i in 1:nrow(outlier_species)) {
    cat(paste("- ", outlier_species$StandardSpeciesName[i], 
              " (", outlier_species$biome_category[i], " biome)\n"))
  }
  cat(paste("\nTotal outliers found:", nrow(outlier_species), "\n"))
} else {
  cat("No outlier species found outside the 95% confidence ellipses.\n")
}

cat("\nExtreme points (top 10% most distant from center):\n")
if(nrow(extreme_points) > 0) {
  for(i in 1:nrow(extreme_points)) {
    cat(paste("- ", extreme_points$StandardSpeciesName[i], 
              " (distance:", round(extreme_points$distance_from_center[i], 3), ")\n"))
  }
}

# Diagnostic: Check what's happening with the ellipse calculations
diagnose_outliers <- function(nmds_data, confidence_level = 0.95) {
  
  cat("=== OUTLIER DIAGNOSTICS ===\n")
  
  for(biome in unique(nmds_data$biome_category)) {
    cat(paste("\n", toupper(biome), "BIOME:\n"))
    
    # Get data for this biome
    biome_data <- nmds_data[nmds_data$biome_category == biome, ]
    cat(paste("  Number of points:", nrow(biome_data), "\n"))
    
    if(nrow(biome_data) < 3) {
      cat("  Too few points for ellipse calculation\n")
      next
    }
    
    # Calculate ellipse parameters
    center <- c(mean(biome_data$MDS1), mean(biome_data$MDS2))
    cov_matrix <- cov(biome_data[, c("MDS1", "MDS2")])
    
    # Calculate Mahalanobis distances
    mahal_dist <- mahalanobis(biome_data[, c("MDS1", "MDS2")], 
                              center = center, 
                              cov = cov_matrix)
    
    # Critical value for confidence ellipse
    critical_value <- qchisq(confidence_level, df = 2)
    cat(paste("  Critical value (95%):", round(critical_value, 3), "\n"))
    
    # Show all distances
    biome_data$mahal_distance <- mahal_dist
    biome_data$is_outlier <- mahal_dist > critical_value
    
    cat("  Mahalanobis distances:\n")
    for(i in 1:nrow(biome_data)) {
      status <- if(biome_data$is_outlier[i]) "OUTLIER" else "normal"
      cat(paste("    ", biome_data$StandardSpeciesName[i], ":", 
                round(biome_data$mahal_distance[i], 3), "(", status, ")\n"))
    }
  }
}

# Run diagnostics
diagnose_outliers(nmds_plot_data)

# Alternative approach: Use multiple confidence levels
identify_outliers_multiple_levels <- function(nmds_data) {
  
  all_potential_outliers <- data.frame()
  confidence_levels <- c(0.95, 0.90, 0.85, 0.80)  # Try different thresholds
  
  for(conf_level in confidence_levels) {
    outliers <- identify_outliers(nmds_data, conf_level)
    if(nrow(outliers) > 0) {
      outliers$confidence_level <- conf_level
      all_potential_outliers <- rbind(all_potential_outliers, outliers)
    }
  }
  
  # Remove duplicates but keep the most stringent confidence level
  all_potential_outliers <- all_potential_outliers[
    !duplicated(all_potential_outliers$StandardSpeciesName), ]
  
  return(all_potential_outliers)
}

# Try multiple confidence levels
outliers_multi_level <- identify_outliers_multiple_levels(nmds_plot_data)

# Visual identification approach: Find points that appear outside the visual hulls
identify_visual_outliers <- function(nmds_data) {
  
  # For each biome, find the convex hull and identify points outside
  visual_outliers <- data.frame()
  
  for(biome in unique(nmds_data$biome_category)) {
    biome_data <- nmds_data[nmds_data$biome_category == biome, ]
    
    if(nrow(biome_data) < 3) next
    
    # Find convex hull
    hull_indices <- chull(biome_data$MDS1, biome_data$MDS2)
    hull_points <- biome_data[hull_indices, ]
    
    # Expand hull slightly to account for ellipse vs convex hull differences
    center_x <- mean(biome_data$MDS1)
    center_y <- mean(biome_data$MDS2)
    
    # Points that are far from both the center and the hull
    for(i in 1:nrow(nmds_data)) {
      point <- nmds_data[i, ]
      
      # Distance from this biome's center
      dist_from_center <- sqrt((point$MDS1 - center_x)^2 + (point$MDS2 - center_y)^2)
      
      # If point is far from center and from the other biome, it might be a visual outlier
      if(dist_from_center > quantile(sqrt((biome_data$MDS1 - center_x)^2 + 
                                          (biome_data$MDS2 - center_y)^2), 0.75)) {
        visual_outliers <- rbind(visual_outliers, point)
      }
    }
  }
  
  return(visual_outliers[!duplicated(visual_outliers$StandardSpeciesName), ])
}

# Get visual outliers
visual_outliers <- identify_visual_outliers(nmds_plot_data)

# Combine all approaches
all_outliers <- rbind(
  outlier_species[, names(nmds_plot_data)],
  outliers_multi_level[, names(nmds_plot_data)],
  visual_outliers
)
all_outliers <- all_outliers[!duplicated(all_outliers$StandardSpeciesName), ]

# Create comprehensive plot with all potential outliers (NO VECTORS)
nmds_plot_comprehensive_no_vectors <- nmds_plot_with_hulls +
  # Label all potential outliers
  geom_text_repel(data = all_outliers,
                  aes(x = MDS1, y = MDS2, label = StandardSpeciesName),
                  size = 2.5,
                  color = "red",
                  fontface = "italic",
                  box.padding = 0.5,
                  point.padding = 0.3,
                  max.overlaps = Inf,
                  force = 2,  # Increase repulsion force
                  inherit.aes = FALSE)

# Simple manual approach: Just label the most extreme points in each direction
extreme_points_directional <- nmds_plot_data[
  nmds_plot_data$MDS1 == max(nmds_plot_data$MDS1) |  # Rightmost
    nmds_plot_data$MDS1 == min(nmds_plot_data$MDS1) |  # Leftmost
    nmds_plot_data$MDS2 == max(nmds_plot_data$MDS2) |  # Topmost
    nmds_plot_data$MDS2 == min(nmds_plot_data$MDS2) |  # Bottommost
    abs(nmds_plot_data$MDS1) > quantile(abs(nmds_plot_data$MDS1), 0.9) |
    abs(nmds_plot_data$MDS2) > quantile(abs(nmds_plot_data$MDS2), 0.9), ]

# Most permissive approach: Label anything that looks like it might be outside (NO VECTORS)
nmds_plot_permissive_no_vectors <- nmds_plot_with_hulls +
  geom_text_repel(data = extreme_points_directional,
                  aes(x = MDS1, y = MDS2, label = StandardSpeciesName),
                  size = 3.5,
                  color = "black",
                  fontface = "italic",
                  box.padding = 0.5,
                  point.padding = 0.3,
                  max.overlaps = Inf,
                  force = 2,
                  inherit.aes = FALSE)

# Statistical outliers only (NO VECTORS)
nmds_plot_statistical_only <- nmds_plot_with_hulls +
  geom_text_repel(data = outlier_species,
                  aes(x = MDS1, y = MDS2, label = StandardSpeciesName),
                  size = 3.5,
                  color = "black",
                  fontface = "italic",
                  box.padding = 0.5,
                  point.padding = 0.3,
                  max.overlaps = Inf,
                  inherit.aes = FALSE)

cat("\n=== SUMMARY ===\n")
cat("Statistical outliers (95%):", nrow(outlier_species), "\n")
cat("Multi-level outliers:", nrow(outliers_multi_level), "\n") 
cat("Visual outliers:", nrow(visual_outliers), "\n")
cat("All combined outliers:", nrow(all_outliers), "\n")
cat("Extreme directional points:", nrow(extreme_points_directional), "\n")

# Display the plots without vectors
cat("\n=== PLOTS WITHOUT VECTORS ===\n")
print("1. Statistical outliers only:")
nmds_plot_statistical_only

print("2. Comprehensive (all methods combined):")
nmds_plot_comprehensive_no_vectors

print("3. Permissive (extreme points):")
nmds_plot_permissive_no_vectors