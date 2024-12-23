library(tidyr)
library(dplyr)
library(ggplot2)
library(sf)
library(spdep)
library(rnaturalearth)
library(rnaturalearthdata)
library(geosphere)
library(igraph)
library(gridExtra)
library(Matrix)
library(patchwork)
library(rlang)
library(xtable)
library(tools)
library(ggnewscale)
library(cowplot)

setwd('/Users/sophie/Documents/implied_weights/superfunds/')

source('../impliedweights_randeffs.R')

states <- st_read('/Users/sophie/Downloads/tl_2010_us_state10/tl_2010_us_state10.shp')

# Read in preprocessed data
load('data/preprocessed_superfunds.RData')
n <- nrow(buffers)
us_outline <- ne_countries(scale = "medium", country = "United States of America", returnclass = "sf")

# Plot Z on the map
# order buffers by Z
buffer_centroids <- st_centroid(buffers)
buffer_centroids_ordered <- buffer_centroids[order(buffers$Z),]
png('images/treatment.png', width = 1500, height = 800, res = 160)
ggplot() +
  # Add map outline
  geom_sf(data = us_outline, fill = NA, color = "black", linetype = "solid") +
  # Plot centroids with both color and shape mapped to factor(Z)
  geom_sf(data = buffer_centroids_ordered, aes(color = factor(Z), shape = factor(Z)), size = 3) +
  coord_sf(xlim = c(-125, -65), ylim = c(25, 50), expand = FALSE) +
  labs(
    title = "Superfund Sites that were cleaned up and removed from the National Priority List between 2005 and 2014",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  ) +
  # Use the same variable for color and shape scales with unified legend
  scale_color_manual(
    values = c("0" = "red", "1" = "green"),
    name = "Treatment (Cleanup)",
    labels = c("Control", "Treated")
  ) +
  scale_shape_manual(
    values = c("0" = 16, "1" = 17), # Circle for control, triangle for treated
    name = "Treatment (Cleanup)",
    labels = c("Control", "Treated")
  )
dev.off()

# Next, visualize the confounders. 
confounder_names <- c(
  'population_density',
  'percent_hispanic',
  'percent_black',
  'percent_indigenous',
  'percent_asian',
  'median_household_income',
  'median_house_value',
  'percent_poverty',
  'percent_high_school_grad',
  'median_year_built'
)

gs <- list()
for (i in seq_along(confounder_names)) {
  print(confounder_names[i])
  # Order buffer_centroids by increasing confounder value
  buffer_centroids <- buffer_centroids[order(buffers[[confounder_names[i]]]),]
  gs[[i]] <- ggplot() +
    geom_sf(data = us_outline, fill = NA, color = "black", linetype = "solid") +
    geom_sf(data = buffer_centroids, 
            aes(color = .data[[confounder_names[i]]], shape = factor(Z)), size = 2) + 
    coord_sf(xlim = c(-125, -65), ylim = c(25, 50), expand = FALSE) +
    labs(
      title = confounder_names[i],
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      legend.position = "bottom",
      legend.key.width = unit(1.2, "cm")
    ) +
    scale_color_viridis_c()+ 
    scale_shape_manual(values = c("0" = 16, "1" = 17), guide = "none") # Remove shape legend
}

# Save the combined plot
png('images/confounders.png', width = 2000, height = 1500, res = 120)
grid.arrange(grobs = gs, ncol = 4)
dev.off()

# Order by cluster
load('data/preprocessed_superfunds.RData')
all_baseweights <- compute_allbaseweights(Z = Z, adjacency_matrix = adjacency_matrix,
                                          statefactor = clusters,
                                          sig2gam = 0.01,
                                          phi = 0.4,
                                          distmat = dmat/10000)
# Check all base weights sum to 1 in each treatment group
sum(all_baseweights$basePooled[Z == 1])
sum(all_baseweights$basePooled[Z == 0])
sum(all_baseweights$baseRE[Z == 1])
sum(all_baseweights$baseRE[Z == 0])
sum(all_baseweights$baseCAR[Z == 1])
sum(all_baseweights$baseCAR[Z == 0])
sum(all_baseweights$baseFE[Z == 1])
sum(all_baseweights$baseFE[Z == 0])
sum(all_baseweights$baseGP[Z == 1])
sum(all_baseweights$baseGP[Z == 0])
sum(all_baseweights$baseSAR[Z == 1])
sum(all_baseweights$baseSAR[Z == 0])

allmint <- min(0,min(c(all_baseweights$baseCAR[Z == 1], all_baseweights$baseFE[Z == 1], all_baseweights$baseGP[Z == 1])))
allmaxt <- max(c(all_baseweights$baseCAR[Z == 1], all_baseweights$baseFE[Z == 1], all_baseweights$baseGP[Z == 1]))
allminc <- min(0,min(c(all_baseweights$baseCAR[Z == 0], all_baseweights$baseFE[Z == 0], all_baseweights$baseGP[Z == 0])))
allmaxc <- max(c(all_baseweights$baseCAR[Z == 0], all_baseweights$baseFE[Z == 0], all_baseweights$baseGP[Z == 0]))

# Merge the base weights with the data
buffer_centroids <- st_centroid(buffers)
buffers_merged <- cbind(buffer_centroids, all_baseweights)

# Plot baseCAR
buffers_merged <- buffers_merged[order(buffers_merged$baseCAR),]
gCAR <- ggplot() +
  geom_sf(data = states, fill = "black", color = "white", linetype = "solid") +
  geom_sf(data = buffers_merged[buffers_merged$Z == 0,], aes(color = baseCAR, shape = 'Control'), size = 1.5) +
  labs(
    shape = "Treatment",
    color = "Base weights (control)"
  ) +
  scale_color_gradient2(
    low = "green", mid = "white", high = "red", midpoint = 0, 
    limits = c(allminc, allmaxc)) + #limits = c(0, max(buffers_merged$baseCAR[buffers_merged$Z == 0]))) + 
  new_scale_color() +
  geom_sf(data = buffers_merged[buffers_merged$Z == 1,], aes(color = baseCAR, shape = 'Treated'), size = 1.5) +
  scale_color_gradient2(
    low = "red", mid = "white", high = "green", midpoint = 0, 
    limits = c(allmint, allmaxt)) + #c(0, max(buffers_merged$baseCAR[buffers_merged$Z == 1]))) + 
  coord_sf(xlim = c(-125, -65), ylim = c(25, 50), expand = FALSE) +
  labs(
    title = "Conditional Autoregressive",
    x = "Longitude",
    y = "Latitude",
    shape = "Treatment",
    color = "Base weights (treated)"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

# Plot base FE
buffers_merged <- buffers_merged[order(buffers_merged$baseFE),]
gRE <- ggplot() +
  geom_sf(data = states, fill = "black", color = "white", linetype = "solid") +
  geom_sf(data = buffers_merged[buffers_merged$Z == 0,], aes(color = baseFE, shape = 'Control'), size = 1.5) +
  labs(
    shape = "Treatment",
    color = "Base weights (control)"
  ) +
  scale_color_gradient2(
    low = "green", mid = "white", high = "red", midpoint = 0, 
    limits = c(allminc, allmaxc)) + #limits = c(0, max(buffers_merged$baseFE[buffers_merged$Z == 0]))) + 
  new_scale_color() +
  geom_sf(data = buffers_merged[buffers_merged$Z == 1,], aes(color = baseFE, shape = 'Treated'), size = 1.5) +
  scale_color_gradient2(
    low = "red", mid = "white", high = "green", midpoint = 0, 
    limits = c(allmint, allmaxt)) + #limits = c(0, max(buffers_merged$baseFE[buffers_merged$Z == 1]))) + 
  coord_sf(xlim = c(-125, -65), ylim = c(25, 50), expand = FALSE) +
  labs(
    title = "Random Effects",
    x = "Longitude",
    y = "Latitude",
    shape = "Treatment",
    color = "Base weights (treated)"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = 'bottom'
  )

# Plot base GP
buffers_merged <- buffers_merged[order(buffers_merged$baseGP),]
gGP <- ggplot() +
  geom_sf(data = states, fill = "black", color = "white", linetype = "solid") +
  geom_sf(data = buffers_merged[buffers_merged$Z == 0,], aes(color = baseGP, shape = 'Control'), size = 1.5) +
  labs(
    shape = "Treatment",
    color = "Base weights (control)"
  ) +
  scale_color_gradient2(
    low = "green", mid = "white", high = "red", midpoint = 0, 
    limits = c(allminc, allmaxc)) +#limits = c(0, 
     #          max(buffers_merged$baseGP[buffers_merged$Z == 0]))) + 
  new_scale_color() +
  geom_sf(data = buffers_merged[buffers_merged$Z == 1,], aes(color = baseGP, shape = 'Treated'), size = 1.5) +
  scale_color_gradient2(
    low = "red", mid = "white", high = "green", midpoint = 0, 
    limits = c(allmint, allmaxt)) +  #limits = c(0,
               #max(buffers_merged$baseGP[buffers_merged$Z == 1]))) + 
  coord_sf(xlim = c(-125, -65), ylim = c(25, 50), expand = FALSE) +
  labs(
    title = "Gaussian Process",
    x = "Longitude",
    y = "Latitude",
    shape = "Treatment",
    color = "Base weights (treated)"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = 'bottom'
  )

gCAR <- gCAR + theme(legend.position = "right")
gRE <- gRE + theme(legend.position = "right")
gGP <- gGP + theme(legend.position = "right")

# Combine the plots with patchwork
combined_plot <- (gRE + gCAR + gGP) +
  plot_layout(ncol = 1, guides = "collect")

png('images/baseweights.png', width = 1800, height = 2200, res = 200)
print(combined_plot)
#grid.arrange(gRE, gCAR, gGP, ncol = 1)
dev.off()

# Create Data Characteristics Table
# Data characteristics Table
covs <- X[,-c(1,12)]
summary_table <- data.frame(
  Mean = round(apply(covs, 2, mean), 3),
  SD = round(apply(covs, 2, sd), 3)
)

# Modify row names: Replace underscores with spaces and capitalize each word
rownames(summary_table) <- gsub("_", " ", rownames(summary_table)) # Replace underscores
rownames(summary_table) <- toTitleCase(rownames(summary_table))    # Capitalize words

# Convert the data frame to an xtable object
xtable_summary <- xtable(summary_table, caption = "Characteristics of the Data")

# Print the table with options for formatting
print(xtable_summary, type = "latex", caption.placement = "top", digits = 3)
mean(Z)