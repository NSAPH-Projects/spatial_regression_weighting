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

# toSci <- function(x, digits = 1) {
#   # Convert to scientific notation string (e.g., "1.9e+02")
#   s <- formatC(x, format = "e", digits = digits)
#   # Split into mantissa and exponent
#   parts <- strsplit(s, "e")[[1]]
#   mantissa <- as.numeric(parts[1])
#   exponent <- as.numeric(parts[2])
#   # Return a plotmath expression
#   bquote(.(mantissa) ~ "x" ~ 10^.(exponent))
# }

set.seed(112)

setwd('/Users/sophie/Documents/implied_weights/superfunds/')

source('../impliedweights_randeffs.R')

states <- st_read('/Users/sophie/Downloads/tl_2010_us_state10/tl_2010_us_state10.shp')
counties <- st_read('/Users/sophie/Downloads/tl_2018_us_county/tl_2018_us_county.shp')
counties <- counties %>% filter(STATEFP == '34')

# Read in preprocessed data
load('data/preprocessed_superfunds_nj.RData')
n <- nrow(buffers)
buffer_centroids <- st_centroid(buffers)

out <- compute_allSigmainv(adjacency_matrix = adjacency_matrix,
                           statefactor = clusters,
                           sig2gam = 1,
                           sig2eps = 1,
                           phi = 0.9,
                           distmat = dmat/10000)
# RE
# Create three vectors with diff smoothness
mapping <- data.frame(cluster = unique(clusters), 
                      cluster_numeric = sample(1:length(unique(clusters)), replace = F))
# Replace clusters with their corresponding numeric values
clusters_numeric <- mapping$cluster_numeric[match(clusters, mapping$cluster)]

# Now add noise to clusters_numeric
# x1 <- clusters_numeric #+ rnorm(length(clusters), sd = 1)
# x2 <- clusters_numeric + rnorm(length(clusters), sd = 15)
# x3 <- clusters_numeric + rnorm(length(clusters), sd = 30)
E <- eigen(out$Sinvs$RE)
kmax <- min(which(round(E$values,3) == 0))
x1 <- E$vectors[,(kmax-5):(kmax)] %*% rnorm(6)
x2 <- E$vectors[,(kmax-15):(kmax)] %*% rnorm(16)
x3 <- E$vectors[,(kmax-20):(kmax)] %*% rnorm(21)
eigvals_S <- 1/E$values[E$values > 0.000001] # 1 over the nonzero eigenvalues
lambdamax <- max(eigvals_S)
s1 <- smoothness(x = x1, Sigmainv = out$Sinvs$RE, lambdamax = lambdamax)
s2 <- smoothness(x = x2, Sigmainv = out$Sinvs$RE, lambdamax = lambdamax)
s3 <- smoothness(x = x3, Sigmainv = out$Sinvs$RE, lambdamax = lambdamax)
s1;s2;s3

# Compute cluster sizes and means
nk <- sapply(unique_clusters, function(k) sum(clusters == k))
xbark <- sapply(unique_clusters, function(k) mean(x1[clusters == k]))
K <- length(xbark)
pairwise_diff <- outer(xbark, xbark, FUN = function(a, b) (a - b)^2)
1/(sum(pairwise_diff) / (2 * K))
s1

nk <- sapply(unique_clusters, function(k) sum(clusters == k))
xbark <- sapply(unique_clusters, function(k) mean(x2[clusters == k]))
K <- length(xbark)
pairwise_diff <- outer(xbark, xbark, FUN = function(a, b) (a - b)^2)
1/(sum(pairwise_diff) / (2 * K))
s2

nk <- sapply(unique_clusters, function(k) sum(clusters == k))
xbark <- sapply(unique_clusters, function(k) mean(x3[clusters == k]))
K <- length(xbark)
pairwise_diff <- outer(xbark, xbark, FUN = function(a, b) (a - b)^2)
1/(sum(pairwise_diff) / (2 * K))
s3


# Plot x1, x2, x3 on map
buffers_merged <- cbind(buffer_centroids, x1, x2, x3)
gx1 <- ggplot() +
  geom_sf(data = counties, fill = NA, color = "darkgray", linetype = "solid", size = 0.5) +
  geom_sf(data = buffers_merged, 
          aes(color = x1), size = 1.5) +
  scale_color_viridis_c() +
  labs(
    title = bquote(xi (x * "," ~ S) == .(round(s1,2))),
    x = "Longitude",
    y = "Latitude",
    color = "x1"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )
gx2 <- ggplot() +
  geom_sf(data = counties, fill = NA, color = "darkgray", linetype = "solid", size = 0.5) +
  geom_sf(data = buffers_merged, 
          aes(color = x2), size = 1.5) +
  scale_color_viridis_c() +
  labs(
    title = bquote(xi (x * "," ~ S) == .(round(s2,2))),
    x = "Longitude",
    y = "Latitude",
    color = "x2"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )
gx3 <- ggplot() +
  geom_sf(data = counties, fill = NA, color = "darkgray", linetype = "solid", size = 0.5) +
  geom_sf(data = buffers_merged, 
          aes(color = x3), size = 1.5) +
  scale_color_viridis_c() +
  labs(
    title = bquote(xi (x * "," ~ S) == .(round(s3,2))),
    x = "Longitude",
    y = "Latitude",
    color = "x3"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )
png('images/smoothness_RE_nj.png', width = 1000, height = 500, res = 150)
# Arrange the plots in a grid
plots <- plot_grid(
  gx1, gx2, gx3,
  ncol = 3,
  align = "hv"
)

# Add a title above the grid
final_plot <- plot_grid(
  ggdraw() + draw_label("Random Effects", fontface = "bold", size = 16, hjust = 0.5),
  plots,
  ncol = 1,
  rel_heights = c(0.1, 1) # Allocate space for the title
)
print(final_plot)
dev.off()

# Sigmainvs$RE <- as.matrix(Sigmainvs$RE)
# Y <- buffers_outcome$Y
# # Impute the two missing values with the mean
# Y[which(is.na(Y))] <- mean(Y, na.rm = TRUE)
# cal_deltas <- calibrate_deltas(Y = Y, # for now for demonstration of sensitivity 
#                                Z = Z, 
#                                X = X[,colnames(X) != "geometry" & colnames(X) != "percent_indigenous"], 
#                                Sigmainv = Sigmainvs$RE)
# p_min1 <- cal_deltas[which(cal_deltas$Covariate == 'percent_black'),]
# p_min2 <- cal_deltas[which(cal_deltas$Covariate == 'population_density'),]
# p_max1 <- cal_deltas[which(cal_deltas$Covariate == 'percent_hispanic'),]
# p_max2 <- cal_deltas[which(cal_deltas$Covariate == 'median_household_income'),]
# delta1s <- seq(0.9*min(cal_deltas$delta1), 1.1*max(cal_deltas$delta1), length.out = 100)
# delta2s <- seq(0.9*min(cal_deltas$delta2), 1.1*max(cal_deltas$delta2), length.out = 100)
# bias_bound <- abs_bias_bound(delta1 = delta1s, 
#                              delta2 = delta2s,
#                              Z = Z, 
#                              X = X[,colnames(X) != "geometry" & colnames(X) != "percent_indigenous"], 
#                              Sigmainv = Sigmainvs$RE)
# 
# # Generate the contour plot with annotations
# annotation_points <- data.frame(
#   Covariatename = c(p_min1$Covariate, p_min2$Covariate, p_max2$Covariate),
#   x = c(p_min1$delta1, p_min2$delta1, p_max2$delta1),
#   y = c(p_min1$delta2, p_min2$delta2, p_max2$delta2)
# )
# # Generate the contour plot with annotations
# annotation_points$Covariatename <- gsub("_", " ", annotation_points$Covariatename)
# annotation_points$y_label <- annotation_points$y + 0.02 # Adjust this value as needed
# annotation_points$x_label <- annotation_points$x + 0.02  # Adjust this value as needed
# 
# ggplot(bias_bound, aexi (x = delta1, y = delta2)) +
#   geom_contour(aes(z = bias_bound), colour = "black") +
#   metR::geom_text_contour(aes(z = bias_bound), stroke = 0.15, skip = 0, color = 'red', size = 5) +  # Center labels
#   geom_point(data = annotation_points, aexi (x = x, y = y), color = "blue", size = 3) +  # Add points
#   geom_text(data = annotation_points, aexi (x = x_label, y = y_label, label = Covariatename),
#             color = "blue", size = 5, fontface = "bold") +  # Use shifted y position
#   labs(
#     title = "Absolute Conditional Bias",
#     x = expression("Smoothness of Unmeasured Confounder (" * delta[1] * ")"),  # Greek delta1
#     y = expression("Strength of Unmeasured Confounder - Outcome Relationship (" * delta[2] * ")"),  # Greek delta2
#     fill = expression(E(hat(tau)[GLS] ~ "|" ~ Z * "," ~ X * "," ~ U) - tau)
#   ) +
#   coord_cartesian(
#     xlim = range(bias_bound$delta1), 
#     ylim = range(bias_bound$delta2)
#   ) +  # Adjust axis limits
#   theme_minimal(base_size = 14)

# CAR
# Create three vectors with diff smoothness
E <- eigen(out$Sinvs$CAR)
kmax <- min(which(round(E$values,3) == 0))
x1 <- E$vectors[,(kmax-10):(kmax)] %*% rnorm(11)
x2 <- E$vectors[,(kmax-30):(kmax)] %*% rnorm(31)
x3 <- E$vectors[,(kmax-100):(kmax)] %*% rnorm(101)
eigvals_S <- 1/E$values[E$values > 0.000001] # 1 over the nonzero eigenvalues
lambdamax <- max(eigvals_S)
s1 <- smoothness(x = x1, Sigmainv = out$Sinvs$CAR, lambdamax = lambdamax)
s2 <- smoothness(x = x2, Sigmainv = out$Sinvs$CAR, lambdamax = lambdamax)
s3 <- smoothness(x = x3, Sigmainv = out$Sinvs$CAR, lambdamax = lambdamax)
s1;s2;s3
# Compute 1/2 sum_{ij} A_{ij}(x_i - x_j)^2
sumtot = 0
for (i in 1:n){
  for (j in 1:n){
    sumtot = sumtot + (1/2)*adjacency_matrix[i,j] * (x1[i] - x1[j])^2
  }
}
1/sumtot
s1

sumtot = 0
for (i in 1:n){
  for (j in 1:n){
    sumtot = sumtot + (1/2)*adjacency_matrix[i,j] * (x2[i] - x2[j])^2
  }
}
1/sumtot
s2

sumtot = 0
for (i in 1:n){
  for (j in 1:n){
    sumtot = sumtot + (1/2)*adjacency_matrix[i,j] * (x3[i] - x3[j])^2
  }
}
1/sumtot
s3

buffers_merged <- cbind(buffer_centroids, x1, x2, x3)
gx1 <- ggplot() +
  geom_sf(data = counties, fill = NA, color = "darkgray", linetype = "solid", size = 0.5) +
  geom_sf(data = buffers_merged, 
          aes(color = x1), size = 1.5) +
  scale_color_viridis_c() +
  labs(
    title = bquote(xi (x * "," ~ S) == .(round(s1,2))),
    x = "Longitude",
    y = "Latitude",
    color = "x1"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )
gx2 <- ggplot() +
  geom_sf(data = counties, fill = NA, color = "darkgray", linetype = "solid", size = 0.5) +
  geom_sf(data = buffers_merged, 
          aes(color = x2), size = 1.5) +
  scale_color_viridis_c() +
  labs(
    title = bquote(xi (x * "," ~ S) == .(round(s2,2))),
    x = "Longitude",
    y = "Latitude",
    color = "x2"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )
gx3 <- ggplot() +
  geom_sf(data = counties, fill = NA, color = "darkgray", linetype = "solid", size = 0.5) +
  geom_sf(data = buffers_merged, 
          aes(color = x3), size = 1.5) +
  scale_color_viridis_c() +
  labs(
    title = bquote(xi (x * "," ~ S) == .(round(s3,2))),
    x = "Longitude",
    y = "Latitude",
    color = "x3"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )
png('images/smoothness_CAR_nj.png', width = 1000, height = 500, res = 150)
# Arrange the plots in a grid
plots <- plot_grid(
  gx1, gx2, gx3,
  ncol = 3,
  align = "hv"
)

# Add a title above the grid
final_plot <- plot_grid(
  ggdraw() + draw_label("CAR", fontface = "bold", size = 16, hjust = 0.5),
  plots,
  ncol = 1,
  rel_heights = c(0.1, 1) # Allocate space for the title
)
print(final_plot)
dev.off()

# cal_deltas <- calibrate_deltas(Y = Y, # for now for demonstration of sensitivity 
#                                Z = Z, 
#                                X = X[,colnames(X) != "geometry" & colnames(X) != "percent_indigenous"], 
#                                Sigmainv = Sigmainvs$CAR)
# p_min1 <- cal_deltas[which(cal_deltas$Covariate == 'percent_black'),]
# p_min2 <- cal_deltas[which(cal_deltas$Covariate == 'population_density'),]
# p_max1 <- cal_deltas[which(cal_deltas$Covariate == 'percent_hispanic'),]
# p_max2 <- cal_deltas[which(cal_deltas$Covariate == 'median_household_income'),]
# delta1s <- seq(0.9*min(cal_deltas$delta1), 1.1*max(cal_deltas$delta1), length.out = 100)
# delta2s <- seq(0.9*min(cal_deltas$delta2), 1.1*max(cal_deltas$delta2), length.out = 100)
# bias_bound <- abs_bias_bound(delta1 = delta1s, 
#                              delta2 = delta2s,
#                              Z = Z, 
#                              X = X[,colnames(X) != "geometry" & colnames(X) != "percent_indigenous"], 
#                              Sigmainv = Sigmainvs$CAR)
# 
# # Generate the contour plot with annotations
# annotation_points <- data.frame(
#   Covariatename = c(p_min1$Covariate, p_min2$Covariate, p_max2$Covariate),
#   x = c(p_min1$delta1, p_min2$delta1, p_max2$delta1),
#   y = c(p_min1$delta2, p_min2$delta2, p_max2$delta2)
# )
# # Generate the contour plot with annotations
# annotation_points$Covariatename <- gsub("_", " ", annotation_points$Covariatename)
# annotation_points$y_label <- annotation_points$y + 0.02 # Adjust this value as needed
# annotation_points$x_label <- annotation_points$x + 0.02  # Adjust this value as needed
# 
# ggplot(bias_bound, aexi (x = delta1, y = delta2)) +
#   geom_contour(aes(z = bias_bound), colour = "black") +
#   metR::geom_text_contour(aes(z = bias_bound), stroke = 0.15, skip = 0, color = 'red', size = 5) +  # Center labels
#   geom_point(data = annotation_points, aexi (x = x, y = y), color = "blue", size = 3) +  # Add points
#   geom_text(data = annotation_points, aexi (x = x_label, y = y_label, label = Covariatename),
#             color = "blue", size = 5, fontface = "bold") +  # Use shifted y position
#   labs(
#     title = "Absolute Conditional Bias",
#     x = expression("Smoothness of Unmeasured Confounder (" * delta[1] * ")"),  # Greek delta1
#     y = expression("Strength of Unmeasured Confounder - Outcome Relationship (" * delta[2] * ")"),  # Greek delta2
#     fill = expression(E(hat(tau)[GLS] ~ "|" ~ Z * "," ~ X * "," ~ U) - tau)
#   ) +
#   coord_cartesian(
#     xlim = range(bias_bound$delta1), 
#     ylim = range(bias_bound$delta2)
#   ) +  # Adjust axis limits
#   theme_minimal(base_size = 14)


# GP
# Create three vectors with diff smoothness
E <- eigen(out$Sinvs$GP)
x1 <- E$vectors[,(n-1):(n)] %*% c(0.1,1)
x2 <- E$vectors[,(n-10):(n)] %*% seq(0,1,0.1)
x3 <- E$vectors[,(n-100):(n)] %*% seq(0,1,0.01)
eigvals_S <- 1/E$values # 1 over the nonzero eigenvalues
lambdamax <- max(eigvals_S)
s1 <- smoothness(x = x1, out$Sinvs$GP, lambdamax = lambdamax)
s2 <- smoothness(x = x2, out$Sinvs$GP, lambdamax = lambdamax)
s3 <- smoothness(x = x3, out$Sinvs$GP, lambdamax = lambdamax)
s1;s2;s3

buffers_merged <- cbind(buffer_centroids, x1, x2, x3)
gx1 <- ggplot() +
  geom_sf(data = counties, fill = NA, color = "darkgray", linetype = "solid", size = 0.5) +
  geom_sf(data = buffers_merged, 
          aes(color = x1), size = 1.5) +
  scale_color_viridis_c() +
  labs(
    title = bquote(xi (x * "," ~ S) == .(round(s1,2))),
    x = "Longitude",
    y = "Latitude",
    color = "x1"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )
gx2 <- ggplot() +
  geom_sf(data = counties, fill = NA, color = "darkgray", linetype = "solid", size = 0.5) +
  geom_sf(data = buffers_merged, 
          aes(color = x2), size = 1.5) +
  scale_color_viridis_c() +
  labs(
    title = bquote(xi (x * "," ~ S) == .(round(s2,4))),
    x = "Longitude",
    y = "Latitude",
    color = "x2"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )
gx3 <- ggplot() +
  geom_sf(data = counties, fill = NA, color = "darkgray", linetype = "solid", size = 0.5) +
  geom_sf(data = buffers_merged, 
          aes(color = x3), size = 1.5) +
  scale_color_viridis_c() +
  labs(
    title = bquote(xi (x * "," ~ S) == .(round(s3,4))),
    x = "Longitude",
    y = "Latitude",
    color = "x3"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )
png('images/smoothness_GP_nj.png', width = 1000, height = 500, res = 150)
# Arrange the plots in a grid
plots <- plot_grid(
  gx1, gx2, gx3,
  ncol = 3,
  align = "hv"
)

# Add a title above the grid
final_plot <- plot_grid(
  ggdraw() + draw_label("GP", fontface = "bold", size = 16, hjust = 0.5),
  plots,
  ncol = 1,
  rel_heights = c(0.1, 1) # Allocate space for the title
)
print(final_plot)
dev.off()

# cal_deltas <- calibrate_deltas(Y = Y, # for now for demonstration of sensitivity 
#                                Z = Z, 
#                                X = X[,colnames(X) != "geometry" & colnames(X) != "percent_indigenous"], 
#                                Sigmainv = Sigmainvs$GP)
# p_min1 <- cal_deltas[which(cal_deltas$Covariate == 'percent_black'),]
# p_min2 <- cal_deltas[which(cal_deltas$Covariate == 'population_density'),]
# p_max1 <- cal_deltas[which(cal_deltas$Covariate == 'percent_hispanic'),]
# p_max2 <- cal_deltas[which(cal_deltas$Covariate == 'median_household_income'),]
# delta1s <- seq(0.9*min(cal_deltas$delta1), 1.1*max(cal_deltas$delta1), length.out = 100)
# delta2s <- seq(0.9*min(cal_deltas$delta2), 1.1*max(cal_deltas$delta2), length.out = 100)
# bias_bound <- abs_bias_bound(delta1 = delta1s, 
#                              delta2 = delta2s,
#                              Z = Z, 
#                              X = X[,colnames(X) != "geometry" & colnames(X) != "percent_indigenous"], 
#                              Sigmainv = Sigmainvs$GP)
# 
# # Generate the contour plot with annotations
# annotation_points <- data.frame(
#   Covariatename = c(p_min1$Covariate, p_min2$Covariate, p_max2$Covariate),
#   x = c(p_min1$delta1, p_min2$delta1, p_max2$delta1),
#   y = c(p_min1$delta2, p_min2$delta2, p_max2$delta2)
# )
# # Generate the contour plot with annotations
# annotation_points$Covariatename <- gsub("_", " ", annotation_points$Covariatename)
# annotation_points$y_label <- annotation_points$y + 0.02 # Adjust this value as needed
# annotation_points$x_label <- annotation_points$x + 0.02  # Adjust this value as needed
# 
# ggplot(bias_bound, aexi (x = delta1, y = delta2)) +
#   geom_contour(aes(z = bias_bound), colour = "black") +
#   metR::geom_text_contour(aes(z = bias_bound), stroke = 0.15, skip = 0, color = 'red', size = 5) +  # Center labels
#   geom_point(data = annotation_points, aexi (x = x, y = y), color = "blue", size = 3) +  # Add points
#   geom_text(data = annotation_points, aexi (x = x_label, y = y_label, label = Covariatename),
#             color = "blue", size = 5, fontface = "bold") +  # Use shifted y position
#   labs(
#     title = "Absolute Conditional Bias",
#     x = expression("Smoothness of Unmeasured Confounder (" * delta[1] * ")"),  # Greek delta1
#     y = expression("Strength of Unmeasured Confounder - Outcome Relationship (" * delta[2] * ")"),  # Greek delta2
#     fill = expression(E(hat(tau)[GLS] ~ "|" ~ Z * "," ~ X * "," ~ U) - tau)
#   ) +
#   coord_cartesian(
#     xlim = range(bias_bound$delta1), 
#     ylim = range(bias_bound$delta2)
#   ) +  # Adjust axis limits
#   theme_minimal(base_size = 14)
# 
