library(tidyr)
library(tidyverse)
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
library(MASS)
library(zoo)

source('../funcs.R')

states <- st_read('../data/tl_2010_us_state10/tl_2010_us_state10.shp')

# Read in preprocessed data
load('../data/preprocessed_superfunds.RData')
n <- nrow(buffers)
us_states <- ne_states(country = "United States of America", 
                       returnclass = "sf")
us_outline <- ne_countries(scale = "medium", country = "United States of America", returnclass = "sf")

################# PLOT BINARY TREATMENT ######################
# Plot Z on the map
# order buffers by Z
buffer_centroids <- st_centroid(buffers)
buffer_centroids_ordered <- buffer_centroids[order(buffers$Z),]
g <- ggplot() +
  # Add map outline
  geom_sf(data = us_outline, fill = NA, color = "black", linetype = "solid") +
  # Plot centroids with both color and shape mapped to factor(Z)
  geom_sf(data = buffer_centroids_ordered, aes(color = factor(Z), shape = factor(Z)), size = 3) +
  labs(
    title = "Superfund Sites that were cleaned up and removed from the National Priority List between 2001 and 2015",
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
    values = c("0" = "orange", "1" = "dodgerblue"),
    name = "Treatment (Cleanup)",
    labels = c("Control", "Treated")
  ) +
  scale_shape_manual(
    values = c("0" = 16, "1" = 17), # Circle for control, triangle for treated
    name = "Treatment (Cleanup)",
    labels = c("Control", "Treated")
  )
png('images/treatment.png', width = 1500, height = 900, res = 160)
plot_with_insets(g)
dev.off()

################# CALCULATE IMPLIED WEIGHTS ##################################
# Scale X
X[,2:12] <- scale(X[,2:12])
X <- as.matrix(X)

# GP model
kappa <- 0.1 # 0.2 # as you increase, becomes more extreme
rangec <- 300
phic <- rangec/(2*sqrt(kappa))
S <- geoR::matern(u =  dmat/10000, phi = phic, kappa = kappa)
Sigma <- diag(n) + 10*S
Sigmainv <- solve(Sigma)
GPiw <- impliedweightsgeneral(X = X,
                              Z = Z,
                              Sigmainv = Sigmainv)
sum(GPiw[Z == 1]);sum(GPiw[Z == 0])

# RE model
statefactor <- factor(clusters)
Jks = list()
for (k in unique(statefactor)){
  nk = sum(statefactor == k)
  Jks[[length(Jks) + 1]] = matrix(1, nrow = nk, ncol = nk)
}
S = bdiag(Jks)
Sigma <- diag(n) + 10*S
Sigmainv <- solve(Sigma)
REiw <- impliedweightsgeneral(X = X,
                              Z = Z,
                              Sigmainv = Sigmainv)
sum(REiw[Z == 1]);sum(REiw[Z == 0])

# CAR model
L <- diag(rowSums(adjacency_matrix)) - adjacency_matrix
S <- ginv(as.matrix(L))
Sigma <- diag(n) + 10*S
Sigmainv <- solve(Sigma)
CARiw <- impliedweightsgeneral(X = X,
                               Z = Z,
                               Sigmainv = Sigmainv)
sum(CARiw[Z == 1]);sum(CARiw[Z == 0])

################## PLOT IMPLIED WEIGHTS ##################
# Merge the Implied weights with the data
buffer_centroids <- st_centroid(buffers)
buffers_merged <- cbind(buffer_centroids, REiw, CARiw, GPiw)

# Plot baseCAR
buffers_merged_geo <- st_transform(buffers_merged, crs = 4326)

buffers_merged_geo <- buffers_merged_geo[order(buffers_merged_geo$CARiw),]

pts <- buffers_merged_geo %>%
  st_centroid() %>%
  cbind(st_coordinates(.)) %>%
  st_set_geometry(NULL)
g <- ggplot() +
  scale_shape_manual(values = c("0" = 21, "1" = 24),
                     guide = guide_legend(order = 1)) + 
  
  # Base map with county outlines
  geom_sf(data = states, fill = NA, color = "darkgray", linetype = "solid", size = 0.5) +
  geom_point(
    data     = filter(pts, Z == 0),
    aes(X, Y, fill = CARiw, shape = factor(Z)),
    size     = 1.5, 
    stroke   = 0.25, 
    color    = "black",
    position = position_nudge(x = 0.5, y = 0.5)  # ← correct placement
  ) +
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "orange", midpoint = 0, 
                       breaks = c(min(buffers_merged_geo[buffers_merged_geo$Z == 0,]$CARiw, na.rm = TRUE), 
                                  max(buffers_merged_geo[buffers_merged_geo$Z == 0,]$CARiw, na.rm = TRUE)),
                       labels = function(x) sprintf("%.3f", x),
                       name = "Implied weights (Control)",
                       guide = guide_colorbar(order = 2)) +  
  
  new_scale_fill() +  
  
  # Treated group (Z == 1): Green triangle with black outline
  geom_sf(data = buffers_merged_geo[buffers_merged_geo$Z == 1,], 
          aes(fill = CARiw, shape = factor(Z)), 
          size = 1.5, stroke = 0.25, color = "black") +  
  
  scale_fill_gradient2(low = "orange", mid = "white", high = "dodgerblue", midpoint = 0, 
                       breaks = c(min(buffers_merged_geo[buffers_merged_geo$Z == 1,]$CARiw, na.rm = TRUE), 
                                  max(buffers_merged_geo[buffers_merged_geo$Z == 1,]$CARiw, na.rm = TRUE)),
                       labels = function(x) sprintf("%.3f", x),
                       name = "Implied weights (Treated)",
                       guide = guide_colorbar(order = 3)) +  
  labs(
    title = "Intrinsic Conditional Autoregressive",
    x = "Longitude",
    y = "Latitude",
    shape = "Treatment",
    fill = "Implied weights"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.key.height = unit(0.2, "cm")
  )

gCAR <- plot_with_insets(g)
# Plot  FE
buffers_merged_geo <- buffers_merged_geo[order(buffers_merged_geo$REiw),]
g <- ggplot() +
  
  # Custom shape: 21 (circle), 24 (triangle)
  scale_shape_manual(values = c("0" = 21, "1" = 24),
                     guide = guide_legend(order = 1)) + 
  
  # Base map with county outlines
  geom_sf(data = states, fill = NA, color = "darkgray", linetype = "solid", size = 0.5) +
  
  # Control group (Z == 0): Red circle with black outline
  geom_point(
    data     = filter(pts, Z == 0),
    aes(X, Y, fill = REiw, shape = factor(Z)),
    size     = 1.5, 
    stroke   = 0.25, 
    color    = "black",
    position = position_nudge(x = 0.5, y = 0.5)  # ← correct placement
  ) +
  
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "orange", midpoint = 0, 
                       breaks = c(min(buffers_merged_geo[buffers_merged_geo$Z == 0,]$REiw, na.rm = TRUE), 
                                  max(buffers_merged_geo[buffers_merged_geo$Z == 0,]$REiw, na.rm = TRUE)),
                       labels = function(x) sprintf("%.3f", x),
                       name = "Implied weights (Control)",
                       guide = guide_colorbar(order = 2)) +  
  
  new_scale_fill() +  
  
  # Treated group (Z == 1): Green triangle with black outline
  geom_sf(data = buffers_merged_geo[buffers_merged_geo$Z == 1,], 
          aes(fill = REiw, shape = factor(Z)),
          size = 1.5, stroke = 0.25, color = "black") +  
  
  scale_fill_gradient2(low = "orange", mid = "white", high = "dodgerblue", midpoint = 0, 
                       breaks = c(min(buffers_merged_geo[buffers_merged_geo$Z == 1,]$REiw, na.rm = TRUE), 
                                  max(buffers_merged_geo[buffers_merged_geo$Z == 1,]$REiw, na.rm = TRUE)),
                       labels = function(x) sprintf("%.3f", x),
                       name = "Implied weights (Treated)",
                       guide = guide_colorbar(order = 3)) +  
  labs(
    title = "Random Effects",
    x = "Longitude",
    y = "Latitude",
    shape = "Treatment",
    fill = "Implied weights"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.key.height = unit(0.2, "cm")
  )

gRE <- plot_with_insets(g)
# Plot  GP
buffers_merged_geo <- buffers_merged_geo[order(buffers_merged_geo$GPiw),]
g <- ggplot() +
  # Custom shape: 21 (circle), 24 (triangle)
  scale_shape_manual(values = c("0" = 21, "1" = 24),
                     guide = guide_legend(order = 1)) + 
  
  # Base map with county outlines
  geom_sf(data = states, fill = NA, color = "darkgray", linetype = "solid", size = 0.5) +
  
  # Control group (Z == 0): Red circle with black outline
  geom_point(
    data     = filter(pts, Z == 0),
    aes(X, Y, fill = GPiw, shape = factor(Z)),
    size     = 1.5, 
    stroke   = 0.25, 
    color    = "black",
    position = position_nudge(x = 0.5, y = 0.5)  # ← correct placement
  ) +

  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "orange", midpoint = 0, 
                       breaks = c(min(buffers_merged_geo[buffers_merged_geo$Z == 0,]$GPiw, na.rm = TRUE), 
                                  max(buffers_merged_geo[buffers_merged_geo$Z == 0,]$GPiw, na.rm = TRUE)),
                       labels = function(x) sprintf("%.3f", x),
                       name = "Implied weights (Control)",
                       guide = guide_colorbar(order = 2)) +  
  
  new_scale_fill() +  
  
  # Treated group (Z == 1): Green triangle with black outline
  geom_sf(data = buffers_merged_geo[buffers_merged_geo$Z == 1,],
          aes(fill = GPiw, shape = factor(Z)),
          size = 1.5, stroke = 0.25, color = "black") +
  
  scale_fill_gradient2(low = "orange", mid = "white", high = "dodgerblue", midpoint = 0,
                       breaks = c(min(buffers_merged_geo[buffers_merged_geo$Z == 1,]$GPiw, na.rm = TRUE), 
                                  max(buffers_merged_geo[buffers_merged_geo$Z == 1,]$GPiw, na.rm = TRUE)),
                       labels = function(x) sprintf("%.3f", x),
                       name = "Implied weights (Treated)",
                       guide = guide_colorbar(order = 3)) +
  labs(
    title = "Gaussian Process",
    x = "Longitude",
    y = "Latitude",
    shape = "Treatment",
    fill = "Implied weights"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.key.height = unit(0.2, "cm")
  )
gGP <- plot_with_insets(g)

png('images/impliedweights_us.png', width = 1500, height = 2000, res = 250)
grid.arrange(gRE, gCAR, gGP, ncol = 1)
dev.off()


######################################## UNMEASURED SPATIAL CONFOUNDER EXAMPLES #########################################
buffers_merged <- st_centroid(buffers)

set.seed(100)
kappa <- 0.1 # 0.2 # as you increase, becomes more extreme
rangec <- 300
phic <- rangec/(2*sqrt(kappa))
S <- geoR::matern(u =  dmat/10000, phi = phic, kappa = kappa)
E <- eigen(S)
V <- E$vectors
Sigma <- diag(n) + 10*S
Sigmainv <- solve(Sigma)
GPiw <- impliedweightsgeneral(X = X,
                              Z = Z,
                              Sigmainv = Sigmainv)
sum(GPiw[Z == 1]);sum(GPiw[Z == 0])

# Squared Eigenvector imbalance
Es <- eigen(S)
m = 0
while (m < 0.9 | m > 0.95){
  buffers_merged$GP_u1 <- V[,1:n] %*% rnorm(ncol(V), mean = 0, sd = 1/(1:ncol(V))^4) #runif(ncol(V), min = 0, max = 1)
  buffers_merged$GP_u1 <- buffers_merged$GP_u1 - mean(buffers_merged$GP_u1)
  buffers_merged$GP_u1 <- buffers_merged$GP_u1/norm(buffers_merged$GP_u1, type = '2')
  m <- moran_adhoc(buffers_merged$GP_u1, Wmat = S, coeff = 1/Es$values[2])
}
print(paste0('Moran I is ', m, ' and imbalance is ', sum(GPiw*(2*Z-1)*buffers_merged$GP_u1)))

m = 0
while (m < 0.45 | m > 0.55){
  buffers_merged$GP_u2 <- V[,1:n] %*% rnorm(ncol(V), mean = 0, sd = 1/(1:ncol(V))^2) #runif(ncol(V), min = 0, max = 1)
  buffers_merged$GP_u2 <- buffers_merged$GP_u2 - mean(buffers_merged$GP_u2)
  buffers_merged$GP_u2 <- buffers_merged$GP_u2/norm(buffers_merged$GP_u2, type = '2')
  m <- moran_adhoc(buffers_merged$GP_u2, Wmat = S, coeff = 1/Es$values[2])
}
print(paste0('Moran I is ', m, ' and imbalance is ', sum(GPiw*(2*Z-1)*buffers_merged$GP_u2)))

m = 1
while (m > 0.1){
  buffers_merged$GP_u3 <- V[,1:n] %*% rnorm(ncol(V), mean = 0, sd = sqrt(1/(1:ncol(V)))) #runif(ncol(V), min = 0, max = 1)
  buffers_merged$GP_u3 <- buffers_merged$GP_u3 - mean(buffers_merged$GP_u3)
  buffers_merged$GP_u3 <- buffers_merged$GP_u3/norm(buffers_merged$GP_u3, type = '2')
  m <- moran_adhoc(buffers_merged$GP_u3, Wmat = S, coeff = 1/Es$values[2])
}
print(paste0('Moran I is ', m, ' and imbalance is ', sum(GPiw*(2*Z-1)*buffers_merged$GP_u3)))

# PLOT each on the map
png('images/gp1.png', width = 1300, height = 1000, res = 250)
g <- ggplot() +
  geom_sf(data = states, fill = NA, color = "darkgray", linetype = "solid", size = 0.5) +
  geom_sf(data = buffers_merged, 
          aes(fill = GP_u1, shape = factor(Z)), size = 4, stroke = 0) +
  # drop the shape legend:
  scale_shape_manual(values = c("0" = 21, "1" = 24),
                     guide = "none") +
  # relabel the fill legend:
  scale_fill_viridis_c(name = "U") +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    panel.grid    = element_blank(),
    axis.text     = element_blank(),
    axis.ticks    = element_blank(),
    axis.title    = element_blank(),
    plot.title    = element_text(hjust = 0.5),
    legend.position = "top",
    legend.key.width = unit(1.5, "cm")
  )
plot_with_insets(g)
dev.off()

png('images/gp2.png', width = 1300, height = 1000, res = 250)
g <- ggplot() +
  geom_sf(data = states, fill = NA, color = "darkgray", linetype = "solid", size = 0.5) +
  geom_sf(data = buffers_merged, 
          aes(fill = GP_u2, shape = factor(Z)), size = 4, stroke = 0) +
  # drop the shape legend:
  scale_shape_manual(values = c("0" = 21, "1" = 24),
                     guide = "none") +
  # relabel the fill legend:
  scale_fill_viridis_c(name = "U") +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    panel.grid    = element_blank(),
    axis.text     = element_blank(),
    axis.ticks    = element_blank(),
    axis.title    = element_blank(),
    plot.title    = element_text(hjust = 0.5),
    legend.position = "top",
    legend.key.width = unit(1.5, "cm")
  )
plot_with_insets(g)
dev.off()

png('images/gp3.png', width = 1300, height = 1000, res = 250)
g <- ggplot() +
  geom_sf(data = states, fill = NA, color = "darkgray", linetype = "solid", size = 0.5) +
  geom_sf(data = buffers_merged, 
          aes(fill = GP_u3, shape = factor(Z)), size = 4, stroke = 0) +
  # drop the shape legend:
  scale_shape_manual(values = c("0" = 21, "1" = 24),
                     guide = "none") +
  # relabel the fill legend:
  scale_fill_viridis_c(name = "U") +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    panel.grid    = element_blank(),
    axis.text     = element_blank(),
    axis.ticks    = element_blank(),
    axis.title    = element_blank(),
    plot.title    = element_text(hjust = 0.5),
    legend.position = "top",
    legend.key.width = unit(1.5, "cm")
  )
plot_with_insets(g)
dev.off()

statefactor <- factor(clusters)
Jks = list()
for (k in unique(statefactor)){
  nk = sum(statefactor == k)
  Jks[[length(Jks) + 1]] = matrix(1, nrow = nk, ncol = nk)
}
S = bdiag(Jks)
V <- eigen(S)$vectors
Sigma <- diag(n) + 10*S
Sigmainv <- solve(Sigma)
REiw <- impliedweightsgeneral(X = X,
                              Z = Z,
                              Sigmainv = Sigmainv)

# Squared Eigenvector imbalance
Es <- eigen(S)
m = 0
while (m < 0.9 | m > 0.95) {
  buffers_merged$RE_u1 <- V[,1:n] %*% c(rnorm(51, mean = 0, sd = 1/(1:ncol(V))), rep(0, ncol(V)-51)) #runif(ncol(V), min = 0, max = 1)
  buffers_merged$RE_u1 <- buffers_merged$RE_u1 - mean(buffers_merged$RE_u1)
  buffers_merged$RE_u1 <- buffers_merged$RE_u1/norm(buffers_merged$RE_u1, type = '2')
  m <- moran_adhoc(buffers_merged$RE_u1, Wmat = S, coeff = 1/Es$values[1])
}
print(paste0('Moran I is ', m, ' and imbalance is ', sum(REiw*(2*Z-1)*buffers_merged$RE_u1)))

m = 0
while (m < 0.45 | m > 0.55){
  buffers_merged$RE_u2 <- V[,1:n] %*% c(rnorm(51, mean = 0, sd = 1/(1:ncol(V))), rep(0.01, ncol(V)-51)) #runif(ncol(V), min = 0, max = 1)
  buffers_merged$RE_u2 <- buffers_merged$RE_u2 - mean(buffers_merged$RE_u2)
  buffers_merged$RE_u2 <- buffers_merged$RE_u2/norm(buffers_merged$RE_u2, type = '2')
  m <- moran_adhoc(buffers_merged$RE_u2, Wmat = S, coeff = 1/Es$values[1])
}
print(paste0('Moran I is ', m, ' and imbalance is ', sum(REiw*(2*Z-1)*buffers_merged$RE_u2)))

m = 1
while (m > 0.1){
  buffers_merged$RE_u3 <- V[,1:n] %*% c(rnorm(51, mean = 0, sd = 1/(1:ncol(V))), rep(0.1, ncol(V)-51)) #runif(ncol(V), min = 0, max = 1)
  buffers_merged$RE_u3 <- buffers_merged$RE_u3 - mean(buffers_merged$RE_u3)
  buffers_merged$RE_u3 <- buffers_merged$RE_u3/norm(buffers_merged$RE_u3, type = '2')
  m <- moran_adhoc(buffers_merged$RE_u3, Wmat = S, coeff = 1/Es$values[1])
}
print(paste0('Moran I is ', m, ' and imbalance is ', sum(REiw*(2*Z-1)*buffers_merged$RE_u3)))

# PLOT each on the map
png('images/re1.png', width = 1300, height = 1000, res = 250)
g <- ggplot() +
  geom_sf(data = states, fill = NA, color = "darkgray", linetype = "solid", size = 0.5) +
  geom_sf(data = buffers_merged, 
          aes(fill = RE_u1, shape = factor(Z)), size = 4, stroke = 0) +
  # drop the shape legend:
  scale_shape_manual(values = c("0" = 21, "1" = 24),
                     guide = "none") +
  # relabel the fill legend:
  scale_fill_viridis_c(name = "U") +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    panel.grid    = element_blank(),
    axis.text     = element_blank(),
    axis.ticks    = element_blank(),
    axis.title    = element_blank(),
    plot.title    = element_text(hjust = 0.5),
    legend.position = "top",
    legend.key.width = unit(1.5, "cm")
  )
plot_with_insets(g)
dev.off()

png('images/re2.png', width = 1300, height = 1000, res = 250)
g <- ggplot() +
  geom_sf(data = states, fill = NA, color = "darkgray", linetype = "solid", size = 0.5) +
  geom_sf(data = buffers_merged, 
          aes(fill = RE_u2, shape = factor(Z)), size = 4, stroke = 0) +
  # drop the shape legend:
  scale_shape_manual(values = c("0" = 21, "1" = 24),
                     guide = "none") +
  # relabel the fill legend:
  scale_fill_viridis_c(name = "U") +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    panel.grid    = element_blank(),
    axis.text     = element_blank(),
    axis.ticks    = element_blank(),
    axis.title    = element_blank(),
    plot.title    = element_text(hjust = 0.5),
    legend.position = "top",
    legend.key.width = unit(1.5, "cm")
  )
plot_with_insets(g)
dev.off()

png('images/re3.png', width = 1300, height = 1000, res = 250)
g <- ggplot() +
  geom_sf(data = states, fill = NA, color = "darkgray", linetype = "solid", size = 0.5) +
  geom_sf(data = buffers_merged, 
          aes(fill = RE_u3, shape = factor(Z)), size = 4, stroke = 0) +
  # drop the shape legend:
  scale_shape_manual(values = c("0" = 21, "1" = 24),
                     guide = "none") +
  # relabel the fill legend:
  scale_fill_viridis_c(name = "U") +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    panel.grid    = element_blank(),
    axis.text     = element_blank(),
    axis.ticks    = element_blank(),
    axis.title    = element_blank(),
    plot.title    = element_text(hjust = 0.5),
    legend.position = "top",
    legend.key.width = unit(1.5, "cm")
  )
plot_with_insets(g)
dev.off()

sigma2 <- 1
rho2 <- 10
L <- diag(rowSums(adjacency_matrix)) - adjacency_matrix
S <- ginv(as.matrix(L))
V <- Re(eigen(S)$vectors)
Sigma <- sigma2*diag(n) + rho2*S
Sigmainv <- solve(Sigma)
CARiw <- impliedweightsgeneral(X = X,
                               Z = Z,
                               Sigmainv = Sigmainv)
sum(CARiw[Z == 1]);sum(CARiw[Z == 0])

# Squared Eigenvector imbalance
Es <- eigen(S)
m = 0
while (m < 0.9 | m > 0.95){
  buffers_merged$CAR_u1 <- V[,1:n] %*% rnorm(ncol(V), mean = 0, sd = 1/(1:ncol(V))^2) #runif(ncol(V), min = 0, max = 1)
  buffers_merged$CAR_u1 <- buffers_merged$CAR_u1 - mean(buffers_merged$CAR_u1)
  buffers_merged$CAR_u1 <- buffers_merged$CAR_u1/norm(buffers_merged$CAR_u1, type = '2')
  m <- moran_adhoc(buffers_merged$CAR_u1, Wmat = S, coeff = 1/Re(Es$values[1]))
}
print(paste0('Moran I is ', m, ' and imbalance is ', sum(CARiw*(2*Z-1)*buffers_merged$CAR_u1)))

m = 0
while (m < 0.45 | m > 0.55){
  buffers_merged$CAR_u2 <- V[,1:n] %*% rnorm(ncol(V), mean = 0, sd = 1/(1:ncol(V))) #runif(ncol(V), min = 0, max = 1)
  buffers_merged$CAR_u2 <- buffers_merged$CAR_u2 - mean(buffers_merged$CAR_u2)
  buffers_merged$CAR_u2 <- buffers_merged$CAR_u2/norm(buffers_merged$CAR_u2, type = '2')
  m <- moran_adhoc(buffers_merged$CAR_u2, Wmat = S, coeff = 1/Re(Es$values[1]))
}
print(paste0('Moran I is ', m, ' and imbalance is ', sum(CARiw*(2*Z-1)*buffers_merged$CAR_u2)))

m = 1
while (m > 0.1){
  buffers_merged$CAR_u3 <- V[,1:n] %*% rnorm(ncol(V), mean = 0, sd = sqrt(1/(1:ncol(V)))) #runif(ncol(V), min = 0, max = 1)
  buffers_merged$CAR_u3 <- buffers_merged$CAR_u3 - mean(buffers_merged$CAR_u3)
  buffers_merged$CAR_u3 <- buffers_merged$CAR_u3/norm(buffers_merged$CAR_u3, type = '2')
  m <- moran_adhoc(buffers_merged$CAR_u3, Wmat = S, coeff = 1/Re(Es$values[1]))
}
print(paste0('Moran I is ', m, ' and imbalance is ', sum(CARiw*(2*Z-1)*buffers_merged$CAR_u3)))

png('images/car1.png', width = 1300, height = 1000, res = 250)
g <- ggplot() +
  geom_sf(data = states, fill = NA, color = "darkgray", linetype = "solid", size = 0.5) +
  geom_sf(data = buffers_merged, 
          aes(fill = CAR_u1, shape = factor(Z)), size = 4, stroke = 0) +
  # drop the shape legend:
  scale_shape_manual(values = c("0" = 21, "1" = 24),
                     guide = "none") +
  # relabel the fill legend:
  scale_fill_viridis_c(name = "U") +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    panel.grid    = element_blank(),
    axis.text     = element_blank(),
    axis.ticks    = element_blank(),
    axis.title    = element_blank(),
    plot.title    = element_text(hjust = 0.5),
    legend.position = "top",
    legend.key.width = unit(1.5, "cm")
  )
plot_with_insets(g)
dev.off()

png('images/car2.png', width = 1300, height = 1000, res = 250)
g <- ggplot() +
  geom_sf(data = states, fill = NA, color = "darkgray", linetype = "solid", size = 0.5) +
  geom_sf(data = buffers_merged, 
          aes(fill = CAR_u2, shape = factor(Z)), size = 4, stroke = 0) +
  # drop the shape legend:
  scale_shape_manual(values = c("0" = 21, "1" = 24),
                     guide = "none") +
  # relabel the fill legend:
  scale_fill_viridis_c(name = "U") +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    panel.grid    = element_blank(),
    axis.text     = element_blank(),
    axis.ticks    = element_blank(),
    axis.title    = element_blank(),
    plot.title    = element_text(hjust = 0.5),
    legend.position = "top",
    legend.key.width = unit(1.5, "cm")
  )
plot_with_insets(g)
dev.off()

png('images/car3.png', width = 1300, height = 1000, res = 250)
g <- ggplot() +
  geom_sf(data = states, fill = NA, color = "darkgray", linetype = "solid", size = 0.5) +
  geom_sf(data = buffers_merged, 
          aes(fill = CAR_u3, shape = factor(Z)), size = 4, stroke = 0) +
  # drop the shape legend:
  scale_shape_manual(values = c("0" = 21, "1" = 24),
                     guide = "none") +
  # relabel the fill legend:
  scale_fill_viridis_c(name = "U") +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    panel.grid    = element_blank(),
    axis.text     = element_blank(),
    axis.ticks    = element_blank(),
    axis.title    = element_blank(),
    plot.title    = element_text(hjust = 0.5),
    legend.position = "top",
    legend.key.width = unit(1.5, "cm")
  )
plot_with_insets(g)
dev.off()

