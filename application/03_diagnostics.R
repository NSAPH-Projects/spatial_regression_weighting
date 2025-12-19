library(tidyr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(gridExtra)
library(Matrix)
library(xtable)
library(ggnewscale)
library(cowplot)
library(MASS)

source('../funcs.R')

states <- st_read('../data/tl_2010_us_state10/tl_2010_us_state10.shp')

# Read in preprocessed data
load('../data/preprocessed_superfunds.RData')
n <- nrow(buffers)
us_states <- ne_states(country = "United States of America", 
                       returnclass = "sf")
us_outline <- ne_countries(scale = "medium", country = "United States of America", returnclass = "sf")

################# PLOT BINARY Exposure ######################
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
    title = "Superfund Sites that were cleaned up and removed from the National Priority List between 1991 and 2015",
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
    name = "Binary treatment (Remediation)",
    labels = c("Control", "Treated")
  ) +
  scale_shape_manual(
    values = c("0" = 16, "1" = 17), # Circle for Control, triangle for Treated
    name = "Binary treatment (Remediation)",
    labels = c("Control", "Treated")
  )
png('images/treatment.png', width = 1500, height = 900, res = 160)
plot_with_insets(g) 
dev.off()

g <- ggplot() +
  # Add map outline
  geom_sf(data = us_outline, fill = NA, color = "black", linetype = "solid") +
  # Plot centroids with both color and shape mapped to factor(Z)
  geom_sf(data = buffers, aes(color = factor(Z), fill = factor(Z)), size = 0.5) +
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
    name = "Binary treatment (Cleanup)",
    labels = c("Control", "Treated")
  ) +
  scale_fill_manual(
    values = c("0" = "orange", "1" = "dodgerblue"),
    name = "Binary treatment (Cleanup)",
    labels = c("Control", "Treated")
  )
png('images/treatment_with_buffers.png', width = 4000, height = 2400, res = 400)
plot_with_insets(g) 
dev.off()

################# CALCULATE IMPLIED WEIGHTS ##################################
# Scale X
X[,2:ncol(X)] <- scale(X[,2:ncol(X)])
X <- as.matrix(X)

# GP model
set.seed(100)
kappa <- 10 
rangec <- 500
phic <- rangec/(2*sqrt(kappa))
S <- geoR::matern(u =  dmat, phi = phic, kappa = kappa)
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
    size = 1.25, 
    alpha = 0.7,
    stroke   = 0.25, 
    color    = "black",
    position = position_nudge(x = -0.3, y = 0)  # ← correct placement
  ) +
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "orange", midpoint = 0, 
                       breaks = c(-0.003, 0.005),
                       limits = c(-0.0031, 0.0051),
                       labels = function(x) sprintf("%.3f", x),
                       name = "Implied weights (Control)",
                       guide = guide_colorbar(order = 2)) +  
  
  new_scale_fill() +  
  
  # Treated group (Z == 1): Green triangle with black outline
  geom_sf(data = buffers_merged_geo[buffers_merged_geo$Z == 1,], 
          aes(fill = CARiw, shape = factor(Z)), 
          size = 1.25, alpha = 0.7, stroke = 0.25, color = "black") +  
  
  scale_fill_gradient2(low = "orange", mid = "white", high = "dodgerblue", midpoint = 0, 
                       breaks = c(0, 0.008),
                       limits = c(-0.0001, 0.0081),
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
    legend.position = "left",
    legend.key.height = unit(0.2, "cm")
  )

inset_CAR <- g + 
  coord_sf(xlim = c(-86,-80), ylim = c(25.25,30.5)) +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    plot.background = element_rect(color = "black", linewidth = 1)
  )

gCAR_main <- plot_with_insets(g,xshift = 0.25)

# Now use ggdraw with extended xlim 
gCAR <- ggdraw(xlim = c(0, 1.1), ylim = c(0, 1)) +
  draw_plot(gCAR_main, x = 0,   y = 0, width = 1.1,   height = 1.1) +
  draw_plot(inset_CAR,  x = 0.825, y = 0.12, width = 0.35, height = 0.2, scale = 1.9)

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
    size = 1.25, 
    alpha = 0.7,
    stroke   = 0.25, 
    color    = "black",
    position = position_nudge(x = -0.3, y = 0)  # ← correct placement
  ) +
  
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "orange", midpoint = 0, 
                       breaks = c(-0.003, 0.005),
                       limits = c(-0.0031, 0.0051),
                       labels = function(x) sprintf("%.3f", x),
                       name = "Implied weights (Control)",
                       guide = guide_colorbar(order = 2)) +  
  
  new_scale_fill() +  
  
  # Treated group (Z == 1): Green triangle with black outline
  geom_sf(data = buffers_merged_geo[buffers_merged_geo$Z == 1,], 
          aes(fill = REiw, shape = factor(Z)),
          size = 1.25, alpha = 0.7, stroke = 0.25, color = "black") +  
  
  scale_fill_gradient2(low = "orange", mid = "white", high = "dodgerblue", midpoint = 0, 
                       breaks = c(0, 0.008),
                       limits = c(-0.0001, 0.0081),
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
    legend.position = "left",
    legend.key.height = unit(0.2, "cm")
  )

inset_RE <- g + 
  coord_sf(xlim = c(-86,-80), ylim = c(25.25,30.5)) +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    plot.background = element_rect(color = "black", linewidth = 1)
  )

gRE_main <- plot_with_insets(g,xshift = 0.25)

# Now use ggdraw with extended xlim 
gRE <- ggdraw(xlim = c(0, 1.1), ylim = c(0, 1)) +
  draw_plot(gRE_main, x = 0,   y = 0, width = 1.1,   height = 1.1) +
  draw_plot(inset_RE,  x = 0.825, y = 0.12, width = 0.35, height = 0.2, scale = 1.9)

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
    size = 1.25, 
    alpha = 0.7,
    stroke   = 0.25, 
    color    = "black",
    position = position_nudge(x = -0.3, y = 0)
  ) +

  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "orange", midpoint = 0, 
                       breaks = c(-0.003, 0.005),
                       limits = c(-0.0031, 0.0051),
                       labels = function(x) sprintf("%.3f", x),
                       name = "Implied weights (Control)",
                       guide = guide_colorbar(order = 2)) +  
  
  new_scale_fill() +  
  
  # Treated group (Z == 1): Green triangle with black outline
  geom_sf(data = buffers_merged_geo[buffers_merged_geo$Z == 1,],
          aes(fill = GPiw, shape = factor(Z)),
          size = 1.25, alpha = 0.7, stroke = 0.25, color = "black") +
  
  scale_fill_gradient2(low = "orange", mid = "white", high = "dodgerblue", midpoint = 0,
                       breaks = c(0, 0.008),
                       limits = c(-0.0001, 0.0081),
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
    legend.position = "left",
    legend.key.height = unit(0.2, "cm")
  )
inset_GP <- g + 
  coord_sf(xlim = c(-86,-80), ylim = c(25.25,30.5)) +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    plot.background = element_rect(color = "black", linewidth = 1)
  )

gGP_main <- plot_with_insets(g, xshift = 0.25)

# Now use ggdraw with extended xlim 
gGP <- ggdraw(xlim = c(0, 1.1), ylim = c(0, 1)) +
  draw_plot(gGP_main, x = 0,   y = 0, width = 1.1,   height = 1.1) +
  draw_plot(inset_GP,  x = 0.825, y = 0.12, width = 0.35, height = 0.2, scale = 1.9)

#png('images/impliedweights_us.png', width = 1500, height = 2000, res = 230)
png('images/impliedweights_us_dec18.png', width = 2300, height = 3150, res = 300) 
grid.arrange(gRE, gCAR, gGP, ncol = 1)
dev.off()


######################################## UNMEASURED SPATIAL CONFOUNDER EXAMPLES #########################################
buffers_merged <- st_centroid(buffers)

# GP 
set.seed(100)
kappa <- 10 
rangec <- 500
phic <- rangec/(2*sqrt(kappa))
S <- geoR::matern(u =  dmat, phi = phic, kappa = kappa)
E <- eigen(S)
V <- E$vectors
moran_adhoc(V[,1:n] %*% rnorm(ncol(V), mean = 0, sd = 1/(1:ncol(V))^5), Wmat = S, coeff = 1/E$values[1])
moran_adhoc(V[,1:n] %*% rnorm(ncol(V), mean = 0, sd = 1/(1:ncol(V))), Wmat = S, coeff = 1/E$values[1])
moran_adhoc(rnorm(ncol(V), mean = 0, sd =1), Wmat = S, coeff = 1/E$values[1])
Sigma <- diag(n) + 10*S
Sigmainv <- solve(Sigma)
GPiw <- impliedweightsgeneral(X = X,
                              Z = Z,
                              Sigmainv = Sigmainv)
sum(GPiw[Z == 1]);sum(GPiw[Z == 0])

# Squared Eigenvector imbalance
#Es <- eigen(S)
m <- 0
while (m < 0.7 | m > 0.95){
  buffers_merged$GP_u1 <- V[,1:n] %*% rnorm(ncol(V), mean = 0, sd = 1/(1:ncol(V))^5) #runif(ncol(V), min = 0, max = 1)
  buffers_merged$GP_u1 <- buffers_merged$GP_u1 - mean(buffers_merged$GP_u1)
  buffers_merged$GP_u1 <- buffers_merged$GP_u1/norm(buffers_merged$GP_u1, type = '2')
  m <- moran_adhoc(buffers_merged$GP_u1, Wmat = S, coeff = 1/E$values[1])
}
print(paste0('Moran I is ', m, ' and imbalance is ', sum(GPiw*(2*Z-1)*buffers_merged$GP_u1)))

m = 0
while (m < 0.45 | m > 0.55){
  buffers_merged$GP_u2 <- V[,1:n] %*% rnorm(ncol(V), mean = 0, sd = 1/1:ncol(V)) #runif(ncol(V), min = 0, max = 1)
  buffers_merged$GP_u2 <- buffers_merged$GP_u2 - mean(buffers_merged$GP_u2)
  buffers_merged$GP_u2 <- buffers_merged$GP_u2/norm(buffers_merged$GP_u2, type = '2')
  m <- moran_adhoc(buffers_merged$GP_u2, Wmat = S, coeff = 1/E$values[1])
}
print(paste0('Moran I is ', m, ' and imbalance is ', sum(GPiw*(2*Z-1)*buffers_merged$GP_u2)))

m = 1
while (m > 0.1){
  buffers_merged$GP_u3 <- V[,1:n] %*% rnorm(ncol(V), mean = 0, sd = sqrt(1/(1:ncol(V))))
  buffers_merged$GP_u3 <- buffers_merged$GP_u3 - mean(buffers_merged$GP_u3)
  buffers_merged$GP_u3 <- buffers_merged$GP_u3/norm(buffers_merged$GP_u3, type = '2')
  m <- moran_adhoc(buffers_merged$GP_u3, Wmat = S, coeff = 1/E$values[1])
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
S <- bdiag(Jks)
E <- eigen(S)
V <- E$vectors
Sigma <- diag(n) + 10*S
Sigmainv <- solve(Sigma)
REiw <- impliedweightsgeneral(X = X,
                              Z = Z,
                              Sigmainv = Sigmainv)

# Squared Eigenvector imbalance
m = 0
while (m < 0.9 | m > 0.95) {
  buffers_merged$RE_u1 <- V[,1:n] %*% c(rnorm(51, mean = 0, sd = 1/(1:ncol(V))), rep(0, ncol(V)-51)) #runif(ncol(V), min = 0, max = 1)
  buffers_merged$RE_u1 <- buffers_merged$RE_u1 - mean(buffers_merged$RE_u1)
  buffers_merged$RE_u1 <- buffers_merged$RE_u1/norm(buffers_merged$RE_u1, type = '2')
  m <- moran_adhoc(buffers_merged$RE_u1, Wmat = S, coeff = 1/E$values[1])
}
print(paste0('Moran I is ', m, ' and imbalance is ', sum(REiw*(2*Z-1)*buffers_merged$RE_u1)))

m = 0
while (m < 0.45 | m > 0.55){
  buffers_merged$RE_u2 <- V[,1:n] %*% c(rnorm(51, mean = 0, sd = 1/(1:ncol(V))), rep(0.01, ncol(V)-51)) #runif(ncol(V), min = 0, max = 1)
  buffers_merged$RE_u2 <- buffers_merged$RE_u2 - mean(buffers_merged$RE_u2)
  buffers_merged$RE_u2 <- buffers_merged$RE_u2/norm(buffers_merged$RE_u2, type = '2')
  m <- moran_adhoc(buffers_merged$RE_u2, Wmat = S, coeff = 1/E$values[1])
}
print(paste0('Moran I is ', m, ' and imbalance is ', sum(REiw*(2*Z-1)*buffers_merged$RE_u2)))

m = 1
while (m > 0.1){
  buffers_merged$RE_u3 <- V[,1:n] %*% c(rnorm(51, mean = 0, sd = 1/(1:ncol(V))), rep(0.1, ncol(V)-51)) #runif(ncol(V), min = 0, max = 1)
  buffers_merged$RE_u3 <- buffers_merged$RE_u3 - mean(buffers_merged$RE_u3)
  buffers_merged$RE_u3 <- buffers_merged$RE_u3/norm(buffers_merged$RE_u3, type = '2')
  m <- moran_adhoc(buffers_merged$RE_u3, Wmat = S, coeff = 1/E$values[1])
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
E <- eigen(S)
V <- Re(E$vectors)
Sigma <- sigma2*diag(n) + rho2*S
Sigmainv <- solve(Sigma)
CARiw <- impliedweightsgeneral(X = X,
                               Z = Z,
                               Sigmainv = Sigmainv)
sum(CARiw[Z == 1]);sum(CARiw[Z == 0])

# Squared Eigenvector imbalance
m = 0
while (m < 0.9 | m > 0.95){
  buffers_merged$CAR_u1 <- V[,1:n] %*% rnorm(ncol(V), mean = 0, sd = 1/(1:ncol(V))^2) #runif(ncol(V), min = 0, max = 1)
  buffers_merged$CAR_u1 <- buffers_merged$CAR_u1 - mean(buffers_merged$CAR_u1)
  buffers_merged$CAR_u1 <- buffers_merged$CAR_u1/norm(buffers_merged$CAR_u1, type = '2')
  m <- moran_adhoc(buffers_merged$CAR_u1, Wmat = S, coeff = 1/Re(E$values[1]))
}
print(paste0('Moran I is ', m, ' and imbalance is ', sum(CARiw*(2*Z-1)*buffers_merged$CAR_u1)))

m = 0
while (m < 0.45 | m > 0.55){
  buffers_merged$CAR_u2 <- V[,1:n] %*% rnorm(ncol(V), mean = 0, sd = 1/(1:ncol(V))) #runif(ncol(V), min = 0, max = 1)
  buffers_merged$CAR_u2 <- buffers_merged$CAR_u2 - mean(buffers_merged$CAR_u2)
  buffers_merged$CAR_u2 <- buffers_merged$CAR_u2/norm(buffers_merged$CAR_u2, type = '2')
  m <- moran_adhoc(buffers_merged$CAR_u2, Wmat = S, coeff = 1/Re(E$values[1]))
}
print(paste0('Moran I is ', m, ' and imbalance is ', sum(CARiw*(2*Z-1)*buffers_merged$CAR_u2)))

m = 1
while (m > 0.1){
  buffers_merged$CAR_u3 <- V[,1:n] %*% rnorm(ncol(V), mean = 0, sd = sqrt(1/(1:ncol(V)))) #runif(ncol(V), min = 0, max = 1)
  buffers_merged$CAR_u3 <- buffers_merged$CAR_u3 - mean(buffers_merged$CAR_u3)
  buffers_merged$CAR_u3 <- buffers_merged$CAR_u3/norm(buffers_merged$CAR_u3, type = '2')
  m <- moran_adhoc(buffers_merged$CAR_u3, Wmat = S, coeff = 1/Re(E$values[1]))
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

################ PLOT OF PROPERTIES OF THE WEIGHTS ##############
# Calculate the proportion of Control neighbors each Treated unit has and the number of Treated neighbors that each Control unit has
buffers_merged_geo <- st_transform(buffers_merged, crs = 4326)

prop_neighbors_opposite <- rep(NA, nrow(buffers_merged_geo))
for (i in 1:nrow(buffers_merged_geo)) {
  if (buffers_merged_geo$Z[i] == 1) {
    prop_neighbors_opposite[i] <- sum(adjacency_matrix[i,buffers_merged_geo$Z == 0])/sum(adjacency_matrix[i,])
  } else {
    prop_neighbors_opposite[i] <- sum(adjacency_matrix[i,buffers_merged_geo$Z == 1])/sum(adjacency_matrix[i,])
  }
}
g1 <- ggplot(buffers_merged_geo, aes(x = prop_neighbors_opposite, y = abs(CARiw))) +
  geom_point() +
  labs(title = "Conditional Autoregressive",
       x = "Proportion of neighbors with the opposite treatment status",
       y = "absolute implied weight") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) 

# Calculate the proportion of Control neighbors each Treated unit has and the number of Treated neighbors that each Control unit has
prop_distanceneighbors_opposite <- rep(NA, nrow(buffers_merged_geo))
for (i in 1:nrow(buffers_merged_geo)) {
  if (buffers_merged_geo$Z[i] == 1) {
    prop_distanceneighbors_opposite[i] <- sum(dmat[i,buffers_merged_geo$Z == 0] < 200)/sum(dmat[i,] < 200) #/1000
  } else {
    prop_distanceneighbors_opposite[i] <- sum(dmat[i,buffers_merged_geo$Z == 1] < 200)/sum(dmat[i,] < 200) #/1000
  }
}

g2 <- ggplot(buffers_merged_geo, 
             aes(x = prop_distanceneighbors_opposite, y = abs(GPiw))) + # 
  geom_point() +
  labs(title = "Gaussian Process",
       x = "Proportion of sites within 200km with the opposite treatment status",
       y = "absolute implied weight") +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5))

# Calculate the proportion of Control units in a cluster for each Treated unit and the proportion of Treated units in a cluster for each Control unit
prop_cluster_opposite <- rep(NA, nrow(buffers_merged_geo))
for (i in 1:nrow(buffers_merged_geo)) {
  if (buffers_merged_geo$Z[i] == 1) {
    prop_cluster_opposite[i] <- length(which(clusters == clusters[i] & buffers_merged_geo$Z == 0))/length(which(clusters == clusters[i]))
  } else {
    prop_cluster_opposite[i] <- length(which(clusters == clusters[i] & buffers_merged_geo$Z == 1))/length(which(clusters == clusters[i]))
  }
}
g3 <- ggplot(buffers_merged_geo, aes(x = prop_cluster_opposite, y = abs(REiw))) +
  geom_point() +
  labs(title = "Random Effects",
       x = "Proportion of units in the same state with the opposite treatment status",
       y = "absolute implied weight") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

png('images/weights_properties.png', width = 700, height = 1500, res = 130)
grid.arrange(g3,g1,g2, ncol = 1)
dev.off()

################## Data characteristics table ##################
# Create a summary table of data characteristics, mean and sd of each column in X
data_characteristics <- data.frame(
  Variable = colnames(X),
  Mean = sapply(1:ncol(X), function(i) mean(X[,i])),
  SD = sapply(1:ncol(X), function(i) sd(X[,i]))
)
# Print the table with xtable for latex without row numbers
print(xtable(data_characteristics, digits = 3), 
      type = "latex", include.rownames = FALSE)
