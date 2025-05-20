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
load('data/preprocessed_superfunds.RData')
n <- nrow(buffers)
us_states <- ne_states(country = "United States of America", returnclass = "sf")

################# CALCULATE IMPLIED WEIGHTS ##################################
# Scale X
X <- st_drop_geometry(X)
X[,2:11] <- scale(X[,2:11])
X <- X[,-12]
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

################## PLOT WEIGHTS AND BASE WEIGHTS SIDE BY SIDE ##################
# allmint <- min(0,min(c(CARiw[Z == 1], 
#                        REiw[Z == 1],
#                        GPiw[Z == 1])))
# allmaxt <- max(c(CARiw[Z == 1], 
#                  REiw[Z == 1],
#                  GPiw[Z == 1]))
# allminc <- min(0,min(c(CARiw[Z == 0], 
#                        REiw[Z == 0],
#                        GPiw[Z == 0])))
# allmaxc <- max(c(CARiw[Z == 0], 
#                  REiw[Z == 0],
#                  GPiw[Z == 0]))

# Merge the Implied weights with the data
buffer_centroids <- st_centroid(buffers)
buffers_merged <- cbind(buffer_centroids, REiw, CARiw, GPiw)

# Plot baseCAR
buffers_merged_geo <- st_transform(buffers_merged, crs = 4326)

buffers_merged_geo <- buffers_merged_geo[order(buffers_merged_geo$CARiw),]
gCAR <- ggplot() +
  scale_shape_manual(values = c("0" = 21, "1" = 24),
                     guide = guide_legend(order = 1)) + 
  
  # Base map with county outlines
  #geom_sf(data = counties, fill = NA, color = "darkgray", linetype = "solid", size = 0.5) +
  geom_sf(data = states, fill = NA, color = "darkgray", linetype = "solid", size = 0.5) +
  # Control group (Z == 0): Red circle with black outline
  # Custom shape: 21 (circle), 24 (triangle)
  geom_sf(data = buffers_merged_geo[buffers_merged_geo$Z == 0,], 
          aes(fill = CARiw, shape = factor(Z)), 
          size = 2, stroke = 0.25, color = "black") +  # Black outline
  
  scale_fill_gradient2(low = "green", mid = "white", high = "red", midpoint = 0, 
                       breaks = c(min(buffers_merged_geo[buffers_merged_geo$Z == 0,]$CARiw, na.rm = TRUE), 
                                  max(buffers_merged_geo[buffers_merged_geo$Z == 0,]$CARiw, na.rm = TRUE)),
                       labels = function(x) sprintf("%.3f", x),
                       #limits = c(allminc, allmaxc), 
                       name = "Implied weights (Control)",
                       guide = guide_colorbar(order = 2)) +  
  
  new_scale_fill() +  
  
  # Treated group (Z == 1): Green triangle with black outline
  geom_sf(data = buffers_merged_geo[buffers_merged_geo$Z == 1,], 
          aes(fill = CARiw, shape = factor(Z)), 
          size = 2, stroke = 0.25, color = "black") +  
  
  scale_fill_gradient2(low = "red", mid = "white", high = "green", midpoint = 0, 
                       breaks = c(min(buffers_merged_geo[buffers_merged_geo$Z == 1,]$CARiw, na.rm = TRUE), 
                                  max(buffers_merged_geo[buffers_merged_geo$Z == 1,]$CARiw, na.rm = TRUE)),
                       labels = function(x) sprintf("%.3f", x),
                       #limits = c(allmint, allmaxt), 
                       name = "Implied weights (Treated)",
                       guide = guide_colorbar(order = 3)) +  
  coord_sf(xlim = c(-125, -65), ylim = c(25, 50)) +
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

# Plot  FE
buffers_merged_geo <- buffers_merged_geo[order(buffers_merged_geo$REiw),]
gRE <- ggplot() +
  
  # Custom shape: 21 (circle), 24 (triangle)
  scale_shape_manual(values = c("0" = 21, "1" = 24),
                     guide = guide_legend(order = 1)) + 
  
  # Base map with county outlines
  geom_sf(data = states, fill = NA, color = "darkgray", linetype = "solid", size = 0.5) +
  
  # Control group (Z == 0): Red circle with black outline
  geom_sf(data = buffers_merged_geo[buffers_merged_geo$Z == 0,], 
          aes(fill = REiw, shape = factor(Z)), 
          size = 2, stroke = 0.25, color = "black") +  # Black outline

  
  scale_fill_gradient2(low = "green", mid = "white", high = "red", midpoint = 0, 
                       breaks = c(min(buffers_merged_geo[buffers_merged_geo$Z == 0,]$REiw, na.rm = TRUE), 
                                  max(buffers_merged_geo[buffers_merged_geo$Z == 0,]$REiw, na.rm = TRUE)),
                       labels = function(x) sprintf("%.3f", x),
                       #limits = c(allminc, allmaxc), 
                       name = "Implied weights (Control)",
                       guide = guide_colorbar(order = 2)) +  
  
  new_scale_fill() +  
  
  # Treated group (Z == 1): Green triangle with black outline
  geom_sf(data = buffers_merged_geo[buffers_merged_geo$Z == 1,], 
          aes(fill = REiw, shape = factor(Z)), 
          size = 2, stroke = 0.25, color = "black") +  
  
  scale_fill_gradient2(low = "red", mid = "white", high = "green", midpoint = 0, 
                       breaks = c(min(buffers_merged_geo[buffers_merged_geo$Z == 1,]$REiw, na.rm = TRUE), 
                                  max(buffers_merged_geo[buffers_merged_geo$Z == 1,]$REiw, na.rm = TRUE)),
                       labels = function(x) sprintf("%.3f", x),
                       #limits = c(allmint, allmaxt), 
                       name = "Implied weights (Treated)",
                       guide = guide_colorbar(order = 3)) +  
  coord_sf(xlim = c(-125, -65), ylim = c(25, 50), expand = FALSE) +
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

# Plot  GP
buffers_merged_geo <- buffers_merged_geo[order(buffers_merged_geo$GPiw),]
gGP <- ggplot() +
  # Custom shape: 21 (circle), 24 (triangle)
  scale_shape_manual(values = c("0" = 21, "1" = 24),
                     guide = guide_legend(order = 1)) + 
  
  # Base map with county outlines
  geom_sf(data = states, fill = NA, color = "darkgray", linetype = "solid", size = 0.5) +
  
  # Control group (Z == 0): Red circle with black outline
  geom_sf(data = buffers_merged_geo[buffers_merged_geo$Z == 0,], 
          aes(fill = GPiw, shape = factor(Z)), 
          size = 2, stroke = 0.25, color = "black") +  # Black outline
  

  scale_fill_gradient2(low = "green", mid = "white", high = "red", midpoint = 0, 
                       breaks = c(min(buffers_merged_geo[buffers_merged_geo$Z == 0,]$GPiw, na.rm = TRUE), 
                                  max(buffers_merged_geo[buffers_merged_geo$Z == 0,]$GPiw, na.rm = TRUE)),
                       labels = function(x) sprintf("%.3f", x),
                       #limits = c(allminc, allmaxc), 
                       name = "Implied weights (Control)",
                       guide = guide_colorbar(order = 2)) +  
  
  new_scale_fill() +  
  
  # Treated group (Z == 1): Green triangle with black outline
  geom_sf(data = buffers_merged_geo[buffers_merged_geo$Z == 1,],
          aes(fill = GPiw, shape = factor(Z)),
          size = 2, stroke = 0.25, color = "black") +
  
  scale_fill_gradient2(low = "red", mid = "white", high = "green", midpoint = 0,
                       breaks = c(min(buffers_merged_geo[buffers_merged_geo$Z == 1,]$GPiw, na.rm = TRUE), 
                                  max(buffers_merged_geo[buffers_merged_geo$Z == 1,]$GPiw, na.rm = TRUE)),
                       labels = function(x) sprintf("%.3f", x),
                       #limits = c(allmint, allmaxt), 
                       name = "Implied weights (Treated)",
                       guide = guide_colorbar(order = 3)) +
  coord_sf(xlim = c(-125, -65), ylim = c(25, 50), expand = FALSE) +
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

# gCAR <- gCAR + theme(legend.position = "bottom")
# gRE <- gRE + theme(legend.position = "bottom")
# gGP <- gGP + theme(legend.position = "bottom")
# 
# # Combine the plots with patchwork
# combined_plot <- (gRE + gCAR + gGP) +
#   plot_layout(ncol = 1, guides = "collect") &
#   theme(legend.position = "right",
#         legend.direction = "vertical")

png('images/impliedweights_us.png', width = 1500, height = 2000, res = 200)
#print(combined_plot)
grid.arrange(gRE, gCAR, gGP, ncol = 1)
dev.off()


######################################## UNMEASURED SPATIAL CONFOUNDER EXAMPLES #########################################
buffers_merged <- st_centroid(buffers)

set.seed(111)
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
buffers_merged$GP_u1 <- V[,1:n] %*% rnorm(ncol(V), mean = 0, sd = 1/(1:ncol(V))^4) #runif(ncol(V), min = 0, max = 1)
buffers_merged$GP_u1 <- buffers_merged$GP_u1 - mean(buffers_merged$GP_u1)
buffers_merged$GP_u1 <- buffers_merged$GP_u1/norm(buffers_merged$GP_u1, type = '2')
sum(GPiw*(2*Z-1)*buffers_merged$GP_u1)
moran_adhoc(buffers_merged$GP_u1, Wmat = S, coeff = 1/eigen(S)$values[2])

buffers_merged$GP_u2 <- V[,1:n] %*% rnorm(ncol(V), mean = 0, sd = 1/(1:ncol(V))^1) #runif(ncol(V), min = 0, max = 1)
buffers_merged$GP_u2 <- buffers_merged$GP_u2 - mean(buffers_merged$GP_u2)
buffers_merged$GP_u2 <- buffers_merged$GP_u2/norm(buffers_merged$GP_u2, type = '2')
sum(GPiw*(2*Z-1)*buffers_merged$GP_u2)
moran_adhoc(buffers_merged$GP_u2, Wmat = S, coeff = 1/eigen(S)$values[2])

buffers_merged$GP_u3 <- V[,1:n] %*% rnorm(ncol(V), mean = 0, sd = sqrt(1/(1:ncol(V)))) #runif(ncol(V), min = 0, max = 1)
buffers_merged$GP_u3 <- buffers_merged$GP_u3 - mean(buffers_merged$GP_u3)
buffers_merged$GP_u3 <- buffers_merged$GP_u3/norm(buffers_merged$GP_u3, type = '2')
sum(GPiw*(2*Z-1)*buffers_merged$GP_u3)
moran_adhoc(buffers_merged$GP_u3, Wmat = S, coeff = 1/eigen(S)$values[2])


# PLOT each on the map
png('images/gp1.png', width = 1300, height = 1000, res = 250)
ggplot() +
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
  coord_sf(xlim = c(-125, -65), ylim = c(25, 50)) +
  theme(
    panel.grid    = element_blank(),
    axis.text     = element_blank(),
    axis.ticks    = element_blank(),
    axis.title    = element_blank(),
    plot.title    = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.key.width = unit(1.5, "cm")
  )
dev.off()

png('images/gp2.png', width = 1300, height = 1000, res = 250)
ggplot() +
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
  coord_sf(xlim = c(-125, -65), ylim = c(25, 50)) +
  theme(
    panel.grid    = element_blank(),
    axis.text     = element_blank(),
    axis.ticks    = element_blank(),
    axis.title    = element_blank(),
    plot.title    = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.key.width = unit(1.5, "cm")
  )
dev.off()

png('images/gp3.png', width = 1300, height = 1000, res = 250)
ggplot() +
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
  coord_sf(xlim = c(-125, -65), ylim = c(25, 50)) +
  theme(
    panel.grid    = element_blank(),
    axis.text     = element_blank(),
    axis.ticks    = element_blank(),
    axis.title    = element_blank(),
    plot.title    = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.key.width = unit(1.5, "cm")
  )
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
buffers_merged$RE_u1 <- V[,1:n] %*% c(rnorm(51, mean = 0, sd = 1/(1:ncol(V))), rep(0, ncol(V)-51)) #runif(ncol(V), min = 0, max = 1)
buffers_merged$RE_u1 <- buffers_merged$RE_u1 - mean(buffers_merged$RE_u1)
buffers_merged$RE_u1 <- buffers_merged$RE_u1/norm(buffers_merged$RE_u1, type = '2')
sum(REiw*(2*Z-1)*buffers_merged$RE_u1)
moran_adhoc(buffers_merged$RE_u1, Wmat = S, coeff = 1/eigen(S)$values[1])

buffers_merged$RE_u2 <- V[,1:n] %*% c(rnorm(51, mean = 0, sd = 1/(1:ncol(V))), rep(0.01, ncol(V)-51)) #runif(ncol(V), min = 0, max = 1)
buffers_merged$RE_u2 <- buffers_merged$RE_u2 - mean(buffers_merged$RE_u2)
buffers_merged$RE_u2 <- buffers_merged$RE_u2/norm(buffers_merged$RE_u2, type = '2')
sum(REiw*(2*Z-1)*buffers_merged$RE_u2)
moran_adhoc(buffers_merged$RE_u2, Wmat = S, coeff = 1/eigen(S)$values[1])

buffers_merged$RE_u3 <- V[,1:n] %*% c(rnorm(51, mean = 0, sd = 1/(1:ncol(V))), rep(0.1, ncol(V)-51)) #runif(ncol(V), min = 0, max = 1)
buffers_merged$RE_u3 <- buffers_merged$RE_u3 - mean(buffers_merged$RE_u3)
buffers_merged$RE_u3 <- buffers_merged$RE_u3/norm(buffers_merged$RE_u3, type = '2')
sum(REiw*(2*Z-1)*buffers_merged$RE_u3)
moran_adhoc(buffers_merged$RE_u3, Wmat = S, coeff = 1/eigen(S)$values[1])

# PLOT each on the map
png('images/re1.png', width = 1300, height = 1000, res = 250)
ggplot() +
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
  coord_sf(xlim = c(-125, -65), ylim = c(25, 50)) +
  theme(
    panel.grid    = element_blank(),
    axis.text     = element_blank(),
    axis.ticks    = element_blank(),
    axis.title    = element_blank(),
    plot.title    = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.key.width = unit(1.5, "cm")
  )
dev.off()

png('images/re2.png', width = 1300, height = 1000, res = 250)
ggplot() +
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
  coord_sf(xlim = c(-125, -65), ylim = c(25, 50)) +
  theme(
    panel.grid    = element_blank(),
    axis.text     = element_blank(),
    axis.ticks    = element_blank(),
    axis.title    = element_blank(),
    plot.title    = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.key.width = unit(1.5, "cm")
  )
dev.off()

png('images/re3.png', width = 1300, height = 1000, res = 250)
ggplot() +
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
  coord_sf(xlim = c(-125, -65), ylim = c(25, 50)) +
  theme(
    panel.grid    = element_blank(),
    axis.text     = element_blank(),
    axis.ticks    = element_blank(),
    axis.title    = element_blank(),
    plot.title    = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.key.width = unit(1.5, "cm")
  )
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
buffers_merged$CAR_u1 <- V[,1:n] %*% rnorm(ncol(V), mean = 0, sd = 1/(1:ncol(V))^2) #runif(ncol(V), min = 0, max = 1)
buffers_merged$CAR_u1 <- buffers_merged$CAR_u1 - mean(buffers_merged$CAR_u1)
buffers_merged$CAR_u1 <- buffers_merged$CAR_u1/norm(buffers_merged$CAR_u1, type = '2')
sum(CARiw*(2*Z-1)*buffers_merged$CAR_u1)
moran_adhoc(buffers_merged$CAR_u1, Wmat = S, coeff = 1/eigen(S)$values[1])

buffers_merged$CAR_u2 <- V[,1:n] %*% rnorm(ncol(V), mean = 0, sd = 1/(1:ncol(V)))  #runif(ncol(V), min = 0, max = 1)
buffers_merged$CAR_u2 <- buffers_merged$CAR_u2 - mean(buffers_merged$CAR_u2)
buffers_merged$CAR_u2 <- buffers_merged$CAR_u2/norm(buffers_merged$CAR_u2, type = '2')
sum(CARiw*(2*Z-1)*buffers_merged$CAR_u2)
moran_adhoc(buffers_merged$CAR_u2, Wmat = S, coeff = 1/eigen(S)$values[1])

buffers_merged$CAR_u3 <- V[,1:n] %*% rnorm(ncol(V), mean = 0, sd = sqrt(1/(1:ncol(V)))) #runif(ncol(V), min = 0, max = 1)
buffers_merged$CAR_u3 <- buffers_merged$CAR_u3 - mean(buffers_merged$CAR_u3)
buffers_merged$CAR_u3 <- buffers_merged$CAR_u3/norm(buffers_merged$CAR_u3, type = '2')
sum(CARiw*(2*Z-1)*buffers_merged$CAR_u3)
moran_adhoc(buffers_merged$CAR_u3, Wmat = S, coeff = 1/eigen(S)$values[1])

png('images/car1.png', width = 1300, height = 1000, res = 250)
ggplot() +
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
  coord_sf(xlim = c(-125, -65), ylim = c(25, 50)) +
  theme(
    panel.grid    = element_blank(),
    axis.text     = element_blank(),
    axis.ticks    = element_blank(),
    axis.title    = element_blank(),
    plot.title    = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.key.width = unit(1.5, "cm")
  )
dev.off()

png('images/car2.png', width = 1300, height = 1000, res = 250)
ggplot() +
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
  coord_sf(xlim = c(-125, -65), ylim = c(25, 50)) +
  theme(
    panel.grid    = element_blank(),
    axis.text     = element_blank(),
    axis.ticks    = element_blank(),
    axis.title    = element_blank(),
    plot.title    = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.key.width = unit(1.5, "cm")
  )
dev.off()

png('images/car3.png', width = 1300, height = 1000, res = 250)
ggplot() +
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
  coord_sf(xlim = c(-125, -65), ylim = c(25, 50)) +
  theme(
    panel.grid    = element_blank(),
    axis.text     = element_blank(),
    axis.ticks    = element_blank(),
    axis.title    = element_blank(),
    plot.title    = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.key.width = unit(1.5, "cm")
  )
dev.off()

############################ Fit methods ###############################
mun <-  st_read('shapefiles/Municipality_Low_Birthweight_Rate_2016_2020.shp')

mun$mun_id <- 1:nrow(mun)
# For each buffer, assign to it a Y value that is average percent low bw
# Transform municipality geometry to the buffer's CRS
buffers <- st_transform(buffers, st_crs(mun))

# Compute intersections between buffers and municipalities
intersections <- st_intersection(buffers, mun)

# Compute area of each intersection
intersections <- intersections %>% 
  mutate(int_area = as.numeric(st_area(.)))

# Assume each municipality has a total area field (mun_area)
# If not, compute it:
mun <- mun %>% mutate(mun_area = as.numeric(st_area(.)))

# Join municipality total area to intersections (if not already present)
# Here we assume an identifier 'mun_id' exists in both data frames
intersections <- intersections %>% 
  left_join(mun %>% st_set_geometry(NULL) %>% select(mun_id, mun_area), by = "mun_id")

# Calculate the area fraction for each intersection
intersections <- intersections %>% 
  mutate(area_fraction = int_area / mun_area)

# Estimate births within each intersection
intersections <- intersections %>% 
  mutate(
    Nm_Brth_int = Nm_Brth * area_fraction,
    Lw_Brth_int = Lw_Brth * area_fraction
  )

# Aggregate by buffer (assuming buffer identifier S_EPA_I)
buffers_outcome <- intersections %>%
  group_by(S_EPA_I) %>%
  summarise(
    COUNTY          = first(COUNTY),
    total_Lw_Brth   = sum(Lw_Brth_int, na.rm = TRUE),
    total_Nm_Brth   = sum(Nm_Brth_int, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    Y = total_Lw_Brth / total_Nm_Brth
  )


summary(buffers_outcome$Y) # this looks super normal
hist(buffers_outcome$Y)
# 
# ############## VERY low birth weight #########
# mun_vlow <-  st_read('shapefiles/Municipality_VLow_Birthweight_Rate_2016_2020.shp')
# 
# mun_vlow$mun_id <- 1:nrow(mun_vlow)
# # For each buffer, assign to it a Y value that is average percent low bw
# # Transform municipality geometry to the buffer's CRS
# buffers <- st_transform(buffers, st_crs(mun_vlow))
# 
# # Compute intersections between buffers and municipalities
# intersections <- st_intersection(buffers, mun_vlow)
# 
# # Compute area of each intersection
# intersections <- intersections %>% 
#   mutate(int_area = as.numeric(st_area(.)))
# 
# # Assume each municipality has a total area field (mun_area)
# # If not, compute it:
# mun_vlow <- mun_vlow %>% mutate(mun_area = as.numeric(st_area(.)))
# 
# # Join municipality total area to intersections (if not already present)
# # Here we assume an identifier 'mun_id' exists in both data frames
# intersections <- intersections %>% 
#   left_join(mun_vlow %>% st_set_geometry(NULL) %>% select(mun_id, mun_area), by = "mun_id")
# 
# # Calculate the area fraction for each intersection
# intersections <- intersections %>% 
#   mutate(area_fraction = int_area / mun_area)
# 
# # Estimate births within each intersection
# intersections <- intersections %>% 
#   mutate(
#     Nm_Brth_int = Nm_Brth * area_fraction,
#     VLw_Brth_int = VLw_Brt * area_fraction
#   )
# 
# # Aggregate by buffer (assuming buffer identifier S_EPA_I)
# buffers_outcome <- intersections %>%
#   group_by(S_EPA_I) %>%
#   summarise(
#     COUNTY          = first(COUNTY),
#     total_VLw_Brth   = sum(VLw_Brth_int, na.rm = TRUE),
#     total_Nm_Brth   = sum(Nm_Brth_int, na.rm = TRUE),
#     .groups = 'drop'
#   ) %>%
#   mutate(
#     Ylow = total_VLw_Brth / total_Nm_Brth
#   )
# 
# 
# summary(buffers_outcome$Ylow) # this looks super normal
# hist(buffers_outcome$Ylow)
# 
# # Merge with buffers
buffers_outcome <- left_join(buffers,
                             st_drop_geometry(buffers_outcome %>%
                                                select(S_EPA_I,Y,COUNTY)),
                             by = c('S_EPA_I' = 'S_EPA_I'))
summary(buffers_outcome$Y)

#################################### What about baseline birth weight ####################################
mun_baseline <-  st_read('shapefiles/Municipality_Low_Birthweight_Rate_2001_2005.shp')

mun_baseline$mun_id <- 1:nrow(mun_baseline)
# For each buffer, assign to it a Y value that is average percent low bw
# Transform municipality geometry to the buffer's CRS
buffers <- st_transform(buffers, st_crs(mun_baseline))

# Compute intersections between buffers and municipalities
intersections <- st_intersection(buffers, mun_baseline)

# Compute area of each intersection
intersections <- intersections %>% 
  mutate(int_area = as.numeric(st_area(.)))

# Assume each municipality has a total area field (mun_area)
# If not, compute it:
mun_baseline <- mun_baseline %>% mutate(mun_area = as.numeric(st_area(.)))

# Join municipality total area to intersections (if not already present)
# Here we assume an identifier 'mun_id' exists in both data frames
intersections <- intersections %>% 
  left_join(mun_baseline %>% st_set_geometry(NULL) %>% select(mun_id, mun_area), by = "mun_id")

# Calculate the area fraction for each intersection
intersections <- intersections %>% 
  mutate(area_fraction = int_area / mun_area)

# Estimate births within each intersection
intersections <- intersections %>% 
  mutate(
    Nm_Brth_int = Nm_Brth * area_fraction,
    Lw_Brth_int = Lw_Brth * area_fraction
  )

# Aggregate by buffer (assuming buffer identifier S_EPA_I)
buffers_outcome_baseline <- intersections %>%
  group_by(S_EPA_I) %>%
  summarise(
    COUNTY          = first(COUNTY),
    total_Lw_Brth   = sum(Lw_Brth_int, na.rm = TRUE),
    total_Nm_Brth   = sum(Nm_Brth_int, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    Ybaseline = total_Lw_Brth / total_Nm_Brth
  )


summary(buffers_outcome_baseline$Ybaseline) # this looks super normal
# Merge buffers_outcome_baseline with buffers_outcome using S_EPA_I
buffers_outcome <- left_join(buffers_outcome, 
                             st_drop_geometry(dplyr::select(buffers_outcome_baseline, S_EPA_I, Ybaseline)), 
                             by = 'S_EPA_I')
hist(buffers_outcome$Y)
hist(buffers_outcome$Ybaseline)
cor(buffers_outcome$Y, buffers_outcome$Ybaseline, use = 'complete.obs') # 0.41

####################### MERGE WITH DESIGN MATRIX  #######################
# Subset to only buffers with Y
buffers_outcome <- buffers_outcome[!is.na(buffers_outcome$Y),]

n <- nrow(buffers_outcome)
X <- buffers_outcome[c('population_density',
               'percent_hispanic',
               'percent_black',
               'percent_indigenous',
               'percent_asian',
               'median_household_income',
               'median_house_value',
               'percent_poverty',
               'percent_high_school_grad',
               'median_year_built',
               'Ybaseline')]
X <- cbind.data.frame(Intercept = 1, X) # Add intercept

X <- st_drop_geometry(X)
X[,2:12] <- scale(X[,2:12])
X <- X[,-13]
X <- as.matrix(X)

dmat <- distm(cbind(buffers_outcome$Longitd, 
                    buffers_outcome$Latitud), fun = distHaversine)
# Convert distances to kilometers
dmat <- dmat / 1000
n <- nrow(dmat)

# Create an empty sparse matrix for 2-nearest neighbors
nn_matrix <- Matrix(0, n, n, sparse = TRUE)

# Loop through each row to find 2-nearest neighbors
for (i in 1:n) {
  # Get indices of the two smallest distances (excluding the diagonal)
  nearest_indices <- order(dmat[i, ], decreasing = FALSE)[2:6]
  
  # Set these indices to 1 in the adjacency matrix
  nn_matrix[i, nearest_indices] <- 1
  nn_matrix[nearest_indices, i] <- 1
}

# Create adjacency matrix: two points are adjacent if they are within 10k of each other
adjacency_matrix <- nn_matrix  #dmat < 50 #
buffers_outcome$cluster[buffers_outcome$cluster == 9] = 40
buffers_outcome$cluster[buffers_outcome$cluster == 36] = 40
clusters <- buffers_outcome$cluster

kappa <- 0.1 # 0.2 # as you increase, becomes more extreme
rangec <- 300
phic <- rangec/(2*sqrt(kappa))
S <- geoR::matern(u =  dmat/10000, phi = phic, kappa = kappa)
E <- eigen(S)
Egp <- E$values
Vgp <- E$vectors
#V <- E$vectors[,2:6] %*% diag(E$values[2:6])/median(E$values[2:6])
Sigma <- diag(n) + 10*S
Sigmainvgp <- solve(Sigma)
# add RE eigens to V
statefactor <- factor(clusters)
Jks <- list()
for (k in unique(statefactor)){
  nk <- sum(statefactor == k)
  Jks[[length(Jks) + 1]] <- matrix(1, nrow = nk, ncol = nk)
}
S <- bdiag(Jks)
E <- eigen(S)
Ere <- E$values
Vre <- E$vectors
Sigma <- diag(n) + 10*S
Sigmainvre <- solve(Sigma)
# add CAR eigens to V
L <- diag(rowSums(adjacency_matrix)) - adjacency_matrix
S <- ginv(as.matrix(L))
E <- eigen(S)
Ecar <- Re(E$values)
Vcar <- Re(E$vectors)
Sigma <- diag(n) + 10*S
Sigmainvcar <- solve(Sigma)

tols <- calculate_optimal_tolerances(X = X, 
                                     Z = buffers_outcome$Z, 
                                     Sigmainvre = Sigmainvre,
                                     Sigmainvcar = Sigmainvcar,
                                     Sigmainvgp = Sigmainvgp,
                                     Vre = Vre,
                                     Ere = Ere,
                                     Vcar = Vcar,
                                     Ecar = Ecar,
                                     Vgp = Vgp,
                                     Egp = Egp)
OLS_att <- fit_method(Y = 100*buffers_outcome$Y,
            X,
            Z = buffers_outcome$Z,
            Vre = Vre,
            Vcar = Vcar,
            Vgp = Vgp,
            Sigmainvre = Sigmainvre,
            Sigmainvcar = Sigmainvcar,
            Sigmainvgp = Sigmainvgp,
            tols = NULL,
            method = 'OLS',
            neigen = 50)

RE_att <- fit_method(Y = 100*buffers_outcome$Y,
                      X,
                      Z = buffers_outcome$Z,
                      Vre = Vre,
                      Vcar = Vcar,
                      Vgp = Vgp,
                      Sigmainvre = Sigmainvre,
                      Sigmainvcar = Sigmainvcar,
                      Sigmainvgp = Sigmainvgp,
                      tols = NULL,
                      method = 'RE',
                      neigen = 50)

CAR_att <- fit_method(Y = 100*buffers_outcome$Y,
                     X,
                     Z = buffers_outcome$Z,
                     Vre = Vre,
                     Vcar = Vcar,
                     Vgp = Vgp,
                     Sigmainvre = Sigmainvre,
                     Sigmainvcar = Sigmainvcar,
                     Sigmainvgp = Sigmainvgp,
                     tols = NULL,
                     method = 'CAR',
                     neigen = 50)

GP_att <- fit_method(Y = 100*buffers_outcome$Y,
                      X,
                      Z = buffers_outcome$Z,
                      Vre = Vre,
                      Vcar = Vcar,
                      Vgp = Vgp,
                      Sigmainvre = Sigmainvre,
                      Sigmainvcar = Sigmainvcar,
                      Sigmainvgp = Sigmainvgp,
                      tols = NULL,
                      method = 'GP',
                      neigen = 50)

outcomemods <- lm(100*buffers_outcome$Y ~ buffers_outcome$Ybaseline + buffers_outcome$Z)
outcomemods <- lm(100*buffers_outcome$Y ~ buffers_outcome$Ybaseline + buffers_outcome$Z)

# Augmented balancing weights
outcomemods <- lm(100*buffers_outcome$Y ~ buffers_outcome$Ybaseline + buffers_outcome$Z)

SW_att <- fit_method(Y = 100*buffers_outcome$Y,
                     X,
                     Z = buffers_outcome$Z,
                     Vre = Vre,
                     Vcar = Vcar,
                     Vgp = Vgp,
                     Sigmainvre = Sigmainvre,
                     Sigmainvcar = Sigmainvcar,
                     Sigmainvgp = Sigmainvgp,
                     tols = tols,
                     method = 'SW',
                     neigen = 20)
SW_att

print(paste('OLS_att = ', round(OLS_att,3), 
            'RE_att = ', round(RE_att,3),
            'CAR_att = ', round(CAR_att,3),
            'GP_att = ', round(GP_att,3),
            'SW_att = ', round(SW_att,3)))
print(mean(buffers_outcome$Y))
print(mean(buffers_outcome$Y[buffers_outcome$Z == 1]))

OLS_att_boot <- boot_func(Y = 100*buffers_outcome$Y,
                         X,
                         Z = buffers_outcome$Z,
                         Vre = Vre,
                         Vcar = Vcar,
                         Vgp = Vgp,
                         Sigmainvre = Sigmainvre,
                         Sigmainvcar = Sigmainvcar,
                         Sigmainvgp = Sigmainvgp,
                         tols = tols,
                         method = 'OLS',
                         neigen = 25)
# Print confidence interval for OLS
OLS_lower <- OLS_att - 1.96*sd(OLS_att_boot)
OLS_upper <- OLS_att + 1.96*sd(OLS_att_boot)
print(paste('OLS_att_boot = ', round(OLS_att,3), 
            'CI = [', round(OLS_lower,3), ', ', round(OLS_upper,3), ']'))


RE_att_boot <- boot_func(Y = 100*buffers_outcome$Y,
                          X,
                          Z = buffers_outcome$Z,
                          Vre = Vre,
                          Vcar = Vcar,
                          Vgp = Vgp,
                          Sigmainvre = Sigmainvre,
                          Sigmainvcar = Sigmainvcar,
                          Sigmainvgp = Sigmainvgp,
                          tols = tols,
                          method = 'RE',
                          neigen = 25)
# Print confidence interval for RE
RE_lower <- RE_att - 1.96*sd(RE_att_boot)
RE_upper <- RE_att + 1.96*sd(RE_att_boot)
print(paste('RE_att_boot = ', round(RE_att,3), 
            'CI = [', round(RE_lower,3), ', ', round(RE_upper,3), ']'))

CAR_att_boot <- boot_func(Y = 100*buffers_outcome$Y,
                          X,
                          Z = buffers_outcome$Z,
                          Vre = Vre,
                          Vcar = Vcar,
                          Vgp = Vgp,
                          Sigmainvre = Sigmainvre,
                          Sigmainvcar = Sigmainvcar,
                          Sigmainvgp = Sigmainvgp,
                          tols = tols,
                          method = 'CAR',
                          neigen = 25)
# Print confidence interval for CAR
CAR_lower <- CAR_att - 1.96*sd(CAR_att_boot)
CAR_upper <- CAR_att + 1.96*sd(CAR_att_boot)
print(paste('CAR_att_boot = ', round(CAR_att,3), 
            'CI = [', round(CAR_lower,3), ', ', round(CAR_upper,3), ']'))

GP_att_boot <- boot_func(Y = 100*buffers_outcome$Y,
                          X,
                          Z = buffers_outcome$Z,
                          Vre = Vre,
                          Vcar = Vcar,
                          Vgp = Vgp,
                          Sigmainvre = Sigmainvre,
                          Sigmainvcar = Sigmainvcar,
                          Sigmainvgp = Sigmainvgp,
                          tols = tols,
                          method = 'GP',
                          neigen = 25)
# Print confidence interval for GP
GP_lower <- GP_att - 1.96*sd(GP_att_boot)
GP_upper <- GP_att + 1.96*sd(GP_att_boot)
print(paste('GP_att_boot = ', round(GP_att,3), 
            'CI = [', round(GP_lower,3), ', ', round(GP_upper,3), ']'))

# Print confidence interval for SW
SW_att_boot <- boot_func(Y = 100*buffers_outcome$Y,
                          X[,-12],
                          Z = buffers_outcome$Z,
                          Vre = Vre,
                          Vcar = Vcar,
                          Vgp = Vgp,
                          Sigmainvre = Sigmainvre,
                          Sigmainvcar = Sigmainvcar,
                          Sigmainvgp = Sigmainvgp,
                          tols = tols,
                          method = 'SW',
                          neigen = 25)
# Print confidence interval for SW
SW_lower <- SW_att - 1.96*sd(SW_att_boot)
SW_upper <- SW_att + 1.96*sd(SW_att_boot)
print(paste('SW_att_boot = ', round(SW_att,3), 
            'CI = [', round(SW_lower,3), ', ', round(SW_upper,3), ']'))

summary_df <- tibble(
  Method   = c("OLS","RE","CAR","GP","SW"),
  Estimate = c(OLS_att, RE_att, CAR_att, GP_att, SW_att),
  lower    = c(OLS_lower, RE_lower, CAR_lower, GP_lower, SW_lower),
  upper    = c(OLS_upper, RE_upper, CAR_upper, GP_upper, SW_upper)
)
# Mutate summary_df so estimates are reported as this
#mean(100*buffers_outcome$Y[buffers_outcome$Z == 1])/(100*mean(buffers_outcome$Y[buffers_outcome$Z == 1])-SW_att)
summary_df <- summary_df %>%
  mutate(Estimate =mean(100*buffers_outcome$Y[buffers_outcome$Z == 1])/(mean(100*buffers_outcome$Y[buffers_outcome$Z == 1]) - Estimate),
         lower = mean(100*buffers_outcome$Y[buffers_outcome$Z == 1])/(mean(100*buffers_outcome$Y[buffers_outcome$Z == 1]) - lower),
         upper = mean(100*buffers_outcome$Y[buffers_outcome$Z == 1])/(mean(100*buffers_outcome$Y[buffers_outcome$Z == 1]) - upper))

# turn Method into a factor with the desired order
summary_df <- summary_df %>%
  mutate(Method = factor(Method, levels = c("OLS","RE","CAR","GP","SW")))

# plot
png("images/results_superfund_att.png", 
    width = 2000, 
    height = 700, 
    res = 270)
ggplot(summary_df, aes(x = Method, y = Estimate)) +
  # 95% CI error bars
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  # point estimates
  geom_point(size = 3, shape = 18, color = "blue") +
  # horizontal zero line
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)) + 
  labs(
    title = "Effect of Superfund cleanup on percent of low birth weight births in NJ and PA",
    x     = "Method",
    y     = "Estimate of the ATT"
  )
dev.off()

# Table of results instead
print(xtable(summary_df, digits = 3), include.rownames = F)

############################ SENSITIVITY ANALYSIS VARYING NONLINEARITY #########################
# Sensitivity 1: SW with algorithm
t_ind <- buffers_outcome$Z
bal_cov <- cbind.data.frame(X,
                            Vre[,1:25],# %*% diag(sqrt(pmax(Ere[1:25],0))), # scale eigens by sqrt eigenvalues
                            Vcar[,1:25],# %*% diag(sqrt(pmax(Ecar[1:25],0))), 
                            Vgp[,1:25])# %*% diag(sqrt(pmax(Egp[1:25],0)))) 
colnames(bal_cov) <- paste0('X', 1:ncol(bal_cov))
data_frame <- as.data.frame(cbind(t_ind, bal_cov))
t_ind <- "t_ind"
bal <- list()
bal$bal_gri <- c(0.1,0.2,0.5) # grid of tuning parameters, 0.1 was smallest of defaults that produced a solution
bal$bal_cov <- colnames(bal_cov)[-1]
bal$bal_alg = T # tuning algorithm in Wang and Zubizarreta (2020) used for automatically selecting the degree of approximate covariates balance.
bal$bal_sam = 1000
bal$bal_std <- 'group'
sbwatttun_object = sbw(dat = data_frame, ind = t_ind, bal = bal, 
                       sol = list(sol_nam = "quadprog"), 
                       par = list(par_est = "att", par_tar = NULL))
SW_alg_att <- sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 1]*100*buffers_outcome$Y[buffers_outcome$Z == 1]) - 
  sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 0]*100*buffers_outcome$Y[buffers_outcome$Z == 0])
# Print causal risk ratio estimate
mean(100*buffers_outcome$Y[buffers_outcome$Z == 1])/(mean(100*buffers_outcome$Y[buffers_outcome$Z == 1]) - SW_alg_att) # 0.94

# Sensitivity 2: Nonlinear only X, algorithm
t_ind <- buffers_outcome$Z
vars   <- colnames(X[,-1])
sq     <- paste0("I(", vars, "^2)", collapse = " + ")
full_f <- as.formula(paste("~ (.)^2 +", sq))
Xquad <- model.matrix(full_f, data = as.data.frame(X[,-1]))
bal_cov <- cbind.data.frame(Xquad,
                            Vre[,1:25],
                            Vcar[,1:25], 
                            Vgp[,1:25]) 
data_frame <- as.data.frame(cbind(t_ind, bal_cov))
t_ind <- "t_ind"
bal <- list()
bal$bal_gri <-  c(0.2,0.5,1,2) # grid of tuning parameters, 0.2 was smallest of defaults that produced a solution
bal$bal_cov <- colnames(bal_cov)[-1]
bal$bal_alg = T # tuning algorithm in Wang and Zubizarreta (2020) used for automatically selecting the degree of approximate covariates balance.
bal$bal_sam = 1000
bal$bal_std <- 'group'
sbwatttun_object = sbw(dat = data_frame, ind = t_ind, bal = bal, 
                       sol = list(sol_nam = "quadprog"), 
                       par = list(par_est = "att", par_tar = NULL))
SW_nonlin_att <- sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 1]*100*buffers_outcome$Y[buffers_outcome$Z == 1]) - 
  sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 0]*100*buffers_outcome$Y[buffers_outcome$Z == 0])
# Print causal risk ratio estimate
mean(100*buffers_outcome$Y[buffers_outcome$Z == 1])/(mean(100*buffers_outcome$Y[buffers_outcome$Z == 1]) - SW_nonlin_att) # 0.94 

# Sensitivity 3: Nonlinear only X, manual
t_ind <- buffers_outcome$Z
vars   <- colnames(X[,-1])
sq     <- paste0("I(", vars, "^2)", collapse = " + ")
full_f <- as.formula(paste("~ (.)^2 +", sq))
Xquad <- model.matrix(full_f, data = as.data.frame(X[,-1]))
bal_cov <- cbind.data.frame(Xquad,
                            Vre[,1:25],
                            Vcar[,1:25], 
                            Vgp[,1:25]) 
data_frame <- as.data.frame(cbind(t_ind, bal_cov))
t_ind <- "t_ind"
bal <- list()
bal$bal_cov <- colnames(bal_cov)[-1]
bal$bal_std <- 'manual'
bal$bal_alg <- F
bal$bal_tol <- c(rep(0, ncol(X)-1), # exact balance on linear terms
                 rep(0.5,ncol(bal_cov)-25*3-ncol(X)), # 0.5 calipers for quadratic terms
                 tols$RE[1:25], tols$CAR[1:25], tols$GP[1:25]) # close to exact balance for eigen
sbwatttun_object = sbw(dat = data_frame, ind = t_ind, bal = bal, 
                       sol = list(sol_nam = "quadprog"), 
                       par = list(par_est = "att", par_tar = NULL))
SW_nonlin_att <- sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 1]*100*buffers_outcome$Y[buffers_outcome$Z == 1]) -
  sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 0]*100*buffers_outcome$Y[buffers_outcome$Z == 0])
# Print causal risk ratio estimate
mean(100*buffers_outcome$Y[buffers_outcome$Z == 1])/(mean(100*buffers_outcome$Y[buffers_outcome$Z == 1]) - SW_nonlin_att) # 0.94

# Sensitivity 4: Nonlinear X and eigen, algorithm
t_ind <- buffers_outcome$Z
bal_cov <- cbind.data.frame(X,
                            Vre[,1:25],
                            Vcar[,1:25], 
                            Vgp[,1:25]) 
colnames(bal_cov) <- paste0('X', 1:ncol(bal_cov))
vars   <- names(bal_cov)
sq     <- paste0("I(", vars, "^2)", collapse = " + ")
full_f <- as.formula(paste("~ (.)^2 +", sq))
bal_cov <- model.matrix(full_f, data = bal_cov)
data_frame <- as.data.frame(cbind(t_ind, bal_cov))
t_ind <- "t_ind"
bal <- list()
bal$bal_gri <-  c(1) # grid of tuning parameters 
bal$bal_tol <- bal$bal_gri     # ← explicitly give sbw() something to work with
bal$bal_cov <- colnames(bal_cov)[-1]
bal$bal_alg = T # tuning algorithm in Wang and Zubizarreta (2020) used for automatically selecting the degree of approximate covariates balance.
bal$bal_sam = 100
bal$bal_std <- 'group'
sbwatttun_object = sbw(dat = data_frame, ind = t_ind, bal = bal, 
                       sol = list(sol_nam = "quadprog"), 
                       par = list(par_est = "att", par_tar = NULL))
SW_nonlin_att <- sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 1]*100*buffers_outcome$Y[buffers_outcome$Z == 1]) - 
  sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 0]*100*buffers_outcome$Y[buffers_outcome$Z == 0])
# Print causal risk ratio estimate
mean(100*buffers_outcome$Y[buffers_outcome$Z == 1])/(mean(100*buffers_outcome$Y[buffers_outcome$Z == 1]) - SW_nonlin_att) # 0.94 

# Sensitivity 5: Nonlinear X and eigen, manual
t_ind <- buffers_outcome$Z
bal_cov <- cbind.data.frame(X[,-c(2,7,8,11)],
                            Vre[,1:25],
                            Vcar[,1:25], 
                            Vgp[,1:25]) 
colnames(bal_cov) <- paste0('X', 1:ncol(bal_cov))
vars   <- names(bal_cov)
sq     <- paste0("I(", vars, "^2)", collapse = " + ")
full_f <- as.formula(paste("~ (.)^2 +", sq))
bal_cov <- model.matrix(full_f, data = bal_cov)
data_frame <- as.data.frame(cbind(t_ind, bal_cov))
t_ind <- "t_ind"
bal <- list()
bal$bal_cov <- colnames(bal_cov)[-1]
bal$bal_std <- 'manual'
bal$bal_alg <- F
bal$bal_tol <- c(rep(0.1, ncol(X)-1), # exact balance on linear terms
                 tols$RE[1:25], tols$CAR[1:25], tols$GP[1:25], # close to exact balance for eigen
                 rep(0.5,ncol(bal_cov)-25*3-ncol(X))) # 0.5 calipers for all quadratic terms
sbwatttun_object = sbw(dat = data_frame, ind = t_ind, bal = bal, 
                       sol = list(sol_nam = "quadprog"), 
                       par = list(par_est = "att", par_tar = NULL))
SW_nonlin_att <- sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 1]*100*buffers_outcome$Y[buffers_outcome$Z == 1]) -
  sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 0]*100*buffers_outcome$Y[buffers_outcome$Z == 0])
# Print causal risk ratio estimate
mean(100*buffers_outcome$Y[buffers_outcome$Z == 1])/(mean(100*buffers_outcome$Y[buffers_outcome$Z == 1]) - SW_nonlin_att) # 0.94


t_ind <- buffers_outcome$Z
bal_cov <- cbind.data.frame(X)#,
                            #Vre[,1:25],
                            #Vcar[,1:25], 
                            #Vgp[,1:25]) 
colnames(bal_cov) <- paste0('X', 1:ncol(bal_cov))
vars   <- names(bal_cov)
sq     <- paste0("I(", vars, "^2)", collapse = " + ")
full_f <- as.formula(paste("~ (.)^2 +", sq))
bal_cov <- model.matrix(full_f, data = bal_cov)
data_frame <- as.data.frame(cbind(t_ind, bal_cov))
t_ind <- "t_ind"
bal <- list()
bal$bal_cov <- colnames(bal_cov)[-1]
bal$bal_std <- 'manual'
bal$bal_alg <- F
bal$bal_tol <- c(rep(0.1, ncol(X)-1),#, # exact balance on linear terms
                 #tols$RE[1:25], tols$CAR[1:25], tols$GP[1:25], # close to exact balance for eigen
                 rep(0.5,ncol(bal_cov)-ncol(X))) # 0.5 calipers for all quadratic terms
sbwatttun_object = sbw(dat = data_frame, ind = t_ind, bal = bal, 
                       sol = list(sol_nam = "quadprog"), 
                       par = list(par_est = "att", par_tar = NULL))
SW_nonlin_att <- sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 1]*100*buffers_outcome$Y[buffers_outcome$Z == 1]) -
  sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 0]*100*buffers_outcome$Y[buffers_outcome$Z == 0])
# Print causal risk ratio estimate
mean(100*buffers_outcome$Y[buffers_outcome$Z == 1])/(mean(100*buffers_outcome$Y[buffers_outcome$Z == 1]) - SW_nonlin_att) # 0.94

########################## SENSITIVITY ANALYSIS VARYING NUMBER of EIGEN ##########################
neigens <- 1:50
SW_att_sens <- rep(NA, length(neigens))
for (i in 1:length(neigens)) {
  SW_att_sens[i] <- fit_method(Y = 100*buffers_outcome$Y,
                               X,
                               Z = buffers_outcome$Z,
                               Vre = Vre,
                               Vcar = Vcar,
                               Vgp = Vgp,
                               Sigmainvre = Sigmainvre,
                               Sigmainvcar = Sigmainvcar,
                               Sigmainvgp = Sigmainvgp,
                               tols = tols,
                               method = 'SW',
                               neigen = neigens[i])
}

sens_df <- tibble(
  neigen = neigens,
  att    = SW_att_sens
) %>%
  mutate(
    att_smooth = rollmean(att, k = 10, fill = NA, align = "center")
  )

# 2) Gather the two curves into long form, renaming 'att' to "SW estimate"
line_df <- sens_df %>%
  select(neigen, att, att_smooth) %>%
  rename(
    `SW estimate`                = att,
    `Rolling mean (window = 10)` = att_smooth
  ) %>%
  pivot_longer(
    cols      = -neigen,
    names_to  = "label",
    values_to = "value"
  )

# 3) Data frame for the single horizontal OLS line
hline_df <- tibble(
  label = "OLS estimate",
  y     = OLS_att
)

png("images/sensitivity_analysis.png", 
    width = 1500, 
    height = 1000, 
    res = 200)

# 4) Combined plot with color and linetype legends
ggplot() +
  # both curves
  geom_line(
    data = line_df,
    aes(x = neigen, y = value, color = label, linetype = label),
    size = 1
  ) +
  # single h‐line for OLS
  geom_hline(
    data = hline_df,
    aes(yintercept = y, color = label, linetype = label),
    size = 0.8
  ) +
  # manual color scale
  scale_color_manual(
    name   = NULL,
    values = c(
      "SW estimate"                = "grey50",
      "Rolling mean (window = 10)" = "red",
      "OLS estimate"               = "darkgreen"
    )
  ) +
  # manual linetype scale
  scale_linetype_manual(
    name   = NULL,
    values = c(
      "SW estimate"                = "solid",
      "Rolling mean (window = 10)" = "dashed",
      "OLS estimate"               = "dotted"
    )
  ) +
  labs(
    title = "Sensitivity Analysis for the Spatial Weighting Method",
    x     = "Number of Eigenvectors",
    y     = "ATT Estimate"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor= element_blank(),
    legend.position = "bottom"
  )
dev.off()

