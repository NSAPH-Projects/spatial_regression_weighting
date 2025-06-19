library(sf)
library(dplyr)
library(utils)
library(gridExtra)
library(ggplot2)
library(sbw)
library(sf)
library(MASS)
source('../funcs.R')

# Import geographies
states <- st_read('../data/tl_2010_us_state10/tl_2010_us_state10.shp')

set.seed(1234)
# Read in X and Z data
load('../data/preprocessed_superfunds.RData')
lat <- buffers$Latitud
long <- buffers$Longitd
n <- nrow(buffers)
X <- st_drop_geometry(X)
X[,2:11] <- scale(X[,2:11])
X <- X[,-12]
X <- as.matrix(X)
kappa <- 0.1 # 0.2 # as you increase, becomes more extreme
rangec <- 300
phic <- rangec/(2*sqrt(kappa))
S <- geoR::matern(u =  dmat/10000, phi = phic, kappa = kappa)
E <- eigen(S)
Egp <- E$values
Vgp <- E$vectors
Wgp <- exp(-dmat/500)
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
Wre <- bdiag(Jks)
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
E <- eigen(adjacency_matrix)
Wcar <- E$vectors %*% diag(E$values^100) %*% t(E$vectors)
Sigma <- diag(n) + 10*S
Sigmainvcar <- solve(Sigma)

# Save simulation data
simlist <- list(
  'X' = X,
  'Z' = Z,
  #'V' = V,
  'Vre' = Vre,
  'Ere' = Ere,
  'Vcar' = Vcar,
  'Ecar' = Ecar,
  'Vgp' = Vgp,
  'Egp' = Egp,
  'Sigmainvre' = Sigmainvre,
  'Sigmainvcar' = Sigmainvcar,
  'Sigmainvgp' = Sigmainvgp,
  'Wre' = Wre,
  'Wcar' = Wcar,
  'Wgp' = Wgp,
  'lat' = lat,
  'long' = long
)
save(simlist, 
     file = 'sim.RData')

####################### Plot one generation of U from 3 confounding mechanisms. #########################
buffers_merged <- st_centroid(buffers)

# RE
buffers_merged$Ure <- as.numeric(gen_U(Z = Z, 
           Wre = Wre, 
           Wcar = Wcar, 
           Wgp = Wgp,
           smoothing = 'clusterbased'
           ))

png('images/Ure.png', width = 1300, height = 1000, res = 250)
ggplot() +
  geom_sf(data = states, fill = NA, color = "darkgray", linetype = "solid", size = 0.5) +
  geom_sf(data = buffers_merged, 
          aes(fill = Ure, shape = factor(Z)), size = 3, stroke = 0) +
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

buffers_merged$Ucar <- as.numeric(gen_U(Z = Z, 
                            Wre = Wre, 
                            Wcar = Wcar, 
                            Wgp = Wgp,
                            smoothing = 'adjacencybased'
))

# CAR
png('images/Ucar.png', width = 1300, height = 1000, res = 250)
ggplot() +
  geom_sf(data = states, fill = NA, color = "darkgray", linetype = "solid", size = 0.5) +
  geom_sf(data = buffers_merged, 
          aes(fill = buffers_merged$Ucar, shape = factor(Z)), size = 3, stroke = 0) +
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

# GP
buffers_merged$Ugp <- as.numeric(gen_U(Z = Z, 
                            Wre = Wre, 
                            Wcar = Wcar, 
                            Wgp = Wgp,
                            smoothing = 'distancebased'
))
# Order buffers_merged by the value of Ugp
buffers_merged <- buffers_merged[order(buffers_merged$Ugp), ]
png('images/Ugp.png', width = 1300, height = 1000, res = 250)
ggplot() +
  geom_sf(data = states, fill = NA, color = "darkgray", linetype = "solid", size = 0.5) +
  geom_sf(data = buffers_merged, 
          aes(fill = buffers_merged$Ugp, shape = factor(Z)), size = 3, stroke = 0) +
  # drop the shape legend:
  scale_shape_manual(values = c("0" = 21, "1" = 24),
                     guide = "none") +
  # relabel the fill legend:
  scale_fill_viridis_c(name = "U", limits = c(0.03, 0.16)) +
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
