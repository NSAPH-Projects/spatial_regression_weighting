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

set.seed(111)

setwd('/Users/sophie/Documents/implied_weights/superfunds/')

source('../impliedweights_randeffs.R')

states <- st_read('/Users/sophie/Downloads/tl_2010_us_state10/tl_2010_us_state10.shp')

# Read in preprocessed data
load('data/preprocessed_superfunds.RData')
n <- nrow(buffers)
buffer_centroids <- st_centroid(buffers)

Sigmainvs <- compute_allSigmainv(adjacency_matrix = adjacency_matrix,
                                 statefactor = clusters,
                                 sig2gam = 1,
                                 sig2eps = 1,
                                 phi = 0.9,
                                 distmat = dmat/10000)
# RE
# Create three vectors with diff smoothness
x1 <- clusters + rnorm(length(clusters), sd = 0.01)
x2 <- clusters + rnorm(length(clusters), sd = 5)
x3 <- clusters + rnorm(length(clusters), sd = 50)
s1 <- smoothness(x = x1, Sigmainv = Sigmainvs$RE)
s2 <- smoothness(x = x2, Sigmainv = Sigmainvs$RE)
s3 <- smoothness(x = x3, Sigmainv = Sigmainvs$RE)

# Plot x1, x2, x3 on map
buffers_merged <- cbind(buffer_centroids, x1, x2, x3)
gx1 <- ggplot() +
  geom_sf(data = states, fill = "black", color = "white", linetype = "solid") +
  geom_sf(data = buffers_merged, 
          aes(color = x1), size = 1.5) +
  scale_color_viridis_c() +
  coord_sf(xlim = c(-125, -65), ylim = c(25, 50), expand = FALSE) +
  labs(
    title = bquote(s(x * "," ~ Sigma) == .(round(s1,2))),
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
  geom_sf(data = states, fill = "black", color = "white", linetype = "solid") +
  geom_sf(data = buffers_merged, 
          aes(color = x2), size = 1.5) +
  scale_color_viridis_c() +
  coord_sf(xlim = c(-125, -65), ylim = c(25, 50), expand = FALSE) +
  labs(
    title = bquote(s(x * "," ~ Sigma) == .(round(s2,2))),
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
  geom_sf(data = states, fill = "black", color = "white", linetype = "solid") +
  geom_sf(data = buffers_merged, 
          aes(color = x3), size = 1.5) +
  scale_color_viridis_c() +
  coord_sf(xlim = c(-125, -65), ylim = c(25, 50), expand = FALSE) +
  labs(
    title = bquote(s(x * "," ~ Sigma) == .(round(s3,2))),
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
png('images/smoothness_RE.png', width = 2000, height = 500, res = 150)
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

# CAR
# Create three vectors with diff smoothness
E <- eigen(Sigmainvs$CAR)
x1 <- E$vectors[,(n-30):(n)] %*% rnorm(31)
x2 <- E$vectors[,(n-500):(n)] %*% rnorm(501)
x3 <- E$vectors[,(n-1000):(n)] %*% rnorm(1001)
s1 <- smoothness(x = x1, Sigmainv = Sigmainvs$CAR)
s2 <- smoothness(x = x2, Sigmainv = Sigmainvs$CAR)
s3 <- smoothness(x = x3, Sigmainv = Sigmainvs$CAR)

buffers_merged <- cbind(buffer_centroids, x1, x2, x3)
gx1 <- ggplot() +
  geom_sf(data = states, fill = "black", color = "white", linetype = "solid") +
  geom_sf(data = buffers_merged, 
          aes(color = x1), size = 1.5) +
  scale_color_viridis_c() +
  coord_sf(xlim = c(-125, -65), ylim = c(25, 50), expand = FALSE) +
  labs(
    title = bquote(s(x * "," ~ Sigma) == .(round(smoothness(x = x1, Sigmainv = Sigmainvs$CAR),2))),
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
  geom_sf(data = states, fill = "black", color = "white", linetype = "solid") +
  geom_sf(data = buffers_merged, 
          aes(color = x2), size = 1.5) +
  scale_color_viridis_c() +
  coord_sf(xlim = c(-125, -65), ylim = c(25, 50), expand = FALSE) +
  labs(
    title = bquote(s(x * "," ~ Sigma) == .(round(smoothness(x = x2, Sigmainv = Sigmainvs$CAR),2))),
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
  geom_sf(data = states, fill = "black", color = "white", linetype = "solid") +
  geom_sf(data = buffers_merged, 
          aes(color = x3), size = 1.5) +
  scale_color_viridis_c() +
  coord_sf(xlim = c(-125, -65), ylim = c(25, 50), expand = FALSE) +
  labs(
    title = bquote(s(x * "," ~ Sigma) == .(round(smoothness(x = x3, Sigmainv = Sigmainvs$CAR),2))),
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
png('images/smoothness_CAR.png', width = 2000, height = 500, res = 150)
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

# GP
# Create three vectors with diff smoothness
E <- eigen(Sigmainvs$GP)
x1 <- E$vectors[,(n-10):(n)] %*% rnorm(11)
x2 <- E$vectors[,(n-50):(n)] %*% rnorm(51)
x3 <- E$vectors[,(n-500):(n)] %*% rnorm(501)
s1 <- smoothness(x = x1, Sigmainv = Sigmainvs$GP)
s2 <- smoothness(x = x2, Sigmainvs$GP)
s3 <- smoothness(x = x3, Sigmainvs$GP)

buffers_merged <- cbind(buffer_centroids, x1, x2, x3)
gx1 <- ggplot() +
  geom_sf(data = states, fill = "black", color = "white", linetype = "solid") +
  geom_sf(data = buffers_merged, 
          aes(color = x1), size = 1.5) +
  scale_color_viridis_c() +
  coord_sf(xlim = c(-125, -65), ylim = c(25, 50), expand = FALSE) +
  labs(
    title = bquote(s(x * "," ~ Sigma) == .(round(smoothness(x = x1, Sigmainv = Sigmainvs$GP),2))),
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
  geom_sf(data = states, fill = "black", color = "white", linetype = "solid") +
  geom_sf(data = buffers_merged, 
          aes(color = x2), size = 1.5) +
  scale_color_viridis_c() +
  coord_sf(xlim = c(-125, -65), ylim = c(25, 50), expand = FALSE) +
  labs(
    title = bquote(s(x * "," ~ Sigma) == .(round(smoothness(x = x2, Sigmainv = Sigmainvs$GP),2))),
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
  geom_sf(data = states, fill = "black", color = "white", linetype = "solid") +
  geom_sf(data = buffers_merged, 
          aes(color = x3), size = 1.5) +
  scale_color_viridis_c() +
  coord_sf(xlim = c(-125, -65), ylim = c(25, 50), expand = FALSE) +
  labs(
    title = bquote(s(x * "," ~ Sigma) == .(round(smoothness(x = x3, Sigmainv = Sigmainvs$GP),2))),
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
png('images/smoothness_GP.png', width = 2000, height = 500, res = 150)
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

