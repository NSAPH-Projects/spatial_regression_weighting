library(ggplot2)
library(sf)
library(spdep)
library(tidyr)
library(dplyr)
library(scales)
library(stringr)

source('../funcs.R')
states <- st_read('../data/tl_2010_us_state10/tl_2010_us_state10.shp')

# Read in preprocessed data
load('../data/preprocessed_superfunds.RData')
n <- nrow(buffers)
us_outline <- ne_countries(scale = "medium", country = "United States of America", returnclass = "sf")

buffer_centroids <- st_centroid(buffers)

X[,2:12] <- scale(X[,2:12])
X <- as.matrix(X)

# RE model
clusters <- substr(buffers$S_EPA_I, 1, 2)
statefactor <- factor(clusters)

# Cluster matrix of indicators
clustermat <- model.matrix(~clusters-1)
statenames <- gsub('clusters', '', colnames(clustermat))

# Goal is to plot imbalance for each of the 50 eigen as a function of rho^2
rho2s <- seq(0, 1, by = 0.01)
imbalancemat <- cbind.data.frame('rho2'= rho2s, 
                                 matrix(NA, nrow = length(rho2), ncol = ncol(clustermat)))
colnames(imbalancemat)[-1] <- statenames
imbalancemat$dispersion <- NA
for (i in 1:length(rho2s)){
  rho2 <- rho2s[i]
  Jks = list()
  for (k in unique(statefactor)){
    nk = sum(statefactor == k)
    Jks[[length(Jks) + 1]] = matrix(1, nrow = nk, ncol = nk)
  }
  S = bdiag(Jks)
  Sigma <- diag(n) + rho2*S
  Sigmainv <- solve(Sigma)
  REiw <- impliedweightsgeneral(X = X,
                                Z = Z,
                                Sigmainv = Sigmainv)
  imbalancemat[i,2:(ncol(imbalancemat)-1)] <- abs(REiw[Z == 1] %*% clustermat[Z == 1,] - REiw[Z == 0] %*% clustermat[Z == 0,])
  imbalancemat$dispersion[i] <- compute_dispersion(REiw)
}

df_long <- imbalancemat %>%
  pivot_longer(
    cols      = 2:(ncol(imbalancemat)-1),
    names_to  = "state",
    values_to = "imbalance"
  )

# 2) Compute ranges, scale, and offset
rng_imb      <- range(df_long$imbalance, na.rm = TRUE)
rng_disp <- range(imbalancemat$dispersion, na.rm = TRUE)
scale_factor <- diff(rng_imb)/diff(rng_disp)

# for shifting axes
disp_offset <- -rng_disp[1]          # shift dispersion min to 0
shift_amount <- rng_imb[1]           # then lift that to imbalance min

png('images/re_fe.png', width = 1500, height = 1000, res = 170)
ggplot() +
  # state curves
  geom_line(
    data = df_long,
    aes(x = rho2, y = imbalance, color = state),
    size = 0.8
  ) +
  geom_line(
    data = imbalancemat,
    aes(
      x = rho2,
      y = (dispersion + disp_offset) * scale_factor + shift_amount,
      color = "dispersion"
    ),
    linetype = "dashed",
    size     = 1
  ) +
  # Legend: all states + ESS
  scale_color_manual(
    name   = NULL,
    # order legend so ESS comes first
    breaks = c("dispersion", statenames),
    values = c(
      dispersion = "black",
      # then your states, in the same order
      setNames(
        hue_pal()(length(statenames)),
        statenames
      )
    )
  )+
  # Dual axes: left = imbalance, right = ESS (undo offset & scale)
  scale_y_continuous(
    name     = "Absolute imbalance",
    sec.axis = sec_axis(
      trans = ~ (. - shift_amount) / scale_factor - disp_offset,
      name  = "Weight Dispersion"
    )
  ) +
  # Labels and title with rho^2
  labs(
    x     = expression(rho^2),
    title = expression(
      paste(
        "Imbalance on state-level indicators and weight dispersion vs. ", 
        rho^2
      )
    )
  ) +
  theme_minimal() +
  theme(
    axis.title.y.left  = element_text(vjust = 2),
    axis.title.y.right = element_text(vjust = 2)
  )
dev.off()