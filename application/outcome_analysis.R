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
library(SuperLearner)
library(sbw)

source('../funcs.R')

states <- st_read('../data/tl_2010_us_state10/tl_2010_us_state10.shp')

# Read in preprocessed data
load('../data/preprocessed_superfunds.RData')
n <- nrow(buffers)
us_states <- ne_states(country = "United States of America", returnclass = "sf")

############################ Fit methods ###############################
mun <-  st_read('../shapefiles/Municipality_Birth_Statistics_2016_2020.shp')

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
  left_join(mun %>% st_set_geometry(NULL) %>% dplyr::select(mun_id, mun_area), by = "mun_id")

# Calculate the area fraction for each intersection
intersections <- intersections %>% 
  mutate(area_fraction = int_area / mun_area)

# Estimate births within each intersection, assuming uniform distribution
intersections <- intersections %>% 
  mutate(
    Lw_Brth_int = Lw_Brth * area_fraction,
    Nm_Br_L_int = Nm_Br_L * area_fraction,
    Ttl_Brt_int = Ttl_Brt * area_fraction,
    Nm_Br_B_int = Nm_Br_B * area_fraction,
    Infnt_M_int = Infnt_M * area_fraction,
    Nm_Br_I_int = Nm_Br_I * area_fraction,
    Preterm_int = Preterm * area_fraction,
    Nm_Br_P_int = Nm_Br_P * area_fraction,
    Avg_G_A_int = Avg_G_A * area_fraction,
    Nm_Br_A_int = Nm_Br_A * area_fraction
  )

# Aggregate by buffer (assuming buffer identifier S_EPA_I)
buffers_outcome <- intersections %>%
  group_by(S_EPA_I) %>%
  summarise(
    COUNTY          = first(COUNTY),
    total_Lw_Brth_int = sum(Lw_Brth_int, na.rm = TRUE),
    total_Nm_Br_L_int = sum(Nm_Br_L_int, na.rm = TRUE),
    total_Ttl_Brt_int = sum(Ttl_Brt_int, na.rm = TRUE),
    total_Nm_Br_B_int = sum(Nm_Br_B_int, na.rm = TRUE),
    total_Infnt_M_int = sum(Infnt_M_int, na.rm = TRUE),
    total_Nm_Br_I_int = sum(Nm_Br_I_int, na.rm = TRUE),
    total_Preterm_int = sum(Preterm_int, na.rm = TRUE),
    total_Nm_Br_P_int = sum(Nm_Br_P_int, na.rm = TRUE),
    total_Avg_G_A_int = sum(Avg_G_A_int, na.rm = TRUE),
    total_Nm_Br_A_int = sum(Nm_Br_A_int, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    Y_lbw = 100*total_Lw_Brth_int/total_Nm_Br_L_int,
    Y_avgbw = (1/453.592)*total_Ttl_Brt_int/total_Nm_Br_B_int,
    Y_infmort = 100*total_Infnt_M_int / total_Nm_Br_I_int,
    Y_preterm = 100*total_Preterm_int / total_Nm_Br_P_int,
    Y_avggest = total_Avg_G_A_int / total_Nm_Br_A_int
  )

summary(buffers_outcome$Y_lbw)
hist(buffers_outcome$Y_lbw)
summary(buffers_outcome$Y_avgbw)
hist(buffers_outcome$Y_avgbw)
summary(buffers_outcome$Y_infmort)
hist(buffers_outcome$Y_infmort)
summary(buffers_outcome$Y_preterm)
hist(buffers_outcome$Y_preterm)
summary(buffers_outcome$Y_avggest)
hist(buffers_outcome$Y_avggest)

buffers_outcome <- left_join(buffers,
                             st_drop_geometry(buffers_outcome %>%
                                                dplyr::select(S_EPA_I,
                                                       Y_lbw,
                                                       Y_avgbw, 
                                                       Y_infmort, 
                                                       Y_preterm,
                                                       Y_avggest,
                                                       COUNTY)),
                             by = c('S_EPA_I' = 'S_EPA_I'))


####################### MERGE WITH DESIGN MATRIX  #######################
# Subset to only buffers with outcome data (NJ)
buffers_outcome <- buffers_outcome[!is.na(buffers_outcome$Y_lbw),]

n <- nrow(buffers_outcome)
lat <- buffers_outcome$Latitud
long <- buffers_outcome$Longitd
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
                       'Sit_Scr')]

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

# Create an empty sparse matrix for 5-nearest neighbors
nn_matrix <- Matrix(0, n, n, sparse = TRUE)

# Loop through each row to find 5-nearest neighbors
for (i in 1:n) {
  # Get indices of the five smallest distances (excluding the diagonal)
  nearest_indices <- order(dmat[i, ], decreasing = FALSE)[2:6]
  
  # Set these indices to 1 in the adjacency matrix
  nn_matrix[i, nearest_indices] <- 1
  nn_matrix[nearest_indices, i] <- 1
}

# Create adjacency matrix
adjacency_matrix <- nn_matrix 
buffers_outcome$cluster <- 32 # assign border PA sites to NJ (ignore state in this analysis)
clusters <- buffers_outcome$cluster

rho2 <- 10
kappa <- 10 
rangec <- 500
phic <- rangec/(2*sqrt(kappa))
S <- geoR::matern(u =  dmat, phi = phic, kappa = kappa)
E <- eigen(S)
Egp <- E$values
Vgp <- E$vectors
Sigma <- diag(n) + rho2*S
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
Sigma <- diag(n) + rho2*S
Sigmainvre <- solve(Sigma)
# add CAR eigens to V
L <- diag(rowSums(adjacency_matrix)) - adjacency_matrix
S <- ginv(as.matrix(L))
E <- eigen(S)
Ecar <- Re(E$values)
Vcar <- Re(E$vectors)
Sigma <- diag(n) + rho2*S
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
# Ignore eigenvectors from RE - just one state here
tols$RE <- rep(10^(10), n)

# Fit comparison methods for each outcome
outcomes <- cbind.data.frame(Y_lbw = buffers_outcome$Y_lbw,
             Y_avgbw = buffers_outcome$Y_avgbw,
             Y_infmort = buffers_outcome$Y_infmort,
             Y_preterm = buffers_outcome$Y_preterm,
             Y_avggest = buffers_outcome$Y_avggest)
methods <- c('OLS', 'RE', 'CAR', 'GP', 'spatialcoord')

# Create a dataframe to store estimates and CIs for each method
ests <- data.frame(matrix(NA, nrow = ncol(outcomes), ncol = length(methods)))
colnames(ests) <- methods
rownames(ests) <- c('Y_lbw', 'Y_avgbw', 'Y_infmort', 'Y_preterm', 'Y_avggest')
cis_upper <- data.frame(matrix(NA, nrow = ncol(outcomes), ncol = length(methods)))
colnames(cis_upper) <- methods
rownames(cis_upper) <- c('Y_lbw', 'Y_avgbw', 'Y_infmort', 'Y_preterm', 'Y_avggest')
cis_lower <- data.frame(matrix(NA, nrow = ncol(outcomes), ncol = length(methods)))
colnames(cis_lower) <- methods
rownames(cis_lower) <- c('Y_lbw', 'Y_avgbw', 'Y_infmort', 'Y_preterm', 'Y_avggest')

for (i in 1:ncol(outcomes)){
  for (j in 1:length(methods)){
    method <- methods[j]
    print(method)
    Y <- outcomes[,i]
    
    # Fit the model
    fit <- fit_method(Y = Y,
                      X,
                      Z = buffers_outcome$Z,
                      Vre = Vre,
                      Vcar = Vcar,
                      Vgp = Vgp,
                      Sigmainvre = Sigmainvre,
                      Sigmainvcar = Sigmainvcar,
                      Sigmainvgp = Sigmainvgp,
                      tols = tols,
                      method = method,
                      lat = lat,
                      long = long,
                      neigen = 25,
                      boot = T)
    
    # Store the estimates
    print(c(fit$est,fit$boot_sd))
    ATT_est <- fit$est
    ATT_upper <- fit$est + 1.96*fit$boot_sd
    ATT_lower <- fit$est - 1.96*fit$boot_sd
    ests[i,j] <- mean(Y[buffers_outcome$Z == 1])/(mean(Y[buffers_outcome$Z == 1]) - ATT_est)
    cis_upper[i,j] <- mean(Y[buffers_outcome$Z == 1])/(mean(Y[buffers_outcome$Z == 1]) - ATT_upper)
    cis_lower[i,j] <- mean(Y[buffers_outcome$Z == 1])/(mean(Y[buffers_outcome$Z == 1]) - ATT_lower)
  }
}

SW_att_ests <- data.frame('neigen' = 1:20,
                          'SW_att_lbw' = NA,
                          'SW_att_avgbw' = NA,
                          'SW_att_infmort' = NA,
                          'SW_att_preterm' = NA,
                          'SW_att_gestage' = NA)
SW_att_ci_lowers <- data.frame('neigen' = 1:20,
                              'SW_lower_lbw' = NA,
                              'SW_lower_avgbw' = NA,
                              'SW_lower_infmort' = NA,
                              'SW_lower_preterm' = NA,
                              'SW_lower_gestage' = NA)
SW_att_ci_uppers <- data.frame('neigen' = 1:20,
                              'SW_upper_lbw' = NA,
                              'SW_upper_avgbw' = NA,
                              'SW_upper_infmort' = NA,
                              'SW_upper_preterm' = NA,
                              'SW_upper_gestage' = NA)

for (i in 1:nrow(SW_att_ests)){
  print(paste('neigen', i, 'of', nrow(SW_att_ests)))
  neigen <- SW_att_ests$neigen[i]
  SW_att_fit_lbw <- fit_method(Y = buffers_outcome$Y_lbw,
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
                               neigen = neigen, 
                               boot = T)
  SW_att_lbw <- SW_att_fit_lbw$est
  SW_lower_lbw <- SW_att_fit_lbw$est - 1.96*SW_att_fit_lbw$boot_sd
  SW_upper_lbw <- SW_att_fit_lbw$est + 1.96*SW_att_fit_lbw$boot_sd
  SW_att_ests$SW_att_lbw[i] = mean(buffers_outcome$Y_lbw[buffers_outcome$Z == 1])/(mean(buffers_outcome$Y_lbw[buffers_outcome$Z == 1]) - SW_att_lbw)
  SW_att_ci_lowers$SW_lower_lbw[i] = mean(buffers_outcome$Y_lbw[buffers_outcome$Z == 1])/(mean(buffers_outcome$Y_lbw[buffers_outcome$Z == 1]) - SW_lower_lbw)
  SW_att_ci_uppers$SW_upper_lbw[i] = mean(buffers_outcome$Y_lbw[buffers_outcome$Z == 1])/(mean(buffers_outcome$Y_lbw[buffers_outcome$Z == 1]) - SW_upper_lbw)
  
  # Repeat for average birth weight
  SW_att_fit_avgbw <- fit_method(Y = buffers_outcome$Y_avgbw,
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
                                 neigen = neigen, 
                                 boot = T)
  SW_att_avgbw <- SW_att_fit_avgbw$est
  SW_lower_avgbw <- SW_att_fit_avgbw$est - 1.96*SW_att_fit_avgbw$boot_sd
  SW_upper_avgbw <- SW_att_fit_avgbw$est + 1.96*SW_att_fit_avgbw$boot_sd
  SW_att_ests$SW_att_avgbw[i] = mean(buffers_outcome$Y_avgbw[buffers_outcome$Z == 1])/(mean(buffers_outcome$Y_avgbw[buffers_outcome$Z == 1]) - SW_att_avgbw)
  SW_att_ci_lowers$SW_lower_avgbw[i] = mean(buffers_outcome$Y_avgbw[buffers_outcome$Z == 1])/(mean(buffers_outcome$Y_avgbw[buffers_outcome$Z == 1]) - SW_lower_avgbw)
  SW_att_ci_uppers$SW_upper_avgbw[i] = mean(buffers_outcome$Y_avgbw[buffers_outcome$Z == 1])/(mean(buffers_outcome$Y_avgbw[buffers_outcome$Z == 1]) - SW_upper_avgbw)
  
  # Repeat for infant mortality 
  SW_att_fit_infmort <- fit_method(Y = buffers_outcome$Y_infmort,
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
                                   neigen = neigen, 
                                   boot = T)
  SW_att_infmort <- SW_att_fit_infmort$est
  SW_lower_infmort <- SW_att_fit_infmort$est - 1.96*SW_att_fit_infmort$boot_sd
  SW_upper_infmort <- SW_att_fit_infmort$est + 1.96*SW_att_fit_infmort$boot_sd
  SW_att_ests$SW_att_infmort[i] = mean(buffers_outcome$Y_infmort[buffers_outcome$Z == 1])/(mean(buffers_outcome$Y_infmort[buffers_outcome$Z == 1]) - SW_att_infmort)
  SW_att_ci_lowers$SW_lower_infmort[i] = mean(buffers_outcome$Y_infmort[buffers_outcome$Z == 1])/(mean(buffers_outcome$Y_infmort[buffers_outcome$Z == 1]) - SW_lower_infmort)
  SW_att_ci_uppers$SW_upper_infmort[i] = mean(buffers_outcome$Y_infmort[buffers_outcome$Z == 1])/(mean(buffers_outcome$Y_infmort[buffers_outcome$Z == 1]) - SW_upper_infmort)
  
  # Repeat for preterm birth
  SW_att_fit_preterm <- fit_method(Y = buffers_outcome$Y_preterm,
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
                                   neigen = neigen, 
                                   boot = T)
  SW_att_preterm <- SW_att_fit_preterm$est
  SW_lower_preterm <- SW_att_fit_preterm$est - 1.96*SW_att_fit_preterm$boot_sd
  SW_upper_preterm <- SW_att_fit_preterm$est + 1.96*SW_att_fit_preterm$boot_sd
  SW_att_ests$SW_att_preterm[i] = mean(buffers_outcome$Y_preterm[buffers_outcome$Z == 1])/(mean(buffers_outcome$Y_preterm[buffers_outcome$Z == 1]) - SW_att_preterm)
  SW_att_ci_lowers$SW_lower_preterm[i] = mean(buffers_outcome$Y_preterm[buffers_outcome$Z == 1])/(mean(buffers_outcome$Y_preterm[buffers_outcome$Z == 1]) - SW_lower_preterm)
  SW_att_ci_uppers$SW_upper_preterm[i] = mean(buffers_outcome$Y_preterm[buffers_outcome$Z == 1])/(mean(buffers_outcome$Y_preterm[buffers_outcome$Z == 1]) - SW_upper_preterm)
  
  # Repeat for gestational age
  SW_att_fit_gestage <- fit_method(Y = buffers_outcome$Y_avggest,
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
                                   neigen = neigen, 
                                   boot = T)
  SW_att_gestage <- SW_att_fit_gestage$est
  SW_lower_gestage <- SW_att_fit_gestage$est - 1.96*SW_att_fit_gestage$boot_sd
  SW_upper_gestage <- SW_att_fit_gestage$est + 1.96*SW_att_fit_gestage$boot_sd
  SW_att_ests$SW_att_gestage[i] = mean(buffers_outcome$Y_avggest[buffers_outcome$Z == 1])/(mean(buffers_outcome$Y_avggest[buffers_outcome$Z == 1]) - SW_att_gestage)
  SW_att_ci_lowers$SW_lower_gestage[i] = mean(buffers_outcome$Y_avggest[buffers_outcome$Z == 1])/(mean(buffers_outcome$Y_avggest[buffers_outcome$Z == 1]) - SW_lower_gestage)
  SW_att_ci_uppers$SW_upper_gestage[i] = mean(buffers_outcome$Y_avggest[buffers_outcome$Z == 1])/(mean(buffers_outcome$Y_avggest[buffers_outcome$Z == 1]) - SW_upper_gestage)
}

# Plot estimates and cis vs neigen
ests_long <- SW_att_ests %>%
  pivot_longer(-neigen, names_to = "outcome", values_to = "att") %>%
  mutate(outcome = str_remove(outcome, "^SW_att_"))

lowers_long <- SW_att_ci_lowers %>%
  pivot_longer(-neigen, names_to = "outcome", values_to = "ci_lower") %>%
  mutate(outcome = str_remove(outcome, "^SW_lower_"))

uppers_long <- SW_att_ci_uppers %>%
  pivot_longer(-neigen, names_to = "outcome", values_to = "ci_upper") %>%
  mutate(outcome = str_remove(outcome, "^SW_upper_"))

## 2. Combine---------------------------------
plot_df <- ests_long %>%
  left_join(lowers_long, by = c("neigen", "outcome")) %>%
  left_join(uppers_long, by = c("neigen", "outcome")) %>%
  arrange(outcome, neigen) %>%
  group_by(outcome) %>%
  mutate(att_smooth = rollmean(att, k = 5, fill = NA, align = "center")) %>%
  ungroup()

## 3. Facet labels -----------------------------------------------------------
label_map <- c(
  avgbw   = "Average Birth Weight",
  gestage = "Gestional Age",        # <- requested spelling
  infmort = "Infant Mortality Rate",
  preterm = "Preterm Birth Rate",
  lbw     = "Rate of Low Birth Weight"
)

OLS_ests <- c(ests['Y_lbw', 'OLS'], 
              ests['Y_avgbw', 'OLS'],
              ests['Y_infmort', 'OLS'],
              ests['Y_preterm', 'OLS'],
              ests['Y_avggest', 'OLS'])
hline_df <- tibble(
  outcome = c("lbw", "avgbw", "infmort", "preterm", "gestage"),
  y       = OLS_ests
)


png('images/spatial_weighting_estimates.png',
    width = 2500, height = 2250, res = 300)
ggplot(plot_df, aes(x = neigen, y = att)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
              fill = "grey80", alpha = 0.4) +
  geom_line(                                   # spatial-weighting curve
    aes(color = "Spatial weighting estimate",
        linetype = "Spatial weighting estimate"),
    size = 1
  ) +
  geom_line(                                   # rolling mean
    aes(y = att_smooth,
        color = "Rolling average", #  (k = 5)
        linetype = "Rolling average"), #  (k = 5)
    size = 0.9
  ) +
  geom_hline(                                  # RR = 1 reference
    aes(yintercept = 1,
        color = "Reference (RR = 1)",
        linetype = "Reference (RR = 1)"),
    size = 0.8
  ) +
  geom_hline(                                  # OLS reference lines
    data = hline_df,
    aes(yintercept = y,
        color = "OLS estimate",
        linetype = "OLS estimate"),
    size = 0.5
  ) +
  facet_wrap(
    ~ outcome,
    ncol     = 2,
    scales   = "free_y",
    labeller = labeller(outcome = label_map)
  ) +
  scale_color_manual(                          # legend styling
    name = NULL,
    values = c(
      "Spatial weighting estimate"   = "steelblue",
      "Rolling average"       = "firebrick", #  (k = 5)
      "OLS estimate"                  = "darkgreen",
      "Reference (RR = 1)"            = "red"
    )
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c(
      "Spatial weighting estimate"   = "solid",
      "Rolling average"       = "dashed", #  (k = 5)
      "OLS estimate"                  = "solid",
      "Reference (RR = 1)"            = "dotted"
    )
  ) +
  labs(
    x     = "Number of Hidden Covariates",
    y     = "Causal Risk Ratio Estimate",
    title = "Spatial Weighting Estimates of the Causal Risk Ratio with 95% CIs"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"      # <- turn the legend on
  )
dev.off()

# PRINT ESTIMATES

ests_t      <- t(ests)
cis_lower_t <- t(cis_lower)
cis_upper_t <- t(cis_upper)

# Subset only the outcomes you want to include
keep <- c("Y_lbw", "Y_avgbw", "Y_infmort", "Y_preterm", "Y_avggest")

# Format each cell as: est (lower, upper)
formatted <- matrix(nrow = nrow(ests_t), ncol = length(keep))
rownames(formatted) <- rownames(ests_t)
colnames(formatted) <- gsub("^Y_", "", keep)  # cleaner column names

for (i in seq_along(keep)) {
  y <- keep[i]
  formatted[, i] <- sprintf(
    "%.2f (%.2f, %.2f)",
    ests_t[, y],
    cis_lower_t[, y],
    cis_upper_t[, y]
  )
}

# Convert to data frame for xtable
formatted_df <- as.data.frame(formatted)

# Print as LaTeX table
print(xtable(formatted_df, align = "lccccc", digits = 3),
      include.rownames = TRUE,
      sanitize.text.function = identity)
SW_att_ests[SW_att_ests$neigen == 9,-1]
SW_att_ci_lowers[SW_att_ci_lowers$neigen == 9,-1]
SW_att_ci_uppers[SW_att_ci_uppers$neigen == 9,-1]

# FINAL SW estimate: use neigen = 9, seems to be the point of stability
SW_att_fit <- fit_method(
  Y = buffers_outcome$Y_lbw, # outcome doesn't matter - just looking at weights
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
  neigen = 9
)

w_df <- data.frame(
  weight = SW_att_fit$weights,
  Z      = factor(buffers_outcome$Z,
                  levels = c(0, 1),
                  labels = c("Control (Z = 0)", "Treated (Z = 1)"))
)

# 2.  Effective sample size (overall) ------------------------------------
ess <- compute_ess(SW_att_fit$weights)          # a single number
title_txt <- sprintf("ESS = %.1f", ess) # Histogram of spatial weights  (

# 3.  Histogram coloured by treatment group ------------------------------
png('images/spatial_weights_histogram.png',
    width = 1200,
    height = 1500,
    res = 350)
ggplot(w_df, aes(x = weight, fill = Z)) +
  geom_histogram(alpha = 0.65,
                 position = "identity",
                 bins = 30,
                 colour = "grey30") +
  scale_fill_manual(values = c("orange", "dodgerblue")) +
  labs(x = "Spatial weight",
       y = "Count",
       fill = NULL,
       title = title_txt) +
  theme_minimal() +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5))
dev.off()

bal <- balance_table(
  w = SW_att_fit$weights,
  X = cbind(
    st_drop_geometry(buffers_outcome[
      c("population_density", "percent_hispanic", "percent_black",
        "percent_indigenous", "percent_asian",
        "median_household_income", "median_house_value",
        "percent_poverty", "percent_high_school_grad",
        "median_year_built",
        "Sit_Scr")
    ]),
    Vcar[, 1:9],
    Vgp[, 1:9]
  ),
  Z = buffers_outcome$Z
)

bal_round <- bal %>% 
  mutate(across(where(is.numeric), ~ round(.x, 3)))

# final print
print(xtable(format(bal_round, nsmall = 2, scientific = FALSE), right = TRUE), include.rownames = F)

# NOW: map of weights! YAY! 
SWiw <- SW_att_fit$weights
buffer_centroids <- st_centroid(buffers_outcome)
buffers_merged <- cbind(buffer_centroids, SWiw)

buffers_merged_geo <- st_transform(buffers_merged, crs = 4326)

buffers_merged_geo <- buffers_merged_geo[order(buffers_merged_geo$SWiw),]
const <- unique(buffers_merged_geo$SWiw[buffers_merged_geo$Z == 1])
gSW <- ggplot() +
  scale_shape_manual(values = c("0" = 21, "1" = 24), guide  = guide_legend(order = 1)) +
  geom_sf(
    data   = states[states$STATEFP10 %in% c("34"), ],
    fill   = NA,
    colour = "darkgray",
    size   = .5
  ) +
  geom_sf(
    data  = subset(buffers_merged_geo, Z == 0),
    aes(fill = SWiw, shape = factor(Z)),
    size   = 2,
    stroke = .25,
    colour = "black"
  ) +
  scale_fill_gradient2(
    low        = "dodgerblue",
    mid        = "white",
    high       = "orange",
    midpoint   = 0,
    labels     = function(x)
      sprintf("%.3f", x),
    name       = "Control weights",
    guide      = guide_colorbar(order = 2, barheight = unit(1.0, "cm"))   # <<--- taller
  ) +
  
  new_scale_fill() +
  geom_sf(
    data  = subset(buffers_merged_geo, Z == 1),
    aes(fill = SWiw, shape = factor(Z)),
    size   = 2,
    stroke = .25,
    colour = "black"
  ) +
  scale_fill_gradient2(
    low        = "dodgerblue",
    mid = "dodgerblue",
    high = "dodgerblue",
    midpoint   = const,
    limits     = c(const - 1e-4, const + 1e-4),
    breaks     = const,
    labels     = sprintf("%.3f", const),
    name       = "Treated weights",
    guide      = guide_colorbar(order = 3, barheight = unit(0.25, "cm"))  # <<--- shorter
  ) +
  labs(shape = "Treatment") + # title = "Spatial Weights", 
  theme_minimal() +
  theme(
    panel.grid   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    axis.title   = element_blank(),
    plot.title   = element_text(hjust = .5),
    legend.position = "right"
  )
png('images/spatial_weights_map.png',
    width = 1000,
    height = 1000,
    res = 250)
print(gSW)
dev.off()


################################### PLOT THE HIDDEN COVARIATES ####################################
gs <- list()
for (i in 1:9){
  vcari <- Vcar[, i]
  buffers_merged <- cbind(buffer_centroids, vcari)
  buffers_merged_geo <- st_transform(buffers_merged, crs = 4326)
  buffers_merged_geo <- buffers_merged_geo[order(buffers_merged_geo$vcari),]
  
  gs[[i]] <- ggplot() +
    scale_shape_manual(values = c("0" = 21, "1" = 24), guide  = guide_legend(order = 1)) +
    geom_sf(
      data   = states[states$STATEFP10 %in% c("34"), ],
      fill   = NA,
      colour = "darkgray",
      size   = .5
    ) +
    geom_sf(
      data  = buffers_merged_geo,
      aes(fill = vcari, shape = factor(Z)),
      size   = 2,
      stroke = 0,
      colour = "black"
    ) +
    # virids
    scale_fill_viridis_c(
      labels     = function(x)
        sprintf("%.3f", x),
      name       = paste0("Vcar", i),
      guide      = guide_colorbar(order = 2, barheight = unit(1.0, "cm"))   # <<--- taller
    ) +
    theme_minimal() +
    # title it with latex \bm{v}_i^{CAR}
    labs(title = bquote(bold(V)[.(i)]^{CAR})) +
    theme(
      panel.grid   = element_blank(),
      axis.text    = element_blank(),
      axis.ticks   = element_blank(),
      axis.title   = element_blank(),
      plot.title   = element_text(hjust = .5),
      legend.position = "none"
    )
}

for (i in 1:9){
  vgpi <- Vgp[, i]
  buffers_merged <- cbind(buffer_centroids, vgpi)
  buffers_merged_geo <- st_transform(buffers_merged, crs = 4326)
  buffers_merged_geo <- buffers_merged_geo[order(buffers_merged_geo$vgpi),]
  
  gs[[9+i]] <- ggplot() +
    scale_shape_manual(values = c("0" = 21, "1" = 24), guide  = guide_legend(order = 1)) +
    geom_sf(
      data   = states[states$STATEFP10 %in% c("34"), ],
      fill   = NA,
      colour = "darkgray",
      size   = .5
    ) +
    geom_sf(
      data  = buffers_merged_geo,
      aes(fill = vgpi, shape = factor(Z)),
      size   = 2,
      stroke = 0,
      colour = "black"
    ) +
    # virids
    scale_fill_viridis_c(
      labels     = function(x)
        sprintf("%.3f", x),
      name       = paste0("Vgp", i),
      guide      = guide_colorbar(order = 2, barheight = unit(1.0, "cm"))   # <<--- taller
    ) +
    theme_minimal() +
    # title it with latex \bm{v}_i^{gp}
    labs(title = bquote(bold(V)[.(i)]^{GP})) +
    theme(
      panel.grid   = element_blank(),
      axis.text    = element_blank(),
      axis.ticks   = element_blank(),
      axis.title   = element_blank(),
      plot.title   = element_text(hjust = .5),
      legend.position = "none"
    )
}
png('images/hidden_covariates.png',
    width = 2000,
    height = 600,
    res = 170)
grid.arrange(grobs = gs, nrow = 2)
dev.off()

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
                            #Vre[,1:9],# %*% diag(sqrt(pmax(Ere[1:9,0))), # scale eigens by sqrt eigenvalues
                            Vcar[,1:9],# %*% diag(sqrt(pmax(Ecar[1:9],0))), 
                            Vgp[,1:9])# %*% diag(sqrt(pmax(Egp[1:9],0)))) 
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
SW_alg_att_lbw <- sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 1]*buffers_outcome$Y_lbw[buffers_outcome$Z == 1]) - 
  sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 0]*buffers_outcome$Y_lbw[buffers_outcome$Z == 0])
mean(buffers_outcome$Y_lbw[buffers_outcome$Z == 1])/(mean(buffers_outcome$Y_lbw[buffers_outcome$Z == 1]) - SW_alg_att_lbw ) # 0.94
SW_alg_att_infmort <- sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 1]*buffers_outcome$Y_infmort[buffers_outcome$Z == 1]) - 
  sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 0]*buffers_outcome$Y_infmort[buffers_outcome$Z == 0])
mean(buffers_outcome$Y_infmort[buffers_outcome$Z == 1])/(mean(buffers_outcome$Y_infmort[buffers_outcome$Z == 1]) - SW_alg_att_infmort ) # 0.94

# Sensitivity 2: Nonlinear only X, manual
t_ind <- buffers_outcome$Z
vars   <- colnames(X[,-1])
sq     <- paste0("I(", vars, "^2)", collapse = " + ")
full_f <- as.formula(paste("~ (.)^2 +", sq))
Xquad <- model.matrix(full_f, data = as.data.frame(X[,-1]))
bal_cov <- cbind.data.frame(Xquad,
                            #Vre[,1:9],
                            Vcar[,1:9], 
                            Vgp[,1:9]) 
data_frame <- as.data.frame(cbind(t_ind, bal_cov))
t_ind <- "t_ind"
bal <- list()
bal$bal_cov <- colnames(bal_cov)[-1]
bal$bal_std <- 'manual'
bal$bal_alg <- F
bal$bal_tol <- c(rep(0, ncol(X)-1), # exact balance on linear terms
                 rep(0.5,ncol(bal_cov)-9*2-ncol(X)), # 0.5 calipers for quadratic terms
                 #tols$RE[1:9], 
                 tols$CAR[1:9], tols$GP[1:9]) # close to exact balance for eigen
sbwatttun_object = sbw(dat = data_frame, ind = t_ind, bal = bal, 
                       sol = list(sol_nam = "quadprog"), 
                       par = list(par_est = "att", par_tar = NULL))
SW_nonlin_att_lbw <- sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 1]*buffers_outcome$Y_lbw[buffers_outcome$Z == 1]) -
  sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 0]*buffers_outcome$Y_lbw[buffers_outcome$Z == 0])
# Print causal risk ratio estimate
mean(buffers_outcome$Y_lbw[buffers_outcome$Z == 1])/(mean(buffers_outcome$Y_lbw[buffers_outcome$Z == 1]) - SW_nonlin_att_lbw) # 0.94
SW_nonlin_att_infmort <- sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 1]*buffers_outcome$Y_infmort[buffers_outcome$Z == 1]) - 
  sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 0]*buffers_outcome$Y_infmort[buffers_outcome$Z == 0])
mean(buffers_outcome$Y_infmort[buffers_outcome$Z == 1])/(mean(buffers_outcome$Y_infmort[buffers_outcome$Z == 1]) - SW_nonlin_att_infmort) # 0.94

# Sensitivity 3: Nonlinear only X, algorithm
t_ind <- buffers_outcome$Z
vars   <- colnames(X[,-1])
sq     <- paste0("I(", vars, "^2)", collapse = " + ")
full_f <- as.formula(paste("~ (.)^2 +", sq))
Xquad <- model.matrix(full_f, data = as.data.frame(X[,-1]))
bal_cov <- cbind.data.frame(Xquad,
                            #Vre[,1:9],
                            Vcar[,1:9], 
                            Vgp[,1:9]) 
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
SW_nonlin_att_lbw <- sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 1]*100*buffers_outcome$Y_lbw[buffers_outcome$Z == 1]) - 
  sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 0]*100*buffers_outcome$Y_lbw[buffers_outcome$Z == 0])
# Print causal risk ratio estimate
mean(100*buffers_outcome$Y_lbw[buffers_outcome$Z == 1])/(mean(100*buffers_outcome$Y_lbw[buffers_outcome$Z == 1]) - SW_nonlin_att_lbw) # 0.94 
SW_nonlin_att_infmort <- sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 1]*buffers_outcome$Y_infmort[buffers_outcome$Z == 1]) - 
  sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 0]*buffers_outcome$Y_infmort[buffers_outcome$Z == 0])
mean(buffers_outcome$Y_infmort[buffers_outcome$Z == 1])/(mean(buffers_outcome$Y_infmort[buffers_outcome$Z == 1]) - SW_nonlin_att_infmort) # 0.94

# Sensitivity 4: Nonlinear X and eigen, manual
t_ind <- buffers_outcome$Z
bal_cov <- cbind.data.frame(X,
                            #Vre[,1:9],
                            Vcar[,1:9], 
                            Vgp[,1:9]) 
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
bal$bal_tol <- c(rep(0, ncol(X)-1), # exact balance on linear terms
                 #tols$RE[1:9], 
                 tols$CAR[1:9], tols$GP[1:9], # close to exact balance for eigen
                 rep(0.5,ncol(bal_cov)-9*2-ncol(X))) # 0.5 calipers for all quadratic terms
sbwatttun_object = sbw(dat = data_frame, ind = t_ind, bal = bal, 
                       sol = list(sol_nam = "quadprog"), 
                       par = list(par_est = "att", par_tar = NULL))
SW_nonlin_att_lbw <- sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 1]*buffers_outcome$Y_lbw[buffers_outcome$Z == 1]) -
  sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 0]*buffers_outcome$Y_lbw[buffers_outcome$Z == 0])
# Print causal risk ratio estimate
mean(buffers_outcome$Y_lbw[buffers_outcome$Z == 1])/(mean(buffers_outcome$Y_lbw[buffers_outcome$Z == 1]) - SW_nonlin_att_lbw) # 0.94
SW_nonlin_att_infmort <- sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 1]*buffers_outcome$Y_infmort[buffers_outcome$Z == 1]) - 
  sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 0]*buffers_outcome$Y_infmort[buffers_outcome$Z == 0])
mean(buffers_outcome$Y_infmort[buffers_outcome$Z == 1])/(mean(buffers_outcome$Y_infmort[buffers_outcome$Z == 1]) - SW_nonlin_att_infmort) # 0.94

# Sensitivity 5: Nonlinear X and eigen, algorithm
t_ind <- buffers_outcome$Z
bal_cov <- cbind.data.frame(X,
                            #Vre[,1:9],
                            Vcar[,1:9], 
                            Vgp[,1:9]) 
colnames(bal_cov) <- paste0('X', 1:ncol(bal_cov))
vars   <- names(bal_cov)
sq     <- paste0("I(", vars, "^2)", collapse = " + ")
full_f <- as.formula(paste("~ (.)^2 +", sq))
bal_cov <- model.matrix(full_f, data = bal_cov)
data_frame <- as.data.frame(cbind(t_ind, bal_cov))
t_ind <- "t_ind"
bal <- list()
bal$bal_gri <-  c(0.56,0.75) # grid of tuning parameters, 0.2 was smallest of defaults that produced a solution
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


############################ Subgroup analysis stratifying on site score ############################
# Top 1/2 of sites
qmid <- quantile(buffers_outcome$Sit_Scr, 0.5)
buffers_outcome_high <- subset(buffers_outcome, Sit_Scr>qmid)
# Top 1/3 is site scores > 50
n <- nrow(buffers_outcome_high)
X <- buffers_outcome_high[c('population_density',
                       'percent_hispanic',
                       'percent_black',
                       'percent_indigenous',
                       'percent_asian',
                       'median_household_income',
                       'median_house_value',
                       'percent_poverty',
                       'percent_high_school_grad',
                       'median_year_built')]
X <- cbind.data.frame(Intercept = 1, X) # Add intercept

X <- st_drop_geometry(X)
X[,2:11] <- scale(X[,2:11])
X <- X[,-12]
X <- as.matrix(X)

dmat <- distm(cbind(buffers_outcome_high$Longitd, 
                    buffers_outcome_high$Latitud), fun = distHaversine)
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
buffers_outcome_high$cluster[buffers_outcome_high$cluster == 9] = 40
buffers_outcome_high$cluster[buffers_outcome_high$cluster == 36] = 40
clusters <- buffers_outcome_high$cluster

rho2 <- 10
kappa <- 10 
rangec <- 500
phic <- rangec/(2*sqrt(kappa))
S <- geoR::matern(u =  dmat, phi = phic, kappa = kappa)
E <- eigen(S)
Egp <- E$values
Vgp <- E$vectors
Sigma <- diag(n) + rho2*S
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
Sigma <- diag(n) + rho2*S
Sigmainvre <- solve(Sigma)
# add CAR eigens to V
L <- diag(rowSums(adjacency_matrix)) - adjacency_matrix
S <- ginv(as.matrix(L))
E <- eigen(S)
Ecar <- Re(E$values)
Vcar <- Re(E$vectors)
Sigma <- diag(n) + rho2*S
Sigmainvcar <- solve(Sigma)

tols <- list('RE' = rep(0.01,n), 
             'CAR' = rep(0.01,n), 
             'GP'  = rep(0.01,n))
SW_att_fit <- fit_method(Y = 100*buffers_outcome_high$Y,
                         X,
                         Z = buffers_outcome_high$Z,
                         Vre = Vre,
                         Vcar = Vcar,
                         Vgp = Vgp,
                         Sigmainvre = Sigmainvre,
                         Sigmainvcar = Sigmainvcar,
                         Sigmainvgp = Sigmainvgp,
                         tols = tols,
                         method = 'SW',
                         neigen = 5)
SW_att <- SW_att_fit$est
mean(100*buffers_outcome_high$Y[buffers_outcome_high$Z == 1])/(mean(100*buffers_outcome_high$Y[buffers_outcome_high$Z == 1]) - SW_att)

# LOWER
qmid <- quantile(buffers_outcome$Sit_Scr, 0.5)
buffers_outcome_low <- subset(buffers_outcome, Sit_Scr<=qmid)
# Top 1/3 is site scores > 50
n <- nrow(buffers_outcome_low)
X <- buffers_outcome_low[c('population_density',
                            'percent_hispanic',
                            'percent_black',
                            'percent_indigenous',
                            'percent_asian',
                            'median_household_income',
                            'median_house_value',
                            'percent_poverty',
                            'percent_high_school_grad',
                            'median_year_built')]
X <- cbind.data.frame(Intercept = 1, X) # Add intercept

X <- st_drop_geometry(X)
X[,2:11] <- scale(X[,2:11])
X <- X[,-12]
X <- as.matrix(X)

dmat <- distm(cbind(buffers_outcome_low$Longitd, 
                    buffers_outcome_low$Latitud), fun = distHaversine)
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
buffers_outcome_low$cluster[buffers_outcome_low$cluster == 9] = 40
buffers_outcome_low$cluster[buffers_outcome_low$cluster == 36] = 40
clusters <- buffers_outcome_low$cluster

rho2 <- 10
kappa <- 10 
rangec <- 500
phic <- rangec/(2*sqrt(kappa))
S <- geoR::matern(u =  dmat, phi = phic, kappa = kappa)
E <- eigen(S)
Egp <- E$values
Vgp <- E$vectors
Sigma <- diag(n) + rho2*S
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
Sigma <- diag(n) + rho2*S
Sigmainvre <- solve(Sigma)
# add CAR eigens to V
L <- diag(rowSums(adjacency_matrix)) - adjacency_matrix
S <- ginv(as.matrix(L))
E <- eigen(S)
Ecar <- Re(E$values)
Vcar <- Re(E$vectors)
Sigma <- diag(n) + rho2*S
Sigmainvcar <- solve(Sigma)

tols <- list('RE' = rep(0.01,n), 
             'CAR' = rep(0.01,n), 
             'GP'  = rep(0.01,n))
SW_att_fit <- fit_method(Y = 100*buffers_outcome_low$Y,
                         X,
                         Z = buffers_outcome_low$Z,
                         Vre = Vre,
                         Vcar = Vcar,
                         Vgp = Vgp,
                         Sigmainvre = Sigmainvre,
                         Sigmainvcar = Sigmainvcar,
                         Sigmainvgp = Sigmainvgp,
                         tols = tols,
                         method = 'SW',
                         neigen = 5)
SW_att <- SW_att_fit$est
mean(100*buffers_outcome_low$Y[buffers_outcome_low$Z == 1])/(mean(100*buffers_outcome_low$Y[buffers_outcome_low$Z == 1]) - SW_att)

# This suggestions that the reduction in low birth weight rate is actually... lower among the 
# more dangerous, contaminated sites. So cleanup does not as effectively reduce low birth weight compared to 
# the less contaminated sites.. this could be a result of some other factor that is not accounted by the model

buffers_outcome_nj <- subset(buffers_outcome, cluster == '32')
# Top 1/3 is site scores > 50
n <- nrow(buffers_outcome_nj)
X <- buffers_outcome_nj[c('population_density',
                            'percent_hispanic',
                            'percent_black',
                            'percent_indigenous',
                            'percent_asian',
                            'median_household_income',
                            'median_house_value',
                            'percent_poverty',
                            'percent_high_school_grad',
                            'median_year_built')]
X <- cbind.data.frame(Intercept = 1, X) # Add intercept

X <- st_drop_geometry(X)
X[,2:11] <- scale(X[,2:11])
X <- X[,-12]
X <- as.matrix(X)

dmat <- distm(cbind(buffers_outcome_nj$Longitd, 
                    buffers_outcome_nj$Latitud), fun = distHaversine)
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
buffers_outcome_nj$cluster[buffers_outcome_nj$cluster == 9] = 40
buffers_outcome_nj$cluster[buffers_outcome_nj$cluster == 36] = 40
clusters <- buffers_outcome_nj$cluster

rho2 <- 10
kappa <- 10 
rangec <- 500
phic <- rangec/(2*sqrt(kappa))
S <- geoR::matern(u =  dmat, phi = phic, kappa = kappa)
E <- eigen(S)
Egp <- E$values
Vgp <- E$vectors
Sigma <- diag(n) + rho2*S
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
Sigma <- diag(n) + rho2*S
Sigmainvre <- solve(Sigma)
# add CAR eigens to V
L <- diag(rowSums(adjacency_matrix)) - adjacency_matrix
S <- ginv(as.matrix(L))
E <- eigen(S)
Ecar <- Re(E$values)
Vcar <- Re(E$vectors)
Sigma <- diag(n) + rho2*S
Sigmainvcar <- solve(Sigma)

tols <- list('RE' = rep(0.01,n), 
             'CAR' = rep(0.01,n), 
             'GP'  = rep(0.01,n))
SW_att_fit <- fit_method(Y = 100*buffers_outcome_nj$Y,
                         X,
                         Z = buffers_outcome_nj$Z,
                         Vre = Vre,
                         Vcar = Vcar,
                         Vgp = Vgp,
                         Sigmainvre = Sigmainvre,
                         Sigmainvcar = Sigmainvcar,
                         Sigmainvgp = Sigmainvgp,
                         tols = tols,
                         method = 'SW',
                         neigen = 5)
SW_att <- SW_att_fit$est
mean(100*buffers_outcome_nj$Y[buffers_outcome_nj$Z == 1])/(mean(100*buffers_outcome_nj$Y[buffers_outcome_nj$Z == 1]) - SW_att)


SW_att_boot <- boot_func(Y = 100*buffers_outcome_nj$Y,
                         X,
                         Z = buffers_outcome_nj$Z,
                         Vre = Vre,
                         Vcar = Vcar,
                         Vgp = Vgp,
                         Sigmainvre = Sigmainvre,
                         Sigmainvcar = Sigmainvcar,
                         Sigmainvgp = Sigmainvgp,
                         tols = tols,
                         method = 'SW',
                         neigen = 5)
SW_lower <- SW_att - 1.96*sd(SW_att_boot)
SW_upper <- SW_att + 1.96*sd(SW_att_boot)
print(paste('SW_att_boot = ', round(SW_att,3), 
            'CI = [', round(SW_lower,3), ', ', round(SW_upper,3), ']'))
mean(100*buffers_outcome_nj$Y[buffers_outcome_nj$Z == 1])/(mean(100*buffers_outcome_nj$Y[buffers_outcome_nj$Z == 1]) - SW_lower)
mean(100*buffers_outcome_nj$Y[buffers_outcome_nj$Z == 1])/(mean(100*buffers_outcome_nj$Y[buffers_outcome_nj$Z == 1]) - SW_att)
mean(100*buffers_outcome_nj$Y[buffers_outcome_nj$Z == 1])/(mean(100*buffers_outcome_nj$Y[buffers_outcome_nj$Z == 1]) - SW_upper)

