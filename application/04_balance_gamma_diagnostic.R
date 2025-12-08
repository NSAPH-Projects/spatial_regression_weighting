library(tidyr)
library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(Matrix)
library(patchwork)
library(rlang)
library(cowplot)
library(nloptr)

source('../funcs.R')
set.seed(111)

states <- st_read('../data/tl_2010_us_state10/tl_2010_us_state10.shp')

# Read in preprocessed data
load('../data/preprocessed_superfunds.RData')
us_states <- ne_states(country = "United States of America", returnclass = "sf")

X[,2:ncol(X)] <- scale(X[,2:ncol(X)])
X <- as.matrix(X)
n <- nrow(X)

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
sum(REiw[Z == 1]);sum(REiw[Z == 0])
l_vec <- REiw*(2*Z-1) # known vector l

# Compute eigen-decomposition to get eigenvalues and eigenvectors of S
eig_S <- eigen(S)
lambda_min <- min(eig_S$values)
lambda_max <- max(eig_S$values)
v_min <- eig_S$vectors[, which.min(eig_S$values)]
v_max <- eig_S$vectors[, which.max(eig_S$values)]

# Choose a value k in the feasible range, i.e. between lambda_min and lambda_max.
cat("Feasible range for k is roughly from", lambda_min, "to", lambda_max, "\n")

k_vals <- seq(lambda_min+0.1, lambda_max, length.out = 20)
# Remove the last k_vals
k_vals <- k_vals[k_vals < lambda_max]

dfre <- data.frame(
  k_vals = k_vals,
  moran = k_vals / lambda_max,
  abs_imbalance = rep(NA, length(k_vals))
)
# Loop over each k value.

k_vals <- sort(dfre$k_vals)
start_time <- Sys.time()

for (i in seq_along(k_vals)) {
  k_val <- dfre$k_vals[i]
  
  # Construct a feasible initial guess U0.
  alpha <- sqrt((lambda_max - k_val) / (lambda_max - lambda_min))
  beta  <- sqrt((k_val - lambda_min) / (lambda_max - lambda_min))
  U0 <- alpha * v_min + beta * v_max + rnorm(n, sd = 0.1)
  U0 <- U0 - mean(U0)
  U0 <- U0/norm(U0, type = '2')
  
  # Diagnostic: check that U0 satisfies the constraints.
  cat("Initial guess check for k =", k_val, ":\n")
  cat("  U0^T U0 =", sum(U0^2), "\n")
  cat("  U0^T S U0 =", as.numeric(t(U0) %*% S %*% U0), "\n")
  
  eval_f <- function(U) {
    f <- - (sum(l_vec * U))^2
    return(f)
  }
  
  eval_grad_f <- function(U) {
    grad <- -2 * sum(l_vec * U) * l_vec
    return(grad)
  }
  
  # Define the equality constraints and their Jacobian.
  eval_g_eq <- function(U) {
    constr1 <- sum(U^2) - 1
    constr2 <- as.numeric(t(U) %*% S %*% U) - k_val
    return(c(constr1, constr2))
  }
  
  eval_jac_g_eq <- function(U) {
    n <- length(U)
    jac1 <- as.vector(2 * U)
    jac2 <- as.vector(2 * (S %*% U))
    matrix(c(jac1, jac2), nrow = 2, byrow = TRUE)
  }
  
  opts <- list("algorithm" = "NLOPT_LD_SLSQP",
               "xtol_rel" = 1.0e-5,
               "maxeval"  = 200)
  
  result <- nloptr(x0 = U0,
                   eval_f = eval_f,
                   eval_grad_f = eval_grad_f,
                   lb = rep(-Inf, length(U0)),
                   ub = rep(Inf, length(U0)),
                   eval_g_eq = eval_g_eq,
                   eval_jac_g_eq = eval_jac_g_eq,
                   opts = opts)
  U_opt <- result$solution
  opt_val <- (sum(l_vec * U_opt))^2
  dfre$abs_imbalance[i] <- sqrt(opt_val)
  
  # Print diagnostic information.
  cat("\nResults for k =", k_val, "\n")
  cat("Optimal U:\n")
  cat("Optimal value (l^T U)^2 =", opt_val, "\n")
  cat("Constraint check:\n")
  cat("  U^T U =", sum(U_opt^2), "\n")
  cat("  U^T S U =", as.numeric(t(U_opt) %*% S %*% U_opt), " (should equal k =", k_val, ")\n\n")
}
print(Sys.time() - start_time)
save(dfre, file = 'misc_data/dfre.RData')

# Create a data frame with absgamma values
gamma_df <- data.frame(abs_gamma = seq(0.01, 1, length.out = 20))

# Perform a cross join to expand the data frame
df_expanded <- merge(dfre, gamma_df, by = NULL)
df_expanded$abs_cond_bias <- df_expanded$abs_gamma*df_expanded$abs_imbalance

plot1 <- ggplot(df_expanded, aes(x = moran, y = abs_gamma, fill = abs_cond_bias)) +
  geom_tile() +
  labs(
    title = "Random Effects Model",
    x = "Moran's I",
    y = expression("|" * gamma * "|"),
    fill = expression(max ~ abs(E(hat(tau)[GLS] ~ "|" ~ (bold(X) * ", " * bold(Z) * ", " * bold(U))) - tau))
  ) +
  scale_fill_viridis_c() +
  theme_minimal() + # increase legend key width
  theme(
    legend.key.width = unit(1.5, "cm")
  )

# CAR
sigma2 <- 1
rho2 <- 10
L <- diag(rowSums(adjacency_matrix)) - adjacency_matrix
S <- ginv(as.matrix(L))
S[abs(S) < 1e-10] <- 0
Sigma <- sigma2*diag(n) + rho2*S
Sigmainv <- solve(Sigma)
CARiw <- impliedweightsgeneral(X = X,
                               Z = Z,
                               Sigmainv = Sigmainv)
sum(CARiw[Z == 1]);sum(CARiw[Z == 0])
l_vec <- CARiw*(2*Z-1)

# Compute eigen-decomposition to get eigenvalues and eigenvectors of S
eig_S <- eigen(S)
eig_S$values <- Re(eig_S$values)
eig_S$vectors <- Re(eig_S$vectors)
# Make eig_S$vectors sparse (lots of numbers < 10^(-15))
eig_S$vectors[abs(eig_S$vectors) < 1e-10] <- 0
lambda_min <- min(eig_S$values)
lambda_max <- max(eig_S$values)
v_min <- eig_S$vectors[, which.min(eig_S$values)]
v_max <- eig_S$vectors[, which.max(eig_S$values)]

# Choose a value k in the feasible range, i.e. between lambda_min and lambda_max.
k_vals <- seq(lambda_min+0.1, lambda_max, length.out = 20)
# Remove the last k_vals
k_vals <- k_vals[k_vals < lambda_max]

dfcar <- data.frame(
  k_vals = k_vals,
  moran = k_vals / lambda_max,
  abs_imbalance = rep(NA, length(k_vals))
)
# Loop over each k value.
for (i in seq_along(k_vals)) {
  k_val <- dfcar$k_vals[i]
  
  # Construct a feasible initial guess U0.
  alpha <- sqrt((lambda_max - k_val) / (lambda_max - lambda_min))
  beta  <- sqrt((k_val - lambda_min) / (lambda_max - lambda_min))
  U0 <- alpha * v_min + beta * v_max + rnorm(n, sd = 0.1)
  U0 <- U0 - mean(U0)
  U0 <- U0/norm(U0, type = '2')
  
  # Diagnostic: check that U0 satisfies the constraints.
  cat("Initial guess check for k =", k_val, ":\n")
  cat("  U0^T U0 =", sum(U0^2), "\n")
  cat("  U0^T S U0 =", as.numeric(t(U0) %*% S %*% U0), "\n")
  
  eval_f <- function(U) {
    f <- - (sum(l_vec * U))^2
    return(f)
  }
  
  eval_grad_f <- function(U) {
    grad <- -2 * sum(l_vec * U) * l_vec
    return(grad)
  }
  
  # Define the equality constraints and their Jacobian.
  eval_g_eq <- function(U) {
    constr1 <- sum(U^2) - 1
    constr2 <- as.numeric(t(U) %*% S %*% U) - k_val
    return(c(constr1, constr2))
  }
  
  eval_jac_g_eq <- function(U) {
    n <- length(U)
    jac1 <- as.vector(2 * U)
    jac2 <- as.vector(2 * (S %*% U))
    matrix(c(jac1, jac2), nrow = 2, byrow = TRUE)
  }
  
  opts <- list("algorithm" = "NLOPT_LD_SLSQP",
               "xtol_rel" = 1.0e-5,
               "maxeval"  = 200)
  
  result <- nloptr(x0 = U0,
                   eval_f = eval_f,
                   eval_grad_f = eval_grad_f,
                   lb = rep(-Inf, length(U0)),
                   ub = rep(Inf, length(U0)),
                   eval_g_eq = eval_g_eq,
                   eval_jac_g_eq = eval_jac_g_eq,
                   opts = opts)
  U_opt <- result$solution
  opt_val <- (sum(l_vec * U_opt))^2
  dfcar$abs_imbalance[i] <- sqrt(opt_val)
  
  # Print diagnostic information.
  cat("\nResults for k =", k_val, "\n")
  cat("Optimal U:\n")
  cat("Optimal value (l^T U)^2 =", opt_val, "\n")
  cat("Constraint check:\n")
  cat("  U^T U =", sum(U_opt^2), "\n")
  cat("  U^T S U =", as.numeric(t(U_opt) %*% S %*% U_opt), " (should equal k =", k_val, ")\n\n")
}
save(dfcar, file = 'misc_data/dfcar.RData')

gamma_df <- data.frame(abs_gamma = seq(0.01, 1, length.out = 20))

# Perform a cross join to expand the data frame
df_expanded <- merge(dfcar, gamma_df, by = NULL)
df_expanded$abs_cond_bias <- df_expanded$abs_gamma*df_expanded$abs_imbalance

plot2 <- ggplot(df_expanded, aes(x = moran, y = abs_gamma, fill = abs_cond_bias)) +
  geom_tile() +
  labs(
    title = "Intrinsic Conditional Autoregressive Model",
    x = "Moran's I",
    y = expression("|" * gamma * "|"),
    fill = expression(max ~ abs(E(hat(tau)[GLS] ~ "|" ~ (bold(X) * ", " * bold(Z) * ", " * bold(U))) - tau))
  ) +
  scale_fill_viridis_c() +
  theme_minimal() + # increase legend key width
  theme(
    legend.key.width = unit(1.5, "cm")
  )

# GP
kappa <- 10 
rangec <- 500
phic <- rangec/(2*sqrt(kappa))
S <- geoR::matern(u =  dmat, phi = phic, kappa = kappa)
V <- eigen(S)$vectors
Sigma <- diag(n) + 10*S
Sigmainv <- solve(Sigma)
GPiw <- impliedweightsgeneral(X = X,
                              Z = Z,
                              Sigmainv = Sigmainv)
sum(GPiw[Z == 1]);sum(GPiw[Z == 0])
l_vec <- GPiw*(2*Z-1) # known vector l
# Compute the eigen-decomposition of S.
eig_S <- eigen(S)
lambda_min <- min(eig_S$values)
lambda_max <- max(eig_S$values)
v_min <- eig_S$vectors[, which.min(eig_S$values)]
v_max <- eig_S$vectors[, which.max(eig_S$values)]

cat("Feasible range for k is roughly from", lambda_min, "to", lambda_max, "\n")

k_vals <- seq(lambda_min+0.1, lambda_max, length.out = 20)
# Remove the last k_vals
k_vals <- k_vals[k_vals < lambda_max]

dfgp <- data.frame(
  k_vals = k_vals,
  moran = k_vals / lambda_max,
  abs_imbalance = rep(NA, length(k_vals))
)
# Loop over each k value.
start_time <- Sys.time()
for (i in seq_along(dfgp$k_vals)) {
  k_val <- dfgp$k_vals[i]
  
  # Construct a feasible initial guess U0.
  alpha <- sqrt((lambda_max - k_val) / (lambda_max - lambda_min))
  beta  <- sqrt((k_val - lambda_min) / (lambda_max - lambda_min))
  U0 <- alpha * v_min + beta * v_max + rnorm(n, sd = 0.1)
  U0 <- U0 - mean(U0)
  U0 <- U0/norm(U0, type = '2')
  
  # Diagnostic: check that U0 satisfies the constraints.
  cat("Initial guess check for k =", k_val, ":\n")
  cat("  U0^T U0 =", sum(U0^2), "\n")
  cat("  U0^T S U0 =", as.numeric(t(U0) %*% S %*% U0), "\n")
  
  eval_f <- function(U) {
    f <- - (sum(l_vec * U))^2
    return(f)
  }
  
  eval_grad_f <- function(U) {
    grad <- -2 * sum(l_vec * U) * l_vec
    return(grad)
  }
  
  # Define the equality constraints and their Jacobian.
  eval_g_eq <- function(U) {
    constr1 <- sum(U^2) - 1
    constr2 <- as.numeric(t(U) %*% S %*% U) - k_val
    return(c(constr1, constr2))
  }
  
  eval_jac_g_eq <- function(U) {
    n <- length(U)
    jac1 <- as.vector(2 * U)
    jac2 <- as.vector(2 * (S %*% U))
    matrix(c(jac1, jac2), nrow = 2, byrow = TRUE)
  }
  
  opts <- list("algorithm" = "NLOPT_LD_SLSQP",
               "xtol_rel" = 1.0e-5,
               "maxeval"  = 200)
  
  result <- nloptr(x0 = U0,
                   eval_f = eval_f,
                   eval_grad_f = eval_grad_f,
                   lb = rep(-Inf, length(U0)),
                   ub = rep(Inf, length(U0)),
                   eval_g_eq = eval_g_eq,
                   eval_jac_g_eq = eval_jac_g_eq,
                   opts = opts)
  U_opt <- result$solution
  opt_val <- (sum(l_vec * U_opt))^2
  dfgp$abs_imbalance[i] <- sqrt(opt_val)
  
  # Print diagnostic information.
  cat("\nResults for k =", k_val, "\n")
  cat("Optimal U:\n")
  cat("Optimal value (l^T U)^2 =", opt_val, "\n")
  cat("Constraint check:\n")
  cat("  U^T U =", sum(U_opt^2), "\n")
  cat("  U^T S U =", as.numeric(t(U_opt) %*% S %*% U_opt), " (should equal k =", k_val, ")\n\n")
}
print(Sys.time() - start_time)
save(dfgp, file = 'misc_data/dfgp.RData')

gamma_df <- data.frame(abs_gamma = seq(0.01, 1, length.out = 20))

# Perform a cross join to expand the data frame
df_expanded <- merge(dfgp, gamma_df, by = NULL)
df_expanded$abs_cond_bias <- df_expanded$abs_gamma*df_expanded$abs_imbalance

plot3 <- ggplot(df_expanded, aes(x = moran, y = abs_gamma, fill = abs_cond_bias)) +
  geom_tile() +
  labs(
    title = "Gaussian Process Model",
    x = "Moran's I",
    y = expression("|" * gamma * "|"),
    fill = expression(max ~ abs(E(hat(tau)[GLS] ~ "|" ~ (bold(X) * ", " * bold(Z) * ", " * bold(U))) - tau))
  ) +
  scale_fill_viridis_c() +
  theme_minimal() + # increase legend key width
  theme(
    legend.key.width = unit(1.5, "cm")
  )

combined_plot <- plot1 + plot2 + plot3 +
  plot_layout(guides = 'collect') &
  theme(legend.position = "bottom")
# Save
png('images/abs_cond_bias.png', width = 2500, 
    height = 1100, res = 200)
combined_plot
dev.off()
