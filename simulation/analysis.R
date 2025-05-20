library(ggplot2)
library(sf)
library(tidyr)
library(dplyr)
library(utils)
library(xtable)
library(Matrix)
library(tidyverse)
library(parallel)
source('../funcs.R')
load('sim.RData')

# read in results files
csvs <- list.files('results_May1/', pattern = '.csv')

# Create storage for metrics 
analysisdf <- data.frame(
  smoothing = character(length(csvs)),
  outcomemod = character(length(csvs)),
  method = character(length(csvs)),
  bias = numeric(length(csvs)),
  RMSE = numeric(length(csvs)),
  se = numeric(length(csvs))
)

# Extract components from filenames
smoothing <- gsub("(_.*$)", "", csvs)
outcomemod <- gsub("^(.*?_)(.*?)(_.*$)", "\\2", csvs)
method <- gsub(".*_([^\\.]+)\\.csv", "\\1", csvs)

# Precompute true estimand for each outcome model and confounding mechanism
# ATTs <- data.frame(expand.grid(smoothing = c('clusterbased', 'adjacencybased', 'distancebased'),
#                                  outcomemod = c('linear', 'linearinteraction', 'nonlinearinteraction')))
# ATTs$theta <- NA
# ATTs$smoothing <- as.character(ATTs$smoothing)
# ATTs$outcomemod <- as.character(ATTs$outcomemod)
# 
# ATTs$tau <- unlist(mclapply(1:nrow(ATTs), function(i) {
#   compute_ATT_MC(X = X,
#                  Z = Z,
#                  beta = c(-0.44,0.46,-0.69,-1.45,0.57,-1.02,-0.02,-0.94,1.10,-0.48,-0.71),
#                  Wre = simlist$Wre,
#                  Wcar = simlist$Wcar,
#                  Wgp = simlist$Wgp,
#                  outcomemod = ATTs$outcomemod[i],
#                 smoothing = ATTs$smoothing[i]
#                 )
# }, mc.cores = 2))  # Adjust the number of cores
# ATTs
# save(ATTs, file = 'results_May1/ATTs.RData')

load('results_May1/ATTs.RData')
                
# Loop through results to calculate metrics and create plots.
for (i in 1:length(csvs)){
  filename <- csvs[i]
  print(filename)
  
  analysisdf$smoothing[i] <- smoothing[i]
  analysisdf$outcomemod[i] <- outcomemod[i]
  analysisdf$method[i] <- method[i]
  df_temp <- read.csv(file.path('results_May1/', filename))
  tauests <- df_temp 
  # Convert muests to a vector, it's just a single column
  tauests <- as.vector(as.matrix(tauests))
  
  # Compute true truncated exposure estimate
  tau <- ATTs[ATTs$smoothing == smoothing[i] & 
                      ATTs$outcomemod == outcomemod[i],]$tau
  df_temp$tau <- tau
  
  # Save metrics in analysisdf
  analysisdf$bias[i] <- mean(tauests, na.rm = T) - tau
  analysisdf$RMSE[i] <- sqrt(mean((tauests - tau)^2, na.rm = T))
  analysisdf$se[i] <- sd(tauests, na.rm = T)
}

# Format the numeric columns (bias, RMSE, se) in scientific notation with 3 decimals
analysisdf_bias <- analysisdf %>%
  mutate(
    bias = round(bias*100, 3),
    RMSE = round(RMSE*100, 3),
    se = round(se*100, 3)
  )

# Pivot wider from the original analysisdf
analysisdf_bias <- analysisdf_bias[, c(1:4)] %>% 
  pivot_wider(names_from = method, values_from = bias)
analysisdf_bias <- analysisdf_bias[, c("smoothing", "outcomemod", "OLS", "RE", "CAR", "GP", "SW"
                                       )]

# Print using xtable and prevent xtable from reformatting the already-formatted text
print(xtable(analysisdf_bias), 
      include.rownames = FALSE, sanitize.text.function = identity)

# Do the same with RMSE
analysisdf_RMSE <- analysisdf %>%
  mutate(
    bias = round(bias*100, 3),
    RMSE = round(RMSE*100, 3),
    se = round(se*100, 3)
  )

analysisdf_RMSE <- analysisdf_RMSE[, c(1:3, 5)] %>% 
  pivot_wider(names_from = method, values_from = RMSE
)
analysisdf_RMSE <- analysisdf_RMSE[, c("smoothing", "outcomemod", "OLS", "RE", "CAR", "GP", "SW"
                                       )]

print(xtable(analysisdf_RMSE), include.rownames = F, sanitize.text.function = identity)

# Create facet_wrap boxplots with ggplot2

folder <- "results_May1"

# List all CSV files in that folder (with full paths)
files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)

read_estimates <- function(i) {
  
  filename <- csvs[i]
  print(filename)
  
  smoothing <- smoothing[i]
  outcomemod <- outcomemod[i]
  method <- method[i]
  dat <- read.csv(file.path('results_May1/', filename))

  # If the CSV doesn't have a header and just one column, name it "estimate"
  if (!"estimate" %in% colnames(dat)) {
    names(dat)[1] <- "estimate"
  }
  # Add the new columns
  dat <- dat %>%
    mutate(smoothing = smoothing,
           outcomemod = outcomemod,
           method = method)
  return(dat)
}

# Read all files and combine into one data frame
df <- map_dfr(1:length(csvs), read_estimates)

desired_order <- c("OLS", 
                  "RE", 
                  "CAR", 
                  "GP", 
                  "SW"
                  )
df$method <- factor(df$method, levels = desired_order)

# Ensure smoothing and outcomemod are factors in both data frames with the same levels:
df <- df %>% 
  mutate(smoothing = factor(smoothing, levels = c('clusterbased', 'adjacencybased', 'distancebased')),
          outcomemod = factor(outcomemod, levels = c("linear", "linearinteraction", "nonlinearinteraction")))
ATTs <- ATTs %>% 
  mutate(smoothing = factor(smoothing, levels = c('clusterbased', 'adjacencybased', 'distancebased')),
         outcomemod = factor(outcomemod, levels = c("linear", "linearinteraction", "nonlinearinteraction")))

# Create a boxplot of results
df2 <- df %>% 
  mutate(
    smoothing = factor(
      smoothing,
      levels = c("clusterbased", "adjacencybased", "distancebased"),
      labels = c(
        "U: cluster-based", 
        "U: adjacency-based", 
        "U: distance-based"
      )
    ),
    outcomemod = factor(
      outcomemod,
      levels = c("linear", "linearinteraction", "nonlinearinteraction"),
      labels = c(
        "Outcome: linear",
        "Outcome: linear + interaction",
        "Outcome: non-linear + interaction"
      )
    )
  )

ATTs2 <- ATTs %>% 
  mutate(
    smoothing = factor(
      smoothing,
      levels = c("clusterbased", "adjacencybased", "distancebased"),
      labels = c(
        "U: cluster-based", 
        "U: adjacency-based", 
        "U: distance-based"
      )
    ),
    outcomemod = factor(
      outcomemod,
      levels = c("linear", "linearinteraction", "nonlinearinteraction"),
      labels = c(
        "Outcome: linear",
        "Outcome: linear + interaction",
        "Outcome: non-linear + interaction"
      )
    )
  )

png("images/boxplot.png", width = 2880, height = 2400, res = 350)
ggplot(df2, aes(x = method, y = estimate, fill = method)) +
  geom_boxplot(alpha = 0.5) +
  facet_grid(
    outcomemod ~ smoothing,
    scales = "free_y"
  ) +
  geom_hline(data = ATTs2, aes(yintercept = tau),
             color = "red", linetype = "dashed", size = 1) +
  labs(x = NULL, 
       y = "Estimate of the ATT", 
       fill = "Method") +
  theme_bw() 
dev.off()
