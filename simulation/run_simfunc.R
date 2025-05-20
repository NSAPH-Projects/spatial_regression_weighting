# Load in all utility functions
source('../funcs.R')
library(sbw)
library(sf)
library(MASS)

# Load in simulation data,
# a list of... 
load('sim.RData')

set.seed(123)
args <- commandArgs(trailingOnly = TRUE)
nsims <- as.integer(args[1])
smoothing <- args[2]
outcomemod <- args[3]

print(c(nsims, smoothing, outcomemod))

simfunc(nsims = nsims,
        X = simlist$X,
        Z = simlist$Z,
        Vre = simlist$Vre,
        Ere = simlist$Ere,
        Vcar = simlist$Vcar,
        Ecar = simlist$Ecar,
        Vgp = simlist$Vgp,
        Egp = simlist$Egp,
        Sigmainvre = simlist$Sigmainvre,
        Sigmainvcar = simlist$Sigmainvcar,
        Sigmainvgp = simlist$Sigmainvgp,
        Wre = simlist$Wre,
        Wcar = simlist$Wcar,
        Wgp = simlist$Wgp,
        smoothing = smoothing,
        outcomemod = outcomemod
        )


  


