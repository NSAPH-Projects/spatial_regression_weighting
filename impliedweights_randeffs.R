library(Matrix)

impliedweightsfast = function(X, Xks, Zks, nks, mks, Xbarks, Xtbarks, sig2gam, sig2eps){
  # Xks is the list of nk x p matrices of covariates including intercept
  # Zks is the list of nk x 1 matrices of treatment indicators
  # nks is the length k vector of cluster sizes
  # mks is the length k vector of number of treated units in each cluster
  # Xbarks is the list of 1 x p matrices of cluster means of covariates
  # Xtbarks is the list of 1 x p matrices of cluster means of treated covariates
  # sig2gam is the variance of random effects
  # sig2eps is the variance of the error term
  # returns the implied weights for the observations
  
  nks = table(clusters)
  rhoks = sig2gam/(sig2eps + nks*sig2gam)
  avgsum = compute_treated_avgsum_fast(nks, mks, rhoks, Xbarks, Xtbarks)
  Sinv = computeSinvfast(X, nks, rhoks, Xbarks)
  # to compute denom, sum mj^2*pj over clusters
  denom = sum(Z)
  for (k in 1:K){
    denom = denom - mks[k]^2*rhoks[k]
  }
  denom = denom - avgsum%*%Sinv%*%t(avgsum) # 1x1
  weights = lapply(1:K, function(k) {
    num = Zks[[k]] - mks[k]*rhoks[k]*matrix(1, nrow = nks[k], ncol = 1) - 
      (Xks[[k]] - nks[k]*rhoks[k]*matrix(1, nrow = nks[k], ncol = 1) %*% Xbarks[[k]]) %*% Sinv %*% t(avgsum)
    return((2*Zks[[k]]-1)*num/as.numeric(denom)) # negative for untreated
  })
  return(weights)
}


impliedweights = function(X, Z, clusters, sig2gam, sig2eps){
  # X is the n x p matrix of covariates including intercept
  # Z is the n x 1 matrix of treatment indicators
  # clusters is the n x 1 matrix of cluster indicators
  # sig2gam is the variance of random effects
  # sig2eps is the variance of the error term
  # returns the implied weights for the observations
  
  avgsum = compute_treated_avgsum(X, Z, clusters, sig2gam, sig2eps)
  K = length(unique(clusters))
  n = nrow(X)
  Sinv = computeSinv(X, Z, clusters, sig2gam, sig2eps)
  # to compute denom, sum mj*pj over clusters
  denom = sum(Z)
  for (k in 1:K){
    dat = extractcluster(X, Z, clusters, sig2gam, sig2eps, k)
    mk = dat$mk
    rhok = dat$rhok
    denom = denom - mk^2*rhok
  }
  denom = denom - avgsum%*%Sinv%*%t(avgsum) # 1x1
  
  weights = list()
  for (k in 1:K){
    dat = extractcluster(X, Z, clusters, sig2gam, sig2eps, k)
    nk = dat$nk
    mk = dat$mk
    rhok = dat$rhok
    Xbark = dat$Xbar
    Xtbark = dat$Xtbar
    Zk = dat$Zk
    Xk = dat$Xk
    # compute the implied weights for cluster k
    num = Zk - mk*rhok*matrix(1, nrow = nk, ncol = 1) - 
      (Xk - nk*rhok*matrix(1, nrow = nk, ncol = 1) %*% Xbark)%*% Sinv %*% t(avgsum)
    weights[[k]] = (2*Zk-1)*num/as.numeric(denom) # negative for untreated
  }
  return(weights)
}

impliedweightsslow = function(X, Z, clusters, sig2gam, sig2eps){
  # compute number of units in each cluster efficiently
  nks = table(clusters)
  sigs = list()
  for (k in 1:length(unique(clusters))){
    sigs[[k]] = sig2eps*diag(1, nks[k]) + sig2gam*matrix(1, nks[k], nks[k])
  }
  Sigma = bdiag(sigs)
  Sinv = solve(Sigma)
  n = nrow(X)
  XSinvX <- solve(t(X) %*% Sinv %*% X)
  lnum = (diag(1,n) - Sinv %*% X %*% XSinvX %*% t(X)) %*% Sinv %*% Z
  ldenom = t(Z) %*% Sinv %*% (diag(1,n) - X %*% XSinvX %*% t(X) %*% Sinv) %*% Z
  l = lnum/as.numeric(ldenom)
  return(as.numeric((2*Z-1)*l))
}

# Create a function that extracts the data for cluster k
extractcluster = function(X, Z, clusters, sig2gam, sig2eps, k){
  # X is the n x p matrix of covariates including intercept
  # Z is the n x 1 matrix of treatment indicators
  # clusters is the n x 1 matrix of cluster indicators
  # sig2gam is the variance of random effects
  # sig2eps is the variance of the error term
  # k is the cluster number
  
  p = ncol(X)
  nk = length(clusters[clusters == k])
  Zk = matrix(Z[clusters == k], ncol = 1, nrow = nk)
  mk = sum(Zk)#length(Z[Z == 1 & clusters == k]) # num treated in cluster k
  rhok = sig2gam/(sig2eps + nk*sig2gam)
  Xk = matrix(X[clusters == k,], ncol = p, nrow = nk)
  Xbar = matrix(colMeans(Xk), ncol = p, nrow = 1)
  if (sum(Zk) > 0){
    Xktreated = matrix(Xk[Zk == 1,], nrow = sum(Zk), ncol = p)
    Xtbar = matrix(colMeans(Xktreated), ncol = p, nrow = 1)
  } 
  else { # to accomodate clusters with no treated units
    Xtbar = matrix(0, nrow = 1, ncol = p)
  }

  # return a list with mk, rhok, Xk, Zk, nk, Xbar, Xtbar
  return(list(mk = mk, rhok = rhok, Xk = Xk, Zk = Zk, nk = nk, Xbar = Xbar, Xtbar = Xtbar))
}

computeSinv = function(X, Z, clusters, sig2gam, sig2eps){
  # X is the n x p matrix of covariates including intercept
  # Z is the n x 1 matrix of treatment indicators
  # clusters is the n x 1 matrix of cluster indicators
  # sig2gam is the variance of random effects
  # sig2eps is the variance of the error term
  # returns the Sinv matrix
  
  p = ncol(X)
  S = t(X)%*%X
  K = length(unique(clusters))
  for (k in 1:K){
    dat = extractcluster(X, Z, clusters, sig2gam, sig2eps, k)
    nk = dat$nk
    rhok = dat$rhok
    Xbark = dat$Xbar
    S = S - nk^2*rhok*t(Xbark)%*%Xbark
  }
  # return the inverse of S
  return(pd.solve(S))
}

computeSinvfast = function(X, nks, rhoks, Xbarks){
  # Multiply each Xbarsum by the corresponding nk^2*rhok
  Xbarsum = lapply(1:length(nks), function(i) nks[i]^2*rhoks[i]*t(Xbarks[[i]])%*%Xbarks[[i]])
  # Sum the Xbarsums
  S = Reduce("+", Xbarsum)
  S = t(X) %*% X - S
  return(pd.solve(S))
}

compute_treated_avgsum = function(X, Z, clusters, sig2gam, sig2eps){
  # X is the n x p matrix of covariates including intercept
  # Z is the n x 1 matrix of treatment indicators
  # clusters is the n x 1 matrix of cluster indicators
  # sig2gam is the variance of random effects
  # sig2eps is the variance of the error term

  p = ncol(X)
  K = length(unique(clusters))
  out = matrix(0, nrow = 1, ncol = p)
  for (k in 1:K){
    dat = extractcluster(X, Z, clusters, sig2gam, sig2eps, k)
    nk = dat$nk
    mk = dat$mk
    rhok = dat$rhok
    Xbark = dat$Xbar
    Xtbark = dat$Xtbar
    out = out + mk*(Xtbark - nk*rhok*Xbark)
  }
  return(out)
}

compute_treated_avgsum_fast = function(nks, mks, rhoks, Xbarks, Xtbarks){
  # Multiply each Xbarsum by the corresponding mk*(Xtbark - nk*rhok*Xbark)
  Xbarsum = lapply(1:length(Xbarks), function(i) mks[i]*(Xtbarks[[i]] - nks[i]*rhoks[i]*Xbarks[[i]]))
  # Sum the Xbarsums
  out = Reduce("+", Xbarsum)
  return(out)
}

fixedeffects = function(X,Z){
  # X is the n x p matrix of covariates including intercept
  # Z is the n x 1 matrix of treatment indicators
  # returns the fixed effects implied weights
  
  nt = sum(Z)
  nc = sum(1-Z)
  n = nrow(X)
  p = ncol(X)
  Xt = matrix(X[Z == 1,], nrow = nt, ncol = p)
  Xbart = t(matrix(colMeans(Xt), nrow = 1, ncol = p))
  Xc = matrix(X[Z == 0,], nrow = nc, ncol = p)
  Xbarc = t(matrix(colMeans(Xc), nrow = 1, ncol = p))
  Xbar = t(matrix(colMeans(X), nrow = 1, ncol = p))
  # St is sample variance of treated, or sum of vecxi*vecxi' for treated- nt*Xbart%*%t(Xbart)
  St = matrix(0, nrow = p, ncol = p)
  for (i in which(Z == 1)){
    St = St + (X[i,]-Xbart)%*%t(X[i,]-Xbart)
  }
  St = St - nt*Xbart%*%t(Xbart)
  # Sc is sample variance of control, or sum of vecxi*vecxi' for control- nc*Xbarc%*%t(Xbarc)
  Sc = matrix(0, nrow = p, ncol = p)
  for (i in which(Z == 0)){
    Sc = Sc + (X[i,]-Xbarc)%*%t(X[i,]-Xbarc)
  }
  # S is the sum of the two sample variances
  S = St + Sc
  # Sinv is the inverse of S
  Sinv = solve(S) 
  weights = rep(NA, n)
  for (i in 1:n){
    if (Z[i] == 1){
      weights[i] = (1/nt) + (n/nc)*t(X[i,]-Xbart) %*% Sinv %*% (Xbar-Xbart)
    }
    else{
      weights[i] = (1/nc) + (n/nt)*t(X[i,]-Xbarc) %*% Sinv %*% (Xbar-Xbarc)
    }
  }
  return(weights)
}

empirical_bayes = function(Y, X, Z, clusters, sig2gam, sig2eps){
  # Y is the n x 1 matrix of outcomes
  # X is the n x p matrix of covariates including intercept
  # Z is the n x 1 matrix of treatment indicators
  # clusters is the n x 1 matrix of cluster indicators
  # sig2gam is the variance of random effects
  # sig2eps is the variance of the error term
  # returns the empirical bayes estimate of the random effects
  
  p = ncol(X)
  K = length(unique(clusters))
  gam_bayes = rep(NA, K)
  betahat = solve(t(X)%*%X)%*%t(X)%*%Y
  for (k in 1:K){
    dat = extractcluster(X, Z, clusters, sig2gam, sig2eps, k)
    nk = dat$nk
    rhok = dat$rhok
    Xk = dat$Xk
    Yk = Y[clusters == k]
    gam_bayes[k] = sig2gam*matrix(1, nrow = 1, ncol = nk) %*% solve(sig2gam*matrix(1, nrow = nk, ncol = nk) + sig2eps*diag(nk)) %*% (Yk - Xk%*%betahat)
  }
  return(gam_bayes)
}

car = function(X, Z, phi, A, sig2eps){
  # X is the n x p matrix of covariates including intercept
  # Z is the n x 1 matrix of treatment indicators
  # phi controls the amount of spatial correlation
  # A is the n x n adjacency matrix
  # sig2eps is the variance of the error term
  # returns implied weights of CAR model
  
  # Compute degree matrix D from A
  D = diag(apply(A, 1, sum))
  # Compute CAR variance
  #Sigmacar = sig2eps*solve(D - phi*A)
  Sinv = 1/sig2eps*(D-phi*A)#solve(Sigmacar)
  n = nrow(X)
  lnum = (diag(1,n) - Sinv %*% X %*% solve(t(X) %*% Sinv %*% X) %*% t(X)) %*% Sinv %*% Z
  ldenom = t(Z) %*% Sinv %*% (diag(1,n) - X %*% solve(t(X) %*% Sinv %*% X) %*% t(X) %*% Sinv) %*% Z
  l = lnum/as.numeric(ldenom)
  return(as.numeric((2*Z-1)*l))
}

sar = function(X, Z, phi, A, sig2eps){
  # X is the n x p matrix of covariates including intercept
  # Z is the n x 1 matrix of treatment indicators
  # phi controls the amount of spatial correlation
  # A is the n x n adjacency matrix
  # sig2eps is the variance of the error term
  # returns implied weights of SAR model
  
  # Compute degree matrix D from A
  D = diag(apply(A, 1, sum))
  Dinv = solve(D)
  n = nrow(X)
  mult = solve(diag(1,n) - phi*Dinv %*% A)
  # Compute SAR variance
  Sigmasar = sig2eps*mult%*% Dinv%*%t(mult)
  Sinv = solve(Sigmasar)
  lnum = (diag(1,n) - Sinv %*% X %*% solve(t(X) %*% Sinv %*% X) %*% t(X)) %*% Sinv %*% Z
  ldenom = t(Z) %*% Sinv %*% (diag(1,n) - X %*% solve(t(X) %*% Sinv %*% X) %*% t(X) %*% Sinv) %*% Z
  l = lnum/as.numeric(ldenom)
  return(as.numeric((2*Z-1)*l))
}

impliedweightsgeneral = function(X, Z, Sigmainv, returnL = F){
  # X is the n x p matrix of covariates including intercept
  # Z is the n x 1 matrix of treatment indicators
  # Sigmainv is the n x n precision matrix
  # returns implied weights of general model
  
  n = nrow(X)
  lnum = (diag(1,n) - Sigmainv %*% X %*% solve(t(X) %*% Sigmainv %*% X) %*% t(X)) %*% Sigmainv %*% Z
  ldenom = t(Z) %*% Sigmainv %*% (diag(1,n) - X %*% solve(t(X) %*% Sigmainv %*% X) %*% t(X) %*% Sigmainv) %*% Z
  l = lnum/as.numeric(ldenom)
  if (returnL){
    return(as.numeric(l))
  }
  return(as.numeric((2*Z-1)*as.numeric(l))) 
}

target = function(w, Z, X, Q, transform = F){
  # w is the n x 1 matrix of weights
  # Z is the n x 1 matrix of treatment indicators
  # X is the n x p matrix of covariates including intercept
  # Q is the n x n square root of precision matrix 
  # transform is a boolean indicating whether we're looking at the transformed or untransformed domain
  # returns the target covariate vector
  
  n = nrow(X)
  if (!transform){
    # Control target is sum of wi Xi in control group (Zi = 0)
    ctrltarget = t((1-Z)*w) %*% X
    treattarget = t(Z*w) %*% X
    return(list(ctrltarget = ctrltarget, treattarget = treattarget))
  }
  if(transform){
    return(NULL)
  }
}

compute_ess = function(w){
  # w is the n x 1 matrix of weights
  # returns the effective sample size
  
  return(sum(abs(w))^2/sum(w^2))
}

disp = function(w, wtilde, scaleweights){
  return(sum((w - wtilde)^2/scaleweights))
}

# Function that performs quadratic programming problem on control(treated) units in section 6
qp = function(N,
              Nc,
              baseweights = rep(1/N, Nc),
              scaleweights = baseweights,
              X, # does this include the intercept? 
              Xstar, # does this include the intercept? 
              caliper = 0,
              sumconstraint = 1) {
  # N is the number of units
  # Nc is the number of control(treated) units
  # baseweights are the Nc x 1 matrix of base weights in control group
  # scaleweights are the Nc x 1 matrix of scale weights in control group
  # X is the Nc x p matrix of covariates in the control group with the intercept
  # Xstar is a p x 1 matrix of covariate averages in the target pop with the intercept
  # caliper is the caliper for the matching, i.e. delta
  # sumconstraint is the sum constraint for the weights
  print(head(baseweights))
  wtilde = baseweights/sum(baseweights)
  #print(wtilde)
  #print(dim(scaleweights))
  # out = constrOptim(theta = baseweights, function(x) disp(x, wtilde, scaleweights), 
  #                   grad = 
  Dmat = 2*diag(as.numeric(1/scaleweights))
  dvec = Dmat %*% wtilde
  Amat = as.matrix(X)
  bvec = Xstar
  #print(Dmat[1:5,1:5])
  #print(dvec[1:5])
  #print(Amat[1:5,1:5])
  out = solve.QP(Dmat, dvec, Amat, bvec, meq = length(bvec))
  return(out$solution)
}

qp_sbw <- function(covs, target, baseweights){
  # covs is a n x p matrix of user-specified covariates to balance 
  # target is a p x 1 matrix of target covariate values
  # baseweights is a n x 1 vector of base weights
  
  n = length(baseweights)
  stopifnot(nrow(covs) == n)
  stopifnot(ncol(covs) == length(target))
  b0 = c(1, target, rep(0,n))
  A = cbind(1, covs, diag(x = 1, nrow = n, ncol = n))
  D = 2*diag(x = 1, nrow = n, ncol = n)
  d = 2*baseweights
  out = solve.QP(Dmat = D, dvec = d, Amat = A, bvec = b0, meq = length(target) + 1)
  return(out$solution)
}

compute_target <- function(w, X){
  # w is the n x 1 matrix of weights
  # X is the n x p matrix of covariates including intercept
  # returns the target covariate vector
  return(t(w) %*% X)
}

extractdata = function(X, Z){
  # X is the n x p matrix of covariates including intercept
  # Z is the n x 1 matrix of treatment indicators
  # returns 
  
  nt = sum(Z)
  nc = sum(1-Z)
  n = nrow(X)
  p = ncol(X)
  Xt = matrix(X[Z == 1,], nrow = nt, ncol = p)
  Xbart = t(matrix(colMeans(Xt), nrow = 1, ncol = p))
  Xc = matrix(X[Z == 0,], nrow = nc, ncol = p)
  Xbarc = t(matrix(colMeans(Xc), nrow = 1, ncol = p))
  Xbar = t(matrix(colMeans(X), nrow = 1, ncol = p))
  # St is sample variance of treated, or sum of vecxi*vecxi' for treated- nt*Xbart%*%t(Xbart)
  St = matrix(0, nrow = p, ncol = p)
  for (i in which(Z == 1)){
    St = St + (X[i,]-Xbart)%*%t(X[i,]-Xbart)
  }
  St = St - nt*Xbart%*%t(Xbart)
  # Sc is sample variance of control, or sum of vecxi*vecxi' for control- nc*Xbarc%*%t(Xbarc)
  Sc = matrix(0, nrow = p, ncol = p)
  for (i in which(Z == 0)){
    Sc = Sc + (X[i,]-Xbarc)%*%t(X[i,]-Xbarc)
  }
  # S is the sum of the two sample variances
  S = St + Sc
  # Sinv is the inverse of S
  Sinv = solve(S) 
  Xstar = Sc %*% Sinv %*% Xbart + St %*% Sinv %*% Xbarc
  
  return(list(n = n, nc = nc, nt = nt, Xt = Xt, Xbart = Xbart, Xc = Xc, Xbarc = Xbarc, Xbar = Xbar, St = St, Sc = Sc, S = S, Sinv = Sinv, Xstar = Xstar))
}

computeXstar = function(X, Z, scaleweights){
  # X is the n x p matrix of covariates including intercept
  # Z is the n x 1 matrix of treatment indicators
  # scaleweights is the n x 1 matrix of scale weights
  
  nt = sum(Z)
  nc = sum(1-Z)
  n = nrow(X)
  p = ncol(X)
  Xt = matrix(X[Z == 1,], nrow = nt, ncol = p)
  Xsbart = t(Xt) %*% scaleweights[Z == 1]/sum(scaleweights[Z == 1])
  Xc = matrix(X[Z == 0,], nrow = nc, ncol = p)
  Xsbarc = t(Xc) %*% scaleweights[Z == 0]/sum(scaleweights[Z == 0])
  # St is sample variance of treated, or sum of vecxi*vecxi' for treated- nt*Xbart%*%t(Xbart)
  Sst = matrix(0, nrow = p, ncol = p)
  for (i in which(Z == 1)){
    Sst = Sst + scaleweights[i]*(X[i,]-Xsbart)%*%t(X[i,]-Xsbart)
  }
  Sst = Sst - nt*Xbart%*%t(Xbart)
  # Sc is sample variance of control, or sum of vecxi*vecxi' for control- nc*Xbarc%*%t(Xbarc)
  Ssc = matrix(0, nrow = p, ncol = p)
  for (i in which(Z == 0)){
    Ssc = Ssc + scaleweights[i]*(X[i,]-Xsbarc)%*%t(X[i,]-Xsbarc)
  }
  # S is the sum of the two sample variances
  Ss = Sst + Ssc
  # Sinv is the inverse of S
  Ssinv = solve(Ss) 
  Xstar = Ssc %*% Ssinv %*% Xsbart + Sst %*% Ssinv %*% Xsbarc
  
  return(Xstar = Xstar)
}

contamination = function(Z1, X1state, Z2){
  mod = lm(Z1 ~ X1state)
  # Create a new column in mapdat that has name val_contam
  contam_weights = rep(NA, length(Z1))
  xtildetilde1 = mod$residuals
  for (i in 1:length(Z1)){
    st = grid_sf2$state[i]
    grid_sf2$lam12nonzero[i] = mean(xtildetilde1*Z2[grid_sf2$state == st])/mean(xtildetilde1^2)
  }
}

plotfunc = function(df, names, stateboundaries = NULL, labels=NULL){
  # df is a sf dataframe including columns names
  # names is a vector of strings with the names of the columns to be plotted
  # labels is a vector of strings to title the ggplots
  
  # returns a list of ggplots
  
  if (is.null(labels)){
    labels = names
  }
  K = length(names)
  gs = list()
  for (k in 1:K){
    # extract the column with names[k] from df
    var = df[[names[k]]]
    qs = round(quantile(var, probs = c(0.1, 0.3, 0.5, 0.7, 0.9)),2)
    # Turn qs into a vector of strings
    qschar = as.character(qs)
    
    gs[[k]] = ggplot(df) +
      #xlim(-125, -65) + 
      #ylim(25, 50) +
      geom_sf(aes_string(fill = names[k]), color=NA, size = 0.005) +
      geom_sf(data = stateboundaries, fill = NA, color = "black", size = 10) +
      scale_fill_gradient2(low = "pink", 
                           mid = "white", 
                           high = "lightblue", 
                           midpoint = 0,
                           breaks = qs,
                           labels = qschar,
                           limits = c(min(var), max(var) + 0.001),
                           na.value = "white") +
      theme_minimal() +
      theme(#plot.title = element_text(size = 24 * 2,hjust = 0.5),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            line = element_blank(),
            axis.title = element_blank(),
            legend.position = "bottom",
            legend.direction = "horizontal", 
            legend.text.align = 0.75,
            legend.key.width = unit(50, "points"),
            panel.grid.major = element_line(colour = "transparent")) + 
      ggtitle(labels[k])
  }
  return(gs)
}

compute_impliedweights = function(X, Z, clusters, adj, 
                                  sig2gam = 0.5, 
                                  sig2eps = 1,
                                  phi = 0.7){
  # X is the n x p matrix of covariates including intercept but WITHOUT cluster indicators
  # Z is the n x 1 matrix of treatment (can be binary or multivalued)
  # clusters is the n x 1 matrix of cluster assignments
  # adj is the adjacency matrix of the graph
  # sig2gam is the variance of the random intercepts for re model
  # sig2eps is the variance of the random errors for re model
  # phi is the spatial autocorrelation parameter for CAR/SAR
  
  # If Z is binary, then we are in the binary treatment case
  if (length(unique(Z)) == 2){
    # Initialize empty data frame with 5 columns titled 
    # Pooled, RE, FE, CAR, SAR
    df = data.frame(Pooled = rep(NA, nrow(X)),
                    RE = rep(NA, nrow(X)),
                    FE = rep(NA, nrow(X)),
                    CAR = rep(NA, nrow(X)),
                    SAR = rep(NA, nrow(X)))
    df$Pooled = fixedeffects(X = X, Z = Z)
    state.dummy = model.matrix(~factor(clusters)-1)
    Xstate = cbind(X, state.dummy[,-1]) # already intercept in X
    Xstate = as.matrix(Xstate)
    df$FE = fixedeffects(X = Xstate, Z = Z)
    df$RE = impliedweightsslow(
      X = X,
      Z = Z,
      clusters = clusters,
      sig2gam = sig2gam,
      sig2eps = sig2eps
    )
    df$CAR = car(
      X = X,
      Z = Z,
      phi = phi,
      A = adj,
      sig2eps = sig2eps
    )
    df$SAR = sar(
      X = X,
      Z = Z,
      phi = phi,
      A = adj,
      sig2eps = sig2eps
    )
    return(df)
  }
  else{
    # Create indicator variables for each level of Z
    multitreatment.dummy = model.matrix(~factor(Z)-1)
    multitreatment = multitreatment.dummy[, -1]
    dfs = list()
    for (i in 1:ncol(multitreatment)){
      Znew = multitreatment[,i]
      Xnew = cbind(X, multitreatment[,-i])
      dfs[[i]] = compute_impliedweights(X = Xnew,
                                        Z = Znew, 
                                        clusters = clusters,
                                        adj = adj,
                                        sig2gam = sig2gam,
                                        sig2eps = sig2eps,
                                        phi = phi)
    }
    return(dfs)
  }
}

# Compute weights using spatial precision 
check_property3 = function(Q,X,Z){
  
  # Projection matrix P_{QX}
  QX <- Q %*% X
  P_QX <- QX %*% solve(t(QX) %*% QX) %*% t(QX)
  
  n = nrow(Q)
  eigen_Q <- eigen(Q)
  eigenvectors_Q <- eigen_Q$vectors
  eigenvalues_Q <- eigen_Q$values
  # Identify eigenvectors of Q orthogonal to QX
  orthogonal_eigenvectors <- list()
  # save their corresponding eigenvalues
  orthogonal_eigenvalues <- list()
  for (i in 1:ncol(eigenvectors_Q)) {
    eigvec <- eigenvectors_Q[, i]
    # Check if the eigenvector is orthogonal to the column space of QX
    if (all(abs(t(X) %*% eigvec) < 1e-6)) {
      orthogonal_eigenvectors[[length(orthogonal_eigenvectors) + 1]] <- eigvec
      orthogonal_eigenvalues[[length(orthogonal_eigenvalues) + 1]] <- eigenvalues_Q[i]
    }
  }
  
  # Convert the list of orthogonal eigenvectors to a matrix
  orthogonal_eigenvectors_matrix <- do.call(cbind, orthogonal_eigenvectors)
  print(ncol(orthogonal_eigenvectors_matrix))
  
  l <- rep(0, n)
  for (k in 1:ncol(orthogonal_eigenvectors_matrix)) {
    vec <- orthogonal_eigenvectors_matrix[, k]
    l <- l + orthogonal_eigenvalues[[k]]^2 * vec %*% t(vec) %*% Z
  }
  c <- t(Z) %*% Q %*% (diag(1,n) - P_QX) %*% Q %*% Z
  return((2*Z-1)*l/c)
}

# Function that computes sig2gam for desired cluster balances
givenbalance_getsig2gam <- function(sig2eps=1, 
                                    X, 
                                    clusters, 
                                    Z, 
                                    vector_to_balance,
                                    balance_thresh = 0.1){
  # sig2eps is the variance of the random errors for RE model
  # X is the n x p matrix of covariates including intercept but WITHOUT cluster indicators
  # clusters is the n x 1 matrix of cluster assignments
  # Z is the n x 1 matrix of binary treatment
  # vector_to_balance is a vector such as a cluster indicator for which we want to achieve balance
  # balance_thresh is the threshold for balance
  
  n = nrow(X)
  p = ncol(X)
  stopifnot(sig2eps > 0)
  stopifnot(length(clusters) == n)
  stopifnot(length(Z) == n)
  stopifnot(length(vector_to_balance) == n)
  
}

computebaseweights <- function(Siginv, Z){
  #Siginv is the inverse of the covariance matrix of errors
  #Z is the n x 1 matrix of treatment
  # if (length(unique(Z)) == 2){
  #   num = (2*Z-1) * (Siginv %*% (2*Z-1))
  #   denom = t(Z) %*% Siginv %*% (2*Z-1)
  # }
  
  #if (length(unique(Z)) == 2){
    nt = sum(Z)
    n = length(Z)
    b = -nt/n
    num = (2*Z-1) * (Siginv %*% (Z + b))
    denom = t(Z) %*% Siginv %*% (Z + b)
  #}
  # else{
  #   num = Siginv %*% Z
  #   denom = t(Z) %*% Siginv %*% Z
  # }
  return(as.numeric(num/as.numeric(denom)))
}

compute_allbaseweights <- function(Z, 
                                   adjacency_matrix, 
                                   statefactor,
                                   sig2gam = 0.5, 
                                   sig2eps = 1,
                                   phi = 0.7,
                                   distmat){
  n = nrow(adjacency_matrix)
  df = data.frame(basePooled = rep(NA, n),
                  baseRE = rep(NA, n),
                  baseFE = rep(NA, n),
                  baseCAR = rep(NA, n),
                  baseSAR = rep(NA, n),
                  baseGP = rep(NA, n))
  
  # Pooled 
  Sigmainv = diag(n)
  df$basePooled = computebaseweights(Siginv = Sigmainv,
                                                 Z = Z)

  # CAR
  Sigmainv = diag(rowSums(adjacency_matrix)) - phi*adjacency_matrix
  df$baseCAR = computebaseweights(Siginv = Sigmainv,
                                              Z = Z)

  # SAR
  Sigmainv = (diag(1,n) - phi*diag(1/rowSums(adjacency_matrix)) %*% adjacency_matrix) %*%
    diag(rowSums(adjacency_matrix)) %*%
    (diag(1,n) - phi*diag(1/rowSums(adjacency_matrix)) %*% adjacency_matrix)
  df$baseSAR = computebaseweights(Siginv = Sigmainv,
                                              Z = Z)

  # RE
  Sigs = list()
  for (k in unique(statefactor)){
    nk = sum(statefactor == k)
    Sigs[[length(Sigs) + 1]] = diag(nk) - (sig2gam/(nk*sig2gam + sig2eps))*matrix(1, nrow = nk, ncol = nk)
  }
  Sigmainv = bdiag(Sigs)
  df$baseRE = computebaseweights(Siginv = Sigmainv,
                                 Z = Z)
  
  # FE
  Sigs = list()
  sig2gam = 10000
  for (k in unique(statefactor)){
    nk = sum(statefactor == k)
    Sigs[[length(Sigs) + 1]] = diag(nk) - (sig2gam/(nk*sig2gam + sig2eps))*matrix(1, nrow = nk, ncol = nk)
  }
  Sigmainv = bdiag(Sigs)
  df$baseFE = computebaseweights(Siginv = Sigmainv,
                                             Z = Z)
  
  # Gaussian Process
  kappa = 0.05 # 0.2 # as you increase, becomes more extreme
  rangec = 300
  phic <- rangec/(2*sqrt(kappa))
  Sigma <- geoR::matern(u = dmat, phi = phic, kappa = kappa)
  Sigmainv <- solve(Sigma)
  df$baseGP = computebaseweights(Siginv = Sigmainv,
                                             Z = Z)
  return(df)
}

spscale <- function(Siginv, x, threshold = 0.05, E = NULL){
  # Siginv is the inverse of the covariance matrix of errors
  # x is the random variable whose scale we want to compute
  # threshold is the threshold for ck
  if (is.null(E)){
    E <- eigen(Siginv)
  }
  V <- E$vectors
  xnormalized <- x/norm(x, type = "2")
  xstar <- t(V) %*% xnormalized
  max_lam <- max(E$values[abs(xstar) > threshold])
  return(1/max_lam)
}

spscale <- function(Siginv, x, E = NULL){
  # Siginv is the inverse of the covariance matrix of errors
  # x is the random variable whose scale we want to compute
  if (is.null(E)){
    E <- eigen(Siginv)
  }
  n <- length(x)
  V <- E$vectors
  xnormalized <- x/norm(x, type = "2")
  xstar <- t(V) %*% xnormalized
  xstar <- rev(xstar)
  # find max {p: xstar[k] < 1/k^p for all k}
  # Try a sequence of ps
  ps <- seq(0.01, 1, by = 0.01)
  pmax <- 0.01
  for (p in ps){
    if (all(xstar < 1/(1:n)^p)){ #  #rev(E$values)
      pmax = p
    }
  }
  return(pmax)
}

spscale <- function(Siginv, x, E = NULL){
  # Siginv is the inverse of the covariance matrix of errors
  # x is the random variable whose scale we want to compute
  if (is.null(E)){
    E <- eigen(Siginv)
  }
  n <- length(x)
  V <- E$vectors
  xnormalized <- x/norm(x, type = "2")
  xstar <- t(V) %*% xnormalized
  return(max(xstar[E$values > 1]))
}

spscale <- function(Siginv, x, E = NULL){
  # Siginv is the inverse of the covariance matrix of errors
  # x is the random variable whose scale we want to compute
  if (is.null(E)){
    E <- eigen(Siginv)
  }
  n <- length(x)
  V <- E$vectors
  xnormalized <- x/norm(x, type = "2")
  xstar <- t(V) %*% xnormalized
  num <- sum(E$values*xstar^2)
  denom <- sum(xstar^2)
  return(num/denom)
}

maximal_imbalance <- function(Siginv, Z, X, s, E = NULL, Q = NULL, threshold = 0.05){
  n <- nrow(Siginv)
  if (is.null(E)){
    E <- eigen(Siginv)
  }
  #print('Took Eigen')
  if (is.null(Q)){
    Q <- E$vectors %*% diag(sqrt(E$values)) %*% t(E$vectors)
  }
  #print('Computed Q')
  QX <- Q %*% X
  PQX <- QX %*% solve(t(QX) %*% QX) %*% t(QX)
  V <- E$vectors
  eigvals <- E$values
  denom <- as.numeric(t(Z) %*% Q %*% (diag(1,n) - PQX) %*% Q %*% Z)
  alphas <- as.numeric(t(Z) %*% Q %*% (diag(1, n) - PQX) %*% V) / denom
  #print('Computed alphas')
  sum1 <- sum((E$values^(3/2)*abs(alphas))[E$values < 1/s])
  sum2 <- sum((E$values^(3/2)*abs(alphas)*threshold)[E$values >= 1/s])
  max_imbalance <- sum1 + sum2
  return(max_imbalance)
}

maximal_imbalance <- function(Siginv, Siginv_spatial, Z, X, s, E = NULL, Q = NULL,
                              Espatial = NULL, Qspatial = NULL, threshold = 0.05){
  n <- nrow(Siginv)
  if (is.null(E)){
    E <- eigen(Siginv)
  }
  if (is.null(Espatial)){
    Espatial <- eigen(Siginv_spatial)
  }
  #print('Took Eigen')
  if (is.null(Q)){
    Q <- E$vectors %*% diag(sqrt(E$values)) %*% t(E$vectors)
  }
  if (is.null(Qspatial)){
    Qspatial <- Espatial$vectors %*% diag(sqrt(Espatial$values)) %*% t(Espatial$vectors)
  }
  #print('Computed Q')
  QX <- Q %*% X
  PQX <- QX %*% solve(t(QX) %*% QX) %*% t(QX)
  Vspatial <- Espatial$vectors
  denom <- as.numeric(t(Z) %*% Q %*% (diag(1,n) - PQX) %*% Q %*% Z)
  alphas <- as.numeric(t(Z) %*% Q %*% (diag(1, n) - PQX) %*% Q %*% Vspatial) / denom
  sum1 <- sum((Espatial$values*abs(alphas))[Espatial$values < 1/s])
  sum2 <- sum((Espatial$values*abs(alphas)*threshold)[Espatial$values >= 1/s])
  print(c(sum1,sum2))
  max_imbalance <- sum1 + sum2
  return(max_imbalance)
}

maximal_imbalance <- function(Siginv, Siginv_spatial, Z, X, p, E = NULL, Q = NULL,
                              Espatial = NULL, Qspatial = NULL){
  n <- nrow(Siginv)
  if (is.null(E)){
    E <- eigen(Siginv)
  }
  if (is.null(Espatial)){
    Espatial <- eigen(Siginv_spatial)
  }
  #print('Took Eigen')
  if (is.null(Q)){
    Q <- E$vectors %*% diag(sqrt(E$values)) %*% t(E$vectors)
  }
  if (is.null(Qspatial)){
    Qspatial <- Espatial$vectors %*% diag(sqrt(Espatial$values)) %*% t(Espatial$vectors)
  }
  #print('Computed Q')
  QX <- Q %*% X
  PQX <- QX %*% solve(t(QX) %*% QX) %*% t(QX)
  Vspatial <- Espatial$vectors
  denom <- as.numeric(t(Z) %*% Q %*% (diag(1,n) - PQX) %*% Q %*% Z)
  alphas <- as.numeric(t(Z) %*% Q %*% (diag(1, n) - PQX) %*% Q %*% Vspatial) / denom
  max_imbalance <-  sum(Espatial$values*abs(alphas)*1/(n:1)^p)
  return(max_imbalance)
}

maximal_imbalance <- function(Siginv, Siginv_spatial, Z, X, p, E = NULL, Q = NULL,
                              Espatial = NULL, Qspatial = NULL, threshold = 100){
  n <- nrow(Siginv)
  if (is.null(E)){
    E <- eigen(Siginv)
  }
  if (is.null(Espatial)){
    Espatial <- eigen(Siginv_spatial)
  }
  #print('Took Eigen')
  if (is.null(Q)){
    Q <- E$vectors %*% diag(sqrt(E$values)) %*% t(E$vectors)
  }
  if (is.null(Qspatial)){
    Qspatial <- Espatial$vectors %*% diag(sqrt(Espatial$values)) %*% t(Espatial$vectors)
  }
  #print('Computed Q')
  QX <- Q %*% X
  PQX <- QX %*% solve(t(QX) %*% QX) %*% t(QX)
  Vspatial <- Espatial$vectors
  denom <- as.numeric(t(Z) %*% Q %*% (diag(1,n) - PQX) %*% Q %*% Z)
  alphas <- as.numeric(t(Z) %*% Q %*% (diag(1, n) - PQX) %*% Q %*% Vspatial) / denom
  #print('Computed alphas')
  #print(summary(alphas))
  #sum1 <- sum((Espatial$values*abs(alphas)))
  sum1 <- sum((Espatial$values*abs(alphas)*(1/(n:1))^p)[(n-threshold):n])
  norm2 <- norm(t(Z) %*% Q %*% (diag(1, n) - PQX) %*% Q %*% Vspatial[,1:(threshold-1)] %*% 
                 diag(Espatial$values[1:(threshold-1)]), type = '2')/denom
  norm22 <- sqrt((n-threshold)*(1/threshold^p))
  sum2 <- norm2*norm22
  max_imbalance <- sum1 + sum2
  #max_imbalance <- sum((Espatial$values*abs(alphas)*(1/E$values)^p))
  return(max_imbalance)
}

maximal_imbalance <- function(Siginv, Siginv_spatial, Z, X, E = NULL, Q = NULL,
                              Espatial = NULL, Qspatial = NULL, nreps = 3000){
  n <- nrow(Siginv)
  if (is.null(E)){
    E <- eigen(Siginv)
  }
  if (is.null(Espatial)){
    Espatial <- eigen(Siginv_spatial)
  }
  #print('Took Eigen')
  if (is.null(Q)){
    Q <- E$vectors %*% diag(sqrt(E$values)) %*% t(E$vectors)
  }
  if (is.null(Qspatial)){
    Qspatial <- Espatial$vectors %*% diag(sqrt(Espatial$values)) %*% t(Espatial$vectors)
  }
  #print('Computed Q')
  QX <- Q %*% X
  PQX <- QX %*% solve(t(QX) %*% QX) %*% t(QX)
  denom <- as.numeric(t(Z) %*% Q %*% (diag(1,n) - PQX) %*% Q %*% Z)
  imbs <- rep(NA, nreps)
  scales <- rep(NA, nreps)
  for (rep in 1:nreps){
    print(rep)
    U <- rnorm(n, sd = 0.1)
    U <- U/norm(U, type = '2')
    scales[rep] <- spscale(Siginv_spatial, U, E = Espatial)
    imbs[rep] <- as.numeric(t(Z) %*% Q %*% (diag(1, n) - PQX) %*% Q %*% U) / denom
  }
  return(data.frame(imbs = imbs, scales = scales))
}

maximal_imbalance <- function(Siginv, Siginv_spatial, Z, X, m, E = NULL, Q = NULL,
                              Espatial = NULL, Qspatial = NULL, threshold = 0.05){
  n <- nrow(Siginv)
  if (is.null(E)){
    E <- eigen(Siginv)
  }
  if (is.null(Espatial)){
    Espatial <- eigen(Siginv_spatial)
  }
  #print('Took Eigen')
  if (is.null(Q)){
    Q <- E$vectors %*% diag(sqrt(E$values)) %*% t(E$vectors)
  }
  if (is.null(Qspatial)){
    Qspatial <- Espatial$vectors %*% diag(sqrt(Espatial$values)) %*% t(Espatial$vectors)
  }
  #print('Computed Q')
  QX <- Q %*% X
  PQX <- QX %*% solve(t(QX) %*% QX) %*% t(QX)
  Vspatial <- Espatial$vectors
  denom <- as.numeric(t(Z) %*% Q %*% (diag(1,n) - PQX) %*% Q %*% Z)
  alphas <- as.numeric(t(Z) %*% Q %*% (diag(1, n) - PQX) %*% Q %*% Vspatial) / denom
  #print('Computed alphas')
  #print(summary(alphas))
  sum1 <- sum((Espatial$values*abs(alphas))[Espatial$values <= 1])
  sum2 <- sum((Espatial$values*abs(alphas)*m)[Espatial$values > 1])
  max_imbalance <- sum1 + sum2
  return(sum1)
}


max_imbalance_moran <- function(Siginv, Z, X, nreps = 1000, nb, lw, E = NULL, Q = NULL){
  n <- nrow(Siginv)
  if (is.null(E)){
    E <- eigen(Siginv)
  }
  if (is.null(Q)){
    Q <- E$vectors %*% diag(sqrt(E$values)) %*% t(E$vectors)
  }
  QX <- Q %*% X
  PQX <- QX %*% solve(t(QX) %*% QX) %*% t(QX)
  V <- E$vectors
  eigvals <- E$values
  denom <- as.numeric(t(Z) %*% Q %*% (diag(1,n) - PQX) %*% Q %*% Z)
  mult <- t(Z) %*% Q %*% (diag(1, n) - PQX) %*% Q
  imbs = rep(NA, nreps)
  morans = rep(NA, nreps)
  for (j in 1:nreps){
    U = rnorm(n, sd = 0.01)
    U = U/norm(U, type = '2') # c sums to 1
    imbs[j] = as.numeric(mult %*% U) / denom
    morans[j] = moran_adhoc(U, Wmat = Siginv)#moran(U, lw, length(nb), Szero(lw))[1]$I
  }
  return(data.frame(imbs = imbs, morans = morans))
}

moran_adhoc <- function(x, Wmat){
  n <- length(x)
  stopifnot(n == nrow(Wmat))
  stopifnot(n == ncol(Wmat))
  xbar <- mean(x)
  #coeff <- n / sum(Wmat)
  coeff <- 1
  
  # Calculate the numerator using matrix operations
  centered_x <- x - xbar
  num <- t(centered_x) %*% Wmat %*% centered_x
  
  # Calculate the denominator
  denom <- sum((x - xbar)^2)
  
  # Return the result
  result <- coeff * num / denom
  return(result)
}

dapsm_ps <- function(treated, 
                     control,
                     w = 0.1,
                     caliper = 0.05){
  # treated is a dataframe of treated units with columns Latitude, Longitude, prop.scores
  # control is a dataframe of control units with columns Latitude, Longitude, prop.scores
  # w is weight used to balance between matching on prop and geography
  # coords is the n x 2 matrix of coordinates
  
  stopifnot(all(c("Latitude", "Longitude", "prop.scores", "id") %in% colnames(treated)))
  stopifnot(all(c("Latitude", "Longitude", "prop.scores", "id") %in% colnames(control)))
  
  dist.mat <- fields::rdist(cbind(treated$Longitude, treated$Latitude),
                           cbind(control$Longitude, control$Latitude))
  stand.dist.mat <- StandDist(dist.mat)
  print(summary(as.numeric(stand.dist.mat)))
          
  ps.diff <- outer(treated$prop.scores, control$prop.scores, "-")
  ps.diff <- abs(ps.diff)
  ps.diff <- StandDist(ps.diff)
  print(summary(as.numeric(ps.diff)))
          
  dapscore <- (1-w)*stand.dist.mat + w*ps.diff
  pairs <- MinDistMatch(dapscore, caliper = caliper)
  return(cbind.data.frame('treated' = treated$id[pairs[,1]],
                          'control' = control$id[pairs[,2]]))
}

compute_dus <- function(U, Q, X){
  # U is the unmeasured confounder
  # Q is the square root of precision
  # X is the design matrix
  # returns the squared mahalanobis distance of U, U^t Q (I-P_QX) QU
  
  QX <- Q %*% X
  PQX <- QX %*% solve(t(QX) %*% QX) %*% t(QX)
  
  return(apply(U, 2, function(u) sqrt(as.numeric(t(u) %*% Q %*% (diag(1, nrow(Q)) - PQX) %*% Q %*% u))))
}

compute_duzs <- function(U, Z, Q, X){
  # U is the unmeasured confounder
  # Z is the treatment
  # Q is the square root of precision
  # X is the design matrix
  # returns the squared mahalanobis distance between U and Z, Z^t Q(I-P_QX) QU
  
  QX <- Q %*% X
  PQX <- QX %*% solve(t(QX) %*% QX) %*% t(QX)
  return(apply(U, 2, function(u) sqrt(as.numeric(t(Z-u) %*% Q %*% (diag(1, nrow(Q)) - PQX) %*% Q %*% (Z-u)))))
}

imbalance_mahalanobis <- function(dus, duzs, dz){
  # dus is the vector of mahalanobis distances of U, U^t Q (I-P_QX) QU 
  # duzs is the vector of mahalanobis distances between U and Z, Z^t Q(I-P_QX) QU
  # dz is the mahalanobis distance of Z, Z^t Q(I-P_QX) QZ
  # returns the imbalance
  
  res <- (-(duzs)^2 + (dus)^2 + (dz)^2) / (2*(dz)^2)
  return(res)
}

correlation_Usize <- function(U, Z, X, Q, SiginvX = NULL){
  # U is the unmeasured confounder
  # Z is the treatment
  # X is the design matrix
  # Q is the square root of precision
  # SiginvX is Q (I-PX) Q
  # returns the absolute correlation between QU and QZ, as well as size of U and Z
  
  if (is.null(SiginvX)){
    QX <- Q %*% X
    PQX <- QX %*% solve(t(QX) %*% QX) %*% t(QX)
    SiginvX <- Q %*% (diag(1, nrow(Q)) - PQX) %*% Q
  }
  
  coruz_num <- apply(U, 2, function(u) as.numeric(t(Z) %*% SiginvX %*% u))
  #print(coruz_num)
  print(as.numeric(t(Z) %*% SiginvX %*% Z))
  print(apply(U,2, function(u) as.numeric(t(u) %*% SiginvX %*% u)))
  coruz_denom <- apply(U,2, function(u) sqrt(as.numeric(t(Z) %*% SiginvX %*% Z) * as.numeric(t(u) %*% SiginvX %*% u)))
  #print(coruz_denom)
  coruz <- coruz_num / coruz_denom
  
  sizeu <- apply(U,2, function(u) sqrt(as.numeric(t(u) %*% SiginvX %*% u)))
  sizez <- sqrt(as.numeric(t(Z) %*% SiginvX %*% Z))
  
  return(list(abscoruz = abs(coruz), 
              sizeu = sizeu, 
              sizez = sizez))
}

abs_imbalance_coruz_size <- function(abscoruz, sizeu, sizez){
  # abscoruz is the absolute correlation between QU and QZ
  # sizeu is the size of U
  # sizez is the size of Z
  # returns the imbalance
  
  res <- abscoruz * sizeu / sizez
  return(res)
}


