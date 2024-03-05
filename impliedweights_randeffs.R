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
  lnum = (diag(1,n) - Sinv %*% X %*% solve(t(X) %*% Sinv %*% X) %*% t(X)) %*% Sinv %*% Z
  ldenom = t(Z) %*% Sinv %*% (diag(1,n) - X %*% solve(t(X) %*% Sinv %*% X) %*% t(X) %*% Sinv) %*% Z
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
  Sigmacar = sig2eps*solve(D - phi*A)
  Sinv = solve(Sigmacar)
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

