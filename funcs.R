gen_U <- function(Z, 
                  noise_sd = 0.1, 
                  clusters = NULL,
                  adjacency_matrix = NULL,
                  dmat = NULL,
                  Wre = NULL,
                  Wcar = NULL,
                  Wgp = NULL,
                  smoothing = c('clusterbased', 
                                'adjacencybased', 
                                'distancebased')){
  # Z is binary treatment
  # noise_sd is the standard deviation of the noise
  # clusters is a vector of cluster membership
  # adjacency_matrix is a matrix of adjacency
  # dmat is a distance matrix
  # smoothing is the type of smoothing to apply to U
  
  smoothing <- match.arg(smoothing)
  if (smoothing == 'clusterbased'){
    #stopifnot(!is.null(clusters))
    stopifnot(!is.null(Wre))
    # the i j entry of W should be I(C_i = C_j)
    # statefactor <- factor(clusters)
    # Jks <- list()
    # for (k in unique(statefactor)){
    #   nk <- sum(statefactor == k)
    #   Jks[[length(Jks) + 1]] <- matrix(1, nrow = nk, ncol = nk)
    # }
    W <- Wre #bdiag(Jks)
  }
  if (smoothing == 'adjacencybased'){
    #stopifnot(!is.null(adjacency_matrix))
    stopifnot(!is.null(Wcar))
    #W <- adjacency_matrix
    W <- Wcar
  }
  if (smoothing == 'distancebased'){
    #stopifnot(!is.null(dmat))
    #W <- dmat < 10 # try correctly specifying these. 
    stopifnot(!is.null(Wgp))
    W <- Wgp
  }
  U <- rnorm(length(Z), mean = Z, sd = noise_sd)
  # Smooth U using W
  #U <- W %*% t(W) %*% U
  U <- W %*% U / rowSums(W)
  return(U)
}

gen_Y <- function(X,
                  Z,
                  U,
                  beta,
                  tau = 1,
                  gamma = -0.5,
                  noise_sd = 0.1,
                  outcomemod = c('linear', 
                                 'linearinteraction', 
                                 'nonlinearinteraction')) {
  # X is the design matrix including intercept
  # Z is the binary treatment vector
  # U is the unobserved confounder  
  # beta is the vector of coefficients for X
  # tau is a coefficient multiplying Z
  # gamma is a coefficient multiplying U
  # noise_sd is the standard deviation of the noise
  # outcomemod is the type of outcome model to use
  
  # Check all dimensions conform
  stopifnot(ncol(X) == length(beta))
  stopifnot(nrow(X) == length(Z))
  stopifnot(nrow(X) == length(U))
  
  outcomemod <- match.arg(outcomemod)
  if (outcomemod == 'linear') {
    Y <- X %*% beta + tau * Z + gamma * U + rnorm(nrow(X), 0, noise_sd)
  } else if (outcomemod == 'linearinteraction') {
    Y <- X %*% beta + tau * Z + gamma * U +  U * Z + rnorm(nrow(X), 0, noise_sd)
  } else if (outcomemod == 'nonlinearinteraction') {
    Y <- X %*% beta + tau * Z + Z*sin(U) - U^2 + U*Z*(X[,2]+1) + 0.3*Z*X[,1]^2 + 
      rnorm(nrow(X), 0, noise_sd)
  }
  Y <- as.numeric(Y)
  return(Y)
}

fit_method <- function(Y,
                       X,
                       Z,
                       #V = NULL,
                       Vre = NULL,
                       Vcar = NULL,
                       Vgp = NULL,
                       Sigmainvre = NULL,
                       Sigmainvcar = NULL,
                       Sigmainvgp = NULL,
                       tols = NULL,
                       Xtol = 0.1,
                       method = c('OLS', 'RE', 'CAR', 'GP', 'SW', 'UW'),
                       neigen = 50) {
  
  # Y is the outcome vector
  # X is the design matrix including intercept
  # Z is the binary treatment vector
  # U is the unobserved confounder  
  # V is the augmenting eigenvector set
  # Sigmainvre is the inverse covariance matrix for RE
  # Sigmainvcar is the inverse covariance matrix for CAR
  # Sigmainvgp is the inverse covariance matrix for GP
  # method is the method fit
  
  design <- cbind(X, Z)
  method <- match.arg(method)
  results <- list()
  
  if ('OLS' == method) {
    # print(class(Y))
    # print(class(X))
    # print(class(Z))
    # print(typeof(Y))
    # print(typeof(X))
    # print(typeof(Z))
    ols_model <- lm(Y ~ X + Z)
    return(coefficients(ols_model)['Z'])
  }
  
  if ('RE' == method) {
    stopifnot(!is.null(Sigmainvre))
    return((solve(t(design) %*% Sigmainvre %*% design) %*% t(design) %*% Sigmainvre %*% Y)['Z',] )
  }
  
  if ('CAR'  == method) {
    stopifnot(!is.null(Sigmainvcar))
    return((solve(t(design) %*% Sigmainvcar %*% design) %*% t(design) %*% Sigmainvcar %*% Y)['Z',])
  }
  
  if ('GP'  == method) {
    stopifnot(!is.null(Sigmainvgp))
    return((solve(t(design) %*% Sigmainvgp %*% design) %*% t(design) %*% Sigmainvgp %*% Y)['Z',])
  }
  
  if ('SW'  == method) {
    stopifnot(!is.null(Vre))
    stopifnot(!is.null(Vcar))
    stopifnot(!is.null(Vgp))
    stopifnot(!is.null(tols))
    t_ind <- Z
    bal_cov <- cbind.data.frame(X,Vre[,1:neigen],Vcar[,1:neigen],Vgp[,1:neigen])
    # D <- model.matrix(~ 0 + X*V, data = bal_cov)  # this gives columns X, V, and X:V
    # # do a QR decomposition
    # qrD <- qr(D)
    # # extract the orthonormal Q
    # Q <- qr.Q(qrD)
    # # label the columns 
    # colnames(Q) <- colnames(D)
    # print(colnames(bal_cov))
    colnames(bal_cov) <- paste0('X', 1:ncol(bal_cov))
    data_frame <- as.data.frame(cbind(t_ind, bal_cov))
    t_ind <- "t_ind"
    bal <- list()
    #bal$bal_gri <- c(1e-4, 1e-3) # grid of tuning parameters 
    bal$bal_cov <- colnames(bal_cov)[-1]
    #bal$bal_alg = T # tuning algorithm in Wang and Zubizarreta (2020) used for automatically selecting the degree of approximate covariates balance.
    #bal$bal_sam = 1000
    bal$bal_std <- 'manual'
    bal$bal_alg <- F
    bal$bal_tol <- c(rep(Xtol, ncol(X)-1), tols$RE[1:neigen], 
                     tols$CAR[1:neigen], tols$GP[1:neigen])
    #print(bal$bal_tol)
    stopifnot(length(bal$bal_tol) == length(bal$bal_cov))
    sbwatttun_object = sbw(dat = data_frame, ind = t_ind, bal = bal, 
                           sol = list(sol_nam = "quadprog"), 
                           par = list(par_est = "att", par_tar = NULL))
    
    return(sum(sbwatttun_object$dat_weights$sbw_weights[Z == 1]*Y[Z == 1]) - 
             sum(sbwatttun_object$dat_weights$sbw_weights[Z == 0]*Y[Z == 0]) )
  }
  
  if ('UW'  == method) {
    t_ind <- Z
    bal_cov <- as.data.frame(X)
    #bal_cov <- model.matrix(~.^2, data = bal_cov)
    colnames(bal_cov) <- paste0('X', 1:ncol(bal_cov))
    data_frame <- as.data.frame(cbind(t_ind, bal_cov))
    t_ind <- "t_ind"
    bal <- list()
    bal$bal_cov <- colnames(bal_cov)[-1]
    print(bal$bal_cov)
    bal$bal_alg = TRUE # tuning algorithm in Wang and Zubizarreta (2020) used for automatically selecting the degree of approximate covariates balance.
    bal$bal_sam = 1000
    sbwatttun_object = sbw(dat = data_frame, ind = t_ind, bal = bal,
                           sol = list(sol_nam = "quadprog"), 
                           par = list(par_est = "att", par_tar = NULL))
    return(sum(sbwatttun_object$dat_weights$sbw_weights[Z == 1]*Y[Z == 1]) - 
             sum(sbwatttun_object$dat_weights$sbw_weights[Z == 0]*Y[Z == 0]))
  }
}

simfunc <- function(nsims = 1000,
                    X,
                    Z,
                    beta = c(-0.44,0.46,-0.69,-1.45,0.57,-1.02,-0.02,-0.94,1.10,-0.48,-0.71), 
                    #V = NULL,
                    Wre = NULL,
                    Wcar = NULL,
                    Wgp = NULL,
                    Sigmainvre = NULL,
                    Sigmainvcar = NULL,
                    Sigmainvgp = NULL,
                    Vre = NULL,
                    Ere = NULL,
                    Vcar = NULL,
                    Ecar = NULL,
                    Vgp = NULL,
                    Egp = NULL,
                    smoothing = c('clusterbased', 'adjacencybased', 'distancebased'),
                    outcomemod = c('linear', 'linearinteraction', 'nonlinearinteraction'),
                    methods = c('OLS', 'RE', 'CAR', 'GP', 'SW', 'UW')) {
  # nsims is the number of simulations
  # X is the design matrix including intercept
  # Z is the binary treatment vector
  # beta is the vector of coefficients on X, generated from a standard normal
  # V is the augmenting eigenvector set
  # Wre is the weight matrix to generate the RE-class of unmeasured conf
  # Wcar is the weight matrix to generate the CAR-class of unmeasured conf
  # Wgp is the weight matrix to generate the GP-class of unmeasured conf
  # Sigmainvre is the inverse covariance matrix to fit the RE model
  # Sigmainvcar is the inverse covariance matrix to fit the CAR model
  # Sigmainvgp is the inverse covariance matrix to fit the GP model
  # smoothing is the type of smoothing to use for U
  # outcomemod is the type of outcome model to use
  
  smoothing <- match.arg(smoothing)
  outcomemod <- match.arg(outcomemod)
  
  tols <- calculate_optimal_tolerances(X = X, 
                                       Z = Z, 
                                       Sigmainvre = Sigmainvre,
                                       Sigmainvcar = Sigmainvcar,
                                       Sigmainvgp = Sigmainvgp,
                                       Vre = Vre,
                                       Ere = Ere,
                                       Vcar = Vcar,
                                       Ecar = Ecar,
                                       Vgp = Vgp,
                                       Egp = Egp)
  
  for (method in methods) {
    # Create filename for csvs containing estimates
    filename <- paste0('results_May1/',
                       smoothing,
                       '_',
                       outcomemod,
                       '_',
                       method,
                       '.csv')
    
    # Create storage for estimates
    tauests <- rep(NA, nsims)
    
    for (sim in 1:nsims) {
      print(sim)
      # Generate unmeasured confounder
      U <- gen_U(
        Z,
        Wre = Wre,
        Wcar = Wcar,
        Wgp = Wgp,
        smoothing = smoothing
      )
      #print(cor(U,Z))
      # Generate outcome
      Y <- gen_Y(X, Z, U, beta = beta, outcomemod = outcomemod)
      tauests[sim] = fit_method(
        Y,
        X,
        Z,
        U,
        #V = V,
        Vre = Vre,
        Vcar = Vcar,
        Vgp = Vgp,
        tols = tols,
        Sigmainvre = Sigmainvre,
        Sigmainvcar = Sigmainvcar,
        Sigmainvgp = Sigmainvgp,
        method = method
      )
      
    }
    df <- tauests
    
    # write results to file
    # Check if file exists
    if (file.exists(filename)) {
      # write new sims to file as new columns
      olddf <- read.csv(filename)
      newdf <- cbind(olddf, tauests)
      write.csv(newdf, filename, row.names = FALSE)
    }
    # if file for estimates does not exist create it and write results
    else{
      write.csv(df, filename, row.names = FALSE)
    }
    
    invisible(filename)
  }
}

compute_ATT_MC <- function(mcreps = 1000, 
                           X, 
                           Z, 
                           beta = c(-0.44,0.46,-0.69,-1.45,0.57,-1.02,-0.02,-0.94,1.10,-0.48,-0.71),  
                           Wre = NULL,
                           Wcar = NULL,
                           Wgp = NULL,
                           smoothing = c('clusterbased', 'adjacencybased', 'distancebased'),
                           outcomemod = c('linear', 'linearinteraction', 'nonlinearinteraction')){
  # mcreps is the number of Monte Carlo repetitions
  # X is the design matrix including intercept
  # Z is the binary treatment vector
  # beta is the vector of coefficients on X, generated from a standard normal
  # Wre is the weight matrix to generate the RE-class of unmeasured conf
  # Wcar is the weight matrix to generate the CAR-class of unmeasured conf
  # Wgp is the weight matrix to generate the GP-class of unmeasured conf
  # smoothing is the type of smoothing to use for U
  # outcomemod is the type of outcome model to use
  
  tauests <- rep(NA, mcreps)
  for (rep in 1:mcreps){
    print(rep)
    # Generate unmeasured confounder
    U <- gen_U(
      Z,
      Wre = Wre,
      Wcar = Wcar,
      Wgp = Wgp,
      smoothing = smoothing
    )
    
    # Generate outcome
    Y1 <- gen_Y(X, rep(1, length(Z)), U, beta = beta, outcomemod = outcomemod)
    Y0 <- gen_Y(X, rep(0, length(Z)), U, beta = beta, outcomemod = outcomemod)
    # Fit the model
    tauests[rep] <- mean(Y1[Z == 1]) - mean(Y0[Z == 1])
  }
  print(summary(tauests))
  return(mean(tauests))
}

impliedweightsgeneral <- function(X, Z, Sigmainv, returnL = F){
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
  return((2*Z-1)*as.numeric(l))
}

calculate_optimal_tolerances <- function(X,
                                         Z,
                                         Sigmainvre = NULL,
                                         Sigmainvcar = NULL,
                                         Sigmainvgp = NULL,
                                         Vre = NULL,
                                         Ere = NULL,
                                         Vcar = NULL,
                                         Ecar = NULL,
                                         Vgp = NULL,
                                         Egp = NULL){
  n <- nrow(X)
  tols <- list()
  # Step 1: Calculate implied weights for RE
  lre <- impliedweightsgeneral(X, Z, Sigmainvre, returnL = T)
  # Step 2: calculate threshold for RE
  S <- Vre %*% diag(Ere) %*% t(Vre)
  lSl <- t(lre) %*% S %*% lre
  tols$RE <- sqrt(c(lSl)/(n*Ere))
  tols$RE[is.na(tols$RE)] <- Inf
  
  # Step 1: Calculate implied weights for CAR
  lcar <- impliedweightsgeneral(X, Z, Sigmainvcar, returnL = T)
  # Step 2: calculate threshold for CAR
  S <- Vcar %*% diag(Ecar) %*% t(Vcar)
  lSl <- t(lcar) %*% S %*% lcar
  tols$CAR <- sqrt(c(lSl)/(n*Ecar))
  tols$CAR[is.na(tols$CAR)] <- Inf
  
  # Step 1: Calculate implied weights for GP
  lgp <- impliedweightsgeneral(X, Z, Sigmainvgp, returnL = T)
  # Step 2: calculate threshold for GP
  S <- Vgp %*% diag(Egp) %*% t(Vgp)
  lSl <- t(lgp) %*% S %*% lgp
  tols$GP <- sqrt(c(lSl)/(n*Egp))
  tols$GP[is.na(tols$GP)] <- Inf
  
  return(tols)
}

# Bootstrap function
boot_func <- function(Y,
                      X,
                      Z,
                      bootreps = 1000,
                      Vre = NULL,
                      Vcar = NULL,
                      Vgp = NULL,
                      Sigmainvre = NULL,
                      Sigmainvcar = NULL,
                      Sigmainvgp = NULL,
                      tols = NULL,
                      method = c('OLS', 'RE', 'CAR', 'GP', 'SW', 'UW'),
                      neigen = 50){
  method <- match.arg(method)
  n <- length(Y)
  idxs <- seq_len(n)
  atts <- numeric(bootreps)
  
  b <- 1
  while(b <= bootreps){
    boot_idxs <- sample(idxs, n, replace = TRUE)
    
    Yb   <- Y[boot_idxs]
    Xb   <- X[boot_idxs, , drop = FALSE]
    Zb   <- Z[boot_idxs]
    Vreb <- if(!is.null(Vre))         Vre[boot_idxs, , drop = FALSE] else NULL
    Vcarb <- if(!is.null(Vcar))       Vcar[boot_idxs, , drop = FALSE] else NULL
    Vgpb <- if(!is.null(Vgp))         Vgp[boot_idxs, , drop = FALSE] else NULL
    Sigmainvreb <- if(!is.null(Sigmainvre)) Sigmainvre[boot_idxs, boot_idxs] else NULL
    Sigmainvcarb <- if(!is.null(Sigmainvcar)) Sigmainvcar[boot_idxs, boot_idxs] else NULL
    Sigmainvgpb <- if(!is.null(Sigmainvgp)) Sigmainvgp[boot_idxs, boot_idxs] else NULL
    
    # try to fit; on error, give a warning and retry (does not advance b)
    res <- tryCatch({
      fit_method(Y           = Yb,
                 X           = Xb,
                 Z           = Zb,
                 Vre         = Vreb,
                 Vcar        = Vcarb,
                 Vgp         = Vgpb,
                 Sigmainvre  = Sigmainvreb,
                 Sigmainvcar = Sigmainvcarb,
                 Sigmainvgp  = Sigmainvgpb,
                 tols        = tols,
                 method      = method,
                 neigen      = neigen)
    }, error = function(e) {
      warning(sprintf("bootstrap sample %d failed: %s — skipping", b, e$message))
      NULL
    })
    
    if (!is.null(res) && length(res) == 1 && !is.na(res)) {
      atts[b] <- res
      b <- b + 1
    }
    # else: res is NULL or NA, so we don't increment b and we try again
    print(b)
  }
  
  return(atts)
}

