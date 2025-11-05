GibbsSampler_mh <- function(df, #functional data
                            a0, #a sigma prior
                            b0, #b sigma prior
                            avec0, #mean for alpha prior (should't be here)
                            Sigma0,  #var/cov for alpha prior
                            V0, #column variance for A prior
                            Aprior, #mean for A prior
                            R, #number of iterations
                            burnin, #burnin lenght
                            nbasis, #number of basis to decompose error component 
                            nu0, #df of Psi prior
                            S0 #scale matrix for Psi prior
){
  #basic quantities
  x <- df[, 1]
  fs <- df[, -1]
  T <- ncol(fs) 
  n <- length(x)
  k <- length(avec0)
  n_acc <- 0
  
  #B and C matrix, fixed matrices
  B <- bs(x, df = length(avec0), degree = 3)
  C <- solve(crossprod(B)) %*% t(B)
  D <- make_fourier(x, nbasis)
  
  #Initialization 
  
  sigma0 <- 1/rgamma(1, a0, b0)
  
  Psi0 <- riwish(nu0, S0)
  Psi0 <- make_posdef(Psi0)
  
  Phi0 <- D %*% Psi0 %*% t(D)
  Phi0 <- make_posdef(Phi0)
  
  ybar <- rowMeans(fs)   #use training data instead of all fs
  BtB <- crossprod(B)
  Bty <- crossprod(B, ybar)
  
  alpha0 <- solve(BtB, Bty)
  
  A0 <- rmatnorm(1, Aprior, V0, C%*%Phi0%*%t(C))    
  AO <- make_stable_from_A(A0)
  
  gamma0 <- rnorm(k, 0, 1)
  Gamma0 <- matrix(0, T, k)
  Gamma0[1,] <- gamma0
  
  for(t in 2:T){
    Gamma0[t, ] <- t(A0) %*% Gamma0[t-1, ]
  }
  
  #Chains creation 
  sigma <- numeric(R)
  atest <- numeric(R)
  btest <- numeric(R)
  alpha <- matrix(0, nrow = k, ncol = R)
  A <- array(NA, dim = c(k, k, R))
  Phi <- array(NA, dim = c(n, n, R))
  Psi <- array(NA, dim = c(2*nbasis, 2*nbasis, R))
  Gamma <- array(NA, dim = c(T, k, R))
  
  sigma[1] <- sigma0
  atest[1] <- a0
  btest[1] <- b0
  alpha[, 1] <- alpha0
  A[, , 1] <- A0
  Phi[, , 1] <- Phi0
  Psi[, , 1] <- Psi0
  Gamma[, , 1] <- Gamma0
  
  #Output structure
  out <- vector("list", R)
  accept_count <- 0
  
  
  for(i in 2:R){
    mus <- numeric(n)
    gs <- matrix(0, n, T)
    
    mus <- B %*% alpha[,i-1]
    mus <- matrix(rep(mus, T), n, T)
    for(t in 1:T){
      gs[,t] <- B %*% Gamma[t, , i-1]
    }
    
    outsigma <- sample_sigma(fs, a0, b0, Phi[, ,i-1], T, n, mus, gs)
    sigma[i] <- outsigma$sigma
    atest[i] <- outsigma$a
    btest[i] <- outsigma$b
    
    alpha[, i] <- sample_alpha(fs, avec0, Sigma0, Phi[, ,i-1], sigma[i], A[, , i-1], k, T, n, B, Gamma[, , i-1])
    
    Gamma[1, , i] <- rnorm(k, 0, 1)
    for(t in 2:T){
      Gamma[t, , i] <- t(A[, , i-1])%*%Gamma[t-1, , i]
    }
    A_draw <- sample_A(fs, Phi[, , i-1], V0, A0, sigma[i], Gamma[, , i], C, mus, k, T) 
    A[, , i] <- as.matrix(make_stable_from_A(A_draw))
    out_Phi <- sample_Phi_fourier_mh(fs, mus, gs, sigma[i], D, Psi[, ,i-1], Phi[, , i-1], nu0,  
                                     S0, n, n_acc)
    Phi[, , i] <- out_Phi$Phi
    Psi[, , i] <- out_Phi$Psi
    n_acc <- out_Phi$n_acc
    
  }
  out <- list(
    sigma = sigma[-c(1:burnin)], 
    alpha = alpha[ ,-c(1:burnin)], 
    A = A[, , -c(1:burnin)],
    Gamma = Gamma[, , -c(1:burnin)],
    Phi = Phi[, , -c(1:burnin)],
    a = atest[-c(1:burnin)], 
    b = btest[-c(1:burnin)], 
    n_acc = n_acc
  )
  return(out)
  
}


GibbsSampler_mtmh <- function(df, #functional data
                              a0, #a sigma prior
                              b0, #b sigma prior
                              avec0, #mean for alpha prior (should't be here)
                              Sigma0,  #var/cov for alpha prior
                              V0, #column variance for A prior
                              Aprior, #mean for A prior
                              R, #number of iterations
                              burnin, #burnin lenght
                              nbasis, #number of basis to decompose error component 
                              nu0, #df of Psi prior
                              S0, #scale matrix for Psi prior
                              m  #number of tries for Psi at every iteration
){
  #basic quantities
  x <- df[, 1]
  fs <- df[, -1]
  T <- ncol(fs) 
  n <- length(x)
  k <- length(avec0)
  n_acc <- 0
  
  #B and C matrix, fixed matrices
  B <- bs(x, df = length(avec0), degree = 3)
  C <- solve(crossprod(B)) %*% t(B)
  D <- make_fourier(x, nbasis)
  
  #Initialization 
  
  sigma0 <- 1/rgamma(1, a0, b0)
  
  Psi0 <- riwish(nu0, S0)
  Psi0 <- make_posdef(Psi0)
  
  Phi0 <- D %*% Psi0 %*% t(D)
  Phi0 <- make_posdef(Phi0)
  
  ybar <- rowMeans(fs)   #use training data instead of all fs
  BtB <- crossprod(B)
  Bty <- crossprod(B, ybar)
  
  alpha0 <- solve(BtB, Bty)
  
  A0 <- rmatnorm(1, Aprior, V0, C%*%Phi0%*%t(C))    
  A0 <- make_stable_from_A(A0)
  
  gamma0 <- rnorm(k, 0, 1)
  Gamma0 <- matrix(0, T, k)
  Gamma0[1,] <- gamma0
  
  for(t in 2:T){
    Gamma0[t, ] <- t(A0) %*% Gamma0[t-1, ]
  }
  
  #Chains creation 
  sigma <- numeric(R)
  atest <- numeric(R)
  btest <- numeric(R)
  alpha <- matrix(0, nrow = k, ncol = R)
  A <- array(NA, dim = c(k, k, R))
  Phi <- array(NA, dim = c(n, n, R))
  Psi <- array(NA, dim = c(2*nbasis, 2*nbasis, R))
  Gamma <- array(NA, dim = c(T, k, R))
  
  sigma[1] <- sigma0
  atest[1] <- a0
  btest[1] <- b0
  alpha[, 1] <- alpha0
  A[, , 1] <- A0
  Phi[, , 1] <- Phi0
  Psi[, , 1] <- Psi0
  Gamma[, , 1] <- Gamma0
  
  #Output structure
  out <- vector("list", R)
  accept_count <- 0
  accept_vec <- numeric(R)
  
  
  for(i in 2:R){
    mus <- numeric(n)
    gs <- matrix(0, n, T)
    
    mus <- B %*% alpha[,i-1]
    mus <- matrix(rep(mus, T), n, T)
    for(t in 1:T){
      gs[,t] <- B %*% Gamma[t, , i-1]
    }
    
    outsigma <- sample_sigma(fs, a0, b0, Phi[, ,i-1], T, n, mus, gs)
    sigma[i] <- outsigma$sigma
    atest[i] <- outsigma$a
    btest[i] <- outsigma$b
    
    alpha[, i] <- sample_alpha(fs, avec0, Sigma0, Phi[, ,i-1], sigma[i], A[, , i-1], k, T, n, B, Gamma[, , i-1])
    
    Gamma[1, , i] <- rnorm(k, 0, 1)
    for(t in 2:T){
      Gamma[t, , i] <- t(A[, , i-1])%*%Gamma[t-1, , i]
    }
    A_draw <- sample_A(fs, Phi[, , i-1], V0, A0, sigma[i], Gamma[, , i], C, mus, k, T) 
    A[, , i] <- as.matrix(make_stable_from_A(A_draw))
    out_Phi <- sample_Phi_fourier_mtmh(fs, mus, gs, sigma[i], D, Psi[, ,i-1], 
                                       nu0, S0, m, n_acc)
    Phi[, , i] <- out_Phi$Phi
    Psi[, , i] <- out_Phi$Psi
    n_acc <- out_Phi$n_acc
    accept_vec[i] <- out_Phi$acc
  }
  out <- list(
    sigma = sigma[-c(1:burnin)], 
    alpha = alpha[ ,-c(1:burnin)], 
    A = A[, , -c(1:burnin)],
    Gamma = Gamma[, , -c(1:burnin)],
    Phi = Phi[, , -c(1:burnin)],
    a = atest[-c(1:burnin)], 
    b = btest[-c(1:burnin)], 
    n_acc = n_acc
  )
  return(out)
  
}


GS_mtmh_p <- function(df, #functional data
                      a0, #a sigma prior
                      b0, #b sigma prior
                      avec0, #mean for alpha prior (should't be here)
                      Sigma0,  #var/cov for alpha prior
                      V0, #column variance for A prior
                      Aprior, #mean for A prior
                      p, #number of autoregressive orders
                      R, #number of iterations
                      burnin, #burnin lenght
                      nbasis, #number of basis to decompose error component 
                      nu0, #df of Psi prior
                      S0, #scale matrix for Psi prior
                      m  #number of tries for Psi at every iteration
){
  #basic quantities
  x <- df[, 1]
  fs <- df[, -1]
  T <- ncol(fs) 
  n <- length(x)
  k <- length(avec0)
  n_acc <- 0
  
  #B and C, fixed matrices
  B <- bs(x, df = length(avec0), degree = 3)
  C <- solve(crossprod(B)) %*% t(B)
  D <- make_fourier(x, nbasis)
  
  #Initialization 
  
  sigma0 <- 1/rgamma(1, a0, b0)
  
  Psi0 <- riwish(nu0, S0)
  Psi0 <- make_posdef(Psi0)
  
  Phi0 <- D %*% Psi0 %*% t(D)
  Phi0 <- make_posdef(Phi0)
  
  ybar <- rowMeans(fs)   #use training data instead of all fs
  BtB <- crossprod(B)
  Bty <- crossprod(B, ybar)
  
  alpha0 <- solve(BtB, Bty)
  
  A0 <- rmatnorm(p, Aprior, V0, C%*%Phi0%*%t(C))    
  AO <- make_stable_from_A(A0)
  
  gamma0 <- rnorm(k, 0, 1)
  Gamma0 <- matrix(0, T, k)
  Gamma0[1,] <- gamma0
  
  for(t in 2:T){
    Gamma0[t, ] <- t(A0) %*% Gamma0[t-1, ]
  }
  
  #Chains creation 
  sigma <- numeric(R)
  atest <- numeric(R)
  btest <- numeric(R)
  alpha <- matrix(0, nrow = k, ncol = R)
  A <- array(NA, dim = c(k, k, p, R))  
  Phi <- array(NA, dim = c(n, n, R))
  Psi <- array(NA, dim = c(2*nbasis, 2*nbasis, R))
  Gamma <- array(NA, dim = c(T, k, R))
  
  sigma[1] <- sigma0
  atest[1] <- a0
  btest[1] <- b0
  alpha[, 1] <- alpha0
  for (r in 1:p) {
    A0 <- rmatnorm(1, Aprior, V0, C %*% Phi0 %*% t(C))
    A[, , r, 1] <- make_stable_from_A(A0)
  }
  Phi[, , 1] <- Phi0
  Psi[, , 1] <- Psi0
  Gamma0 <- matrix(0, T, k)
  Gamma0[1:p, ] <- matrix(rnorm(k * p), nrow = p, ncol = k)
  
  for (t in (p+1):T) {
    Gamma0[t, ] <- 0
    for (r in 1:p) {
      Gamma0[t, ] <- Gamma0[t, ] + t(A[, , r, 1]) %*% Gamma0[t - r, ]
    }
  }  
  
  Gamma[, , 1] <- Gamma0
  
  #Output structure
  out <- vector("list", R)
  accept_count <- 0
  accept_vec <- numeric(R)
  
  
  for(i in 2:R){
    mus <- numeric(n)
    gs <- matrix(0, n, T)
    
    mus <- B %*% alpha[,i-1]
    mus <- matrix(rep(mus, T), n, T)
    for(t in 1:T){
      gs[,t] <- B %*% Gamma[t, , i-1]
    }
    
    outsigma <- sample_sigma(fs, a0, b0, Phi[, ,i-1], T, n, mus, gs)
    sigma[i] <- outsigma$sigma
    atest[i] <- outsigma$a
    btest[i] <- outsigma$b
    
    alpha[, i] <- sample_alpha(fs, avec0, Sigma0, Phi[, ,i-1], sigma[i], A[, , i-1], k, T, n, B, Gamma[, , i-1])
    
    Gamma[1:p, , i] <- matrix(rnorm(k * p), nrow = p, ncol = k)
    for (t in (p+1):T) {
      Gamma[t, , i] <- 0
      for (r in 1:p) {
        Gamma[t, , i] <- Gamma[t, , i] + t(A[, , r, i-1]) %*% Gamma[t - r, , i]
      }
    }
    
    for (r in 1:p) {
      A_draw <- sample_A_multi(fs, Phi[, , i-1], V0, Aprior, sigma[i],
                               Gamma[, , i], C, mus, k, T,
                               A[, , , i-1], lag = r)
      A[, , r, i] <- make_stable_from_A(A_draw)
    }
    
    out_Phi <- sample_Phi_fourier_mtmh(fs, mus, gs, sigma[i], D, Psi[, ,i-1], 
                                       nu0, S0, m, n_acc)
    Phi[, , i] <- out_Phi$Phi
    Psi[, , i] <- out_Phi$Psi
    n_acc <- out_Phi$n_acc
    accept_vec[i] <- out_Phi$acc
  }
  out <- list(
    sigma = sigma[-c(1:burnin)],
    alpha = alpha[, -c(1:burnin)],
    A = A[, , , -c(1:burnin)],
    Gamma = Gamma[, , -c(1:burnin)],
    Phi = Phi[, , -c(1:burnin)],
    a = atest[-c(1:burnin)],
    b = btest[-c(1:burnin)],
    n_acc = n_acc, 
    accept_vec = accept_vec
  )
  return(out)
  
}
