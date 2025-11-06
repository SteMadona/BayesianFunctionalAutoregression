sample_sigma <- function(fs, a0, b0, Phi, T, n, mus, gs){
  Phi_pd <- add_nugget(Phi)
  Phi_inv <- chol2inv(chol(Phi_pd))
  
  resid_sum <- 0
  for(t in 1:T){
    e <- fs[, t] - mus[, t] - gs[, t]          # n x 1
    resid_sum <- resid_sum + as.numeric(t(e) %*% Phi_inv %*% e)
  }
  asigma <- a0 + (n*T)/2
  bsigma <- b0 + 0.5*resid_sum
  if(!is.finite(bsigma) || bsigma <= 0) bsigma <- b0 + 1e-6
  list(sigma = rinvgamma(1, shape = asigma, scale = bsigma), a = asigma, b = bsigma)
}

sample_alpha <- function(fs, avec0, Sigma0, Phi, sigma, k, T, n, B, Gamma){
  
  Gamma <- as.matrix(Gamma)
  Phi_pd <- add_nugget(Phi)
  Phi_inv <- chol2inv(chol(Phi_pd))
  
  SigmaT <- solve(solve(Sigma0) + (T/sigma)*t(B)%*%Phi_inv%*%B)
  SigmaT <- (SigmaT + t(SigmaT))/2 #to enforce simmetry, if it's not the case
  
  gs <- B %*% t(Gamma) 
  
  D <- numeric(k)
  for(i in 1:T){
    D <- D + t(B) %*% Phi_inv %*% (fs[, i] - gs[, i])
  }
  avecT <- SigmaT %*%(solve(Sigma0)%*%avec0 + (1/sigma)*D)
  
  as.vector(mvtnorm::rmvnorm(1, mean = avecT, sigma = SigmaT))
}

sample_A <- function(fs, Phi, V0, A0, sigma, Gamma0, C, mus, k, T){
  
  Mu <- matrix(rep(mus, T), ncol = T)
  
  Y <- t(fs - mus)
  
  Gamma2 <- make_posdef(t(Gamma0)%*%Gamma0) + diag(1e-6, k)
  M <- solve(V0) + (1/sigma)*Gamma2
  VT <- add_nugget(chol2inv(chol(M)))      
  
  AT <- VT %*% (solve(V0)%*%A0 + (1/sigma)*(t(Gamma0)%*%Y%*%t(C))) 
  
  U <- C %*% Phi %*% t(C)
  U <- add_nugget(U)
  U <- (U + t(U)) / 2    #enforce simmetry numerically
  
  rmatnorm(1, AT, VT, U)
}

sample_A_multi <- function(fs, Phi, V0, Aprior, sigma, Gamma, C, mus, k, T, A_all, lag){
  stopifnot(lag >= 1, lag < T)         # at least one obs after the lag
  stopifnot(ncol(Gamma) == k)
  
  Y <- fs - mus  
  
  p <- dim(A_all)[3]
  for (t in (lag+1):T) {
    temp <- 0
    for (j in 1:p) {
      if (j != lag && (t-j) > 0) {
        temp <- temp + t(C) %*% A_all[, , j] %*% Gamma[t - j, ]
      }
    }
    Y[, t] <- Y[, t] - temp
  }
  #removes the contribution of all other autoregressive lags
  
  Y_t <- t(Y[, (lag+1):T])          
  Gamma_tlag <- Gamma[1:(T-lag), ]  
  
  Gamma2 <- crossprod(Gamma_tlag)   
  Gamma2 <- make_posdef(Gamma2) + diag(1e-6, k)
  
  M  <- solve(V0) + (1 / sigma) * Gamma2
  VT <- tryCatch(chol2inv(chol(M)), error = function(e) MASS::ginv(M))
  VT <- make_posdef(VT) + diag(1e-6, k)
  
  AT <- VT %*% (solve(V0) %*% Aprior + (1 / sigma) * (t(Gamma_tlag) %*% Y_t %*% t(C))) 
  
  U <- C %*% Phi %*% t(C)
  U <- make_posdef(U) + diag(1e-6, k)
  
  rmatnorm(1, AT, VT, U)
}



