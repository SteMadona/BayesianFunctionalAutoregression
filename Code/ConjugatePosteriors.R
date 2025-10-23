sample_sigma <- function(fs, a0, b0, Phi, T, n, mus, gs){
  Phi_inv <- solve(make_posdef(Phi))
  resid_sum <- 0
  for(t in 1:T){
    e <- fs[, t] - mus[, t] - gs[, t]          # n x 1
    resid_sum <- resid_sum + as.numeric(t(e) %*% Phi_inv %*% e)
  }
  asigma <- a0 + (n*T)/2
  bsigma <- b0 + 0.5*resid_sum
  if(!is.finite(bsigma) || bsigma <= 0) bsigma <- b0 + 1e-6
  list(sigma = 1/rgamma(1, asigma, bsigma), a = asigma, b = bsigma)
}

sample_alpha <- function(fs, avec0, Sigma0, Phi, sigma, A, k, T, n, B, Gamma){
  
  Gamma <- as.matrix(Gamma)
  Phi_inv <- solve(make_posdef(Phi))
  
  SigmaT <- solve(solve(Sigma0) + (T/sigma)*t(B)%*%solve(Phi)%*%B)
  SigmaT <- (SigmaT + t(SigmaT))/2 #to enforce simmetry, if it's not the case
  
  gs <- B %*% t(Gamma) 
  
  D <- numeric(k)
  Phi_inv <- solve(Phi)
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
  VT <- ginv(M)     #using the pseudo inverse
  #  VT <- chol2inv(chol(nearPD(M, corr = FALSE)$mat))  #should be fine, but high costly computation
  #  VT <- chol2inv(chol(make_posdef(M)))
  #  VT <- solve(solve(V0) + (1/sigma)*(t(Gamma0)%*%(Gamma0)))
  #  VT <- make_posdef(VT)
  #  VT <- (VT + t(VT)) / 2   
  
  AT <- VT %*% (solve(V0)%*%A0 + (1/sigma)*(t(Gamma0)%*%Y%*%t(C))) 
  
  U <- C %*% Phi %*% t(C)
  U <- make_posdef(U)
  U <- (U + t(U)) / 2    #enforce simmetry numerically
  
  rmatnorm(1, AT, VT, U)
}
