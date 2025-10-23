make_posdef <- function(M, tol = 1e-8){
  Msym <- (M + t(M)) / 2  #enforce symmetry
  # first try small jitter
  jitter <- tol * max(1, max(abs(diag(Msym))))
  Mj <- Msym + diag(jitter, nrow(Msym))
  np <- tryCatch(nearPD(Mj, corr = FALSE, doSym = TRUE), error = function(e) NULL)
  if(!is.null(np)){
    Mpd <- as.matrix(np$mat)
  } else {
    # fallback: eigen adjustment
    e <- eigen(Msym, symmetric = TRUE)
    vals <- pmax(e$values, tol)
    Mpd <- e$vectors %*% diag(vals) %*% t(e$vectors)
    Mpd <- (Mpd + t(Mpd))/2
  }
  return(Mpd)
}

make_stable_from_A <- function(A, target_rho = 0.95){
  ev <- eigen(A)$values
  rho <- max(abs(ev))
  if(rho == 0) return(A)
  A_stable <- A * (target_rho / rho)
  return(A_stable)
}

logdet_spd <- function(M){
  d <- determinant(M, logarithm = TRUE)
  as.numeric(d$modulus)  # assume SPD; sign should be +1
}

logsumexp <- function(a) {
  m <- max(a)
  m + log(sum(exp(a - m)))
}

softmax_from_logs <- function(a) {
  a_shift <- a - max(a)
  exp(a_shift) / sum(exp(a_shift))
}

log_multigamma <- function(a, p){
  0.25 * p * (p-1)*log(pi) + sum(lgamma(a + (1 - (1:p))/2))
}

log_diwish <- function(X, nu, Psi){
  p <- nrow(Psi)
  
  term1 <- 0.5*nu*log(det(Psi)) - 0.5*nu*p*log(2)
  term2 <- log_multigamma(nu/2, p)
  term3 <- ((nu + p + 1)/2)*log(det(X))
  term4 <- 0.5*sum(diag(Psi %*% solve(X)))
  
  out <- as.numeric(term1 - term2 - term3 - term4)
  
  return(out)
}

make_fourier <- function(x, K) {
  L <- max(x) - min(x)
  xb <- (x - min(x)) / L  # scala su [0,1]
  B <- matrix(NA, nrow=length(x), ncol=2*K)
  for (k in 1:K) {
    B[, 2*k-1] <- sin(2*pi*k*xb)
    B[, 2*k]   <- cos(2*pi*k*xb)
  }
  return(B) # m x (2K)
}
