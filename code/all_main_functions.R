# ==============================================================================
# 1. CONFIGURATION AND LIBRARIES
# ==============================================================================
library(splines)
library(ggplot2)
library(dplyr)
library(tidyr)
library(mvtnorm)
library(matrixNormal)
library(Matrix)
library(MASS)
library(MCMCpack)
library(tictoc)
library(reshape2)
library(plot3D)

# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

#function for ensuring positive definiteness.

make_posdef <- function(M, tol = 1e-8){
  Msym <- (M + t(M)) / 2    
  e <- eigen(Msym, symmetric = TRUE)
  vals <- pmax(e$values, tol)
  Mpd <- e$vectors %*% diag(vals) %*% t(e$vectors)
  return((Mpd + t(Mpd))/2)
}

#add noise to the diagonal and then normalize the scale of the matrix, in order 
#to avoid inversion problems. 

add_nugget <- function(Phi, rel = 1e-6, tau = NULL, tol = 1e-10) {
  Phi <- make_posdef(Phi, tol = tol)
  n <- nrow(Phi)
  
  # se tau non è fornito, lo scegliamo proporzionale alla scala media della diagonale
  if (is.null(tau)) {
    scale_diag <- mean(diag(Phi))
    if (!is.finite(scale_diag) || scale_diag <= 0) scale_diag <- 1
    tau <- rel * scale_diag
  }
  
  Phi2 <- Phi + diag(tau, n)
  make_posdef(Phi2, tol = tol)
}


# function to create stabilized version of Phi
make_Phi <- function(D, Psi, rel = 1e-6, tau = NULL, tol = 1e-10) {
  Psi_pd <- make_posdef(Psi, tol = tol)
  Phi_raw <- D %*% Psi_pd %*% t(D)
  
  add_nugget(Phi_raw, rel = rel, tau = tau, tol = tol)
}



#creation of fourier basis, used in error matrix decomposition

make_fourier <- function(x, K) {
  L <- max(x) - min(x); xb <- (x - min(x)) / L
  B <- matrix(NA, length(x), 2*K)
  for (k in 1:K) { B[, 2*k-1] <- sin(2*pi*k*xb); B[, 2*k] <- cos(2*pi*k*xb) }
  return(B)
}


eig <- function(Amat){max(abs(eigen(Amat)$value))}


robust_invert <- function(M, jitter = 1e-5){
  if(any(!is.finite(M))) return(diag(1, nrow(M))) 
  
  tryCatch(chol2inv(chol(M)), error = function(e){
    tryCatch(solve(M), error = function(e2){
      M_jit <- M + diag(jitter * max(abs(diag(M))), nrow(M))
      tryCatch(solve(M_jit), error = function(e3){
        MASS::ginv(M) 
      })
    })
  })
}


#add controls in cases with extremely high values
stabilize_list <- function(A_list, target=0.99){
  rho <- get_rho(A_list)
  if(!is.finite(rho) || rho > 1e5) { 
    k <- nrow(A_list[[1]])
    return(lapply(A_list, function(x) diag(0.01, k)))
  }
  if(rho < target) return(A_list)
  factor <- target / rho
  return(lapply(A_list, function(x) x * factor))
}

#computing the root of the characteristic polynomial using companion matrix
#to identify non stationary lists of matrices

get_rho <- function(A_list) {
  p <- length(A_list); k <- nrow(A_list[[1]])
  C_mat <- matrix(0, k*p, k*p)
  for(j in 1:p) C_mat[1:k, ((j-1)*k + 1):(j*k)] <- A_list[[j]]
  if(p > 1) C_mat[cbind((k+1):(k*p), 1:(k*(p-1)))] <- 1 
  
  ev <- tryCatch(eigen(C_mat, only.values=TRUE)$values, error=function(e) 100)
  return(max(Mod(ev)))
}


# computing model's likelihood 

loglik <- function(fs, gs, mu, sigma_sq, Psi, D){
  TT <- ncol(fs)
  n  <- nrow(fs)
  
  Phi <- make_Phi(D, Psi, rel = 1e-6)
  R   <- chol(Phi)                           # Cholesky
  Phi_inv <- chol2inv(R)
  logdet  <- 2 * sum(log(diag(R)))           # log|Phi|
  
  quad_term <- 0
  for(t in 1:TT){
    e <- fs[,t] - mu - gs[,t]
    quad_term <- quad_term + as.numeric(t(e) %*% Phi_inv %*% e)
  }
  
  ll <- -0.5 * (quad_term / sigma_sq) -
    0.5 * TT * logdet -
    0.5 * n * TT * log(sigma_sq)
  
  as.numeric(ll)
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

logsumexp <- function(a) {
  m <- max(a)
  m + log(sum(exp(a - m)))
}

softmax_from_logs <- function(a) {
  a_shift <- a - max(a)
  exp(a_shift) / sum(exp(a_shift))
}




# ==============================================================================
# 3. POSTERIOR SAMPLING FUNCTIONS
# ==============================================================================

rw_step_chol <- function(Psi, rw_scale = 0.05){
  L <- t(chol(Psi))
  
  eps <- matrix(rnorm(length(L), mean = 0, sd = rw_scale), nrow = nrow(L))
  eps[upper.tri(eps)] <- 0
  
  L_star <- L + eps
  
  if (any(diag(L_star) <= 0) || any(!is.finite(diag(L_star)))) {
    return(list(Psi_star = NULL, L = L, L_star = L_star, valid = FALSE))
  }
  
  Psi_star <- make_posdef(L_star %*% t(L_star))
  list(Psi_star = Psi_star, L = L, L_star = L_star, valid = TRUE)
}


sample_Phi_fourier_mtmh <- function(fs, gs, mu, sigma_sq, D, Psi, nu0, S0, m,
                                    nu_q = NULL, s = 1, rw_prob = 1, rw_scale = 0.05) {
  
  p <- nrow(Psi)
  if (is.null(nu_q)) nu_q <- p + 4
  
  # log Jacobiano (a costante additiva vicino che cancella nel ratio)
  
  logJ <- function(L) {
    w <- p:1
    d <- diag(L)
    if (any(d <= 0) || any(!is.finite(d))) return(-Inf)
    sum(w * log(d))
  }
  
  if (runif(1) < rw_prob) {
    # RW step on Cholesky decomp 
    prop <- rw_step_chol(Psi, rw_scale = rw_scale)
    
    if (!isTRUE(prop$valid) || is.null(prop$Psi_star)) {
      return(list(
        Phi = make_posdef(D %*% Psi %*% t(D)),
        Psi = Psi,
        acc = 0
      ))
    }
    
    Psi_star <- prop$Psi_star
    L <- prop$L
    L_star <- prop$L_star
    
    logpost_star <- loglik(fs, gs, mu, sigma_sq, Psi_star, D) + log_diwish(Psi_star, nu0, S0)
    logpost_curr <- loglik(fs, gs, mu, sigma_sq, Psi,      D) + log_diwish(Psi,      nu0, S0)
    
    # Jacobian correction
    log_alpha <- (logpost_star - logpost_curr) + (logJ(L_star) - logJ(L))
    
    if (!is.finite(log_alpha)) {
      return(list(
        Phi = make_posdef(D %*% Psi %*% t(D)),
        Psi = Psi,
        acc = 0
      ))
    }
    
    
    if (runif(1) < min(1, exp(log_alpha))) {
      return(list(
        Phi = add_nugget(D %*% Psi_star %*% t(D), rel = 1e-3),
        Psi = Psi_star,
        acc = 1
      ))
    } else {
      return(list(
        Phi = make_posdef(D %*% Psi %*% t(D)),
        Psi = Psi,
        acc = 0
      ))
    }
  }
  
  # Independence MTMH branch
  S_q <- s * make_posdef(Psi)
  
  log_target <- function(P) {
    loglik(fs, gs, mu, sigma_sq, P, D) + log_diwish(P, nu0, S0)
  }
  
  Psi_list <- replicate(m, make_posdef(riwish(nu_q, S_q)), simplify = FALSE)
  logw_f <- sapply(Psi_list, function(P) log_target(P) - log_diwish(P, nu_q, S_q))
  
  idx <- sample.int(m, 1, prob = exp(logw_f - max(logw_f)))
  Psi_star <- Psi_list[[idx]]
  
  Psi_back <- c(replicate(m - 1, make_posdef(riwish(nu_q, S_q)), simplify = FALSE), list(Psi))
  logw_b <- sapply(Psi_back, function(P) log_target(P) - log_diwish(P, nu_q, S_q))
  
  log_alpha <- logsumexp(logw_f) - logsumexp(logw_b)
  alpha <- min(1, exp(log_alpha))
  
  if (runif(1) < alpha) {
    list(Phi = make_posdef(D %*% Psi_star %*% t(D)),
         Psi = Psi_star,
         acc = 1)
  } else {
    list(Phi = make_posdef(D %*% Psi %*% t(D)),
         Psi = Psi,
         acc = 0)
  }
}



sample_sigma <- function(fs, a0, b0, Phi, TT, n, gs, mu){
  Phi_pd <- add_nugget(Phi)
  Phi_inv <- chol2inv(chol(Phi_pd))
  resid_sum <- 0
  for(t in 1:TT){
    e <- fs[, t] - mu - gs[, t]
    resid_sum <- resid_sum + as.numeric(t(e) %*% Phi_inv %*% e)
  }
  rinvgamma(1, shape = a0 + (n*TT)/2, scale = b0 + 0.5*resid_sum)
}


sample_alpha <- function(fs, avec0, Sigma0, Phi, sigma, k, TT, n, B, gs){
  
  Phi_pd <- add_nugget(Phi)
  Phi_inv <- chol2inv(chol(Phi_pd))
  
  SigmaT <- solve(solve(Sigma0) + (TT/sigma)*t(B)%*%Phi_inv%*%B)
  SigmaT <- (SigmaT + t(SigmaT))/2 #to enforce simmetry, if it's not the case
  
  D <- numeric(k)
  for(i in 1:TT){
    D <- D + t(B) %*% Phi_inv %*% (fs[, i] - gs[, i])
  }
  avecT <- SigmaT %*%(solve(Sigma0)%*%avec0 + (1/sigma)*D)
  
  as.vector(mvtnorm::rmvnorm(1, mean = avecT, sigma = SigmaT))
}

#reconstruct Gamma deterministically
reconstruct_Gamma_det <- function(A_list, Gamma_init, TT) {
  p_ar <- length(A_list)
  k <- ncol(A_list[[1]])
  
  Gamma <- matrix(0, TT, k)
  Gamma[1:p_ar, ] <- Gamma_init
  
  for (t in (p_ar + 1):TT) {
    pred <- rep(0, k)
    for (j in 1:p_ar) {
      pred <- pred + A_list[[j]] %*% Gamma[t - j, ]
    }
    Gamma[t, ] <- pred
  }
  
  Gamma
}

vec_A <- function(A_list) unlist(lapply(A_list, c))


flip_Lambda <- function(L, a, b) {
  L2 <- L
  L2[a, b] <- 1 - L2[a, b]
  L2
}

propose_A_single_entry <- function(A_list, j, a, b, sd = 0.01) {
  A_prop <- A_list
  A_prop[[j]][a, b] <- A_prop[[j]][a, b] + rnorm(1, 0, sd)
  A_prop
}

propose_A_diag_entry <- function(A_list, j, a, sd = 0.01) {
  A_prop <- A_list
  A_prop[[j]][a, a] <- A_prop[[j]][a, a] + rnorm(1, 0, sd)
  A_prop
}

logprior_A_spikeslab <- function(A_list, Lambda_list,
                                 tau0 = 1e-3, tau1 = 0.1,
                                 pi = 0.1,
                                 diag_mean = 0.8, diag_sd = 0.15,
                                 include_Lambda = TRUE) {
  p_ar <- length(A_list)
  k <- nrow(A_list[[1]])
  lp <- 0
  
  for (j in 1:p_ar) {
    A <- A_list[[j]]
    L <- Lambda_list[[j]]
    
    for (a in 1:k) {
      lp <- lp + dnorm(A[a, a], mean = diag_mean, sd = diag_sd, log = TRUE)
    }
    
    for (a in 1:k) for (b in 1:k) if (a != b) {
      sd_ab <- if (L[a, b] == 1) tau1 else tau0
      lp <- lp + dnorm(A[a, b], mean = 0, sd = sd_ab, log = TRUE)
      if (include_Lambda) lp <- lp + dbinom(L[a, b], size = 1, prob = pi, log = TRUE)
    }
  }
  
  as.numeric(lp)
}

init_death_first <- function(k, p_ar,
                             prob_on = 0,          
                             diag_mean = 0.8, diag_sd = 0.05,
                             off_sd = 0.05) {        
  
  Lambda_list <- vector("list", p_ar)
  A_list <- vector("list", p_ar)
  
  for (j in 1:p_ar) {
    L <- matrix(0, k, k)
    A <- matrix(0, k, k)
    
    diag(L) <- 1
    diag(A) <- rnorm(k, mean = diag_mean, sd = diag_sd)
    
    for (a in 1:k) for (b in 1:k) if (a != b) {
      L[a, b] <- rbinom(1, 1, prob_on)
      if (L[a, b] == 1) A[a, b] <- rnorm(1, 0, off_sd) else A[a, b] <- 0
    }
    
    Lambda_list[[j]] <- L
    A_list[[j]] <- A
  }
  
  list(A = A_list, Lambda = Lambda_list)
}


mh_update_A_sparse_deathfirst <- function(fs, B, D, Psi, alpha, sigma_sq,
                                          A_curr, Gamma_init,
                                          Lambda_curr,
                                          tau0 = 1e-3, tau1 = 0.1, pi = 0.1,
                                          diag_mean = 0.8, diag_sd = 0.15,
                                          rho_max = 1.01,
                                          prop_sd_on = 0.02,
                                          prop_sd_diag = 0.03,
                                          n_moves = 120,
                                          diag_update_prob = 0.30,
                                          flip_move_prob = 0.7,     
                                          p_death = 0.90,            
                                          birth_sd = 0.02) {         
  
  TT <- ncol(fs)
  mu <- as.vector(B %*% alpha)
  k <- ncol(A_curr[[1]])
  p_ar <- length(A_curr)
  
  # target log-posterior
  log_target <- function(A_list, Lambda_list) {
    A_eff <- A_list
    for (j in 1:p_ar) {
      L <- Lambda_list[[j]]
      A_eff[[j]][L == 0] <- 0
      diag(A_eff[[j]]) <- diag(A_list[[j]])
    }
    
    rho <- get_rho(A_eff)
    if (!is.finite(rho) || rho >= rho_max) return(-Inf)
    
    Gamma <- reconstruct_Gamma_det(A_eff, Gamma_init, TT)
    if (!all(is.finite(Gamma)) || max(abs(Gamma)) > 1e8) return(-Inf)
    
    gs <- B %*% t(Gamma)
    ll <- tryCatch(loglik(fs, gs, mu, sigma_sq, Psi, D), error = function(e) -Inf)
    if (!is.finite(ll)) return(-Inf)
    
    lp <- logprior_A_spikeslab(A_list, Lambda_list,
                               tau0 = tau0, tau1 = tau1, pi = pi,
                               diag_mean = diag_mean, diag_sd = diag_sd,
                               include_Lambda = TRUE)
    
    ll + lp
  }
  
  A <- A_curr
  Lambda <- Lambda_curr
  logpost <- log_target(A, Lambda)
  
  acc_diag <- 0
  acc_rw   <- 0
  acc_flip <- 0
  
  for (m in 1:n_moves) {
    
    j <- sample.int(p_ar, 1)
    
    #diagonal update
    if (runif(1) < diag_update_prob) {
      a <- sample.int(k, 1)
      A_prop <- propose_A_diag_entry(A, j, a, sd = prop_sd_diag)
      
      lp_prop <- log_target(A_prop, Lambda)
      la <- lp_prop - logpost
      if (is.finite(la) && runif(1) < exp(min(0, la))) {
        A <- A_prop; logpost <- lp_prop
        acc_diag <- acc_diag + 1
      }
      next
    }
    
    check <- F
    while(check == F){
      a <- sample.int(k, 1)
      b <- sample.int(k, 1)
      if (a != b) check = T
    }
    
    #off diagonal update
    do_flip <- (runif(1) < flip_move_prob)
    
    if (do_flip) {
      Lcur <- Lambda[[j]][a, b]
      
      q_fwd <- if (Lcur == 1) p_death else (1 - p_death)
      q_rev <- if (Lcur == 1) (1 - p_death) else p_death
      
      if (q_fwd <= 0) next
      
      Lambda_prop <- Lambda
      Lambda_prop[[j]] <- flip_Lambda(Lambda[[j]], a, b)
      
      A_prop <- A
      if (Lambda_prop[[j]][a, b] == 1) {
        if (A_prop[[j]][a, b] == 0) A_prop[[j]][a, b] <- rnorm(1, 0, birth_sd)
      } else {
        A_prop[[j]][a, b] <- 0
      }
      
      lp_prop <- log_target(A_prop, Lambda_prop)
      la <- (lp_prop - logpost) + (log(q_rev) - log(q_fwd))  # Hastings correction
      
      if (is.finite(la) && runif(1) < exp(min(0, la))) {
        A <- A_prop; Lambda <- Lambda_prop; logpost <- lp_prop
        acc_flip <- acc_flip + 1
      }
      
    } else {
      if (Lambda[[j]][a, b] == 0) next
      
      A_prop <- propose_A_single_entry(A, j, a, b, sd = prop_sd_on)
      lp_prop <- log_target(A_prop, Lambda)
      la <- lp_prop - logpost
      
      if (is.finite(la) && runif(1) < exp(min(0, la))) {
        A <- A_prop; logpost <- lp_prop
        acc_rw <- acc_rw + 1
      }
    }
  }
  
  A_eff <- A
  for (jj in 1:p_ar) {
    A_eff[[jj]][Lambda[[jj]] == 0] <- 0
    diag(A_eff[[jj]]) <- diag(A[[jj]])
  }
  Gamma_out <- reconstruct_Gamma_det(A_eff, Gamma_init, TT)
  
  list(
    A = A_eff,
    A_raw = A,
    Lambda = Lambda,
    Gamma = Gamma_out,
    stats = c(acc_diag = acc_diag, acc_rw = acc_rw, acc_flip = acc_flip, n_moves = n_moves)
  )
}


# ==============================================================================
# 4. GIBBS SAMPLER 
# ==============================================================================

GibbsSampler_p <- function(data_obj, 
                           R, 
                           burnin, 
                           k, 
                           nbasis, 
                           p_ar, 
                           fixed_gamma1=NULL, 
                           m = 50, 
                           prop_sd_Phi_init = 2e-3, 
                           adapt_every = 100, 
                           target_acc = 0.2, 
                           adapt_mult = 1.2,
                           prop_sd_min = 1e-8, 
                           prop_sd_max = 1e-1,
                           prop_sd_diag = 0.005,
                           prop_sd_on = 0.001, 
                           n_moves = 200, 
                           adapt_until = R){
  
  data <- as.matrix(data_obj[, -1])
  x <- data_obj[, 1]
  
  TT <- ncol(data)
  n <- length(x)
  
  fs <- data[, 1:TT]
  
  # Fixed Matrices
  B <- bs(x, df = k, degree = 3, intercept = TRUE) 
  C <- solve(crossprod(B)) %*% t(B)
  D <- make_fourier(x, nbasis)
  
  #Empirical Bayes for alpha
  ybar <- rowMeans(fs)
  BtB <- crossprod(B)
  BtY <- crossprod(B, ybar)
  alpha <- solve(BtB, BtY)
  
  
  # Hyperparameters
  a0 <- 1e3; b0 <- 2
  nu0 <- 2 * nbasis + 2; S0 <- diag(1, 2*nbasis)
  avec0 <- matrix(0, nrow = k, ncol = 1)
  Sigma0 <- diag(k)
  
  
  # Inizialization
  sigma_sq <- 0.001
  Psi <- diag(1, 2*nbasis)
  Phi <- make_Phi(D, Psi, rel = 1e-6)
  A <- list()
  init <- init_death_first(k = k, p_ar = p_ar,
                           prob_on = 0,
                           diag_mean = 0.8, diag_sd = 0.05,
                           off_sd = 0.03)
  
  A <- init$A
  Lambda_list <- init$Lambda
  Gamma <- matrix(0, TT, k)
  prop_sd_Phi <- prop_sd_Phi_init
  
  
  # Containers
  store_alpha <- array(NA, dim = c(k, R))
  store_A <- array(NA, dim=c(k, k, p_ar, R))
  store_sigmasq <- numeric(R)
  store_Psi <- array(NA, dim = c(2*nbasis, 2*nbasis, R))
  store_Phi <- array(NA, dim = c(n, n, R))
  store_Gamma <- array(NA, dim = c(TT, k, R))
  A_acceptance_vector_diag <- numeric(R)
  A_acceptance_vector_on <- numeric(R)
  Phi_acceptance_vector <- numeric(R)
  likelihood <- numeric(R)
  prop_sd_trace_Phi <- numeric(R)
  
  
  cat("  > Start Gibbs (R=", R, ")\n")
  pb <- txtProgressBar(min = 0, max = R, style = 3)
  
  for(i in 1:R){
    
    # 1. Gamma Reconstruction 
    Gamma <- matrix(0, TT, k)
    if(!is.null(fixed_gamma1)) Gamma[1:p_ar, ] <- fixed_gamma1 else Gamma[1:p_ar, ] <- t(C %*% fs[,1:p_ar])
    fixed_gamma1 <- Gamma[1:p_ar, ]
    for(t in (p_ar+1):TT){
      pred <- rep(0, k)
      for(j in 1:p_ar) pred <- pred + A[[j]] %*% Gamma[t-j, ]
      Gamma[t, ] <- pred 
    }
    
    # 2. Sample Sigma
    gs <- B %*% t(Gamma)
    mu <- B %*% alpha
    sigma_sq <- sample_sigma(fs, a0, b0, Phi, TT, n, gs, mu)
    
    # 3. Sample Alpha
    alpha <- sample_alpha(fs, avec0, Sigma0, Phi, sigma_sq, k, TT, n, B, gs)
    store_alpha[, i] <- alpha
    mu <- B %*% alpha
    
    # 4. pointwise MH update A 
    A_res <- mh_update_A_sparse_deathfirst(
      fs, B, D, Psi, alpha, sigma_sq,
      A_curr = A,
      Gamma_init = fixed_gamma1,
      Lambda_curr = Lambda_list,
      tau0 = 1e-3, tau1 = 0.1, pi = 0.2,
      diag_mean = 0.8, diag_sd = 0.15,
      rho_max = 1.02,
      prop_sd_on = prop_sd_on, 
      prop_sd_diag = prop_sd_diag, 
      n_moves = n_moves,
      p_death = 0.7,
      birth_sd = 0.03
    )
    
    A_acceptance_vector_diag[i] <- A_res$stats[1]
    A_acceptance_vector_on[i] <- A_res$stats[2]
    
    if (i %% 200 == 0) print(A_res$stats)
    
    #adaptive MH
    if (i %% adapt_every == 0 && i <= adapt_until) {
      idx <- (i - adapt_every + 1):i
      acc_rate_diag <- mean(A_acceptance_vector_diag[idx])/(n_moves*0.30)
      acc_rate_on <- mean(A_acceptance_vector_on[idx])/(n_moves*0.15)
      
      if (acc_rate_diag > target_acc) {
        prop_sd_diag <- prop_sd_diag * adapt_mult
      } else if (acc_rate_diag < target_acc) {
        prop_sd_diag <- prop_sd_diag / adapt_mult
      }
      
      if (acc_rate_on > target_acc) {
        prop_sd_on <- prop_sd_on * adapt_mult
      } else if (acc_rate_on < target_acc) {
        prop_sd_on <- prop_sd_on / adapt_mult     
      }
      
    }
    
    A <- A_res$A
    Gamma <- A_res$Gamma
    Lambda_list <- A_res$Lambda
    
    gs <- B %*% t(Gamma)
    
    for(j in 1:p_ar) store_A[, , j, i] <- A[[j]]
    store_Gamma[, , i] <- Gamma
  
    # 5. Sample Phi
    phi_res <- sample_Phi_fourier_mtmh(fs, gs, mu, sigma_sq, D, Psi, nu0, S0, m, rw_scale = prop_sd_Phi)
    Phi_acceptance_vector[i] <- phi_res$acc
    Psi <- phi_res$Psi
    
    #scaling Psi (identifiebility in sigma^2Phi variance term)
    scalePsi <- mean(diag(Psi))
    Psi <- Psi / scalePsi
    sigma_sq <- sigma_sq * scalePsi 
    
    store_sigmasq[i] <- sigma_sq
    Phi <- make_Phi(D, Psi, rel = 1e-6)
    store_Psi[, , i] <- Psi
    store_Phi[, , i] <- Phi
    prop_sd_trace_Phi[i] <- prop_sd_Phi
    
    #adaptive MH
    if (i %% adapt_every == 0 && i <= adapt_until) {
      idx <- (i - adapt_every + 1):i
      acc_rate <- mean(Phi_acceptance_vector[idx])
      
      if (acc_rate > target_acc) {
        prop_sd_Phi <- min(prop_sd_Phi * adapt_mult, prop_sd_max)
      } else if (acc_rate < target_acc) {
        prop_sd_Phi <- max(prop_sd_Phi / adapt_mult, prop_sd_min)
      }
      
    }
    
    # 6. Likelihood computing
    likelihood[i] <- loglik(fs, gs, mu, sigma_sq, Psi, D)
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  return(list(alpha = store_alpha[, -c(1:burnin)],
              A = store_A[,,,-c(1:burnin)], 
              Phi = store_Phi[, ,-c(1:burnin)], 
              Psi = store_Psi[, , -c(1:burnin)],
              sigmasq = store_sigmasq[-c(1:burnin)], 
              likelihood = likelihood[-c(1:burnin)],
              Gamma = store_Gamma[, , -c(1:burnin)], 
              A_acceptance_vector_on = A_acceptance_vector_on[-c(1:burnin)],
              A_acceptance_vector_diag = A_acceptance_vector_diag[-c(1:burnin)],
              Phi_acceptance_vector = Phi_acceptance_vector[-c(1:burnin)], 
              prop_sd_trace_Phi = prop_sd_trace_Phi[-c(1:burnin)]
  )) 
}



# ==============================================================================
# 5. SYNTHETIC DATA GENERATION
# ==============================================================================

run_sim <- function(TT, noise, alpha, A_true, Phi, Gamma0_true, x, k, p_ar){
  B <- bs(x, df=k, degree=3, intercept=TRUE)
  Gamma <- matrix(0, TT, k)
  Gamma[1:p_ar, ] <- Gamma0_true
  for(t in (p_ar+1):TT){
    sum <- 0
    for(r in 1:p_ar){
      sum <- sum + A_true[[r]] %*% Gamma[t-r, ]
    }
    Gamma[t, ] <- sum
  }
  
  mu <- B %*% alpha
  Mu <- mu %*% t(c(rep(1, TT)))
  
  noise_mat <- t(mvtnorm::rmvnorm(TT, mean=rep(0, length(x)), sigma=(noise^2)*Phi))
  f <- Mu + (B %*% t(Gamma)) + noise_mat
  
  list(data=data.frame(x=x, f), Gamma0 = Gamma0_true)
}




