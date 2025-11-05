loglik <- function(fs, mus, gs, sigma, D, Psi, nugget = 1e-8){
  n <- nrow(fs); T <- ncol(fs)
  
  fs <- as.matrix(fs)
  mus <- as.matrix(mus)
  gs <- as.matrix(gs)
  
  Phi <- D %*% Psi %*% t(D)
  Phi <- make_posdef(Phi) #+ diag(nugget, n)
  Phi_inv <- solve(Phi)
  ld_Phi <- logdet_spd(Phi)
  
  E <- fs - mus - gs  #n x t, every column is e_t
  
  quad <- sum(vapply(1:T, function(t){
    as.numeric(crossprod(E[, t], Phi_inv %*% E[, t]))
  }, numeric(1)))   #use of vapply to ensure getting a scalar output
  
  ll <- -.5 * (quad/sigma) - .5 * T * (n * log(sigma) + ld_Phi)
  return(as.numeric(ll))
}

sample_Phi_fourier_mh <- function(fs, mus, gs, sigma, D, Psi, Phi, nu0, S0, n, n_acc){
  p <- ncol(D)
  
  Psi_prop <- riwish(nu0, S0)  #takes Phi from the step before as scale matrix
  Psi_prop <- make_posdef(Psi_prop)
  
  lp_cur <- log_diwish(Psi, nu0, S0)
  lp_prop <- log_diwish(Psi_prop, nu0, S0)
  
  Phi_prop <- D %*% Psi_prop %*% t(D)
  Phi_prop <- make_posdef(Phi_prop)
  
  ll_cur <- loglik(fs, mus, gs, sigma, D, Psi) #così facendo viene calcolata due volte,
  #forse conviene salvare quella dell'iterazione precedente
  ll_prop <- loglik(fs, mus, gs, sigma, D, Psi_prop)
  
  log_alpha <- (ll_prop + lp_prop) - (ll_cur + lp_cur)
  if(is.na(log_alpha)) log_alpha <- -Inf
  
  if(log(runif(1)) < log_alpha){
    n_acc <- n_acc + 1
    Phi_out <- Phi_prop
    Psi_out <- Psi_prop
  } else{
    Phi_out <- Phi
    Psi_out <- Psi
    n_acc <- n_acc
  }
  
  return(list(Phi = Phi_out, 
              Psi = Psi_out, 
              n_acc = n_acc))
}


rw_step_chol <- function(Psi, rw_scale = 0.05){
  L <- chol(make_posdef(Psi))
  L <- t(L)
  
  eps <- matrix(rnorm(length(L), mean = 0, sd = rw_scale), nrow = nrow(L))
  eps[upper.tri(eps)] <- 0
  
  L_star <- L + eps
  
  Psi_star <- make_posdef(L_star %*% t(L_star))
  
  return(Psi_star)
}

sample_Phi_fourier_mtmh <- function(fs, mus, gs, sigma, D, Psi, nu0, S0, 
                                    m, n_acc, 
                                    nu_q = NULL, s = 1.0, rw_prob = 0.3, rw_scale = 0.05) {
  
  #Hybrid MTMH for Psi:
  # - with prob (1 - rw_prob): Multiple-Try MH using independence IW(nu_q, s*Psi_current)
  # - with prob rw_prob: single symmetric RW step on Cholesky with scale = rw_scale
  
  p <- nrow(Psi)
  if (is.null(nu_q)) {
    nu_q <- p + 4
  }
  
  if (runif(1) < rw_prob) {
    # --- Symmetric RW branch (no proposal terms in acceptance) ---
    Psi_star <- rw_step_chol(Psi, rw_scale = rw_scale)
    log_alpha <- (loglik(fs, mus, gs, sigma, D, Psi_star) + log_diwish(Psi_star, nu0, S0)) -
                 (loglik(fs, mus, gs, sigma, D, Psi)      + log_diwish(Psi,      nu0, S0))
    if (runif(1) < min(1, exp(log_alpha))) {
      return(list(
        Phi = add_nugget(D %*% Psi_star %*% t(D), rel = 1e-2),
        Psi = Psi_star,
        n_acc = n_acc + 1L,
        acc = 1
      ))
    } else {
      return(list(
        Phi = add_nugget(D %*% Psi %*% t(D), rel = 1e-2),
        Psi = Psi,
        n_acc = n_acc,
        acc = 0
      ))
    }
  }
  
  # --- Independence MTMH branch ---
  S_q <- s * make_posdef(Psi)
  
  log_target <- function(P) {
    loglik(fs, mus, gs, sigma, D, P) + log_diwish(P, nu0, S0)
  }
  
  # forward proposals and weights target/proposal
  Psi_list <- replicate(m, make_posdef(riwish(nu_q, S_q)), simplify = FALSE)
  logw_f <- sapply(Psi_list, function(P) log_target(P) - log_diwish(P, nu_q, S_q))
  
  idx <- sample.int(m, 1, prob = exp(logw_f - max(logw_f)))
  Psi_star <- Psi_list[[idx]]
  
  # backward set: m-1 fresh from q + current
  Psi_back <- c(replicate(m - 1, make_posdef(riwish(nu_q, S_q)), simplify = FALSE), list(Psi))
  logw_b <- sapply(Psi_back, function(P) log_target(P) - log_diwish(P, nu_q, S_q))
  
  log_alpha <- logsumexp(logw_f) - logsumexp(logw_b)
  alpha <- min(1, exp(log_alpha))
  
  if (runif(1) < alpha) {
    list(Phi = add_nugget(D %*% Psi_star %*% t(D), rel = 1e-2),
         Psi = Psi_star,
         n_acc = n_acc + 1L,
         acc = 1)
  } else {
    list(Phi = add_nugget(D %*% Psi %*% t(D), rel = 1e-2),
         Psi = Psi,
         n_acc = n_acc,
         acc = 0)
  }
}



