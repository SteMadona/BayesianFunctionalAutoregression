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

sample_Phi_fourier_mtmh <- function(fs, mus, gs, sigma, D, Psi, nu0, S0, n, m, n_acc){
  # bersaglio: usare SOLO la likelihood (q = prior indipendente)
  log_target <- function(Psi){
    Psi <- make_posdef(Psi)
    as.numeric(loglik(fs, mus, gs, sigma, D, Psi))  # <-- solo ll
  }
  
  ## m proposte forward
  Psi_prop_list <- vector("list", m)
  Phi_prop_list <- vector("list", m)
  logw_f_raw <- numeric(m)
  
  for (i in 1:m) {
    Psi_prop_i <- make_posdef(riwish(nu0, S0))
    Phi_prop_i <- make_posdef(D %*% Psi_prop_i %*% t(D))
    Psi_prop_list[[i]] <- Psi_prop_i
    Phi_prop_list[[i]] <- Phi_prop_i
    logw_f_raw[i] <- log_target(Psi_prop_i)   # log-pesi GREZZI
  }
  
  # scelta dell'indice con softmax dei log-pesi grezzi
  idx_star <- sample.int(m, 1, prob = softmax_from_logs(logw_f_raw))
  Psi_star <- Psi_prop_list[[idx_star]]
  Phi_star <- Phi_prop_list[[idx_star]]
  
  ## m-1 backward + stato corrente come m-esimo
  Psi_back_list <- vector("list", m)
  Phi_back_list <- vector("list", m)
  for (i in 1:(m-1)) {
    Psi_back_list[[i]] <- make_posdef(riwish(nu0, S0))
    Phi_back_list[[i]] <- make_posdef(D %*% Psi_back_list[[i]] %*% t(D))
  }
  Psi_back_list[[m]] <- Psi
  Phi_back_list[[m]] <- make_posdef(D %*% Psi %*% t(D))
  
  logw_b_raw <- vapply(Psi_back_list, log_target, numeric(1))
  
  # log-sum-exp CORRETTO e rapporto FORWARD/BACKWARD
  log_den <- logsumexp(logw_f_raw)  # forward
  log_num <- logsumexp(logw_b_raw)  # backward
  log_alpha <- log_den - log_num     # <-- forward - backward
  alpha <- min(1, exp(log_alpha))
  
  if (runif(1) < alpha) {
    list(Phi = Phi_star, Psi = Psi_star, n_acc = n_acc + 1L)
  } else {
    list(Phi = make_posdef(D %*% Psi %*% t(D)), Psi = Psi, n_acc = n_acc)
  }
}
