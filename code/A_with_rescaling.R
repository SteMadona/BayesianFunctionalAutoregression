rm(list = ls()) ; graphics.off() ; cat("\014")

# ==============================================================================
# 1. CONFIGURAZIONE E LIBRERIE
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
library(gridExtra)

# Crea cartella per salvare i risultati parziali
#if(!dir.exists("results_diagnostic")) dir.create("results_diagnostic")

# ==============================================================================
# 2. FUNZIONI HELPER E CORE (INVARIATE E CORRETTE)
# ==============================================================================

make_posdef <- function(M, tol = 1e-8){
  Msym <- (M + t(M)) / 2
  e <- eigen(Msym, symmetric = TRUE)
  vals <- pmax(e$values, tol)
  Mpd <- e$vectors %*% diag(vals) %*% t(e$vectors)
  return((Mpd + t(Mpd))/2)
}

make_stable_from_A <- function(A, target_rho = 0.99){ # Target rho più alto per non schiacciare troppo
  ev <- eigen(A)$values
  rho <- max(abs(ev))
  if(rho <= target_rho) return(A)
  return(A * (target_rho / rho))
}

add_nugget <- function(Phi, rel = 0.1){
  n <- nrow(Phi)
  Phi <- make_posdef(Phi)
  tau <- rel * (sum(diag(Phi)) / n)
  Phi <- Phi + diag(tau, n)
  Phi / mean(diag(Phi))   # normalize scale
}

# Funzione per generare A a banda (più realistica)
generate_A_banded <- function(k, diag_val=0.75, off_val=0.1) {
  A <- diag(diag_val, k)
  # Aggiungi un po' di struttura fuori diagonale
  for(i in 1:(k-1)) {
    A[i, i+1] <- off_val
    A[i+1, i] <- -off_val # Asimmetria per rendere autovalori complessi
  }
  # Normalizza per stabilità
  make_stable_from_A(A, 0.85)
}

make_fourier <- function(x, K) {
  L <- max(x) - min(x); xb <- (x - min(x)) / L
  B <- matrix(NA, length(x), 2*K)
  for (k in 1:K) { B[, 2*k-1] <- sin(2*pi*k*xb); B[, 2*k] <- cos(2*pi*k*xb) }
  return(B)
}

eig <- function(Amat){max(abs(eigen(Amat)$value))}


# --- Funzioni di Campionamento ---

sample_sigma <- function(fs, a0, b0, Phi, TT, n, gs){
  Phi_pd <- make_posdef(Phi)
  Phi_inv <- chol2inv(chol(Phi_pd))
  resid_sum <- 0
  for(t in 1:TT){
    e <- fs[, t] - gs[, t]
    resid_sum <- resid_sum + as.numeric(t(e) %*% Phi_inv %*% e)
  }
  # Ritorna sigma^2
  rinvgamma(1, shape = a0 + (n*TT)/2, scale = b0 + 0.5*resid_sum)
}

sample_A_stable_CORRETTO <- function(fs, Phi, V0, Aprior, sigma_sq, Gamma, C, k, max_tries=200){
  TT <- ncol(fs)
  U <- C %*% Phi %*% t(C); U <- make_posdef((U + t(U))/2)
  
  # Lagged Regression setup
  F_tilde <- t(C %*% fs)         # Coefficienti osservati
  Y_data <- F_tilde[2:TT, ]      # t = 2..T
  X_data <- Gamma[1:(TT-1), ]    # t = 1..T-1
  
  V0_inv <- chol2inv(chol(make_posdef(V0)))
  XX <- t(X_data) %*% X_data
  XY <- t(X_data) %*% Y_data
  
  Vn_inv <- V0_inv + XX
  Vn <- chol2inv(chol(make_posdef(Vn_inv)))
  An <- Vn %*% (V0_inv %*% Aprior + XY)
  
  # Campionamento
  tries <- 0
  rhos <- numeric(max_tries)
  for(i in 1:max_tries){
    A_star <- rmatnorm(1, An, Vn, sigma_sq*U)
    rhos[i] <- max(Mod(eigen(A_star)$values))
    
    if (max(Mod(eigen(A_star)$values)) < 0.995) return(list( A = A_star, 
                                                             rhos = rhos))
    
  }
  return(list(A = make_stable_from_A(An, target_rho = 0.99), 
              rhos = rhos))
  
}



sample_Phi_mh <- function(fs, gs, sigma_sq, D, Psi_curr, nu0, S0, rw_step=0.05){
  # Log-Posterior function per Psi
  log_post <- function(P) {
    Phi <- make_posdef(D %*% P %*% t(D))
    Phi_inv <- chol2inv(chol(Phi))
    logdet <- determinant(Phi, log=TRUE)$modulus
    
    quad_term <- sum(sapply(1:ncol(fs), function(t) {
      e <- fs[,t] - gs[,t]
      t(e) %*% Phi_inv %*% e
    }))
    
    # Likelihood + Prior (Inverse Wishart)
    ll <- -0.5 * (quad_term/sigma_sq) - 0.5 * ncol(fs) * logdet
    lprior <- -(nu0 + nrow(P) + 1)/2 * log(det(P)) - 0.5 * sum(diag(S0 %*% solve(P)))
    return(ll + lprior)
  }
  
  # Proposal: RW su Cholesky
  L <- t(chol(Psi_curr))
  E <- matrix(rnorm(length(L), 0, rw_step), nrow(L)); E[upper.tri(E)] <- 0
  L_prop <- L + E
  Psi_prop <- make_posdef(L_prop %*% t(L_prop))
  
  lp_curr <- log_post(Psi_curr)
  lp_prop <- log_post(Psi_prop)
  
  if(runif(1) < exp(lp_prop - lp_curr)){
    return(list(Psi=Psi_prop, acc=1))
  } else {
    return(list(Psi=Psi_curr, acc=0))
  }
}

# ==============================================================================
# 3. GIBBS SAMPLER (ENGINE)
# ==============================================================================

GibbsSampler_Engine <- function(data_obj, R=500, burnin=0, k=4, nbasis=3, fixed_gamma1=NULL){
  
  fs <- as.matrix(data_obj$data[, -1])
  x <- data_obj$data[, 1]
  TT <- ncol(fs)
  n <- length(x)
  
  # Matrici fisse
  B <- bs(x, df = k, degree = 3, intercept = TRUE) # Importante: intercept=TRUE per flessibilità
  C <- solve(crossprod(B)) %*% t(B)
  D <- make_fourier(x, nbasis)
  
  # Hyperparametri
  a0 <- 2; b0 <- 1
  nu0 <- 2 * nbasis + 2; S0 <- diag(1, 2*nbasis)
  # Prior lasco per A per lasciar parlare i dati
  Aprior <- diag(0.2, k_dim, k_dim)
  V0 <- diag(2, k) 
  
  # Inizializzazione
  sigma_sq <- 0.01
  Psi <- diag(1, 2*nbasis)
  Phi <- D %*% Psi %*% t(D)
  A <- diag(0.25, k) # Start neutro
  Gamma <- matrix(0, TT, k)
  
  # Containers
  store_A <- array(NA, dim=c(k, k, R))
  store_sigmasq <- numeric(R)
  store_Phi <- array(NA, dim = c(n, n, R))
  Rho_matrix <- matrix(0, R, 200)
  
  cat("  > Start Gibbs (R=", R, ")\n")
  pb <- txtProgressBar(min = 0, max = R, style = 3)
  
  for(i in 1:R){
    
    # 1. Gamma Reconstruction (Deterministico)
    if(!is.null(fixed_gamma1)){
      Gamma[1,] <- fixed_gamma1
    } else {
      # Se non fisso, inizializza col dato osservato proiettato (migliore guess)
      Gamma[1,] <- (C %*% fs)[,1] 
    }
    for(t in 2:TT) Gamma[t,] <- A %*% Gamma[t-1,]
    
    # 2. Sample Sigma
    gs <- B %*% t(Gamma)
    sigma_sq <- sample_sigma(fs, a0, b0, Phi, TT, n, gs)
    store_sigmasq[i] <- sigma_sq
    
    # 3. Sample A
    # Nota: Passiamo sigma_sq stimato. Se il modello è buono, convergerà.
    A_sample <- sample_A_stable_CORRETTO(fs, Phi, V0, Aprior, sigma_sq, Gamma, C, k, max_tries = 200)
    Rho_matrix[i, ] <- A_sample$rhos
    A <- A_sample$A
    store_A[,,i] <- A
    
    # 4. Sample Phi
    # Aggiorniamo gs con la nuova Gamma (che dipende dal nuovo A se non usiamo fixed_gamma nel loop A)
    # In questo schema deterministico, Gamma è rigenerata al passo 1, quindi usiamo quella vecchia per Phi
    phi_res <- sample_Phi_mh(fs, gs, sigma_sq, D, Psi, nu0, S0)
    Psi <- phi_res$Psi
    Phi <- add_nugget(D %*% Psi %*% t(D))
    store_Phi[, , i] <- Phi
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  return(list(A = store_A[,,-c(1:burnin)], 
              Phi = store_Phi[, ,-c(1:burnin)], 
              sigmasq = store_sigmasq[-c(1:burnin)], 
              rhos = Rho_matrix[-c(1:burnin), ]))
}

# ==============================================================================
# 4. LOOP DEGLI SCENARI
# ==============================================================================

# Parametri Globali
#set.seed(999)
set.seed(100)
x_grid <- seq(0, 1, length.out = 100)
k_dim <- 6

Phi_true <- MCMCpack::riwish(105, diag(100)) # Phi fissa per tutti
gamma1_true <- rnorm(k_dim)

reps <- 5

A_list <- list()

for(i in 1:reps){
  A <- diag(sample(c(-1, 1), k_dim, T)*round(runif(k_dim, 0.6, 0.9), 2), k_dim)
  for(j in 1:(k_dim-1)) {
    A[j, j+1] <- round(runif(1, 0.05, 0.25), 1)
    A[j+1, j] <- -round(runif(1, 0.05, 0.25), 1) 
  }
  A[sample(1:k_dim, 1), sample(1:k_dim, 1)] <- round(runif(1, -0.7, 0.7), 2)
  A[sample(1:k_dim, 1), sample(1:k_dim, 1)] <- round(runif(1, -0.7, 0.7), 2)
  A[sample(1:k_dim, 1), sample(1:k_dim, 1)] <- round(runif(1, -0.7, 0.7), 2)
  A[sample(1:k_dim, 1), sample(1:k_dim, 1)] <- round(runif(1, -0.7, 0.7), 2)
  A_list[[i]] <- A
}

sapply(A_list, function(u) max(abs(eigen(u)$values)))
mean(sapply(A_list, function(u) max(abs(eigen(u)$values))))
sum(sapply(A_list, function(u) max(abs(eigen(u)$values))) > 1)

TT <- 100
noise <- 0.03
R <- 300
burnin <- 100

run_sim <- function(TT, noise, A, Phi, g1, x, k){
  B <- bs(x, df=k, degree=3, intercept=TRUE)
  Gamma <- matrix(0, TT, k)
  Gamma[1,] <- g1
  for(t in 2:TT) Gamma[t,] <- A %*% Gamma[t-1,]
  
  # Aggiungi rumore spaziale correlato
  noise_mat <- t(mvtnorm::rmvnorm(TT, mean=rep(0, length(x)), sigma=(noise^2)*Phi))
  f <- (B %*% t(Gamma)) + noise_mat
  
  list(data=data.frame(x=x, f), gamma1=g1)
}

results <- list()
results_big <- list

for(i in 1:reps){
  
  # 1. Genera Dati
  sim_data <- run_sim(TT, noise, A_list[[i]], Phi_true, gamma1_true, x_grid, k_dim)
  
  # 2. Esegui Sampler
  # Usiamo fixed_gamma1 = gamma1_true per isolare SOLO la capacità di stimare A
  tic()
  res <- GibbsSampler_Engine(sim_data, R=R, burnin=burnin, k=k_dim, fixed_gamma1=gamma1_true)
  time_elapsed <- toc()
  
  # 3. Calcola Statistiche
  A_est_mean <- apply(res$A, 1:2, mean)
  A_est_median <- apply(res$A, 1:2, median)
  max_eig_true <- max(abs(eigen(A_list[[i]])$values))
  max_eig_est <- max(abs(eigen(A_est_mean)$values))
  max_eig_med <- max(abs(eigen(A_est_median)$values))
  mse <- mean((A_est_mean - A_list[[i]])^2)
  mse2 <- mean((A_est_median - A_list[[i]])^2)
  rho_median <- median(res$rhos[res$rhos != 0])
  mat_rescaled <- A_est_mean
  diag(mat_rescaled) <- diag(mat_rescaled)*rho_median
  mse_rescaled <- mean((mat_rescaled - A_list[[i]])^2)
  max_eig_rescaled <- max(abs(eigen(mat_rescaled)$values))
  A_diag <- diag(A_list[[i]])
  A_diag_mean <- diag(A_est_mean)
  A_diag_median <- diag(A_est_median)
  A_diag_resc <- diag(mat_rescaled)
  
  cat("Max Eigen (True):", max_eig_true, "\n")
  cat("Max Eigen (Est) :", max_eig_est, "\n")
  cat("Max Eigen (Med) :", max_eig_med, "\n")
  cat("Max Eigen (Res) :", max_eig_rescaled, "\n")
  cat("MSE Matrix Mean :", mse, "\n")
  cat("MSE Matrix Med  :", mse2, "\n")
  cat("MSE Rescaled Mat;", mse_rescaled, "\n")
  cat("diag true:", round(A_diag, 2), "\n")
  cat("diag mean:", round(A_diag_mean, 2), "\n")
  cat("diag med:", round(A_diag_median, 2), "\n")
  cat("diag resc:", round(A_diag_resc, 2), "\n")
  
  # 4. Salva
  results[[i]] <- list(
    true = A_list[[i]],
    est = A_est_mean,
    chain = res$A[1,1, ], # Salva traccia solo di un elemento per leggerezza
    eigen_true = max_eig_true,
    eigen_est = max_eig_est,
    mse = mse
  )
  
  #results_big[[i]] <- res
  #  saveRDS(results, file="results_diagnostic/192011")
}

mean(unlist(A_list))
mean((A_list[[1]] - mean(unlist(A_list)))^2)
mean((A_list[[2]] - mean(unlist(A_list)))^2)
mean((A_list[[3]] - mean(unlist(A_list)))^2)
mean((A_list[[4]] - mean(unlist(A_list)))^2)
mean((A_list[[5]] - mean(unlist(A_list)))^2)

results[[1]]$est
make_stable_from_A <- function(A, target_rho = 0.99){ # Target rho più alto per non schiacciare troppo
  ev <- eigen(A)$values
  rho <- max(abs(ev))
  if(rho <= target_rho) return(A)
  return(A * (target_rho / rho))
}

prova <- results[[1]]$est
diag(prova) <- diag(results[[1]]$est) *(0.848/results[[1]]$eigen_est)
mean((A_list[[1]] - prova)^2)

prova <- results[[2]]$est
diag(prova) <- diag(results[[2]]$est) *(0.831/results[[2]]$eigen_est)
mean((A_list[[2]] - prova)^2)

prova <- results[[3]]$est
diag(prova) <- diag(results[[3]]$est) *(0.782/results[[3]]$eigen_est)
mean((A_list[[3]] - prova)^2)

prova <- results[[4]]$est
diag(prova) <- diag(results[[4]]$est) *(0.856/results[[4]]$eigen_est)
mean((A_list[[4]] - prova)^2)

prova <- results[[5]]$est
diag(prova) <- diag(results[[5]]$est) *(0.858/results[[5]]$eigen_est)
mean((A_list[[5]] - prova)^2)
