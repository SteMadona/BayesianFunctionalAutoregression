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

sample_A_stable_CORRETTO <- function(fs, Phi, V0, Aprior, sigma_sq, Gamma, C, k, max_tries=2000){
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
  repeat{
    tries <- tries + 1
    A_star <- rmatnorm(1, An, Vn, sigma_sq * U)
    
    if (max(Mod(eigen(A_star)$values)) < 0.995) return(A_star)
    
    if (tries > max_tries) {
      # Fallback: Stabilizza la media a posteriori (più robusto dell'ultimo sample)
      return(make_stable_from_A(An, target_rho = 0.99))
    }
  }
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

GibbsSampler_Engine <- function(data_obj, R=2000, burnin=500, k=4, nbasis=3, fixed_gamma1=NULL){
  
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
  Aprior <- diag(0.8, k)
  V0 <- diag(10, k) 
  
  # Inizializzazione
  sigma_sq <- 0.01
  Psi <- diag(1, 2*nbasis)
  Phi <- D %*% Psi %*% t(D)
  A <- diag(0.5, k) # Start neutro
  Gamma <- matrix(0, TT, k)
  
  # Containers
  store_A <- array(NA, dim=c(k, k, R))
  store_sigmasq <- numeric(R)
  store_Phi <- array(NA, dim = c(n, n, R))
  
  
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
    A <- sample_A_stable_CORRETTO(fs, Phi, V0, Aprior, sigma_sq, Gamma, C, k)
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
              sigmasq = store_sigmasq[-c(1:burnin)]))
}

# ==============================================================================
# 4. LOOP DEGLI SCENARI
# ==============================================================================

# Parametri Globali
set.seed(999)
x_grid <- seq(0, 1, length.out = 100)
k_dim <- 8
A_true <- generate_A_banded(k_dim, diag_val = 0.75, off_val = 0.15)
Phi_true <- MCMCpack::riwish(105, diag(100)) # Phi fissa per tutti
gamma1_true <- rnorm(k_dim)

# Definizione Scenari
scenarios <- list(
  #  list(name="Ideal",     TT=200, noise=0.001, R=200, burn=50),
  list(name="HighData",  TT=200, noise=0.03,  R=1000, burn=0)
  #  list(name="Baseline",  TT=50,  noise=0.03,  R=200, burn=50)
)

results <- list()

# Funzione simulazione dati on-the-fly
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

# --- ESECUZIONE ---

for(scen in scenarios){
  cat("\n\n================================================\n")
  cat("RUNNING SCENARIO:", scen$name, "\n")
  cat("TT =", scen$TT, "| Noise =", scen$noise, "\n")
  cat("================================================\n")
  
  # 1. Genera Dati
  sim_data <- run_sim(scen$TT, scen$noise, A_true, Phi_true, gamma1_true, x_grid, k_dim)
  
  # 2. Esegui Sampler
  # Usiamo fixed_gamma1 = gamma1_true per isolare SOLO la capacità di stimare A
  tic()
  res <- GibbsSampler_Engine(sim_data, R=scen$R, burnin=scen$burn, k=k_dim, fixed_gamma1=gamma1_true)
  time_elapsed <- toc()
  
  # 3. Calcola Statistiche
  A_est_mean <- apply(res$A, 1:2, mean)
  max_eig_true <- max(abs(eigen(A_true)$values))
  max_eig_est <- max(abs(eigen(A_est_mean)$values))
  mse <- mean((A_est_mean - A_true)^2)
  
  cat("Max Eigen (True):", max_eig_true, "\n")
  cat("Max Eigen (Est) :", max_eig_est, "\n")
  cat("MSE Matrix      :", mse, "\n")
  
  # 4. Salva
  results[[scen$name]] <- list(
    true = A_true,
    est = A_est_mean,
    chain = res$A[1,1,], # Salva traccia solo di un elemento per leggerezza
    eigen_true = max_eig_true,
    eigen_est = max_eig_est,
    mse = mse
  )
  
  #  saveRDS(results, file="results_diagnostic/192011")
}

# ==============================================================================
# 5. PLOT FINALE COMPARATIVO
# ==============================================================================
par(mfrow=c(2, 3), mar=c(3,3,3,1))

# Riga 1: Heatmaps delle stime
zlims <- range(c(A_true, lapply(results, function(x) x$est)), na.rm=TRUE)
for(n in names(results)){
  image(results[[n]]$est, main=paste(n, "\nEigen:", round(results[[n]]$eigen_est, 2)), zlim=zlims)
}

# Riga 2: Tracce di A[1,1]
for(n in names(results)){
  plot(results[[n]]$chain, type='l', main=paste("Trace A[1,1] -", n), ylab="Val", xlab="Iter")
  abline(h=A_true[1,1], col="red", lwd=2)   
}

par(mfrow = c(1, 2))
image(results[[1]]$est, main = "A est")
image(A_true, main = "A true")

eig(results[[1]]$est)
eig(results[[1]]$mse)

traceplot(log(as.mcmc(res$sigmasq)))
summary(res$sigmasq)
mean(res$sigmasq)

apply(res$A, c(1,2), mean)

View(apply(res$Phi, c(1, 2), mean))
