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

source("code/UsefulFunctions.R")
source("code/ConjugatePosteriors.R")
source("code/GibbsSamplers.R")
source("code/MHforPsi.R")

simulation_test_VARp_selective <- function(x, n_fun, nbasis, noise_sd = 0.5,
                                           alpha, Phi, A_list, k) {
  B <- bs(x, df = nbasis, degree = 3)
  C <- solve(crossprod(B)) %*% t(B)
  
  A_list <- make_stable_from_A_list(A_list)
  
  n <- length(x)
  mu <- as.vector(B %*% alpha)
  Phi <- make_posdef(Phi)
  Sigma_eps <- (noise_sd^2) * Phi
  
  p <- length(A_list)
  
  Gamma <- matrix(0, nrow = n_fun, ncol = k)
  G <- matrix(0, nrow = n, ncol = n_fun)
  
  # random initial states for the first p time points
  for (t in 1:p) {
    Gamma[t, ] <- rnorm(k)
    G[, t] <- as.vector(B %*% Gamma[t, ])
  }
  
  for (t in (p+1):n_fun) {
    Gamma[t, ] <- rep(0, k)
    for (r in 1:p) {
      if (!is.null(A_list[[r]])) {
        Gamma[t, ] <- Gamma[t, ] + t(A_list[[r]]) %*% Gamma[t - r, ]
      }
    }
    G[, t] <- as.vector(B %*% Gamma[t, ])
  }

  functions <- matrix(NA, nrow = n, ncol = n_fun)
  for (t in 1:n_fun) {
    f_t <- mvtnorm::rmvnorm(1, mean = mu + G[, t], sigma = Sigma_eps)
    functions[, t] <- as.numeric(f_t)
  }
  
  data <- data.frame(x = x, functions)
  colnames(data) <- c("x", paste0("t", 1:n_fun))
  
  return(data)

}

x <- seq(0, 1, length.out = 200)
alpha_test1 <- c(-2, 0.5, 1, 0, 1.5, -1, 0, 1.5, -1.5, -1, 1, 0)
Phi_test1 <- riwish(205, diag(200))
k <- 12
p <- 10

A_list <- vector("list", p)
A_list[[1]] <- 0.5 * diag(k)
A_list[[2]] <- -0.3 * diag(k)
A_list[[7]] <- 0.2 * matrix(runif(k * k, -0.1, 0.1), k, k)


df_test3 <- simulation_test_VARp_selective(x, n_fun = 100, nbasis = 12,
                                           noise_sd = 0.1,
                                           alpha = alpha_test1,
                                           Phi = Phi_test1,
                                           A_list = A_list,
                                           k = k)


#PSI with MTMH and mulitple orders --------------------------------------------

tic()
out_psi_mtmh_p <- GS_mtmh_p(df = df_test3, 
                            a0 = 2, 
                            b0 = 1, 
                            avec0 = test1, 
                            Sigma0 = diag(12), 
                            V0 = diag(12), 
                            Aprior = test3, 
                            p = 10, 
                            R = 20,
                            burnin = 0,
                            nbasis = 5,
                            nu0 = 12, 
                            S0 = diag(10), 
                            m = 10
)
toc()

out_psi_mtmh_p$sigma
out_psi_mtmh_p$b
out_psi_mtmh_p$Gamma[ , , 1]

