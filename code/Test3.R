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

source("code/UsefulFunctions.R")
source("code/ConjugatePosteriors.R")
source("code/GibbsSamplers.R")
source("code/MHforPsi.R")

set.seed(100)

simulation_test_VARp_selective <- function(x, n_fun, nbasis, noise_sd = 0.5,
                                           alpha, Phi, A_list, k) {
  B <- bs(x, df = nbasis, degree = 3)
  C <- solve(crossprod(B)) %*% t(B)
  
  A_list <- make_stable_from_A_list(A_list, target_rho = 0.9)
  
  n <- length(x)
  mu <- as.vector(B %*% alpha)
  Phi <- make_posdef(Phi)
  Sigma_eps <- (noise_sd^2) * Phi
  
  p <- length(A_list)
  
  Gamma <- matrix(0, nrow = n_fun, ncol = k)
  G <- matrix(0, nrow = n, ncol = n_fun)
  
  # random initial states for the first p time points
  for (t in 1:p) {
    Gamma[t, ] <- as.vector(mvtnorm::rmvnorm(1, sigma = diag(1e-4, k)))
    G[, t] <- as.vector(B %*% Gamma[t, ])
  }
  
  for (t in (p+1):n_fun) {
    Gamma[t, ] <- rep(0, k)
    for (r in 1:p) {
      if (!is.null(A_list[[r]])) {
        Gamma[t, ] <- Gamma[t, ] + A_list[[r]] %*% Gamma[t - r, ] 
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
  
  list(data = data, 
       gamma = Gamma)
}


#A_test2 <- matrix(rnorm(12*12), 12, 12) + diag(0.5, 12)


x <- seq(0, 1, length.out = 200)
alpha_test1 <- c(-2, 0.5, 1, 0, 1.5, -1)
Phi_test1 <- riwish(205, diag(200))
k <- 6
p <- 3

test1 <- rnorm(k, 0, 1) #alpha0 vector
test3 <- matrix(0, k, k)  #A prior

A_list <- vector("list", p)
A_list[[1]] <- rmatnorm(1, matrix(0, k, k), diag(0.5, k, k), diag(0.5, k, k))
A_list[[2]] <- rmatnorm(1, matrix(0, k, k), diag(0.5, k, k), diag(0.5, k, k))
A_list[[3]] <- matrix(0, k, k)

df_test3 <- simulation_test_VARp_selective(x, n_fun = 100, nbasis = 6,
                                           noise_sd = 0.1,
                                           alpha = alpha_test1,
                                           Phi = Phi_test1,
                                           A_list = A_list,
                                           k = k)

df_test3$gamma

#PSI with MTMH and mulitple orders --------------------------------------------

tic()
out_psi_mtmh_p <- GS_mtmh_p(df = df_test3$data, 
                            a0 = 2, 
                            b0 = 1, 
                            avec0 = test1, 
                            Sigma0 = diag(6), 
                            V0 = diag(6), 
                            Aprior = test3, 
                            p = 3, 
                            R = 1500,
                            burnin = 500,
                            nbasis = 5,
                            nu0 = 12, 
                            S0 = diag(10), 
                            m = 10
)
toc()

out_psi_mtmh_p$sigma
out_psi_mtmh_p$b
out_psi_mtmh_p$Gamma[ , , 20]
out_psi_mtmh_p$A

sigma <- as.mcmc(out_psi_mtmh_p$sigma)
traceplot(log(sigma))

matrix_list <- (apply(out_psi_mtmh_p$A, c(1, 2, 3), mean))
mean(abs(matrix_list[, , 1]))
mean(abs(matrix_list[, , 2]))
mean(abs(matrix_list[, , 3]))

(apply(out_psi_mtmh_p$A, c(1, 2, 3), mean))

View(df_test3$data)
