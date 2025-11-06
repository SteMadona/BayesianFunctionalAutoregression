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

source("code/UsefulFunctions.R")
source("code/ConjugatePosteriors.R")
source("code/GibbsSamplers.R")
source("code/MHforPsi.R")


simulation_test2 <- function(x, n_fun, nbasis, noise_sd = 0.5, alpha, Phi, A, k){
  # spline basis
  B <- bs(x, df = nbasis, degree = 3)
  C <- solve(crossprod(B)) %*% t(B)
  
  n_points <- length(x)
  n <- length(x)
  
  A <- make_stable_from_A(A)
  
  mu <- as.vector(B %*% alpha)
  Phi <- make_posdef(Phi)
  Sigma_eps <- (noise_sd^2) * Phi
  
  Gamma <- matrix(0, n_fun, k)
  Gamma[1, ] <- rnorm(k)
  G <- matrix(0, n, n_fun)
  
  for(t in 2:n_fun){
    Gamma[t, ] <- t(A) %*% Gamma[t-1, ]
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


set.seed(100)

x <- seq(0, 1, length.out = 200)

test1 <- rnorm(12, 0, 1) #alpha0 vector
test2 <-  matrix(rnorm(10*12), 10, 12)  #Gamma0 matrix
test3 <- matrix(rnorm(12*12), 12, 12) + diag(0.5, 12)   #A0 matrix


alpha_test2 <- c(-2, 0.5, 1, 0, 1.5, -1, 0, 1.5, -1.5, -1, 1, 0)
Phi_test2 <- riwish(205, diag(200))
A_test2 <- matrix(rnorm(12*12), 12, 12) + diag(0.5, 12)
df_test2 <- simulation_test2(x, 50, 12, 0.1, alpha_test2, Phi_test2, A_test2, 12)


# PSI with MTMH and p=1 --------------------------------------------------------

tic()
out_psi_mtmh <- GibbsSampler_mtmh(df = df_test2, 
                                  a0 = 2, 
                                  b0 = 1, 
                                  avec0 = test1, 
                                  Sigma0 = diag(12), 
                                  V0 = diag(12), 
                                  Aprior = test3, 
                                  R = 2000,
                                  burnin = 500,
                                  nbasis = 5,
                                  nu0 = 14, 
                                  S0 = diag(10), 
                                  m = 30
)
toc()

#salvare in output anche le matrici fissate per generare le funzioni in modo da 
#poterle confrontare

