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



simulation_test2 <- function(x, n_fun, nbasis, noise_sd = 0.5, alpha, Phi, A, k){
  # spline basis
  B <- bs(x, df = nbasis, degree = 3)
  C <- solve(crossprod(B)) %*% t(B)
  
  n_points <- length(x)
  n <- length(x)
  
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
                                  R = 200,
                                  burnin = 50,
                                  nbasis = 5,
                                  nu0 = 12, 
                                  S0 = diag(10), 
                                  m = 2000
)
toc()

#58855 second elapsed

sigma_trace <- as.mcmc(out_psi_mtmh$sigma)
out_psi_mtmh$n_acc  
traceplot(log(sigma_trace))

alpha_trace_phi_mh <- mcmc(t(out_psi_mtmh$alpha))

View(apply(out_psi_mtmh$A, c(1, 2), mean))
View(A_test2)

image(apply(out_psi_mtmh$Phi, c(1, 2), mean), 
      col = gray(seq(1, 0, length = 256)))
image(Phi_test2, col = gray(seq(1, 0, length = 256)))

View(out_psi_mtmh$Phi[, , 150])
apply(out_psi_mtmh$alpha, 1, mean)

out_psi_mtmh$n_acc
