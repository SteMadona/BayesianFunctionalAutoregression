source("code/UsefulFunctions.R")
source("code/ConjugatePosteriors.R")
source("code/GibbsSamplers.R")
source("code/MHforPsi.R")



simulation_test1 <- function(x, n_fun, nbasis, noise_sd = 0.5, alpha, Phi){
  # spline basis
  B <- bs(x, df = nbasis, degree = 3)
  n_points <- length(x)
  
  n <- length(x)
  
  mu <- as.vector(B %*% alpha)
  
  Phi <- make_posdef(Phi)
  Sigma_eps <- (noise_sd^2) * Phi
  
  functions <- matrix(NA, nrow = n, ncol = n_fun)
  for (t in 1:n_fun) {
    f_t <- mvtnorm::rmvnorm(1, mean = mu, sigma = Sigma_eps)
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
test3 <- matrix(rnorm(12*12), 12, 12)   #A0 matrix


alpha_test1 <- c(-2, 0.5, 1, 0, 1.5, -1, 0, 1.5, -1.5, -1, 1, 0)
Phi_test1 <- riwish(205, diag(200))
df_test1 <- simulation_test1(x, 100, 12, 0.1, alpha_test1, Phi_test1)


# Plot
df_test1_long <- df_test1 %>% 
  pivot_longer(
    cols = starts_with("t"),   
    names_to = "funzione", 
    values_to = "valore"
  )

ggplot(df_test1_long, aes(x = x, y = valore, color = funzione)) +
  geom_line(size = 1) +
  theme_minimal() +
  labs(title = "Synthetic Functions", y = "f(x)", color = "Funzione")



# PSI with MH -------------------------------------------------------------

library(tictoc)
tic()
out_psi_mh <- gibbs_sample_mh(df = df_test1, 
                              a0 = 2, 
                              b0 = 1, 
                              avec0 = test1, 
                              Sigma0 = diag(12), 
                              V0 = diag(12), 
                              Aprior = test3, 
                              R = 2000,
                              burnin = 1000, 
                              nbasis = 5, 
                              nu0 = 12, 
                              S0 = diag(10)
)
toc()

#577.01 sec elapsed

out_psi_mh$n_acc  
sigma_trace_psi_mh <- mcmc(out_psi_mh$sigma)
summary((sigma_trace))
traceplot(log(sigma_trace))

alpha_trace_psi_mh <- mcmc(t(out_psi_mh$alpha))

#for A and Gamma look into ggplots (grid plot or boxplots for every component)

View(apply(out_psi_mh$A, c(1, 2), mean))
max(abs(apply(out_psi_mh$A, c(1, 2), mean)))
min(abs(apply(out_psi_mh$A, c(1, 2), mean)))
mean(abs(apply(out_psi_mh$A, c(1, 2), mean)))

View(apply(out_psi_mh$Phi, c(1, 2), mean))
apply(out_psi_mh$alpha, 1, mean)



# PHI with MH -------------------------------------------------------------

tic()
out_phi_mh <- gibbs_sample_phi_mh(df = df_test1, 
                                  a0 = 2, 
                                  b0 = 1, 
                                  avec0 = test1, 
                                  Sigma0 = diag(12), 
                                  V0 = diag(12), 
                                  Aprior = test3, 
                                  R = 2000,
                                  burnin = 1000, 
                                  nu0 = 210, 
                                  S0 = diag(200)
)
toc()

#749.01 sec elapsed

out_phi_mh$n_acc  
sigma_trace_phi_mh <- mcmc(out_phi_mh$sigma)
summary((sigma_trace))
traceplot(log(sigma_trace))

alpha_trace_phi_mh <- mcmc(t(out_phi_mh$alpha))

#for A and Gamma look into ggplots (grid plot or boxplots for every component)

View(apply(out_phi_mh$A, c(1, 2), mean))
max(abs(apply(out_phi_mh$A, c(1, 2), mean)))
min(abs(apply(out_phi_mh$A, c(1, 2), mean)))
mean(abs(apply(out_phi_mh$A, c(1, 2), mean)))

View(apply(out_phi_mh$Phi, c(1, 2), mean))
apply(out_phi_mh$alpha, 1, mean)



# PSI with MTMH -----------------------------------------------------------

tic()
out_psi_mtmh <- GibbsSampler_mtmh(df = df_test1, 
                                  a0 = 2, 
                                  b0 = 1, 
                                  avec0 = test1, 
                                  Sigma0 = diag(12), 
                                  V0 = diag(12), 
                                  Aprior = test3, 
                                  R = 200,
                                  burnin = 0,
                                  nbasis = 5,
                                  nu0 = 12, 
                                  S0 = diag(10), 
                                  m = 2000
)
toc()


out_psi_mtmh$n_acc  
sigma_trace_phi_mtmh <- mcmc(out_psi_mtmh$sigma)
summary((sigma_trace_phi_mtmh))
traceplot(log(sigma_trace_phi_mtmh))

alpha_trace_phi_mh <- mcmc(t(out_psi_mtmh$alpha))

#for A and Gamma look into ggplots (grid plot or boxplots for every component)

View(apply(out_psi_mtmh$A, c(1, 2), mean))
max(abs(apply(out_psi_mtmh$A, c(1, 2), mean)))
min(abs(apply(out_psi_mtmh$A, c(1, 2), mean)))
mean(abs(apply(out_psi_mtmh$A, c(1, 2), mean)))

View(apply(out_psi_mtmh$Phi, c(1, 2), mean))
apply(out_psi_mtmh$alpha, 1, mean)
View(Phi_test1)

View(out_psi_mtmh$Phi[, ,800])

saveRDS(out_psi_mtmh, file = "results/out_psi_mtmh.rds")

out_prova <- readRDS("results/out_psi_mtmh.rds")


#PSI with MTMH and mulitple orders

tic()
out_psi_mtmh_p <- GS_mtmh_p(df = df_test1, 
                            a0 = 2, 
                            b0 = 1, 
                            avec0 = test1, 
                            Sigma0 = diag(12), 
                            V0 = diag(12), 
                            Aprior = test3, 
                            p = 10, 
                            R = 10,
                            burnin = 0,
                            nbasis = 5,
                            nu0 = 12, 
                            S0 = diag(10), 
                            m = 10
)
toc()

