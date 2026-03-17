set.seed(999)

x_grid <- seq(0, 1, length.out = 50)
k_dim <- 10
p_gen <- 1
reps <- 5
nbasis <- 4

D <- make_fourier(x_grid, nbasis)
Psi_true <- diag(2*nbasis)
Phi_true <- make_Phi(D, Psi_true, rel = 1e-2)

#Gamma0_true <- matrix(rnorm(p_gen*k_dim, 0, 3), p_gen, k_dim)
Gamma0_true <- rnorm(k_dim, 0, 2)

TT <- 50 
noise <- 0.01  
R <- 4e4
burnin <- 3.5e4

set.seed(123)

alpha_list <- list()
A_list <- list()

for(r in 1:10){
  alpha_list[[r]] <- rnorm(k_dim, 0, 2)
  
  A <- diag(runif(k_dim, 0.35, 0.95))
  
  off_diag_positions <- which(row(A) != col(A), arr.ind = TRUE)
  
  selected_positions <- off_diag_positions[
    sample(nrow(off_diag_positions), 10),
  ]
  
  for(i in 1:10){
    A[selected_positions[i, 1], selected_positions[i, 2]] <-
      round(runif(1, -0.75, 0.75), 2)
  }
  
  A_list[[r]] <- A
}


rho_vals <- sapply(A_list, eig)      
rho_vals
A_list <- A_list[rho_vals < 1 & rho_vals > 0.9] 
alpha_list <- alpha_list[rho_vals < 1 & rho_vals > 0.9]

results321 <- list()

for(i in 1:reps){
  
  # 1. Data generation and visualization
  sim_data <- run_sim(TT, noise, alpha_list[[i]], A_list[[i]], Phi_true, Gamma0_true, x_grid, k_dim, p_gen)
  
  x <- sim_data$data$x
  Y <- as.matrix(sim_data$data[ , -1])   
  cols <- adjustcolor(sample(colors(), ncol(Y), replace = TRUE), alpha.f = 0.6)
  ltys <- sample(1:6, ncol(Y), replace = TRUE)
  
  plot(x, Y[,1], type="n",
       xlab = "x", ylab = "meanValue",
       ylim = range(Y))
  
  for(j in 1:ncol(Y)){
    lines(x, Y[,j], col = cols[j], lty = ltys[j], lwd = 1)
  }
  
  lines(x, rowMeans(Y), col = "black", lwd = 4)
  
  
  # 2. Sampler
  tic()
  res <- GibbsSampler_p(sim_data$data, 
                        R=R, 
                        burnin=burnin, 
                        k = k_dim, 
                        nbasis = nbasis, 
                        p_ar = p_gen)
  time_elapsed <- toc()
  
  
  # 3. Partial results
  
  alpha_est_mean <- apply(res$alpha, 1, mean)
  alpha_est_med <- apply(res$alpha, 1, median)
  mse_alpha <- mean((alpha_est_mean - alpha_list[[i]])^2)
  mse_alpha2 <- mean((alpha_est_med - alpha_list[[i]])^2)
  
  A_est <- apply(res$A, c(1, 2), mean)
  mse_A <- mean((A_list[[i]] - A_est)^2)
  A_est_med <- apply(res$A, c(1, 2), median)
  mse_A_med <- mean((A_list[[i]] - A_est_med)^2)
  
  cat("MSE alpha mean :", mse_alpha, "\n", 
      "MSE alpha med :", mse_alpha2, "\n", 
      "MSE A mean:", mse_A, "\n", 
      "MSE A med:", mse_A_med, "\n")
  
  results321[[i]] <- list(
    results = res, 
    true = A_list[[i]], 
    data = sim_data
  )
}




