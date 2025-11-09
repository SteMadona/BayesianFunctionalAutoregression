out <- readRDS("out6_11")

out$n_acc

#sigma
sigma_trace <- as.mcmc(out$sigma)
traceplot(log(sigma_trace))
summary(sigma_trace)

#alpha
alpha_trace_phi_mh <- mcmc(t(out$alpha))
apply(out$alpha, 1, mean)
alpha_test2 <- c(-2, 0.5, 1, 0, 1.5, -1, 0, 1.5, -1.5, -1, 1, 0)

#A
A_post_mean <- apply(out$A, c(1, 2), mean)

A_prova <- apply(out$A, c(1, 2), 
                 function(u) mean(abs(u - A_test)))
A_test <- make_stable_from_A(A_test2)

A_error <- abs(A_post_mean - A_test)

A_melt <- melt(A_error)
ggplot(A_melt, aes(Var2, Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  labs(title = "A Posterior Mean vs A True Value",
       x = NULL, y = NULL, fill = "|error|") +
  coord_fixed() +
  theme_minimal()

#A (p > 1)
plot_A_posteriors_array <- function(A_array, A_test_list, titles = NULL) {
  # Controlli di base
  stopifnot(length(dim(A_array)) == 4)
  
  k <- dim(A_array)[1]
  p <- dim(A_array)[3]
  n_iter <- dim(A_array)[4]
  
  stopifnot(length(A_test_list) == p)
  
  plots <- vector("list", p)
  
  for (i in seq_len(p)) {
    # Estrae i campioni per la matrice A_i (kxkxn_iter)
    A_i <- A_array[, , i, ]
    A_post_mean <- apply(A_i, c(1, 2), mean)
    
    # Matrice "vera"
    A_test <- make_stable_from_A(A_test_list[[i]])
    
    # Errore
    A_error <- abs(A_post_mean - A_test)
    A_melt <- melt(A_error)
    
    title_i <- if (!is.null(titles) && length(titles) >= i) {
      titles[i]
    } else {
      paste0("Matrix A[", i, "] : Posterior Mean vs True")
    }
    
    p_i <- ggplot(A_melt, aes(Var2, Var1, fill = value)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "red") +
      labs(title = title_i,
           x = NULL, y = NULL, fill = "|error|") +
      coord_fixed() +
      theme_minimal()
    
    plots[[i]] <- p_i
  }
  
  # Mostra tutti i plot insieme (massimo 2 per riga)
  do.call(grid.arrange, c(plots, ncol = min(p, 2)))
  
  invisible(plots)
}

plot_A_posteriors_array(out_psi_mtmh_p$A, A_list)



#Phi
Phi_post_mean <- apply(out$Phi, c(1, 2), mean)

Phi_error <- abs(Phi_post_mean - Phi_test2)

Phi_melt <- melt(Phi_error)
ggplot(Phi_melt, aes(Var2, Var1, fill = value)) + 
  geom_tile() + 
  scale_fill_gradient(low = "white", high = "red") + 
  labs(title = "Phi posterior mean vs Phi True Value", 
       x = NULL, y = NULL, fill = "|error|") + 
  coord_fixed() + 
  theme_minimal()



# POSTERIOR RECONSTRUCTION ------------------------------------------------

posterior_reconstruction <- function(out, fs, x, k, t_list = c(5, 10, 15, 20)) {
  # reconstruct curves for selected time points 
  # plots posterior mean and +- sd bands vs observed f_t
  
  library(ggplot2)
  library(splines)
  
  B <- bs(x, df = k, degree = 3)
  
  alpha_draws <- out$alpha
  Gamma_draws <- out$Gamma
  
  S <- ncol(alpha_draws)
  n <- length(x)
  T <- dim(Gamma_draws)[1]
  
  BA <- B %*% alpha_draws
  
  for (i in seq_along(t_list)) {
    t <- t_list[i]
    BGt <- B %*% Gamma_draws[t, , ]
    curves <- BA + BGt
    
    mean_curve <- apply(curves, 1, mean)
    sd_curve <- apply(curves, 1, sd)
    lo <- mean_curve - 2 * sd_curve
    hi <- mean_curve + 2 * sd_curve
    
    d <- data.frame(
      x = x,
      mean = mean_curve,
      lo = lo,
      hi = hi,
      obs = fs[, t]
    )
    
    plot <- ggplot(d, aes(x = x)) +
      # banda di credibilità
      geom_ribbon(aes(ymin = lo, ymax = hi, fill = "95% Credible Band"), alpha = 0.2) +
      # funzione stimata
      geom_line(aes(y = mean, color = "Posterior Mean"), size = 1.2) +
      # funzione reale
      geom_line(aes(y = obs, color = "Observed Function"), size = 0.8, alpha = 0.8) +
      # legende e colori
      scale_color_manual(values = c("Posterior Mean" = "blue",
                                    "Observed Function" = "red")) +
      scale_fill_manual(values = c("95% Credible Band" = "skyblue")) +
      labs(title = paste("Posterior Reconstruction at time t =", t),
           x = "x", y = "f(x, t)",
           color = NULL, fill = NULL) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
    
    print(plot)
  }
}


posterior_reconstruction(out, df_test2[, -1], x, 12)

apply(B %*% out$Gamma[15, ,] + B %*% out$alpha, 2, mean)


#start with a lower number of basis to get similar A matrix to the true one
#frobenius distance for distanza tra la matrice vera e le matrici
#simulate che ti lascia comunque l'identificabilità di A 
#(bisogna normalizzare per la grid)
#mi aspetto che la distanza aumenti all'aumentare del numero di 
#basi scelte


#per p > 1 partire dal caso p = 2 e vedere quando si rompe, 
#(il problema potrebbero essere le matrici che non sono 
#campionate in maniera indipendente)
#posso ridurmi al caso vettoriale per vedere come fnuziona 
#(che significa semplicemente ridurre a 1 il numero di basi che 
#sto utilizzando)

#sono due problemi diversi: l'identificabilità di A e il sampling 
#indipendente delle matrici A