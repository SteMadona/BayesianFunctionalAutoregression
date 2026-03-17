smooth_matrix_bspline <- function(Y, x = seq_len(nrow(Y)),
                                  nbasis = 40,
                                  lambda_grid = 10^seq(-6, 2, length.out = 25),
                                  norder = 4, Lfd = 2) {
  rng <- range(x)
  
  basis <- create.bspline.basis(rangeval = rng, nbasis = nbasis, norder = norder)
  Lfdobj <- int2Lfd(Lfd)
  
  gcv_mean <- sapply(lambda_grid, function(lam) {
    fdParobj <- fdPar(basis, Lfdobj = Lfdobj, lambda = lam)
    sm <- smooth.basis(argvals = x, y = Y, fdParobj = fdParobj)
    mean(sm$gcv, na.rm = TRUE)
  })
  lambda_star <- lambda_grid[which.min(gcv_mean)]
  
  fdPar_star <- fdPar(basis, Lfdobj = Lfdobj, lambda = lambda_star)
  sm_star <- smooth.basis(argvals = x, y = Y, fdParobj = fdPar_star)
  
  Y_smooth <- eval.fd(x, sm_star$fd)
  
  list(Y_smooth = Y_smooth, fd = sm_star$fd, lambda = lambda_star,
       lambda_grid = lambda_grid, gcv_mean = gcv_mean)
}


patient <- read.csv("paziente155.csv", header = F)
n <- nrow(patient)
x <- seq(1, 288)
indexes <- seq(from = 0, to = n, length.out = n/3)

patient <- as.matrix(patient)

patient_aligned <- align_fPCA(patient[indexes, ], x[indexes])

patient_smoothed <- smooth_matrix_bspline(patient_aligned$fn, x_aligned, nbasis = 40)
patient_smooth <- patient_smoothed$Y_smooth

