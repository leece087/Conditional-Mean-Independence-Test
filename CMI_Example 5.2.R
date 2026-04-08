library(EDMeasure)
library(splines)
library(MASS)

matpower <- function(a, alpha) {
  small <- 1e-8
  eig <- eigen(a)
  values <- eig$values
  vectors <- eig$vectors
  vectors <- vectors / t(matrix(sqrt(diag(crossprod(vectors))), nrow(a), nrow(a)))
  keep <- abs(values) > small
  values[keep] <- values[keep]^alpha
  vectors %*% diag(values) %*% t(vectors)
}

fit_with_qr <- function(design, response) {
  qr_design <- qr(design)
  list(
    qr = qr_design,
    fitted = qr.fitted(qr_design, response),
    residuals = qr.resid(qr_design, response)
  )
}

generate_sample <- function(i, n, r = 0) {
  set.seed(i)
  Sigma <- matrix(0.25, 2, 2); 
  diag(Sigma) <- 1
  data <- mvrnorm(n , rep(0, 2), Sigma)
  Z<-data[,1]; X<-data[,2]
  Y=(Z^2-1)+0.1*rnorm(n, 0, sd=1)+r*X
  list(Z = Z, X = X, Y = Y)
}

run_case <- function(p, i, n, BB = 499, r = 0) {
  sample_data <- generate_sample(i = i, n = n, r = r)
  Z <- sample_data$Z
  X <- sample_data$X
  Y <- sample_data$Y
  
  B <- bs(Z, df = p, degree = 3, Boundary.knots = c(-3.5, 3.5))
  B <- sweep(B, 2, colMeans(B), FUN = "-")
  B <- B %*% matpower(cov(B), -1 / 2)
  
  design <- cbind(1, B)
  fit_obj <- fit_with_qr(design, Y)
  fit <- fit_obj$fitted
  res <- fit_obj$residuals
  ZX <- cbind(Z, X)
  
  test <- n * mdd(ZX, res, compute = "C", center = "U")
  
  test_boot <- numeric(BB)
  for (k in seq_len(BB)) {
    set.seed(100 + k)
    Y_boot <- fit + res * rnorm(n)
    res_boot <- qr.resid(fit_obj$qr, Y_boot)
    test_boot[k] <- n * mdd(ZX, res_boot, compute = "C", center = "U")
  }
  
  crit <- as.numeric(quantile(test_boot, c(0.9, 0.95, 0.99)))
  test_mdd <- as.numeric(test > crit)
  names(test_mdd) <- c("alpha_0.10", "alpha_0.05", "alpha_0.01")
  test_mdd
}

results_by_p <- vector("list", 3)
names(results_by_p) <- paste0("p", c(10, 25, 50))

for (i in 1:500) {
  for (p_idx in seq_along(c(10, 25, 50))) {
    p <- c(10, 25, 50)[p_idx]
    if (i == 1) {
      results_by_p[[p_idx]] <- matrix(0, nrow = 500, ncol = 3)
      colnames(results_by_p[[p_idx]]) <- c("alpha_0.10", "alpha_0.05", "alpha_0.01")
    }
    results_by_p[[p_idx]][i, ] <- run_case(p = p, i = i, n=200)
  }
}

sizes <- t(vapply(results_by_p, colMeans, numeric(3)))


for (i in 1:500) {
  for (p_idx in seq_along(c(10, 25, 50))) {
    p <- c(10, 25, 50)[p_idx]
    if (i == 1) {
      results_by_p[[p_idx]] <- matrix(0, nrow = 500, ncol = 3)
      colnames(results_by_p[[p_idx]]) <- c("alpha_0.10", "alpha_0.05", "alpha_0.01")
    }
    results_by_p[[p_idx]][i, ] <- run_case(p = p, i = i, n=200, r=0.5)
  }
}

powers <- t(vapply(results_by_p, colMeans, numeric(3)))
