library(EDMeasure)
library(splines)

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
  Z <- rnorm(n)
  Z2 <- rnorm(n)
  Z3 <- rnorm(n)
  X <- Z + Z3
  Y <- Z + Z2 + r * X
  list(Z = Z, X = X, Y = Y)
}

build_basis <- function(Z, p, basis_type = c("bs", "ns"),  degree = 3,
                        boundary_knots) {
  if (basis_type == "bs") {
    B <- bs(Z, df = p, degree = degree, Boundary.knots = boundary_knots)
  } else {
    B <- ns(Z, df = p, Boundary.knots = boundary_knots)
  }
  B <- sweep(B, 2, colMeans(B), FUN = "-")
  B %*% matpower(cov(B), -1 / 2)
}

select_p <- function(Z, Y, p_method= c("gcv", "mallows"), p_candidates, basis_type = c("bs", "ns"), degree = 3,
                     boundary_knots) {
  p_method <- match.arg(p_method)
  basis_type <- match.arg(basis_type)
  scores <- numeric(length(p_candidates))
  n <- length(Y)
  
  for (j in seq_along(p_candidates)) {
    p_try <- p_candidates[j]
    B_try <- build_basis(Z = Z, p = p_try, basis_type = basis_type,
                         degree = degree, boundary_knots = boundary_knots)
    design_try <- cbind(1, B_try)
    qr_try <- qr(design_try)
    res_try <- qr.resid(qr_try, Y)
    rss <- sum(res_try^2)
    
    if (p_method == "gcv") {
      scores[j] <- rss / (1 - p_try / n)^2
    } else if (p_method == "mallows") {
      sigma2_hat <- rss / (n - p_try - 1)
      scores[j] <- rss + 2 * p_try * sigma2_hat
    }
  }
  
  p_candidates[which.min(scores)]
}


run_case <- function(p=NULL, i, n, BB = 499, r = 0, basis_type = c("bs", "ns"),
                     degree = 3,
                     boundary_knots, p_method = c("fixed", "gcv", "mallows"),
                     p_candidates = 3:15) {
  basis_type <- match.arg(basis_type)
  p_method <- match.arg(p_method)
  p_candidates <- p_candidates[p_candidates < n - 1]
  
  if (length(p_candidates) == 0) {
    stop("No valid p_candidates: need p < n - 1.")
  }
  sample_data <- generate_sample(i = i, n = n, r = r)
  Z <- sample_data$Z
  X <- sample_data$X
  Y <- sample_data$Y
  
  if (p_method == "fixed") {
    if (is.null(p)) stop("Provide p when p_method = 'fixed'.")
    p_used <- p
  } else {
    p_used <- select_p(Z = Z, Y = Y, p_method = p_method, p_candidates = p_candidates,
                       basis_type = basis_type, degree = degree, boundary_knots = boundary_knots)
  }
  
  if (basis_type == "bs") {
    B <- bs(Z, df = p_used, degree = degree, Boundary.knots = boundary_knots)
  } else {
    B <- ns(Z, df = p_used, Boundary.knots = boundary_knots)
  }
  
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
  
  pvalue <-sum(test_boot>test)/BB
  crit <- as.numeric(quantile(test_boot, c(0.9, 0.95, 0.99)))
  test_mdd <- as.numeric(test > crit)
  names(test_mdd) <- c("alpha_0.10", "alpha_0.05", "alpha_0.01")
  test_mdd
  
  list(
    rejection = test_mdd,
    pvalue = pvalue,
    test_stat = test,
    crit = crit,
    p_selected = p_used
  )
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
    results_by_p[[p_idx]][i, ] <- run_case(p = p, i = i, n=200, p_method = "fixed",basis_type = "bs",
                                           degree = 3, boundary_knots=c(-3.5,3.5))$rejection
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
    results_by_p[[p_idx]][i, ] <- run_case(p = p, i = i, n=200, p_method = "fixed",basis_type = "bs", r=0.5,
                                           degree = 3, boundary_knots=c(-3.5,3.5))$rejection
  }
}

powers <- t(vapply(results_by_p, colMeans, numeric(3)))
