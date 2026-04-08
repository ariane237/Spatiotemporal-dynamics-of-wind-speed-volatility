# ============================================================
# Script: 05_sdpd_mean_model.R
# Purpose:
#   Estimate the SDPD mean model with heterogeneous AR(1)
#   coefficients for each station, using multiple spatial
#   weight matrices.
#
# Inputs:
#   - STL + AR(1) processed panel (des_mat)
#   - Spatial weight matrices (distance, kNN, directional)
#
# Outputs:
#   - Residual matrices for each weight matrix
#   - Estimated parameters (rho, lambda, gamma_i)
#   - Master object for downstream volatility modelling
#
# Author: Ariane Meli Chrisko
# ============================================================

# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(data.table)
library(Rsolnp)
library(Matrix)

set.seed(123)

# ------------------------------------------------------------
# 2. Settings
# ------------------------------------------------------------

height    <- "ws10"   # "ws10" or "ws100"
drop_burn <- 8

input_dir  <- "data/processed"
weights_dir <- "results/spatial_weights/objects"
output_dir <- file.path("results/sdpd_mean", height)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 3. Load panel data
# ------------------------------------------------------------

load(file.path(input_dir, paste0("wind_STL_AR1_", height, ".RData")))

des_mat   <- get(paste0("des_", height, "_mat"))
dates_vec <- as.Date(get(paste0("dates_", height)))
stations  <- as.character(get(paste0("stations_", height)))
coords    <- get("coords")

# Convert to SDPD format (N x T)
Y <- t(des_mat)
N <- nrow(Y)

# ------------------------------------------------------------
# 4. Load spatial weight matrices
# ------------------------------------------------------------

W_distance <- readRDS(file.path(weights_dir, "W_distance_band_55km.rds"))$W
W_knn      <- readRDS(file.path(weights_dir, "W_knn_k5.rds"))$W
W_directional <- readRDS(file.path(weights_dir, "W_directional_advection.rds"))$W

# Ensure sparse format
to_sparse <- function(W) {
  if (inherits(W, "Matrix")) W else Matrix(W, sparse = TRUE)
}

W_distance    <- to_sparse(W_distance)
W_knn         <- to_sparse(W_knn)
W_directional <- to_sparse(W_directional)

# Remove diagonal (no self-loops)
diag(W_distance)    <- 0
diag(W_knn)         <- 0
diag(W_directional) <- 0

# ------------------------------------------------------------
# 5. Row-standardisation utilities
# ------------------------------------------------------------

is_rowstd <- function(W, tol = 1e-8) {
  rs <- Matrix::rowSums(W)
  max(abs(rs - 1)) < tol
}

row_standardise <- function(W) {
  rs <- Matrix::rowSums(W)
  rs[rs == 0] <- 1
  W / rs
}

ensure_rowstd <- function(W, name) {
  if (!is_rowstd(W)) {
    message(name, ": applying row-standardisation")
    W <- row_standardise(W)
  }
  W
}

W_distance    <- ensure_rowstd(W_distance, "Distance")
W_knn         <- ensure_rowstd(W_knn, "kNN")
W_directional <- ensure_rowstd(W_directional, "Directional")

# ------------------------------------------------------------
# 6. Sanity checks
# ------------------------------------------------------------

check_W <- function(W, N, name) {
  stopifnot(
    inherits(W, "Matrix"),
    nrow(W) == N,
    ncol(W) == N,
    all(is.finite(W@x))
  )
}

check_W(W_distance, N, "W_distance")
check_W(W_knn, N, "W_knn")
check_W(W_directional, N, "W_directional")

# ------------------------------------------------------------
# 7. Helper functions
# ------------------------------------------------------------

safe_gamma_hat <- function(X, RHS, eps = 1e-10) {
  num <- rowSums(X * RHS, na.rm = TRUE)
  den <- rowSums(X^2, na.rm = TRUE)
  den[den < eps] <- NA
  g <- num / den
  g[!is.finite(g)] <- 0
  pmin(pmax(g, -0.99), 0.99)
}

LL_conc <- function(par, Y, W, drop_burn = 0) {
  
  rho    <- par[1]
  lambda <- par[2]
  
  N  <- nrow(Y)
  TT <- ncol(Y)
  
  S <- Diagonal(N) - rho * W
  logdetS <- as.numeric(determinant(S, logarithm = TRUE)$modulus)
  
  if (!is.finite(logdetS)) return(1e10)
  
  X <- Y[, 1:(TT-1)]
  
  RHS <- sapply(2:TT, function(t) {
    as.numeric(S %*% Y[, t]) - lambda * as.numeric(W %*% Y[, t-1])
  })
  
  if (drop_burn > 0) {
    X   <- X[, (drop_burn+1):ncol(X)]
    RHS <- RHS[, (drop_burn+1):ncol(RHS)]
  }
  
  gamma_hat <- safe_gamma_hat(X, RHS)
  
  E <- RHS - gamma_hat * X
  sigma2 <- mean(E^2)
  
  - (ncol(X) * logdetS - length(E)/2 * log(sigma2))
}

# ------------------------------------------------------------
# 8. Estimation function
# ------------------------------------------------------------

fit_SDPD <- function(Y, W, name) {
  
  sol <- solnp(
    pars = c(0.2, 0.2),
    fun  = LL_conc,
    LB   = c(1e-6, 0),
    UB   = c(0.99, 1),
    Y = Y, W = W, drop_burn = drop_burn
  )
  
  list(
    name = name,
    rho = sol$pars[1],
    lambda = sol$pars[2],
    sol = sol
  )
}

# ------------------------------------------------------------
# 9. Run estimation for all W
# ------------------------------------------------------------

W_list <- list(
  Distance    = W_distance,
  kNN         = W_knn,
  Directional = W_directional
)

fits <- lapply(names(W_list), function(nm) {
  fit_SDPD(Y, W_list[[nm]], nm)
})
names(fits) <- names(W_list)

# ------------------------------------------------------------
# 10. Save outputs
# ------------------------------------------------------------

master <- list(
  panel = height,
  des_mat = des_mat,
  dates = dates_vec,
  stations = stations,
  coords = coords,
  fits = fits
)

saveRDS(master, file.path(output_dir, paste0("SDPD_mean_", height, ".rds")))

cat("SDPD estimation completed and saved.\n")