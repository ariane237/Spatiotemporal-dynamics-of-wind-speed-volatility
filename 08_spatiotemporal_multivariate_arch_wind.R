## =========================================================
## COMPLETE PIPELINE:
## Multivariate ST mean  -> residuals -> vec-spARCH_p volatility
## Train/test split, diagnostics, OOS forecasting, evaluation
## =========================================================

rm(list = ls())

library(Matrix)
library(spdep)
library(zoo)
library(numDeriv)
library(Rsolnp)

## =========================================================
## 0. USER SETTINGS
## =========================================================

tau_mean <- 1
tau_vol  <- 1
eps_log  <- 1e-6
lambda_ewma <- 0.94
use_same_W_for_mean_and_vol <- TRUE

## =========================================================
## 1. HELPER FUNCTIONS
## =========================================================

safe_logsq <- function(x, eps = 1e-6) {
  log(pmax(x^2, eps))
}

vec <- function(x) {
  as.vector(x)
}

combine_two_vars_to_array <- function(mat1, mat2) {
  stopifnot(identical(dim(mat1), dim(mat2)))
  n  <- nrow(mat1)
  TT <- ncol(mat1)
  
  out <- array(NA_real_, dim = c(n, 2, TT))
  out[, 1, ] <- mat1
  out[, 2, ] <- mat2
  out
}

build_4d_intercept <- function(n, p, Tt) {
  array(1, dim = c(n, p, Tt, 1))
}

matrix_to_listw <- function(W) {
  Wm <- as.matrix(W)
  diag(Wm) <- 0
  spdep::mat2listw(Wm, style = "W", zero.policy = TRUE)
}

make_realized_proxies <- function(res_mat, lambda = 0.94) {
  n  <- nrow(res_mat)
  TT <- ncol(res_mat)
  
  RV <- res_mat^2
  
  RV5_sq <- t(apply(RV, 1, function(x) {
    zoo::rollapply(x, width = 5, FUN = mean, fill = NA, align = "right")
  }))
  
  RV5_abs <- t(apply(abs(res_mat), 1, function(x) {
    zoo::rollapply(x, width = 5, FUN = mean, fill = NA, align = "right")
  }))^2
  
  EWMA <- matrix(NA_real_, n, TT)
  EWMA[, 1] <- RV[, 1]
  if (TT >= 2) {
    for (tt in 2:TT) {
      EWMA[, tt] <- lambda * EWMA[, tt - 1] + (1 - lambda) * RV[, tt]
    }
  }
  
  list(
    RV = RV,
    RV5_sq = RV5_sq,
    RV5_abs = RV5_abs,
    EWMA = EWMA
  )
}




#make_realized_proxies <- function(res_mat, lambda = 0.94) {
#  n  <- nrow(res_mat)
#  TT <- ncol(res_mat)

#  RV <- res_mat^2

#  RV5_sq <- t(apply(RV, 1, function(x) {
#    zoo::rollapply(x, width = 5, FUN = mean, fill = NA, align = "right")
#  }))

#  RV5_abs <- t(apply(abs(res_mat), 1, function(x) {
#    zoo::rollapply(x, width = 5, FUN = mean, fill = NA, align = "right")
#  }))^2

#  EWMA <- matrix(NA_real_, n, TT)
#  EWMA[, 1] <- RV[, 1]
#  if (TT >= 2) {
#    for (tt in 2:TT) {
#      EWMA[, tt] <- lambda * EWMA[, tt - 1] + (1 - lambda) * RV[, tt - 1]
#    }
#  }

#  list(
#    RV = RV,
#    RV5_sq = RV5_sq,
#    RV5_abs = RV5_abs,
#    EWMA = EWMA
#  )
#}
evaluate_forecasts <- function(Hhat, proxy10, proxy100, on_log_scale = TRUE) {
  proxy_list <- list(
    RV      = combine_two_vars_to_array(proxy10$RV,      proxy100$RV),
    RV5_sq  = combine_two_vars_to_array(proxy10$RV5_sq,  proxy100$RV5_sq),
    RV5_abs = combine_two_vars_to_array(proxy10$RV5_abs, proxy100$RV5_abs),
    EWMA    = combine_two_vars_to_array(proxy10$EWMA,    proxy100$EWMA)
  )
  
  out <- data.frame()
  
  for (nm in names(proxy_list)) {
    R <- proxy_list[[nm]]
    
    for (j in 1:2) {
      h <- Hhat[, j, ]
      r <- R[, j, ]
      
      ok <- is.finite(h) & is.finite(r) & (h > 0) & (r > 0)
      
      if (on_log_scale) {
        e <- log(h[ok]) - log(r[ok])
        rmsfe <- sqrt(mean(e^2))
        mafe  <- mean(abs(e))
      } else {
        e <- h[ok] - r[ok]
        rmsfe <- sqrt(mean(e^2))
        mafe  <- mean(abs(e))
      }
      
      out <- rbind(out, data.frame(
        variable = c("ws10", "ws100")[j],
        proxy    = nm,
        RMSFE    = rmsfe,
        MAFE     = mafe,
        QLIKE    = mean(log(h[ok]) + r[ok] / h[ok]),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  out
}

ljung_box_summary <- function(U_arr, lag = 10) {
  n  <- dim(U_arr)[1]
  p  <- dim(U_arr)[2]
  
  out <- data.frame(
    variable = character(),
    lb_resid_pass = numeric(),
    lb_sq_pass = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (j in 1:p) {
    pval_res <- rep(NA_real_, n)
    pval_sq  <- rep(NA_real_, n)
    
    for (i in 1:n) {
      x <- as.numeric(U_arr[i, j, ])
      x <- x[is.finite(x)]
      
      if (length(x) > lag + 5) {
        pval_res[i] <- Box.test(x, lag = lag, type = "Ljung-Box")$p.value
        pval_sq[i]  <- Box.test(x^2, lag = lag, type = "Ljung-Box")$p.value
      }
    }
    
    out <- rbind(out, data.frame(
      variable = c("ws10", "ws100")[j],
      lb_resid_pass = mean(pval_res > 0.05, na.rm = TRUE),
      lb_sq_pass    = mean(pval_sq  > 0.05, na.rm = TRUE),
      stringsAsFactors = FALSE
    ))
  }
  
  out
}

moran_summary <- function(U_arr, listw_obj) {
  p  <- dim(U_arr)[2]
  TT <- dim(U_arr)[3]
  
  out <- data.frame(
    variable = character(),
    moran_resid_pass = numeric(),
    moran_sq_pass = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (j in 1:p) {
    pval_res <- rep(NA_real_, TT)
    pval_sq  <- rep(NA_real_, TT)
    
    for (tt in 1:TT) {
      x <- U_arr[, j, tt]
      if (all(is.finite(x))) {
        pval_res[tt] <- spdep::moran.test(x, listw_obj, zero.policy = TRUE)$p.value
        pval_sq[tt]  <- spdep::moran.test(x^2, listw_obj, zero.policy = TRUE)$p.value
      }
    }
    
    out <- rbind(out, data.frame(
      variable = c("ws10", "ws100")[j],
      moran_resid_pass = mean(pval_res > 0.05, na.rm = TRUE),
      moran_sq_pass    = mean(pval_sq  > 0.05, na.rm = TRUE),
      stringsAsFactors = FALSE
    ))
  }
  
  out
}

## =========================================================
## 2. MEAN MODEL FUNCTIONS
## =========================================================

qml_spatiotemporal_multivariate_p <- function(Y, W, X, tau = 1, ...) {
  
  dimY <- dim(Y)
  n    <- dimY[1]
  p    <- dimY[2]
  Tt   <- dimY[3]
  k    <- dim(X)[4]
  
  if (!all(dim(W) == c(n, n))) stop("Dimension of W is wrong")
  
  tau <- sort(unique(as.integer(tau)))
  if (any(tau < 1)) stop("All tau must be >= 1")
  max_tau <- max(tau)
  if (Tt <= max_tau) stop("Need T > max(tau)")
  
  vec <- function(x) as.vector(x)
  L <- length(tau)
  
  idx_Psi  <- 1:(p^2)
  idx_Pi   <- (max(idx_Psi) + 1):(max(idx_Psi) + L * p^2)
  idx_Beta <- (max(idx_Pi)  + 1):(max(idx_Pi)  + k * p)
  idx_sig  <- max(idx_Beta) + 1
  
  unpack_pars <- function(pars) {
    Psi <- matrix(pars[idx_Psi], p, p)
    Pi_vec <- pars[idx_Pi]
    Pi_arr <- array(Pi_vec, dim = c(p, p, L))
    Beta <- matrix(pars[idx_Beta], k, p)
    sig_u <- pars[idx_sig]
    list(Psi = Psi, Pi_arr = Pi_arr, Beta = Beta, sig_u = sig_u)
  }
  
  LogLikelihood <- function(pars, Y, W, X) {
    
    up <- unpack_pars(pars)
    Psi   <- up$Psi
    PiArr <- up$Pi_arr
    Beta  <- up$Beta
    sig_u <- up$sig_u
    
    S <- Matrix::Matrix(diag(n * p) - t(Psi) %x% W)
    
    rc <- suppressWarnings(rcond(as.matrix(S)))
    if (!is.finite(rc) || rc < 1e-8) return(1e12)
    
    detS <- tryCatch(Matrix::determinant(S, logarithm = TRUE), error = function(e) NULL)
    if (is.null(detS)) return(1e12)
    log_det_S <- as.numeric(detS$modulus)
    if (!is.finite(log_det_S)) return(1e12)
    
    sum_eps_2 <- 0
    
    for (tt in (max_tau + 1):Tt) {
      
      constant <- Reduce("+", lapply(1:k, function(x)
        X[, , tt, x] * matrix(Beta[x, ], n, p, byrow = TRUE)
      ))
      
      temporal_term <- matrix(0, n, p)
      for (l in 1:L) {
        temporal_term <- temporal_term + Y[, , tt - tau[l]] %*% PiArr[, , l]
      }
      
      vec_u_t <- S %*% vec(Y[, , tt]) - vec(constant) - vec(temporal_term)
      if (any(!is.finite(vec_u_t))) return(1e12)
      
      sum_eps_2 <- sum_eps_2 + sum(vec_u_t^2)
    }
    
    timept <- Tt - max_tau
    
    LL <- - (timept * n * p) / 2 * log(2 * pi) +
      timept * log_det_S -
      (timept * n * p) / 2 * log(sig_u) -
      (1 / (2 * sig_u)) * sum_eps_2
    
    (-1) * LL
  }
  
  residuals_fun <- function(pars, Y, W, X) {
    
    up <- unpack_pars(pars)
    Psi   <- up$Psi
    PiArr <- up$Pi_arr
    Beta  <- up$Beta
    
    S <- diag(n * p) - t(Psi) %x% W
    U_t <- array(NA_real_, dim = c(n, p, Tt))
    
    for (tt in (max_tau + 1):Tt) {
      
      constant <- Reduce("+", lapply(1:k, function(x)
        X[, , tt, x] * matrix(Beta[x, ], n, p, byrow = TRUE)
      ))
      
      temporal_term <- matrix(0, n, p)
      for (l in 1:L) {
        temporal_term <- temporal_term + Y[, , tt - tau[l]] %*% PiArr[, , l]
      }
      
      vec_u_t <- S %*% vec(Y[, , tt]) - vec(constant) - vec(temporal_term)
      U_t[, , tt] <- matrix(vec_u_t, n, p)
    }
    
    U_t
  }
  
  start_Psi  <- matrix(rep(0.1, p^2), p, p)
  start_Pi   <- array(0.2, dim = c(p, p, L))
  start_Beta <- matrix(rep(0.1, k * p), k, p)
  start_pars <- c(vec(start_Psi), vec(start_Pi), vec(start_Beta), 1)
  
  LB_Psi  <- matrix(rep(-0.8, p^2), p, p)
  LB_Pi   <- array(rep(-0.95, L * p^2), dim = c(p, p, L))
  LB_Beta <- matrix(rep(-1000, k * p), k, p)
  LB      <- c(vec(LB_Psi), vec(LB_Pi), vec(LB_Beta), 1e-5)
  
  UB_Psi  <- matrix(rep(0.95, p^2), p, p)
  UB_Pi   <- array(rep(0.95, L * p^2), dim = c(p, p, L))
  UB_Beta <- matrix(rep(1000, k * p), k, p)
  UB      <- c(vec(UB_Psi), vec(UB_Pi), vec(UB_Beta), 10000)
  
  out <- Rsolnp::solnp(
    pars    = start_pars,
    fun     = LogLikelihood,
    Y = Y, W = W, X = X,
    control = list(trace = 0),
    LB = LB, UB = UB
  )
  
  res <- residuals_fun(out$pars, Y = Y, W = W, X = X)
  
  up_hat <- unpack_pars(out$pars)
  
  theta_hat <- out$pars
  
  H <- tryCatch(
    numDeriv::hessian(
      func = LogLikelihood,
      x    = theta_hat,
      Y = Y, W = W, X = X
    ),
    error = function(e) matrix(NA_real_, length(theta_hat), length(theta_hat))
  )
  
  Hsym <- 0.5 * (H + t(H))
  Hinv <- tryCatch(solve(Hsym), error = function(e) NULL)
  if (is.null(Hinv)) {
    Hinv <- tryCatch(
      solve(Hsym + diag(1e-8, nrow(Hsym))),
      error = function(e) matrix(NA_real_, nrow(Hsym), ncol(Hsym))
    )
  }
  
  list(
    Psi_est   = up_hat$Psi,
    Pi_arr_est = up_hat$Pi_arr,
    Beta_est  = up_hat$Beta,
    sig_u_est = up_hat$sig_u,
    pars_est  = theta_hat,
    objective = out$values[length(out$values)],
    conv      = out$convergence,
    residuals = res,
    H         = H,
    Hsym      = Hsym,
    Hinv      = Hinv,
    tau       = tau,
    n         = n,
    p         = p,
    Tt        = Tt,
    k         = k
  )
}

forecast_mean_one_step <- function(fit, Y_hist, X_next, W) {
  # Y_hist: list of lagged n x p matrices, most recent last
  # X_next: n x p x k array reduced to n x p if k=1 represented in [,,x]
  n <- fit$n
  p <- fit$p
  k <- fit$k
  tau <- fit$tau
  L <- length(tau)
  
  Psi  <- fit$Psi_est
  PiArr <- fit$Pi_arr_est
  Beta <- fit$Beta_est
  
  S <- diag(n * p) - t(Psi) %x% as.matrix(W)
  
  const_mat <- matrix(0, n, p)
  for (x in 1:k) {
    const_mat <- const_mat + X_next[, , x] * matrix(Beta[x, ], n, p, byrow = TRUE)
  }
  
  temporal_term <- matrix(0, n, p)
  for (l in 1:L) {
    temporal_term <- temporal_term + Y_hist[[l]] %*% PiArr[, , l]
  }
  
  rhs <- vec(const_mat + temporal_term)
  yhat_vec <- solve(S, rhs)
  matrix(yhat_vec, nrow = n, ncol = p)
}

recursive_mean_test_residuals <- function(fit, Y_train, Y_test, X_test, W) {
  n <- dim(Y_test)[1]
  p <- dim(Y_test)[2]
  TT_test <- dim(Y_test)[3]
  tau <- fit$tau
  max_tau <- max(tau)
  
  U_test <- array(NA_real_, dim = c(n, p, TT_test))
  
  history_list <- vector("list", length(tau))
  for (l in seq_along(tau)) {
    history_list[[l]] <- Y_train[, , dim(Y_train)[3] - tau[l] + 1]
  }
  
  for (tt in 1:TT_test) {
    X_next <- X_test[, , tt, , drop = FALSE]
    X_next <- array(X_next, dim = c(n, p, dim(X_test)[4]))
    
    yhat <- forecast_mean_one_step(
      fit = fit,
      Y_hist = history_list,
      X_next = X_next,
      W = W
    )
    
    ytrue <- Y_test[, , tt]
    U_test[, , tt] <- ytrue - yhat
    
    if (length(tau) == 1) {
      history_list[[1]] <- ytrue
    } else {
      old_hist <- history_list
      for (l in seq_along(tau)) {
        history_list[[l]] <- if (l == 1) ytrue else old_hist[[l - 1]]
      }
    }
  }
  
  U_test
}

## =========================================================
## 3. SAFER vec-spARCH_p FUNCTION
## =========================================================

qml_vec_spARCH_p_safe <- function(Y, W, errortype = "norm",
                                  n_starts = 5,
                                  trace = FALSE,
                                  eps = 1e-6) {
  # Y is expected to be residual array n x p x T (not already log-squared)
  # model is fitted to Z_t = log(Y_t^2 + eps)
  
  dimY <- dim(Y)
  n    <- dimY[1]
  p    <- dimY[2]
  Tt   <- dimY[3]
  
  if (!all(dim(W) == c(n, n))) stop("Dimension of W is wrong")
  
  if (errortype == "norm") {
    E_logsq_errors <- digamma(1) - log(2)
    sig_u2 <- 4.934
  } else if (errortype == "t") {
    E_logsq_errors <- -0.901
    sig_u2 <- 5.870
  } else {
    stop("errortype must be 'norm' or 't'")
  }
  
  Z <- safe_logsq(Y, eps = eps)
  
  LogLikelihood <- function(pars, Z, W, sig_u2) {
    
    A_tilde <- matrix(rep(pars[1:p], n), n, p, byrow = TRUE)
    Psi     <- matrix(pars[(p + 1):(p^2 + p)], p, p)
    Pi      <- matrix(pars[(p^2 + p + 1):(2 * p^2 + p)], p, p)
    
    S <- diag(n * p) - t(Psi) %x% W
    rc <- suppressWarnings(rcond(S))
    if (!is.finite(rc) || rc < 1e-8) return(1e12)
    
    detS <- tryCatch(determinant(S, logarithm = TRUE), error = function(e) NULL)
    if (is.null(detS)) return(1e12)
    log_det_S <- as.numeric(detS$modulus)
    if (!is.finite(log_det_S)) return(1e12)
    
    sum_eps_2 <- 0
    for (tt in 2:Tt) {
      vec_u_t <- S %*% vec(Z[, , tt]) - vec(A_tilde) -
        (diag(p) %x% Z[, , tt - 1]) %*% vec(Pi)
      
      if (any(!is.finite(vec_u_t))) return(1e12)
      sum_eps_2 <- sum_eps_2 + sum(vec_u_t^2)
    }
    
    LL <- -(Tt - 1) * n * p / 2 * log(2 * pi) -
      (Tt - 1) * n * p / 2 * log(sig_u2) +
      (Tt - 1) * log_det_S -
      (1 / (2 * sig_u2)) * sum_eps_2
    
    (-1) * LL
  }
  
  residuals_fun <- function(pars, Z, W) {
    A_tilde <- matrix(rep(pars[1:p], n), n, p, byrow = TRUE)
    Psi     <- matrix(pars[(p + 1):(p^2 + p)], p, p)
    Pi      <- matrix(pars[(p^2 + p + 1):(2 * p^2 + p)], p, p)
    
    S <- diag(n * p) - t(Psi) %x% W
    U_t <- array(NA_real_, dim = c(n, p, Tt))
    
    for (tt in 2:Tt) {
      vec_u_t <- S %*% vec(Z[, , tt]) - vec(A_tilde) -
        (diag(p) %x% Z[, , tt - 1]) %*% vec(Pi)
      U_t[, , tt] <- matrix(vec_u_t, n, p)
    }
    U_t
  }
  
  make_start <- function() {
    start_a   <- runif(p, 0.5, 1.5)
    start_Psi <- matrix(runif(p^2, -0.15, 0.15), p, p)
    start_Pi  <- matrix(runif(p^2,  0.05, 0.35), p, p)
    c(start_a, vec(start_Psi), vec(start_Pi))
  }
  
  LB <- c(rep(-20, p), rep(-0.8, p^2), rep(-0.95, p^2))
  UB <- c(rep( 20, p), rep( 0.8, p^2), rep( 0.95, p^2))
  
  best_out <- NULL
  best_val <- Inf
  
  for (s in 1:n_starts) {
    start_pars <- make_start()
    
    out <- tryCatch(
      Rsolnp::solnp(
        pars = start_pars,
        fun  = LogLikelihood,
        Z = Z, W = W, sig_u2 = sig_u2,
        LB = LB, UB = UB,
        control = list(trace = trace)
      ),
      error = function(e) NULL
    )
    
    if (!is.null(out)) {
      val <- out$values[length(out$values)]
      if (is.finite(val) && val < best_val) {
        best_val <- val
        best_out <- out
      }
    }
  }
  
  if (is.null(best_out)) stop("All optimisation attempts failed in qml_vec_spARCH_p_safe().")
  
  res <- residuals_fun(best_out$pars, Z = Z, W = W)
  
  A_tilde_est <- matrix(rep(best_out$pars[1:p], n), n, p, byrow = TRUE)
  Psi_est     <- matrix(best_out$pars[(p + 1):(p^2 + p)], p, p)
  Pi_est      <- matrix(best_out$pars[(p^2 + p + 1):(2 * p^2 + p)], p, p)
  
  list(
    A_tilde_est = A_tilde_est,
    A_est       = A_tilde_est - E_logsq_errors,
    Psi_est     = Psi_est,
    Pi_est      = Pi_est,
    objective   = best_val,
    residuals   = res,
    Z_used      = Z,
    errortype   = errortype,
    sig_u2      = sig_u2,
    E_logsq_errors = E_logsq_errors,
    n = n,
    p = p,
    Tt = Tt
  )
}



forecast_vec_spARCH_one_step <- function(fit, Z_prev, W) {
  n <- fit$n
  p <- fit$p
  
  A_tilde <- fit$A_tilde_est
  Psi     <- fit$Psi_est
  Pi      <- fit$Pi_est
  
  S <- diag(n * p) - t(Psi) %x% as.matrix(W)
  rhs <- vec(A_tilde + Z_prev %*% Pi)
  
  zhat_vec <- solve(S, rhs)
  zhat <- matrix(zhat_vec, n, p)
  
  # zhat forecasts log(epsilon^2), not directly log(H)
  log_hhat <- zhat - fit$E_logsq_errors
  hhat <- exp(log_hhat)
  
  list(zhat = zhat, log_hhat = log_hhat, hhat = hhat)
}



## =========================================================
## 4. LOAD DATA
## =========================================================

env10  <- new.env()
env100 <- new.env()

load("wind_STL_AR1_ws10.RData",  envir = env10)
load("wind_STL_AR1_ws100.RData", envir = env100)

des10  <- t(env10$des_ws10_mat)    # station x time
des100 <- t(env100$des_ws100_mat)  # station x time

dates_all <- as.Date(env10$dates_ws10)

stopifnot(identical(env10$stations_ws10, env100$stations_ws100))
stopifnot(identical(env10$dates_ws10, env100$dates_ws100))
stopifnot(identical(dim(des10), dim(des100)))

station_ids <- env10$stations_ws10
n_total <- nrow(des10)
TT      <- ncol(des10)

cat("Data dimensions:\n")
cat("des10 :", dim(des10), "\n")
cat("des100:", dim(des100), "\n")

## =========================================================
## 5. LOAD SPATIAL WEIGHT MATRICES
## =========================================================

distance <- readRDS("W_distance_band_r55km_sparse_meta_ws10.rds")
knn      <- readRDS("W_5NN_k5_sparse_meta_ws10.rds")
direct   <- readRDS("W_directional_combined.rds")

weight_matrices <- list(
  distance = Matrix::Matrix(distance$W, sparse = TRUE),
  knn      = Matrix::Matrix(knn$W, sparse = TRUE),
  direct   = Matrix::Matrix(direct, sparse = TRUE)
)



for (nm in names(weight_matrices)) {
  stopifnot(all(dim(weight_matrices[[nm]]) == c(n_total, n_total)))
}


check_W <- function(W, name) {
  Wm <- as.matrix(W)
  diag(Wm) <- 0
  rs <- rowSums(Wm)
  cat("\n", name, "\n")
  cat("Range row sums:", range(rs), "\n")
  cat("Max diagonal:", max(abs(diag(Wm))), "\n")
  cat("Any NA:", anyNA(Wm), "\n")
}

lapply(names(weight_matrices), function(nm) check_W(weight_matrices[[nm]], nm))

## =========================================================
## 6. TRAIN / TEST SPLIT
## =========================================================

test_start <- max(dates_all) - 364

train_idx <- which(dates_all < test_start)
test_idx  <- which(dates_all >= test_start)

cat("\nTrain range:", as.character(min(dates_all[train_idx])),
    "to", as.character(max(dates_all[train_idx])), "\n")
cat("Test range :", as.character(min(dates_all[test_idx])),
    "to", as.character(max(dates_all[test_idx])), "\n")
cat("Train size :", length(train_idx), "\n")
cat("Test size  :", length(test_idx), "\n")

## =========================================================
## 7. BUILD MULTIVARIATE DESEASONALISED ARRAY
## =========================================================

Y_all <- combine_two_vars_to_array(des10, des100)  # n x 2 x T

Y_train <- Y_all[, , train_idx, drop = FALSE]
Y_test  <- Y_all[, , test_idx,  drop = FALSE]

X_train_mean <- build_4d_intercept(
  n = dim(Y_train)[1],
  p = dim(Y_train)[2],
  Tt = dim(Y_train)[3]
)

X_test_mean <- build_4d_intercept(
  n = dim(Y_test)[1],
  p = dim(Y_test)[2],
  Tt = dim(Y_test)[3]
)

## =========================================================
## 8. FIT MULTIVARIATE SPATIOTEMPORAL MEAN MODEL
## =========================================================
set.seed(123)

fit_mean_all <- list()
mean_diag_all <- list()
mean_test_resid_all <- list()

for (wm in names(weight_matrices)) {
  cat("\n============================\n")
  cat("Fitting mean model for:", wm, "\n")
  cat("============================\n")
  
  Wm <- weight_matrices[[wm]]
  
  fit_mean_all[[wm]] <- qml_spatiotemporal_multivariate_p(
    Y   = Y_train,
    W   = Wm,
    X   = X_train_mean,
    tau = tau_mean
  )
  
  cat("Convergence:", fit_mean_all[[wm]]$conv, "\n")
  cat("Objective  :", fit_mean_all[[wm]]$objective, "\n")
  
  U_train <- fit_mean_all[[wm]]$residuals
  valid_train <- (max(fit_mean_all[[wm]]$tau) + 1):dim(U_train)[3]
  U_train_use <- U_train[, , valid_train, drop = FALSE]
  
  lw <- matrix_to_listw(Wm)
  
  mean_diag_all[[wm]] <- list(
    ljung = ljung_box_summary(U_train_use, lag = 10),
    moran = moran_summary(U_train_use, lw)
  )
  
  U_test <- recursive_mean_test_residuals(
    fit = fit_mean_all[[wm]],
    Y_train = Y_train,
    Y_test  = Y_test,
    X_test  = X_test_mean,
    W = Wm
  )
  
  mean_test_resid_all[[wm]] <- U_test
}

## =========================================================
## 9. FIT vec-spARCH VOLATILITY MODEL ON MEAN RESIDUALS
## =========================================================

fit_vol_all <- list()
vol_diag_all <- list()
fcst_all <- list()
eval_results <- list()

for (wm in names(weight_matrices)) {
  cat("\n============================\n")
  cat("Fitting volatility model for:", wm, "\n")
  cat("============================\n")
  
  W_mean <- weight_matrices[[wm]]
  W_vol  <- if (use_same_W_for_mean_and_vol) weight_matrices[[wm]] else weight_matrices[[wm]]
  
  U_train_full <- fit_mean_all[[wm]]$residuals
  valid_train <- (max(fit_mean_all[[wm]]$tau) + 1):dim(U_train_full)[3]
  U_train <- U_train_full[, , valid_train, drop = FALSE]
  
  fit_vol_all[[wm]] <- qml_vec_spARCH_p_safe(
    Y = U_train,
    W = W_vol,
    errortype = "norm",
    n_starts = 5,
    trace = FALSE,
    eps = eps_log
  )
  
  cat("Vol objective:", fit_vol_all[[wm]]$objective, "\n")
  
  lw <- matrix_to_listw(W_vol)
  vol_diag_all[[wm]] <- list(
    ljung = ljung_box_summary(fit_vol_all[[wm]]$residuals[, , -1, drop = FALSE], lag = 10),
    moran = moran_summary(fit_vol_all[[wm]]$residuals[, , -1, drop = FALSE], lw)
  )
  
  ## --- OOS volatility forecasts ---
  U_test <- mean_test_resid_all[[wm]]
  TT_test <- dim(U_test)[3]
  
  Hhat <- array(NA_real_, dim = c(n_total, 2, TT_test))
  
  ## forecast origin = last train residual
  U_last_train <- U_train[, , dim(U_train)[3]]
  Z_prev <- safe_logsq(U_last_train, eps = eps_log)
  
  for (tt in 1:TT_test) {
    ftt <- forecast_vec_spARCH_one_step(
      fit = fit_vol_all[[wm]],
      Z_prev = Z_prev,
      W = W_vol
    )
    
    Hhat[, , tt] <- ftt$hhat
    
    ## update with realised test residual from mean stage
    Z_prev <- safe_logsq(U_test[, , tt], eps = eps_log)
  }
  
  fcst_all[[wm]] <- Hhat
  
  ## realised proxies from same mean-model test residuals
  proxy10 <- make_realized_proxies(U_test[, 1, ], lambda = lambda_ewma)
  proxy100 <- make_realized_proxies(U_test[, 2, ], lambda = lambda_ewma)
  
  eval_results[[wm]] <- evaluate_forecasts(
    Hhat    = Hhat,
    proxy10 = proxy10,
    proxy100 = proxy100,
    on_log_scale = TRUE
  )
}

## =========================================================
## 10. COMBINE RESULTS
## =========================================================

final_eval <- do.call(rbind, lapply(names(eval_results), function(wm) {
  cbind(weight = wm, eval_results[[wm]])
}))
rownames(final_eval) <- NULL

final_mean_diag <- do.call(rbind, lapply(names(mean_diag_all), function(wm) {
  lj <- mean_diag_all[[wm]]$ljung
  mo <- mean_diag_all[[wm]]$moran
  
  rbind(
    data.frame(
      stage = "mean",
      weight = wm,
      test = "LjungBox",
      variable = lj$variable,
      resid_pass = lj$lb_resid_pass,
      sq_pass = lj$lb_sq_pass,
      stringsAsFactors = FALSE
    ),
    data.frame(
      stage = "mean",
      weight = wm,
      test = "Moran",
      variable = mo$variable,
      resid_pass = mo$moran_resid_pass,
      sq_pass = mo$moran_sq_pass,
      stringsAsFactors = FALSE
    )
  )
}))
rownames(final_mean_diag) <- NULL

final_vol_diag <- do.call(rbind, lapply(names(vol_diag_all), function(wm) {
  lj <- vol_diag_all[[wm]]$ljung
  mo <- vol_diag_all[[wm]]$moran
  
  rbind(
    data.frame(
      stage = "volatility",
      weight = wm,
      test = "LjungBox",
      variable = lj$variable,
      resid_pass = lj$lb_resid_pass,
      sq_pass = lj$lb_sq_pass,
      stringsAsFactors = FALSE
    ),
    data.frame(
      stage = "volatility",
      weight = wm,
      test = "Moran",
      variable = mo$variable,
      resid_pass = mo$moran_resid_pass,
      sq_pass = mo$moran_sq_pass,
      stringsAsFactors = FALSE
    )
  )
}))
rownames(final_vol_diag) <- NULL

build_vol_param_table <- function(fit, weight_name) {
  p <- fit$p
  
  labels <- c(
    paste0("A_tilde_", 1:p),
    paste0("Psi_", rep(1:p, each = p), rep(1:p, times = p)),
    paste0("Pi_", rep(1:p, each = p), rep(1:p, times = p))
  )
  
  estimates <- c(
    fit$A_tilde_est[1, ],
    vec(fit$Psi_est),
    vec(fit$Pi_est)
  )
  
  data.frame(
    weight = weight_name,
    parameter = labels,
    estimate = estimates,
    stringsAsFactors = FALSE
  )
}

param_table_vol <- do.call(rbind, lapply(names(fit_vol_all), function(wm) {
  build_vol_param_table(fit_vol_all[[wm]], wm)
}))
rownames(param_table_vol) <- NULL

## =========================================================
## 11. PRINT KEY OUTPUTS
## =========================================================

cat("\n================ FINAL EVALUATION ================\n")
print(final_eval)

cat("\n================ MEAN DIAGNOSTICS ================\n")
print(final_mean_diag)

cat("\n============= VOLATILITY DIAGNOSTICS =============\n")
print(final_vol_diag)

cat("\n============= VOLATILITY PARAMETERS ==============\n")
print(param_table_vol)

## =========================================================
## 12. SAVE RESULTS
## =========================================================

loglik_vec_spARCH_p <- function(pars, Y, W, sig_u2, eps = 1e-6) {
  dimY <- dim(Y)
  n <- dimY[1]
  p <- dimY[2]
  Tt <- dimY[3]
  
  Z <- safe_logsq(Y, eps = eps)
  
  A_tilde <- matrix(rep(pars[1:p], n), n, p, byrow = TRUE)
  Psi     <- matrix(pars[(p + 1):(p^2 + p)], p, p)
  Pi      <- matrix(pars[(p^2 + p + 1):(2 * p^2 + p)], p, p)
  
  S <- diag(n * p) - t(Psi) %x% W
  
  detS <- tryCatch(determinant(S, logarithm = TRUE), error = function(e) NULL)
  if (is.null(detS)) return(NA_real_)
  log_det_S <- as.numeric(detS$modulus)
  if (!is.finite(log_det_S)) return(NA_real_)
  
  sum_eps_2 <- 0
  for (tt in 2:Tt) {
    vec_u_t <- S %*% vec(Z[, , tt]) -
      vec(A_tilde) -
      (diag(p) %x% Z[, , tt - 1]) %*% vec(Pi)
    
    if (any(!is.finite(vec_u_t))) return(NA_real_)
    sum_eps_2 <- sum_eps_2 + sum(vec_u_t^2)
  }
  
  LL <- -(Tt - 1) * n * p / 2 * log(2 * pi) -
    (Tt - 1) * n * p / 2 * log(sig_u2) +
    (Tt - 1) * log_det_S -
    (1 / (2 * sig_u2)) * sum_eps_2
  
  LL
}


## =========================================================
## 12. VOLATILITY MODEL STANDARD ERRORS
## =========================================================

compute_se_vec_spARCH <- function(fit, Y, W, eps = 1e-6) {
  
  pars <- c(
    fit$A_tilde_est[1, ],
    as.vector(fit$Psi_est),
    as.vector(fit$Pi_est)
  )
  
  sig_u2 <- fit$sig_u2
  
  cat("Computing Hessian...\n")
  
  H <- numDeriv::hessian(
    func = function(p) loglik_vec_spARCH_p(
      p, Y = Y, W = W, sig_u2 = sig_u2, eps = eps
    ),
    x = pars
  )
  
  Hsym <- 0.5 * (H + t(H))
  
  Hinv <- tryCatch(
    solve(-Hsym),
    error = function(e) {
      cat("Hessian inversion failed, adding ridge...\n")
      solve(-Hsym + diag(1e-6, nrow(Hsym)))
    }
  )
  
  se <- sqrt(pmax(diag(Hinv), 0))
  tval <- pars / se
  pval <- 2 * (1 - pnorm(abs(tval)))
  
  list(
    pars = pars,
    se = se,
    tval = tval,
    pval = pval,
    vcov = Hinv,
    H = H,
    Hsym = Hsym
  )
}

se_results <- list()

for (wm in names(fit_vol_all)) {
  
  cat("\n============================\n")
  cat("SE computation for:", wm, "\n")
  
  fit_w <- fit_vol_all[[wm]]
  
  U_train_full <- fit_mean_all[[wm]]$residuals
  valid_train <- (max(fit_mean_all[[wm]]$tau) + 1):dim(U_train_full)[3]
  U_train <- U_train_full[, , valid_train, drop = FALSE]
  
  se_results[[wm]] <- compute_se_vec_spARCH(
    fit = fit_w,
    Y   = U_train,
    W   = weight_matrices[[wm]],
    eps = eps_log
  )
}

build_param_labels_vol <- function(p) {
  labels <- c()
  
  labels <- c(labels, paste0("A_tilde_", 1:p))
  
  for (i in 1:p) {
    for (j in 1:p) {
      labels <- c(labels, paste0("Psi_", i, j))
    }
  }
  
  for (i in 1:p) {
    for (j in 1:p) {
      labels <- c(labels, paste0("Pi_", i, j))
    }
  }
  
  labels
}

param_table_se <- do.call(rbind, lapply(names(se_results), function(wm) {
  
  res <- se_results[[wm]]
  p <- fit_vol_all[[wm]]$p
  
  data.frame(
    weight = wm,
    parameter = build_param_labels_vol(p),
    estimate = res$pars,
    se = res$se,
    t_value = res$tval,
    p_value = res$pval,
    signif = cut(
      res$pval,
      c(-Inf, 0.01, 0.05, 0.1, Inf),
      labels = c("***", "**", "*", "")
    ),
    stringsAsFactors = FALSE
  )
}))

print(param_table_se)

## =========================================================
## 13. MEAN MODEL PARAMETER TABLE
## =========================================================

build_param_labels_mean <- function(fit) {
  
  p <- fit$p
  k <- fit$k
  L <- length(fit$tau)
  
  labels <- c()
  
  for (i in 1:p) {
    for (j in 1:p) {
      labels <- c(labels, paste0("Psi_", i, j))
    }
  }
  
  for (l in 1:L) {
    for (i in 1:p) {
      for (j in 1:p) {
        labels <- c(labels, paste0("Pi(", fit$tau[l], ")_", i, j))
      }
    }
  }
  
  for (x in 1:k) {
    for (j in 1:p) {
      labels <- c(labels, paste0("Beta_", x, "_", j))
    }
  }
  
  labels <- c(labels, "sigma_u")
  
  labels
}

param_table_mean <- do.call(rbind, lapply(names(fit_mean_all), function(wm) {
  
  fit <- fit_mean_all[[wm]]
  pars <- fit$pars_est
  se <- sqrt(pmax(diag(fit$Hinv), 0))
  
  tval <- pars / se
  pval <- 2 * (1 - pnorm(abs(tval)))
  
  data.frame(
    weight = wm,
    parameter = build_param_labels_mean(fit),
    estimate = pars,
    se = se,
    t_value = tval,
    p_value = pval,
    signif = cut(
      pval,
      c(-Inf, 0.01, 0.05, 0.1, Inf),
      labels = c("***", "**", "*", "")
    ),
    stringsAsFactors = FALSE
  )
}))

print(param_table_mean)

## =========================================================
## 14. COMBINED PARAMETER TABLE
## =========================================================

param_table_mean$stage <- "Mean"
param_table_se$stage   <- "Volatility"

param_table_all <- rbind(
  param_table_mean,
  param_table_se
)

print(param_table_all)

## =========================================================
## 15. SAVE RESULTS
## =========================================================

save(
  Y_all, Y_train, Y_test,
  fit_mean_all, mean_diag_all, mean_test_resid_all,
  fit_vol_all, vol_diag_all, fcst_all,
  eval_results, final_eval, final_mean_diag, final_vol_diag,
  param_table_vol, param_table_mean, param_table_se, param_table_all,
  se_results,
  station_ids, dates_all, train_idx, test_idx,
  file = "Cfull_pipeline_multivariate_STmean_vecspARCH.RData"
)