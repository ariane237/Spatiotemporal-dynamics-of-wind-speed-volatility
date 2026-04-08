# ============================================================
# Script: 06_starmagarch_fit_and_forecast.R
# Purpose:
#   Fit STARMAGARCH models using STL + AR(1) residuals,
#   compute temporal and spatial diagnostics, and evaluate
#   recursive one-step-ahead variance forecasts under
#   alternative spatial weight matrices.
#
# Inputs:
#   - STL + AR(1) residual panels
#   - Spatial weight matrices
#   - STARMAGARCH TMB library and helper functions
#
# Outputs:
#   - fitted STARMAGARCH objects with diagnostics
#   - recursive forecast matrices
#   - forecast accuracy summary tables
#
# Author: Ariane Meli Chrisko
# ============================================================

suppressPackageStartupMessages({
  library(parallel)
  library(doParallel)
  library(foreach)
  library(TMB)
  library(spdep)
  library(Matrix)
  library(data.table)
  library(dplyr)
})

set.seed(123)

# ------------------------------------------------------------
# 1. Paths and settings
# ------------------------------------------------------------

height <- "ws10"

input_dir <- "data/processed"
weights_dir <- "results/spatial_weights/objects"
output_dir <- file.path("results/starmagarch", height)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

source("scripts/utils/functions_spgarch_hol.R")
dyn.load(TMB::dynlib("STARMAGARCH"))

# ------------------------------------------------------------
# 2. Load residual panel and define train/test split
# ------------------------------------------------------------

load(file.path(input_dir, paste0("wind_STL_AR1_", height, ".RData")))

res_mat   <- get(paste0("res_", height, "_mat"))
dates_vec <- as.Date(get(paste0("dates_", height)))

D_full <- as.matrix(res_mat)        # T x N
T_total <- nrow(D_full)
N <- ncol(D_full)

stopifnot(length(dates_vec) == T_total)

train_idx <- dates_vec >= as.Date("2016-01-01") & dates_vec <= as.Date("2020-12-31")
test_idx  <- dates_vec >= as.Date("2021-01-01") & dates_vec <= as.Date("2021-12-31")

D_train <- D_full[train_idx, , drop = FALSE]
D_test  <- D_full[test_idx,  , drop = FALSE]

Y_train <- t(D_train)   # N x T_train
Y_test  <- t(D_test)    # N x T_test

train_length <- ncol(Y_train)
test_length  <- ncol(Y_test)

cat(
  "Total T =", T_total,
  "| Train =", train_length,
  "| Test =", test_length,
  "| N =", N, "\n"
)

# ------------------------------------------------------------
# 3. Load spatial weight matrices
# ------------------------------------------------------------

W_distance <- readRDS(file.path(weights_dir, "W_distance_band_55km.rds"))$W
W_knn      <- readRDS(file.path(weights_dir, "W_knn_k5.rds"))$W
W_directional <- readRDS(file.path(weights_dir, "W_directional_advection.rds"))$W

to_dense <- function(W) {
  if (inherits(W, "Matrix")) as.matrix(W) else W
}

weight_matrices <- list(
  distance    = to_dense(W_distance),
  knn         = to_dense(W_knn),
  directional = to_dense(W_directional)
)

# ------------------------------------------------------------
# 4. Fit STARMAGARCH for one weight matrix
# ------------------------------------------------------------

fit_starmagarch_oneW <- function(Y, W, label, lb_lag = 20, moran_alpha = 0.05) {
  
  W_array <- array(W, c(nrow(W), ncol(W), 1))
  
  n <- nrow(Y)
  Tt <- ncol(Y)
  
  init_vec <- pmax(apply(Y, 1, var), 1e-6)
  
  init_par <- list(
    mu    = mean(Y),
    phi   = matrix(0.7,  ncol = 1),
    theta = matrix(0.01, ncol = 1),
    omega = 1,
    alpha = matrix(0.01, ncol = 1),
    beta  = matrix(0.01, ncol = 1)
  )
  
  map  <- parameterlist2maptemplate(init_par)
  fobj <- CreateLikelihood(Y, W_array, init = init_vec, parameters = init_par, map = map)
  fit  <- fitSTARMAGARCH(fobj, Y, print = FALSE)
  
  sigma_train <- tryCatch(
    sigma(fit, newdata = Y),
    error = function(e) sigma(fit)
  )
  
  if (!is.matrix(sigma_train)) {
    stop("sigma() did not return a matrix for ", label)
  }
  
  eps_hat <- Y / pmax(sigma_train, 1e-12)
  
  # Ljung-Box diagnostics by station
  p_lb_res <- apply(eps_hat, 1, function(e) {
    Box.test(e, lag = lb_lag, type = "Ljung-Box")$p.value
  })
  
  p_lb_res_sq <- apply(eps_hat^2, 1, function(e) {
    Box.test(e, lag = lb_lag, type = "Ljung-Box")$p.value
  })
  
  ljung_df <- data.frame(
    Series = seq_len(n),
    P_Residuals = p_lb_res,
    P_SqResiduals = p_lb_res_sq
  )
  
  # Moran diagnostics by time
  lw <- mat2listw(W, style = "W", zero.policy = TRUE)
  
  p_moran_res <- numeric(Tt)
  p_moran_res_sq <- numeric(Tt)
  
  for (tt in seq_len(Tt)) {
    p_moran_res[tt] <- tryCatch(
      moran.test(eps_hat[, tt], listw = lw, zero.policy = TRUE)$p.value,
      error = function(e) NA_real_
    )
    
    p_moran_res_sq[tt] <- tryCatch(
      moran.test(eps_hat[, tt]^2, listw = lw, zero.policy = TRUE)$p.value,
      error = function(e) NA_real_
    )
  }
  
  moran_df <- data.frame(
    Time = seq_len(Tt),
    P_Residuals = p_moran_res,
    P_SqResiduals = p_moran_res_sq
  )
  
  list(
    label = label,
    fit = fit,
    W = W,
    sigma_train = sigma_train,
    eps_hat = eps_hat,
    ljung = list(
      df = ljung_df,
      pass_rate = c(
        residuals = mean(p_lb_res > 0.05, na.rm = TRUE) * 100,
        squared   = mean(p_lb_res_sq > 0.05, na.rm = TRUE) * 100
      )
    ),
    moran = list(
      df = moran_df,
      pass_rate = c(
        residuals = mean(p_moran_res > moran_alpha, na.rm = TRUE) * 100,
        squared   = mean(p_moran_res_sq > moran_alpha, na.rm = TRUE) * 100
      )
    )
  )
}

# ------------------------------------------------------------
# 5. Estimate STARMAGARCH models in parallel
# ------------------------------------------------------------

n_cores <- min(length(weight_matrices), max(1, detectCores() - 1))
cluster_type <- if (.Platform$OS.type == "windows") "PSOCK" else "FORK"

cl <- makeCluster(n_cores, type = cluster_type)
registerDoParallel(cl)

clusterEvalQ(cl, {
  suppressPackageStartupMessages({
    library(TMB)
    library(spdep)
    library(Matrix)
  })
})

clusterExport(
  cl,
  c("fit_starmagarch_oneW", "Y_train", "weight_matrices"),
  envir = environment()
)

fits_list <- foreach(
  w_name = names(weight_matrices),
  .packages = c("TMB", "spdep", "Matrix")
) %dopar% {
  fit_starmagarch_oneW(
    Y = Y_train,
    W = weight_matrices[[w_name]],
    label = w_name,
    lb_lag = 20,
    moran_alpha = 0.05
  )
}

stopCluster(cl)

names(fits_list) <- names(weight_matrices)

saveRDS(
  fits_list,
  file.path(output_dir, paste0(height, "_starmagarch_fits_with_diagnostics.rds"))
)

# ------------------------------------------------------------
# 6. Extract parameter and information-criteria tables
# ------------------------------------------------------------

extract_param_table <- function(fits_list) {
  bind_rows(lapply(names(fits_list), function(w_name) {
    fit_obj <- fits_list[[w_name]]$fit
    mc <- as.data.frame(fit_obj$matcoef)
    
    data.frame(
      Weight = w_name,
      Param = rownames(mc),
      Estimate = as.numeric(mc[, "Estimates"]),
      StdErr = as.numeric(mc[, "SD"]),
      z_value = as.numeric(mc[, "Zscore"]),
      p_value = as.numeric(mc[, "Pvalue"]),
      row.names = NULL
    )
  }))
}

extract_ic_table <- function(fits_list) {
  bind_rows(lapply(names(fits_list), function(w_name) {
    fit_obj <- fits_list[[w_name]]$fit
    data.frame(
      Weight = w_name,
      AIC = as.numeric(fit_obj$aic),
      BIC = as.numeric(fit_obj$bic),
      row.names = NULL
    )
  }))
}

param_table <- extract_param_table(fits_list)
ic_table <- extract_ic_table(fits_list)

fwrite(param_table, file.path(output_dir, paste0(height, "_starmagarch_params.csv")))
fwrite(ic_table, file.path(output_dir, paste0(height, "_starmagarch_ic.csv")))

# ------------------------------------------------------------
# 7. Build realised variance proxies on the test period
# ------------------------------------------------------------

test_rows <- (train_length + 1):(train_length + test_length)

RV_var <- D_full[test_rows, , drop = FALSE]^2

RV5_ms_var <- t(sapply(seq_len(test_length), function(i) {
  idx_end   <- train_length + i
  idx_start <- idx_end - 4
  colMeans(D_full[idx_start:idx_end, , drop = FALSE]^2)
}))

RV5_abs_var <- t(sapply(seq_len(test_length), function(i) {
  idx_end   <- train_length + i
  idx_start <- idx_end - 4
  mabs <- colMeans(abs(D_full[idx_start:idx_end, , drop = FALSE]))
  mabs^2
}))

ewma_var <- function(x, lambda = 0.94) {
  out <- numeric(length(x))
  out[1] <- x[1]^2
  for (tt in 2:length(x)) {
    out[tt] <- lambda * out[tt - 1] + (1 - lambda) * x[tt]^2
  }
  out
}

EWMA_full <- apply(D_full, 2, ewma_var)
EWMA_var <- EWMA_full[test_rows, , drop = FALSE]

log_RV_var      <- log(pmax(RV_var, 1e-12))
log_RV5_ms_var  <- log(pmax(RV5_ms_var, 1e-12))
log_RV5_abs_var <- log(pmax(RV5_abs_var, 1e-12))
log_EWMA_var    <- log(pmax(EWMA_var, 1e-12))

# ------------------------------------------------------------
# 8. One-step-ahead recursive forecasting
# ------------------------------------------------------------

sigma_1step_ahead <- function(fit, Y_hist_NxT) {
  n <- nrow(Y_hist_NxT)
  Y_ext <- cbind(Y_hist_NxT, rep(0, n))
  
  sigma_ext <- tryCatch(
    sigma(fit, newdata = Y_ext),
    error = function(e) {
      stop("sigma(fit, newdata = ...) failed.")
    }
  )
  
  as.numeric(sigma_ext[, ncol(sigma_ext)])
}

eval_metric <- function(h_var, proxy_log_var) {
  stopifnot(all(dim(h_var) == dim(proxy_log_var)))
  
  log_h <- log(pmax(h_var, 1e-12))
  err <- log_h - proxy_log_var
  
  c(
    RMSFE = sqrt(mean(err^2, na.rm = TRUE)),
    MAFE  = mean(abs(err), na.rm = TRUE)
  )
}

forecast_results <- list()

for (w_name in names(weight_matrices)) {
  cat("Recursive forecasting for", w_name, "\n")
  
  fit_obj <- fits_list[[w_name]]$fit
  forecasts <- matrix(NA_real_, nrow = test_length, ncol = N)
  
  for (j in seq_len(test_length)) {
    t_index <- train_length + j
    Y_hist <- t(D_full[1:(t_index - 1), , drop = FALSE])
    sigma_next <- sigma_1step_ahead(fit_obj, Y_hist)
    forecasts[j, ] <- sigma_next^2
  }
  
  metrics <- list(
    RV_var      = eval_metric(forecasts, log_RV_var),
    RV5_ms_var  = eval_metric(forecasts, log_RV5_ms_var),
    RV5_abs_var = eval_metric(forecasts, log_RV5_abs_var),
    EWMA_var    = eval_metric(forecasts, log_EWMA_var)
  )
  
  forecast_results[[w_name]] <- list(
    forecasts = forecasts,
    metrics = metrics
  )
}

# ------------------------------------------------------------
# 9. Forecast summary table
# ------------------------------------------------------------

forecast_summary <- bind_rows(lapply(names(forecast_results), function(w_name) {
  m <- forecast_results[[w_name]]$metrics
  data.frame(
    Weight = w_name,
    RMSFE_RV      = m$RV_var["RMSFE"],
    MAFE_RV       = m$RV_var["MAFE"],
    RMSFE_RV5_ms  = m$RV5_ms_var["RMSFE"],
    MAFE_RV5_ms   = m$RV5_ms_var["MAFE"],
    RMSFE_RV5_abs = m$RV5_abs_var["RMSFE"],
    MAFE_RV5_abs  = m$RV5_abs_var["MAFE"],
    RMSFE_EWMA    = m$EWMA_var["RMSFE"],
    MAFE_EWMA     = m$EWMA_var["MAFE"]
  )
}))

fwrite(
  forecast_summary,
  file.path(output_dir, paste0(height, "_starmagarch_forecast_summary.csv"))
)

# ------------------------------------------------------------
# 10. Diagnostic summary table
# ------------------------------------------------------------

diagnostic_summary <- bind_rows(lapply(names(fits_list), function(w_name) {
  obj <- fits_list[[w_name]]
  
  data.frame(
    Weight = w_name,
    LB_pass_residuals = obj$ljung$pass_rate["residuals"],
    LB_pass_squared   = obj$ljung$pass_rate["squared"],
    Moran_pass_residuals = obj$moran$pass_rate["residuals"],
    Moran_pass_squared   = obj$moran$pass_rate["squared"]
  )
}))

fwrite(
  diagnostic_summary,
  file.path(output_dir, paste0(height, "_starmagarch_diagnostic_summary.csv"))
)

# ------------------------------------------------------------
# 11. Save full output bundle
# ------------------------------------------------------------

saveRDS(
  list(
    fits_list = fits_list,
    param_table = param_table,
    ic_table = ic_table,
    forecast_results = forecast_results,
    forecast_summary = forecast_summary,
    diagnostic_summary = diagnostic_summary
  ),
  file.path(output_dir, paste0(height, "_starmagarch_full_results.rds"))
)

cat("STARMAGARCH pipeline completed.\n")