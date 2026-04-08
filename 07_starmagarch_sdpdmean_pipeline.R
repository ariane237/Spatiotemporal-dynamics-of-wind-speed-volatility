# ============================================================
# Script: 07_starmagarch_sdpdmean_pipeline.R
# Purpose:
#   Fit STARMAGARCH using SDPD-mean residuals,
#   evaluate diagnostics and recursive forecasts.
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
# 1. SETTINGS
# ------------------------------------------------------------

height <- "ws10"
meanW_choice <- "kNN"

input_dir   <- "results/sdpd_mean"
weights_dir <- "results/spatial_weights/objects"
output_dir  <- file.path("results/starmagarch_sdpdmean", height, meanW_choice)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

lb_lag <- 20
moran_alpha <- 0.05

# ------------------------------------------------------------
# 2. Load SDPD residuals
# ------------------------------------------------------------

master_path <- file.path(
  input_dir,
  height,
  paste0(height, "_SDPDNL_heteroAR_MASTER.rds")
)

if (!file.exists(master_path)) stop("Missing SDPD MASTER file.")

master <- readRDS(master_path)

if (!meanW_choice %in% names(master$fits)) {
  stop("Invalid meanW_choice.")
}

E_NT <- master$fits[[meanW_choice]]$residuals_NT_aligned

D_full <- t(E_NT)  # T x N
dates_vec <- as.Date(master$dates_residuals)
stations  <- master$stations

T_total <- nrow(D_full)
N <- ncol(D_full)

# ------------------------------------------------------------
# 3. Train / test split
# ------------------------------------------------------------

train_idx <- dates_vec >= "2016-01-01" & dates_vec <= "2020-12-31"
test_idx  <- dates_vec >= "2021-01-01" & dates_vec <= "2021-12-31"

D_train <- D_full[train_idx, ]
D_test  <- D_full[test_idx, ]

Y_train <- t(D_train)
Y_test  <- t(D_test)

train.l <- ncol(Y_train)
out.l   <- ncol(Y_test)

residuals_all <- rbind(D_train, D_test)

# ------------------------------------------------------------
# 4. Load weight matrices
# ------------------------------------------------------------

loadW <- function(name) {
  readRDS(file.path(weights_dir, paste0(name, "_", height, ".rds")))$W
}

weight_matrices <- list(
  distance    = loadW("W_distance_band_55km"),
  knn         = loadW("W_knn_k5"),
  directional = loadW("W_directional_advection")
)

# ------------------------------------------------------------
# 5. STARMAGARCH fit function
# ------------------------------------------------------------

fit_oneW <- function(Y, W, label) {
  
  W_dense <- if (inherits(W, "Matrix")) as.matrix(W) else W
  W_arr <- array(W_dense, c(nrow(W_dense), ncol(W_dense), 1))
  
  init_par <- list(
    mu = mean(Y),
    phi = matrix(0.7, 1),
    theta = matrix(0.01, 1),
    omega = 1,
    alpha = matrix(0.01, 1),
    beta  = matrix(0.01, 1)
  )
  
  init_vec <- pmax(apply(Y, 1, var), 1e-6)
  
  map <- parameterlist2maptemplate(init_par)
  fobj <- CreateLikelihood(Y, W_arr, init = init_vec, parameters = init_par, map = map)
  fit <- fitSTARMAGARCH(fobj, Y, print = FALSE)
  
  sigma_tr <- sigma(fit)
  eps_hat <- Y / pmax(sigma_tr, 1e-12)
  
  # Ljung-Box
  lb_res <- apply(eps_hat, 1, function(e)
    Box.test(e, lag = lb_lag, type = "Ljung-Box")$p.value)
  
  lb_sq <- apply(eps_hat^2, 1, function(e)
    Box.test(e, lag = lb_lag, type = "Ljung-Box")$p.value)
  
  # Moran
  lw <- mat2listw(W_dense, style = "W", zero.policy = TRUE)
  
  mo_res <- apply(eps_hat, 2, function(x)
    tryCatch(moran.test(x, lw)$p.value, error = function(e) NA))
  
  mo_sq <- apply(eps_hat^2, 2, function(x)
    tryCatch(moran.test(x, lw)$p.value, error = function(e) NA))
  
  list(
    label = label,
    fit = fit,
    ljung = list(
      res = mean(lb_res > 0.05, na.rm = TRUE),
      sq  = mean(lb_sq > 0.05, na.rm = TRUE)
    ),
    moran = list(
      res = mean(mo_res > 0.05, na.rm = TRUE),
      sq  = mean(mo_sq > 0.05, na.rm = TRUE)
    )
  )
}

# ------------------------------------------------------------
# 6. Parallel estimation
# ------------------------------------------------------------

cl <- makeCluster(min(3, detectCores()-1))
registerDoParallel(cl)

fits_list <- foreach(w = names(weight_matrices),
                     .packages = c("TMB","spdep","Matrix")) %dopar% {
                       fit_oneW(Y_train, weight_matrices[[w]], w)
                     }

stopCluster(cl)
names(fits_list) <- names(weight_matrices)

saveRDS(fits_list, file.path(output_dir, "fits.rds"))

# ------------------------------------------------------------
# 7. Forecasting
# ------------------------------------------------------------

sigma_next <- function(fit, Y_hist) {
  n <- nrow(Y_hist)
  Y_ext <- cbind(Y_hist, rep(0, n))
  as.numeric(sigma(fit, newdata = Y_ext)[, ncol(Y_ext)])
}

eval_metric <- function(h, proxy) {
  err <- log(pmax(h,1e-12)) - log(pmax(proxy,1e-12))
  c(RMSFE = sqrt(mean(err^2)), MAFE = mean(abs(err)))
}

results <- list()

for (w in names(weight_matrices)) {
  
  fit <- fits_list[[w]]$fit
  H <- matrix(NA, out.l, N)
  
  for (t in 1:out.l) {
    Y_hist <- t(residuals_all[1:(train.l + t - 1), ])
    H[t, ] <- sigma_next(fit, Y_hist)^2
  }
  
  RV <- residuals_all[(train.l+1):(train.l+out.l), ]^2
  
  results[[w]] <- eval_metric(H, RV)
}

summary_table <- bind_rows(lapply(names(results), function(w) {
  data.frame(Weight = w, t(results[[w]]))
}))

fwrite(summary_table, file.path(output_dir, "forecast_summary.csv"))

# ------------------------------------------------------------
# 8. Diagnostics summary
# ------------------------------------------------------------

diag_table <- bind_rows(lapply(names(fits_list), function(w) {
  data.frame(
    Weight = w,
    LB_res = fits_list[[w]]$ljung$res,
    LB_sq  = fits_list[[w]]$ljung$sq,
    MO_res = fits_list[[w]]$moran$res,
    MO_sq  = fits_list[[w]]$moran$sq
  )
}))

fwrite(diag_table, file.path(output_dir, "diagnostics.csv"))

cat("STARMAGARCH SDPD pipeline completed.\n")