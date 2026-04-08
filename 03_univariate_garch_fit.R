# ============================================================
# Fit station-wise GARCH and EGARCH models
# ============================================================

library(rugarch)
library(data.table)

set.seed(123)

# -----------------------------
# SETTINGS
# -----------------------------
height <- "ws100"
input_dir  <- "data/processed"
output_dir <- file.path("results/univariate_garch", height)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# LOAD DATA
# -----------------------------
load(file.path(input_dir, paste0("wind_STL_AR1_", height, ".RData")))

wind_wide <- get(paste0("res_", height, "_wide"))
wind_wide <- wind_wide[complete.cases(wind_wide)]

dates <- wind_wide$Date
Y <- as.matrix(wind_wide[, -"Date"])

cutoff <- as.Date("2020-12-31")
Y_train <- Y[dates <= cutoff, , drop = FALSE]

station_ids <- colnames(Y_train)

# -----------------------------
# MODEL SPECS
# -----------------------------
spec_garch <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model     = list(armaOrder = c(0,0), include.mean = FALSE)
)

spec_egarch <- ugarchspec(
  variance.model = list(model = "eGARCH", garchOrder = c(1,1)),
  mean.model     = list(armaOrder = c(0,0), include.mean = FALSE)
)

# -----------------------------
# HELPERS
# -----------------------------
get_aic_bic <- function(fit) {
  ic <- infocriteria(fit)
  if (is.matrix(ic)) {
    list(aic = ic["Akaike", 1], bic = ic["Bayes", 1])
  } else {
    list(aic = ic["Akaike"], bic = ic["Bayes"])
  }
}

fit_one <- function(x, spec, name) {
  x <- x[is.finite(x)]
  if (length(x) < 100) return(NULL)
  
  fit <- try(ugarchfit(spec, x, solver = "hybrid"), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  
  pars <- coef(fit)
  ic   <- get_aic_bic(fit)
  z    <- residuals(fit, standardize = TRUE)
  
  list(
    model = name,
    alpha1 = pars["alpha1"],
    beta1  = pars["beta1"],
    gamma1 = if ("gamma1" %in% names(pars)) pars["gamma1"] else NA,
    aic = ic$aic,
    bic = ic$bic,
    lb_resid_10 = Box.test(z, lag=10, type="Ljung")$p.value,
    lb_sq_10    = Box.test(z^2, lag=10, type="Ljung")$p.value
  )
}

# -----------------------------
# RUN
# -----------------------------
results <- rbindlist(lapply(seq_along(station_ids), function(j) {
  
  sid <- station_ids[j]
  x <- Y_train[, j]
  
  rbindlist(list(
    cbind(IDStations = sid, as.data.table(fit_one(x, spec_garch, "sGARCH"))),
    cbind(IDStations = sid, as.data.table(fit_one(x, spec_egarch, "eGARCH")))
  ), fill = TRUE)
}))

# -----------------------------
# SAVE
# -----------------------------
saveRDS(results, file.path(output_dir, "garch_results.rds"))
fwrite(results, file.path(output_dir, "garch_results.csv"))

cat("GARCH estimation done.\n")