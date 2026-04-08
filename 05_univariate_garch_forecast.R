library(rugarch)
library(data.table)

set.seed(123)

height <- "ws100"
input_dir <- "data/processed"

load(file.path(input_dir, paste0("wind_STL_AR1_", height, ".RData")))

wind_wide <- get(paste0("res_", height, "_wide"))
wind_wide <- wind_wide[complete.cases(wind_wide)]

dates <- wind_wide$Date
Y <- as.matrix(wind_wide[, -"Date"])

cutoff <- as.Date("2020-12-31")

Y_train <- Y[dates <= cutoff, ]
Y_test  <- Y[dates > cutoff, ]

spec_garch <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model     = list(armaOrder = c(0,0), include.mean = FALSE)
)

forecast_one <- function(train, test, spec) {
  fit <- ugarchfit(spec, train)
  setfixed(spec) <- as.list(coef(fit))
  
  sapply(1:length(test), function(i) {
    hist <- c(train, test[1:(i-1)])
    tail(sigma(ugarchfilter(spec, hist)), 1)^2
  })
}

H <- apply(Y_train, 2, function(x)
  forecast_one(x, Y_test[,1], spec_garch)
)