# === Libraries ===
library(foreach)
library(doParallel)
library(parallel)
library(dplyr)
library(MASS)
library(rugarch)
library(rmgarch)

# === Load Data ===
residuals <- readRDS("residuals_VAR1.rds")
train.l <- nrow(residuals) - 252
out.l   <- 252

# === Realized Volatility Proxies ===
RV_proxy <- residuals[(train.l + 1):(train.l + out.l), ]^2

RV5_abs <- t(sapply(1:out.l, function(i) {
  colMeans(abs(residuals[(train.l + i - 4):(train.l + i), ]))
}))

RV5_sq <- t(sapply(1:out.l, function(i) {
  sqrt(colMeans(residuals[(train.l + i - 4):(train.l + i), ]^2))
}))

ewma_filter <- function(x, lambda = 0.94) {
  result <- numeric(length(x))
  result[1] <- x[1]^2
  for (t in 2:length(x)) {
    result[t] <- lambda * result[t - 1] + (1 - lambda) * x[t]^2
  }
  result
}

test_resids <- residuals[(train.l + 1):(train.l + out.l), ]
EWMA_vol <- apply(test_resids, 2, ewma_filter)

# === Log-transform Proxies ===
log_RV        <- log(pmax(RV_proxy, 1e-12))
log_RV5_abs   <- log(pmax(RV5_abs^2, 1e-12))
log_RV5_sq    <- log(pmax(RV5_sq^2, 1e-12))
log_EWMA      <- log(pmax(EWMA_vol^2, 1e-12))

# === Error Metrics ===
rmsfe_fn <- function(error_mat) mean(apply(error_mat, 2, function(x) sqrt(mean(x^2, na.rm = TRUE))), na.rm = TRUE)
mafe_fn  <- function(error_mat) mean(apply(error_mat, 2, function(x) mean(abs(x), na.rm = TRUE)), na.rm = TRUE)

# === DCC-GARCH Specification ===
spec_garch <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
  distribution.model = "norm"
)

dcc_spec <- dccspec(
  uspec = multispec(replicate(ncol(residuals), spec_garch)),
  dccOrder = c(1, 1),
  distribution = "mvnorm"
)

# === Fit Once on Training Sample ===
cat("Fitting DCC-GARCH model on training sample...\n")

dcc_fit <- dccfit(
  dcc_spec,
  data = residuals[1:(train.l + out.l), ],
  out.sample = out.l
)

# Log-likelihood and BIC from training period
full_loglik <- likelihood(dcc_fit)
full_bic <- infocriteria(dcc_fit)[1]


# === Forecast All OOS Steps Without Re-estimation ===
cat("Generating forecasts without re-estimation...\n")

dcc_fc <- dccforecast(dcc_fit, n.ahead = 1, n.roll = out.l - 1)

# Extract sigma forecasts correctly

# Forecast with no re-estimation
dcc_fc <- dccforecast(dcc_fit, n.ahead = 1, n.roll = out.l - 1)

# Extract the 3D array: [1, N, out.l]
sigma_array <- sigma(dcc_fc)

# Convert to matrix: [out.l, N]
OOSf_dcc <- t(sigma_array[1, , ])^2  # Square to get variance




# === Forecast Errors ===
fERR_RV        <- log_RV - log(pmax(OOSf_dcc, 1e-12))
fERR_RV5_abs   <- log_RV5_abs - log(pmax(OOSf_dcc, 1e-12))
fERR_RV5_sq    <- log_RV5_sq - log(pmax(OOSf_dcc, 1e-12))
fERR_EWMA      <- log_EWMA - log(pmax(OOSf_dcc, 1e-12))
fERR_logresid  <- log_residuals[(train.l + 1):nrow(log_residuals), ] - log(pmax(OOSf_dcc, 1e-12))

# === Forecast Performance Summary ===
metrics_list <- list(
  LogLikelihood = full_loglik,
  BIC = full_bic,
  RMSFE = rmsfe_fn(fERR_logresid), MAFE = mafe_fn(fERR_logresid),
  RMSFE_RV = rmsfe_fn(fERR_RV), MAFE_RV = mafe_fn(fERR_RV),
  RMSFE_RV5_abs = rmsfe_fn(fERR_RV5_abs), MAFE_RV5_abs = mafe_fn(fERR_RV5_abs),
  RMSFE_RV5_sq = rmsfe_fn(fERR_RV5_sq), MAFE_RV5_sq = mafe_fn(fERR_RV5_sq),
  RMSFE_EWMA = rmsfe_fn(fERR_EWMA), MAFE_EWMA = mafe_fn(fERR_EWMA)
)

cat("\n=== Forecast Performance Metrics ===\n")
print(metrics_list)

# === Save Workspace ===
save.image("dcc_no_reestimation_forecast.RData")

