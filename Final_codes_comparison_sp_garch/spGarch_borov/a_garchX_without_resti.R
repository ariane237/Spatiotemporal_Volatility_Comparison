# === Setup ===
library(foreach)
library(doParallel)
library(parallel)
library(dplyr)
library(MASS)
library(rugarch)

# === Load Data ===
residuals <- readRDS("residuals_VAR1.rds")
train.l <- nrow(residuals) - 252
out.l <- 252

# === Load Weight Matrices ===
weight_matrices <- list(
  W1  = readRDS("W1_eucl.rds"),
  W2  = readRDS("W1_corr.rds"),
  W3  = readRDS("W1_picollo.rds"),
  Wg1 = readRDS("W1_eucl_GF.rds"),
  Wg2 = readRDS("W2_corr_GF.rds"),
  Wg3 = readRDS("W3_picollo_GF.rds"),
  W1_5nn = readRDS("W1_eucl5nn.rds"),
  W2_5nn = readRDS("W2_corr5nn.rds"),
  W3_5nn = readRDS("W1_picollo5nn.rds"),
  Wg = readRDS("Wmat_granger.rds")
)

# === Preprocess Weight Matrices ===
preprocess_weight_matrix <- function(W) {
  W[is.na(W)] <- 0
  W[is.infinite(W)] <- 0
  row_sums <- rowSums(W)
  W[row_sums == 0, ] <- 1e-10
  W / rowSums(W)
}
weight_matrices <- lapply(weight_matrices, preprocess_weight_matrix)

# === Realized Volatility Proxies ===
RV_proxy <- residuals[(train.l + 1):(train.l + out.l), ]^2
RV5_abs <- t(sapply(1:out.l, function(i) colMeans(abs(residuals[(train.l + i - 4):(train.l + i), ]))))
RV5_sq <- t(sapply(1:out.l, function(i) sqrt(colMeans(residuals[(train.l + i - 4):(train.l + i), ]^2))))
EWMA_vol <- apply(residuals[(train.l + 1):(train.l + out.l), ], 2, function(x) {
  lambda <- 0.94
  Reduce(function(acc, val) lambda * acc + (1 - lambda) * val^2, x[-1], init = x[1]^2, accumulate = TRUE)
})

# === Log-transformed proxies
log_RV        <- log(pmax(RV_proxy, 1e-12))
log_RV5_abs   <- log(pmax(RV5_abs^2, 1e-12))
log_RV5_sq    <- log(pmax(RV5_sq^2, 1e-12))
log_EWMA      <- log(pmax(EWMA_vol^2, 1e-12))

# === Error Metrics
rmsfe_fn <- function(error_mat) mean(apply(error_mat, 2, function(x) sqrt(mean(x^2, na.rm = TRUE))))
mafe_fn  <- function(error_mat) mean(apply(error_mat, 2, function(x) mean(abs(x), na.rm = TRUE)))

# === Parallel Setup
n.cores <- min(124, detectCores())
cl <- makeCluster(n.cores, type = "FORK")
registerDoParallel(cl)

results_list <- list()
log_residuals <- log(pmax(residuals^2, 1e-12))

# === Loop over Weight Matrices (NO re-estimation)
for (w_name in names(weight_matrices)) {
  cat("Processing", w_name, "...\n")
  W <- weight_matrices[[w_name]]
  n_assets <- ncol(residuals)

  # === Estimate once per asset
  fitted_models <- vector("list", n_assets)
  for (i in 1:n_assets) {
    X_train <- sapply(2:train.l, function(t) sum(W[i, ] * residuals[t - 1, ]^2))
    spec <- ugarchspec(
      variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
      mean.model = list(armaOrder = c(0, 0), include.mean = FALSE, external.regressors = matrix(X_train, ncol = 1)),
      distribution.model = "norm"
    )
    fitted_models[[i]] <- ugarchfit(spec = spec, data = residuals[2:train.l, i], solver = "hybrid")
  }

  # === Forecast without re-estimation
  OOSf_static <- foreach(j = 1:out.l, .combine = 'rbind', .packages = "rugarch") %dopar% {
    fc_step <- numeric(n_assets)
    for (i in 1:n_assets) {
      # external regressor at forecast time
      X_new <- sum(W[i, ] * residuals[train.l + j - 1, ]^2)
      forecast <- tryCatch({
        ugarchforecast(fitted_models[[i]], n.ahead = 1, external.forecasts = list(mregfor = matrix(X_new, ncol = 1)))
      }, error = function(e) NULL)
      fc_step[i] <- if (!is.null(forecast)) as.numeric(sigma(forecast))^2 else NA
    }
    fc_step
  }

  log_fc <- log(pmax(OOSf_static, 1e-12))

  # === Errors
  fERR_RV <- log_RV - log_fc
  fERR_RV5_abs <- log_RV5_abs - log_fc
  fERR_RV5_sq <- log_RV5_sq - log_fc
  fERR_EWMA <- log_EWMA - log_fc

  metrics <- list(
    RMSFE_RV = rmsfe_fn(fERR_RV), MAFE_RV = mafe_fn(fERR_RV),
    RMSFE_RV5_abs = rmsfe_fn(fERR_RV5_abs), MAFE_RV5_abs = mafe_fn(fERR_RV5_abs),
    RMSFE_RV5_sq = rmsfe_fn(fERR_RV5_sq), MAFE_RV5_sq = mafe_fn(fERR_RV5_sq),
    RMSFE_EWMA = rmsfe_fn(fERR_EWMA), MAFE_EWMA = mafe_fn(fERR_EWMA)
  )

  results_list[[w_name]] <- c(list(OOSf_static = OOSf_static), metrics)
}

stopCluster(cl)

# === Summary Table
proxy_summary <- bind_rows(lapply(names(results_list), function(w) {
  res <- results_list[[w]]
  data.frame(
    Weight = w,
    RMSFE_RV = res$RMSFE_RV, MAFE_RV = res$MAFE_RV,
    RMSFE_RV5_abs = res$RMSFE_RV5_abs, MAFE_RV5_abs = res$MAFE_RV5_abs,
    RMSFE_RV5_sq = res$RMSFE_RV5_sq, MAFE_RV5_sq = res$MAFE_RV5_sq,
    RMSFE_EWMA = res$RMSFE_EWMA, MAFE_EWMA = res$MAFE_EWMA
  )
}))

print(proxy_summary)

save.image("spatial_garchX_no_reestimation.RData")
