# === Libraries ===
library(foreach)
library(doParallel)
library(parallel)
library(MASS)
library(dplyr)
library(Rsolnp)
library(nleqslv)
library(Matrix)

source("functions_spegarch.R")

# === Load Data ===
residuals <- readRDS("residuals_VAR1.rds")
train.l <- nrow(residuals) - 252
out.l   <- 252

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

# === Realized Volatility Proxies ===
RV_proxy <- residuals[(train.l + 1):(train.l + out.l), ]^2
RV5_abs <- t(sapply(1:out.l, function(i) colMeans(abs(residuals[(train.l + i - 4):(train.l + i), ]))))
RV5_sq  <- t(sapply(1:out.l, function(i) sqrt(colMeans(residuals[(train.l + i - 4):(train.l + i), ]^2))))
EWMA_vol <- apply(residuals[(train.l + 1):(train.l + out.l), ], 2, function(x) {
  res <- numeric(length(x)); res[1] <- x[1]^2
  for (t in 2:length(x)) res[t] <- 0.94 * res[t - 1] + 0.06 * x[t]^2
  res
})

log_RV        <- log(pmax(RV_proxy, 1e-12))
log_RV5_abs   <- log(pmax(RV5_abs^2, 1e-12))
log_RV5_sq    <- log(pmax(RV5_sq^2, 1e-12))
log_EWMA      <- log(pmax(EWMA_vol^2, 1e-12))

# === Parallel Setup ===
n.cores <- min(124, detectCores())
cl <- makeCluster(n.cores, type = "FORK")
registerDoParallel(cl)

results_list <- list()

for (w_name in names(weight_matrices)) {
  message("Processing ", w_name)
  W_1 <- weight_matrices[[w_name]]
  W_2 <- W_1
  
  # Estimate once using training data
  Y_train <- t(residuals[1:train.l, ])
  fit <- estimate_spEGARCH(Y = Y_train, W_1 = W_1, W_2 = W_2, n_init = 10)
  
  # Forecast all 252 steps using fixed parameters
  forecasts <- foreach(j = 1:out.l, .combine = rbind, .packages = c("Matrix")) %dopar% {
    t_index_j <- train.l + j
    if (t_index_j > nrow(residuals)) return(rep(NA, ncol(residuals)))
    
    Y_test <- residuals[t_index_j, ]
    Y_lag  <- residuals[t_index_j - 1, ]
    
    pred <- g(
      y         = as.numeric(Y_test),
      alpha     = fit$params[5],
      rhoW_1    = fit$params[1] * W_1,
      lambdaW_2 = fit$params[7] * W_2,
      theta_0   = fit$params[3],
      theta_1   = fit$params[4],
      xi        = 1,
      rho_1     = fit$params[2],
      lambda_1  = fit$params[6],
      eps_lag   = as.numeric(fit$eps[, ncol(fit$eps)]),
      Y_lag     = as.numeric(Y_lag)
    )
    exp(pred$lnh)
  }
  
  # === Evaluation ===
  eval_metric <- function(h, proxy) {
    log_h <- log(pmax(h, 1e-12))
    err <- log_h - proxy
    c(RMSFE = sqrt(mean(err^2, na.rm = TRUE)), MAFE = mean(abs(err), na.rm = TRUE))
  }
  
  results_list[[w_name]] <- list(
    forecasts = forecasts,
    metrics = list(
      RV        = eval_metric(forecasts, as.matrix(log_RV)),
      RV5_abs   = eval_metric(forecasts, log_RV5_abs),
      RV5_sq    = eval_metric(forecasts, log_RV5_sq),
      EWMA      = eval_metric(forecasts, log_EWMA)
    )
  )
}

stopCluster(cl)

# === Summary Table ===
summary_table <- bind_rows(lapply(names(results_list), function(name) {
  m <- results_list[[name]]$metrics
  data.frame(
    Weight = name,
    RMSFE_RV = m$RV["RMSFE"], MAFE_RV = m$RV["MAFE"],
    RMSFE_RV5_abs = m$RV5_abs["RMSFE"], MAFE_RV5_abs = m$RV5_abs["MAFE"],
    RMSFE_RV5_sq  = m$RV5_sq["RMSFE"],  MAFE_RV5_sq  = m$RV5_sq["MAFE"],
    RMSFE_EWMA    = m$EWMA["RMSFE"],    MAFE_EWMA    = m$EWMA["MAFE"]
  )
}))

print(summary_table)
save.image("spEGARCH_noReestimation_results.RData")
