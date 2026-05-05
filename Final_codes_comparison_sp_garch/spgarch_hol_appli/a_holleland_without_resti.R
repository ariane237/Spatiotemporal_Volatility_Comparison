# === Libraries ===
library(foreach)
library(doParallel)
library(parallel)
library(dplyr)
library(TMB)

source("functions_spgarch_hol.R")  # Your Holleland functions

# === Load Data ===
residuals <- readRDS("residuals_VAR1.rds")
train.l <- nrow(residuals) - 252  # training length
out.l   <- 252                    # forecasting length

dyn.load(TMB::dynlib("STARMAGARCH"))

# === Load Weight Matrices ===
weight_matrices <- list(
  W1      = readRDS("W1_eucl.rds"),
  W2      = readRDS("W1_corr.rds"),
  W3      = readRDS("W1_picollo.rds"),
  Wg1     = readRDS("W1_eucl_GF.rds"),
  Wg2     = readRDS("W2_corr_GF.rds"),
  Wg3     = readRDS("W3_picollo_GF.rds"),
  W1_5nn  = readRDS("W1_eucl5nn.rds"),
  W2_5nn  = readRDS("W2_corr5nn.rds"),
  W3_5nn  = readRDS("W1_picollo5nn.rds"),
  Wg      = readRDS("Wmat_granger.rds")
)

# ============================================================
# 2. Realised variance proxies
# ============================================================

eps_test <- residuals[(train.l + 1):(train.l + out.l), , drop = FALSE]

RV_1d <- eps_test^2

RV5_sq <- t(sapply(1:out.l, function(i) {
  idx <- max(1, train.l + i - 4):(train.l + i)
  colMeans(residuals[idx, , drop = FALSE]^2, na.rm = TRUE)
}))

RV5_abs <- t(sapply(1:out.l, function(i) {
  idx <- max(1, train.l + i - 4):(train.l + i)
  colMeans(abs(residuals[idx, , drop = FALSE]), na.rm = TRUE)^2
}))

ewma_filter <- function(x, lambda = 0.94) {
  result <- numeric(length(x))
  result[1] <- x[1]^2
  
  for (tt in 2:length(x)) {
    result[tt] <- lambda * result[tt - 1] + (1 - lambda) * x[tt]^2
  }
  
  result
}

EWMA_var <- apply(eps_test, 2, ewma_filter)

log_RV_1d   <- log(pmax(RV_1d, eps_floor))
log_RV5_sq  <- log(pmax(RV5_sq, eps_floor))
log_RV5_abs <- log(pmax(RV5_abs, eps_floor))
log_EWMA    <- log(pmax(EWMA_var, eps_floor))

# ============================================================
# 3. Evaluation metrics
# ============================================================

RMSFE <- function(e) {
  e <- as.matrix(e)
  e[!is.finite(e)] <- NA_real_
  sqrt(mean(e^2, na.rm = TRUE))
}

MAFE <- function(e) {
  e <- as.matrix(e)
  e[!is.finite(e)] <- NA_real_
  mean(abs(e), na.rm = TRUE)
}

QLIKE <- function(h_hat, rv) {
  h_hat <- as.matrix(h_hat)
  rv    <- as.matrix(rv)
  
  h_hat[!is.finite(h_hat)] <- NA_real_
  rv[!is.finite(rv)]       <- NA_real_
  
  h_hat <- pmax(h_hat, eps_floor)
  rv    <- pmax(rv, eps_floor)
  
  loss <- log(h_hat) + rv / h_hat
  loss[!is.finite(loss)] <- NA_real_
  
  mean(loss, na.rm = TRUE)
}

make_qlike_loss <- function(h_hat, rv) {
  h_hat <- as.matrix(h_hat)
  rv    <- as.matrix(rv)
  
  h_hat[!is.finite(h_hat)] <- NA_real_
  rv[!is.finite(rv)]       <- NA_real_
  
  h_hat <- pmax(h_hat, eps_floor)
  rv    <- pmax(rv, eps_floor)
  
  loss <- log(h_hat) + rv / h_hat
  loss[!is.finite(loss)] <- NA_real_
  
  loss
}

# ============================================================
# 4. Forecast helper functions
# ============================================================

coef_to_parameter_list <- function(fit) {
  
  cf <- fit$coefficients
  
  list(
    mu    = unname(cf["mu"]),
    phi   = matrix(unname(cf[grep("^phi", names(cf))]),   ncol = 1),
    theta = matrix(unname(cf[grep("^theta", names(cf))]), ncol = 1),
    omega = unname(cf["omega"]),
    alpha = matrix(unname(cf[grep("^alpha", names(cf))]), ncol = 1),
    beta  = matrix(unname(cf[grep("^beta", names(cf))]),  ncol = 1)
  )
}

sigma_1step_ahead_fixed <- function(fit, Y_hist_NxT, W_dense) {
  
  Y_hist_NxT <- as.matrix(Y_hist_NxT)
  W_dense    <- as.matrix(W_dense)
  
  n <- nrow(Y_hist_NxT)
  
  if (!all(dim(W_dense) == c(n, n))) {
    stop("W dimension does not match Y_hist_NxT.")
  }
  
  Y_ext <- cbind(Y_hist_NxT, rep(0, n))
  Warr  <- array(W_dense, c(n, n, 1))
  
  pars_list <- coef_to_parameter_list(fit)
  map <- parameterlist2maptemplate(pars_list)
  
  init_vec <- pmax(apply(Y_hist_NxT, 1, var, na.rm = TRUE), 1e-6)
  
  fobj_ext <- CreateLikelihood(
    data       = Y_ext,
    W          = Warr,
    init       = init_vec,
    parameters = pars_list,
    map        = map,
    silent     = TRUE
  )
  
  rep_ext <- fobj_ext$report(unname(fit$coefficients))
  
  # report$sigma is variance. Protect against tiny negative numerical values.
  sig_ext <- sqrt(pmax(rep_ext$sigma, eps_floor))
  
  if (!is.matrix(sig_ext)) {
    stop("sigma report is not a matrix.")
  }
  
  out <- as.numeric(sig_ext[, ncol(sig_ext)])
  out[!is.finite(out)] <- NA_real_
  
  out
}

# ============================================================
# 5. Main loop over weight matrices
# ============================================================

summary_list  <- list()
forecast_list <- list()
error_list    <- list()
loss_list     <- list()
fit_list      <- list()
param_list    <- list()

for (wname in names(weight_matrices)) {
  
  cat("\n====================================================\n")
  cat("Weight matrix:", wname, "\n")
  cat("====================================================\n")
  
  W_raw <- as.matrix(weight_matrices[[wname]])
  W_raw[!is.finite(W_raw)] <- 0
  
  if (!all(dim(W_raw) == c(n, n))) {
    warning("Skipping ", wname, ": W has wrong dimension.")
    next
  }
  
  W_arr <- array(W_raw, c(n, n, 1))
  
  Y_train <- t(residuals[1:train.l, , drop = FALSE])
  
  init_vec <- pmax(apply(Y_train, 1, var, na.rm = TRUE), 1e-6)
  
  init_par <- list(
    mu    = mean(Y_train, na.rm = TRUE),
    phi   = matrix(0.7,  ncol = 1),
    theta = matrix(0.01, ncol = 1),
    omega = 1,
    alpha = matrix(0.01, ncol = 1),
    beta  = matrix(0.01, ncol = 1)
  )
  
  map <- parameterlist2maptemplate(init_par)
  
  fobj <- tryCatch(
    CreateLikelihood(
      data       = Y_train,
      W          = W_arr,
      init       = init_vec,
      parameters = init_par,
      map        = map,
      silent     = TRUE
    ),
    error = function(e) {
      message("Likelihood creation failed for ", wname, ": ", e$message)
      NULL
    }
  )
  
  if (is.null(fobj)) next
  
  fit_fixed <- tryCatch(
    fitSTARMAGARCH(
      f     = fobj,
      data  = Y_train,
      print = FALSE
    ),
    error = function(e) {
      message("Estimation failed for ", wname, ": ", e$message)
      NULL
    }
  )
  
  fit_list[[wname]] <- fit_fixed
  
  if (is.null(fit_fixed)) next
  
  cat("Estimated coefficients:\n")
  print(fit_fixed$coefficients)
  
  param_list[[wname]] <- data.frame(
    Weight = wname,
    mu     = unname(fit_fixed$coefficients["mu"]),
    phi    = unname(fit_fixed$coefficients["phi"]),
    theta  = unname(fit_fixed$coefficients["theta"]),
    omega  = unname(fit_fixed$coefficients["omega"]),
    alpha  = unname(fit_fixed$coefficients["alpha"]),
    beta   = unname(fit_fixed$coefficients["beta"])
  )
  
  sigma_fc <- matrix(
    NA_real_,
    nrow = out.l,
    ncol = n
  )
  
  colnames(sigma_fc) <- colnames(residuals)
  
  for (k in 1:out.l) {
    
    if (k %% 25 == 0 || k == 1) {
      cat("Forecast step", k, "of", out.l, "\n")
    }
    
    t_end <- train.l + k - 1
    Y_hist <- t(residuals[1:t_end, , drop = FALSE])
    
    sig_next <- tryCatch(
      sigma_1step_ahead_fixed(
        fit        = fit_fixed,
        Y_hist_NxT = Y_hist,
        W_dense    = W_raw
      ),
      error = function(e) {
        message("Forecast failed for ", wname, " at step ", k, ": ", e$message)
        rep(NA_real_, n)
      }
    )
    
    sigma_fc[k, ] <- pmax(sig_next^2, eps_floor)
  }
  
  log_sig <- log(pmax(sigma_fc, eps_floor))
  log_sig[!is.finite(log_sig)] <- NA_real_
  
  e_1d   <- log_sig - as.matrix(log_RV_1d)
  e_5abs <- log_sig - as.matrix(log_RV5_abs)
  e_5sq  <- log_sig - as.matrix(log_RV5_sq)
  e_EWMA <- log_sig - as.matrix(log_EWMA)
  
  qlike_1d   <- QLIKE(sigma_fc, RV_1d)
  qlike_5abs <- QLIKE(sigma_fc, RV5_abs)
  qlike_5sq  <- QLIKE(sigma_fc, RV5_sq)
  qlike_EWMA <- QLIKE(sigma_fc, EWMA_var)
  
  summary_list[[wname]] <- data.frame(
    Weight = wname,
    
    RMSFE_1d = RMSFE(e_1d),
    MAFE_1d  = MAFE(e_1d),
    QLIKE_1d = qlike_1d,
    
    RMSFE_R5_abs = RMSFE(e_5abs),
    MAFE_R5_abs  = MAFE(e_5abs),
    QLIKE_R5_abs = qlike_5abs,
    
    RMSFE_R5_sq = RMSFE(e_5sq),
    MAFE_R5_sq  = MAFE(e_5sq),
    QLIKE_R5_sq = qlike_5sq,
    
    RMSFE_EWMA = RMSFE(e_EWMA),
    MAFE_EWMA  = MAFE(e_EWMA),
    QLIKE_EWMA = qlike_EWMA,
    
    NA_forecasts = sum(is.na(sigma_fc)),
    Nonfinite_forecasts = sum(!is.finite(sigma_fc))
  )
  
  forecast_list[[wname]] <- sigma_fc
  
  error_list[[wname]] <- list(
    e_1d   = e_1d,
    e_5abs = e_5abs,
    e_5sq  = e_5sq,
    e_EWMA = e_EWMA
  )
  
  loss_list[[wname]] <- list(
    QLIKE_1d   = make_qlike_loss(sigma_fc, RV_1d),
    QLIKE_5abs = make_qlike_loss(sigma_fc, RV5_abs),
    QLIKE_5sq  = make_qlike_loss(sigma_fc, RV5_sq),
    QLIKE_EWMA = make_qlike_loss(sigma_fc, EWMA_var)
  )
}

# ============================================================
# 6. Final output
# ============================================================

summary_table <- bind_rows(summary_list)
param_table   <- bind_rows(param_list)

cat("\n=== FORECAST PERFORMANCE ===\n")
print(summary_table)

cat("\n=== PARAMETER SUMMARY ===\n")
print(param_table)

save(
  summary_table,
  param_table,
  forecast_list,
  error_list,
  loss_list,
  fit_list,
  RV_1d,
  RV5_abs,
  RV5_sq,
  EWMA_var,
  log_RV_1d,
  log_RV5_abs,
  log_RV5_sq,
  log_EWMA,
  file = "spMAGARCH_Holleland_fixed_1step_forecasts_safe.RData"
)

cat("\nDone. Results saved in: spMAGARCH_Holleland_fixed_1step_forecasts_safe.RData\n")
