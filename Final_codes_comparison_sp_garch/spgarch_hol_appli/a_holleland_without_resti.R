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

# Log-transform proxies
log_RV      <- log(pmax(RV_proxy, 1e-12))
log_RV5_abs <- log(pmax(RV5_abs^2, 1e-12))
log_RV5_sq  <- log(pmax(RV5_sq^2, 1e-12))
log_EWMA    <- log(pmax(EWMA_vol^2, 1e-12))

# === Evaluation Metrics ===
RMSFE <- function(e) sqrt(mean(e^2, na.rm = TRUE))
MAFE  <- function(e) mean(abs(e),   na.rm = TRUE)

summary_list <- list()

for (wname in names(weight_matrices)) {
  cat("\n▶ Weight matrix:", wname, "\n")
  W_raw <- weight_matrices[[wname]]
  W_arr <- if (is.matrix(W_raw)) array(W_raw, c(dim(W_raw)[1], dim(W_raw)[1], 1)) else W_raw
  
  sigma_fc <- matrix(NA_real_, nrow = out.l, ncol = ncol(residuals))
  
  # === Estimate once on training data ===
  Y_train <- t(residuals[1:train.l, , drop = FALSE])
  init_vec <- pmax(apply(Y_train, 1, var), 1e-6)
  init_par <- list(
    mu    = mean(Y_train),
    phi   = matrix(c(0.7), ncol = 1),
    theta = matrix(c(0.01), ncol = 1),
    omega = 1,
    alpha = matrix(c(0.01), ncol = 1),
    beta  = matrix(c(0.01), ncol = 1)
  )
  map <- parameterlist2maptemplate(init_par)
  fobj <- CreateLikelihood(Y_train, W_arr, init = init_vec, parameters = init_par, map = map)
  fit_fixed <- fitSTARMAGARCH(fobj, Y_train, print = FALSE)
  
  # === Forecast using fixed parameters ===
  for (k in 1:out.l) {
    start_idx <- (k - 1) + 1
    t_end     <- (k - 1) + train.l
    Y_in      <- t(residuals[start_idx:t_end, , drop = FALSE])
    
    sig <- sigma(fit_fixed, newdata = Y_in)
    if (!is.matrix(sig)) stop("sigma() returned non-matrix")
    sigma_fc[k, ] <- as.numeric(sig[, ncol(sig)])^2
  }
  
  log_sig <- log(pmax(sigma_fc, 1e-12))
  rownames(log_RV) <- NULL
  
  e1  <- log_sig - as.matrix(log_RV)
  e5a <- log_sig - log_RV5_abs
  e5s <- log_sig - log_RV5_sq
  eE  <- log_sig - log_EWMA
  
  summary_list[[wname]] <- data.frame(
    Weight = wname,
    RMSFE_1d   = RMSFE(e1),  MAFE_1d   = MAFE(e1),
    RMSFE_R5_a = RMSFE(e5a), MAFE_R5_a = MAFE(e5a),
    RMSFE_R5_s = RMSFE(e5s), MAFE_R5_s = MAFE(e5s),
    RMSFE_EW   = RMSFE(eE),  MAFE_EW   = MAFE(eE)
  )
}

# === Final output ===
summary_table <- bind_rows(summary_list)
print(summary_table)
save.image("spMAGARCH_Holleland_summary_noreestimate.RData")
