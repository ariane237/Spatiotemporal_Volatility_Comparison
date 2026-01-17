library(foreach)
library(doParallel)
library(parallel)
library(dplyr)
library(MASS)


source("functions_dysparch.R")


# === Load Data ===
residuals <- readRDS("residuals_VAR1.rds")
train.l <- nrow(residuals) - 252
out.l   <- 252

# get log sqaured residuals
log_residuals <- residuals
for (i in 1:ncol(residuals)) {
  log_residuals[,i] <- ifelse(residuals[,i]==0, log(min(residuals[residuals[,i]!=0,i]^2)), log(residuals[,i]^2)) 
}


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

# Modified to include training data for rolling computation
RV5_abs <- t(sapply(1:out.l, function(i) {
  colMeans(abs(residuals[(train.l + i - 4):(train.l + i), ]))
}))

RV5_sq <- t(sapply(1:out.l, function(i) {
  sqrt(colMeans(residuals[(train.l + i - 4):(train.l + i), ]^2))
}))


# EWMA
ewma_filter <- function(x, lambda = 0.74) {
  result <- numeric(length(x))
  result[1] <- x[1]^2
  for (t in 2:length(x)) {
    result[t] <- lambda * result[t - 1] + (1 - lambda) * x[t]^2
  }
  result
}

test_resids <- residuals[(train.l + 1):(train.l + out.l), ]
EWMA_vol <- apply(test_resids, 2, ewma_filter)

# === Log-transform Proxies (once only) ===
log_RV        <- log(pmax(RV_proxy, 1e-12))
log_RV5_abs   <- log(pmax(RV5_abs^2, 1e-12))
log_RV5_sq    <- log(pmax(RV5_sq^2, 1e-12))
log_EWMA      <- log(pmax(EWMA_vol^2, 1e-12))

# === Define evaluation functions ===
rmsfe_fn <- function(error_mat) mean(apply(error_mat, 2, function(x) sqrt(mean(x^2, na.rm = TRUE))))
mafe_fn  <- function(error_mat) mean(apply(error_mat, 2, function(x) mean(abs(x), na.rm = TRUE)))


# setting up the parallel computing

n.cores <- min(124, detectCores())
cl <- makeCluster(n.cores, type = "FORK")
registerDoParallel(cl)

# Store results for different weight matrices
results_list <- list()
n <- ncol(log_residuals)

for (w_name in names(weight_matrices)) {
  Wmat <- weight_matrices[[w_name]]
  cat("\n>> Weight matrix:", w_name, "\n")
  
  # === Estimate model only ONCE on training data ===
  y_train <- t(log_residuals[1:train.l, ])
  model <- GMM_SDPD_2SLS_ARCH_ind_timelags(Y = y_train, X = NULL, W = Wmat,
                                           info = list(ksy = 10, ksx = 10, stl = 0, tl = 1, ted = 0))
  
  mu_0_hat <- apply(
    y_train[, -ncol(y_train)] - model$theta[1,1] * Wmat %*% y_train[, -ncol(y_train)] -
      array(model$theta[2:(n+1), 1], dim = dim(y_train[, -1])) * y_train[, -1],
    1, mean
  )
  
  # === Forecast all steps without re-estimation ===
  OOSf_netw <- foreach(j = 1:out.l, .combine = 'rbind', .packages = c("MASS")) %dopar% {
    lagged_resid <- log_residuals[train.l + j - 1, ]
    forecast <- t(solve(diag(n) - model$theta[1,1] * Wmat) %*%
                    t(model$theta[2:(n+1), 1] * lagged_resid + mu_0_hat))
    forecast
  }
  
  # === Compute forecast errors
  fERR_RV        <- log_RV - OOSf_netw
  fERR_RV5_abs   <- log_RV5_abs - OOSf_netw
  fERR_RV5_sq    <- log_RV5_sq - OOSf_netw
  fERR_EWMA      <- log_EWMA - OOSf_netw
  fERR_logresid  <- log_residuals[(train.l + 1):nrow(log_residuals), ] - OOSf_netw
  
  # === Compute metrics
  metrics_list <- list(
    RMSFE = rmsfe_fn(fERR_logresid),
    MAFE =  mafe_fn(fERR_logresid),
    
    RMSFE_RV = rmsfe_fn(fERR_RV),
    MAFE_RV  = mafe_fn(fERR_RV),
    
    RMSFE_RV5_abs = rmsfe_fn(fERR_RV5_abs),
    MAFE_RV5_abs  = mafe_fn(fERR_RV5_abs),
    
    RMSFE_RV5_sq = rmsfe_fn(fERR_RV5_sq),
    MAFE_RV5_sq  = mafe_fn(fERR_RV5_sq),
    
    RMSFE_EWMA = rmsfe_fn(fERR_EWMA),
    MAFE_EWMA  = mafe_fn(fERR_EWMA)
  )
  
  results_list[[w_name]] <- c(
    list("OOSf_netw" = OOSf_netw, "fERR_netw" = fERR_logresid),
    metrics_list
  )
}

stopCluster(cl)






proxy_summary <- lapply(names(results_list), function(w) {
  res <- results_list[[w]]
  data.frame(
    Weight = w,
    RMSFE = res$RMSFE,
    MAFE = res$MAFE,
    RMSFE_RV = res$RMSFE_RV,
    MAFE_RV = res$MAFE_RV,
    RMSFE_RV5_abs = res$RMSFE_RV5_abs,
    MAFE_RV5_abs = res$MAFE_RV5_abs,
    RMSFE_RV5_sq = res$RMSFE_RV5_sq,
    MAFE_RV5_sq = res$MAFE_RV5_sq,
    RMSFE_EWMA = res$RMSFE_EWMA,
    MAFE_EWMA = res$MAFE_EWMA
  )
}) %>% bind_rows()

print(proxy_summary)





save.image("dysparch_no_reestimation.RData")