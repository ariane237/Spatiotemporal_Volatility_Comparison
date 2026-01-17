# === Libraries ===
library(MASS)
library(dplyr)

source("functions_dysparch.R")  # includes GMM_SDPD_2SLS_ARCH_ind_timelags

# === Load Data ===
residuals <- readRDS("residuals_VAR1.rds")
train.l <- nrow(residuals) - 252
resid_train <- residuals[1:train.l, ]

# === Preprocess Log-Squared Residuals (Handle zeros, NA, Inf) ===
log_resid_train <- log(pmax(resid_train^2, 1e-12))
log_resid_train[!is.finite(log_resid_train)] <- NA  # Remove -Inf and NaN
log_resid_train <- na.omit(log_resid_train)  # Drop rows with NA

# === Load Weight Matrices ===
weight_matrices <- list(
  W1     = readRDS("W1_eucl.rds"),
  W2     = readRDS("W1_corr.rds"),
  W3     = readRDS("W1_picollo.rds"),
  Wg1    = readRDS("W1_eucl_GF.rds"),
  Wg2    = readRDS("W2_corr_GF.rds"),
  Wg3    = readRDS("W3_picollo_GF.rds"),
  W1_5nn = readRDS("W1_eucl5nn.rds"),
  W2_5nn = readRDS("W2_corr5nn.rds"),
  W3_5nn = readRDS("W1_picollo5nn.rds"),
  Wg     = readRDS("Wmat_granger.rds")
)

# === Log-likelihood and BIC computation function ===
compute_loglik_bic_DST <- function(Y_log, W) {
  tryCatch({
    n <- ncol(Y_log)
    T <- nrow(Y_log)

    # Transpose to match expected input: matrix (n x T)
    Y_t <- t(Y_log)

    # Fit the model
    fit <- GMM_SDPD_2SLS_ARCH_ind_timelags(
      Y = Y_t, X = NULL, W = W,
      info = list(ksy = 10, ksx = 10, stl = 0, tl = 1, ted = 0)
    )

    e <- fit$e
    sigma2 <- fit$sigma2

    # Compute log-likelihood
    ll <- sum(dnorm(as.vector(e), mean = 0, sd = sqrt(sigma2), log = TRUE), na.rm = TRUE)

    # BIC
    p_total <- length(fit$theta)
    n_obs <- length(e)
    BIC <- -2 * ll + p_total * log(n_obs)

    return(data.frame(LogLik = ll, BIC = BIC))
  }, error = function(e) {
    cat("⚠️ Error for weight matrix. Reason:", conditionMessage(e), "\n")
    return(data.frame(LogLik = NA, BIC = NA))
  })
}

# === Loop over all matrices ===
results_list <- lapply(names(weight_matrices), function(wname) {
  cat("Processing:", wname, "\n")
  W <- weight_matrices[[wname]]
  res <- compute_loglik_bic_DST(log_resid_train, W)
  cbind(Weight = wname, res)
})

# === Combine and display results
dstarch_loglik_bic_table <- bind_rows(results_list)
print(dstarch_loglik_bic_table)


# === Save
save.image("DSTARCH_loglik_bic_otto.Rdata")
