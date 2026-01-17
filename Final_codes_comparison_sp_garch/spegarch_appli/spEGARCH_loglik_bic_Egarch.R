# === Libraries ===
library(Matrix)
library(Rsolnp)
library(nleqslv)
library(dplyr)

# === Load Functions ===
source("functions_spegarch.R")  # Must contain: estimate_spEGARCH, g, likelihood, etc.

# === Load Data ===
residuals <- readRDS("residuals_VAR1.rds")
train.l <- nrow(residuals) - 252
resid_train <- residuals[1:train.l, ]

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
  Wg     = readRDS("Wmat_granger.rds")
)

# === Compute log-likelihood and BIC ===
compute_loglik_bic <- function(Y, W) {
  # Fit the model
  fit <- estimate_spEGARCH(Y = t(Y), W_1 = W, W_2 = W, n_init = 20)
  
  # Recompute log-likelihood using final parameters
  ll <- -likelihood(fit$params, list(Y = t(Y), W_1 = W, W_2 = W))  # negative because function returns -loglik
  
  # Number of parameters: 7
  p <- length(fit$params)
  
  # Number of observations (T × N)
  n_obs <- length(Y)
  
  bic <- -2 * ll + p * log(n_obs)
  
  return(data.frame(LogLik = ll, BIC = bic))
}

# === Run for each matrix ===
results_list <- lapply(names(weight_matrices), function(wname) {
  cat("Computing for:", wname, "\n")
  W <- weight_matrices[[wname]]
  result <- compute_loglik_bic(resid_train, W)
  cbind(Weight = wname, result)
})

# === Combine and Display ===
loglik_bic_table <- bind_rows(results_list)
print(loglik_bic_table)

# === Save Output ===

save.image("spEGARCH_loglik_bic_Egarch.Rdata")

