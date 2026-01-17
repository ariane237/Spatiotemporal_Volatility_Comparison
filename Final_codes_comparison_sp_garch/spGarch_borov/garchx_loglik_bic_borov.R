# === Libraries ===
library(rugarch)
library(dplyr)

# === Load Data ===
residuals <- readRDS("residuals_VAR1.rds")
train.l <- nrow(residuals) - 252
train_data <- residuals[1:train.l, ]

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

# === Function to compute log-likelihood and BIC ===
compute_loglik_bic_garchx <- function(Y, W) {
  T <- nrow(Y)
  N <- ncol(Y)
  total_loglik <- 0
  failed <- 0

  for (i in 1:N) {
    # Spatial lag of squared returns for asset i
    X_ext <- sapply(2:T, function(t) {
      sum(W[i, ] * Y[t - 1, ]^2)
    })

    y_i <- Y[2:T, i]
    x_i <- matrix(X_ext, ncol = 1)

    spec <- ugarchspec(
      variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
      mean.model = list(armaOrder = c(0, 0), include.mean = FALSE,
                        external.regressors = x_i),
      distribution.model = "norm"
    )

    fit <- tryCatch(
      ugarchfit(spec = spec, data = y_i, solver = "hybrid", fit.control = list(stationarity = 1)),
      error = function(e) NULL
    )

    if (!is.null(fit)) {
      total_loglik <- total_loglik + likelihood(fit)
    } else {
      failed <- failed + 1
    }
  }

  n_obs <- (T - 1) * N
  p_per_asset <- 3  # omega, alpha, beta + 1 external regressor = 4?
  p_total <- (N - failed) * p_per_asset

  BIC <- -2 * total_loglik + p_total * log(n_obs)
  return(data.frame(LogLik = total_loglik, BIC = BIC, FailedAssets = failed))
}

# === Run for each matrix ===
loglik_bic_results <- lapply(names(weight_matrices), function(wname) {
  message("Estimating for ", wname)
  W <- weight_matrices[[wname]]
  res <- compute_loglik_bic_garchx(train_data, W)
  cbind(Weight = wname, res)
})

# === Combine results
loglik_bic_table <- bind_rows(loglik_bic_results)
print(loglik_bic_table)

# ===  Save
save.image("garchx_loglik_bic_borov.Rdata")
