# === Libraries ===
library(TMB)
library(dplyr)

# === Load Data ===
residuals <- readRDS("residuals_VAR1.rds")
train.l <- nrow(residuals) - 252
resid_train <- residuals[1:train.l, ]

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

# === Load and Compile TMB Model ===
model_name <- "STARMAGARCH"
dyn.load(dynlib(model_name))

# === Initialize container for results ===
loglik_bic_list <- list()

for (wname in names(weight_matrices)) {
  cat("Estimating for:", wname, "\n")
  
  W_raw <- weight_matrices[[wname]]
  W_arr <- if (is.matrix(W_raw)) array(W_raw, c(ncol(resid_train), ncol(resid_train), 1)) else W_raw
  Y_in <- t(resid_train)  # shape: assets × time
  init_vec <- pmax(apply(Y_in, 1, var), 1e-6)

  init_par <- list(
    mu    = mean(Y_in),
    phi   = matrix(0.7, ncol = 1),
    theta = matrix(0.01, ncol = 1),
    omega = 1,
    alpha = matrix(0.01, ncol = 1),
    beta  = matrix(0.01, ncol = 1)
  )
  
  map <- parameterlist2maptemplate(init_par)
  
  # Error-handled estimation
  result <- tryCatch({
    fobj <- CreateLikelihood(Y_in, W_arr, init = init_vec, parameters = init_par, map = map)
    fit  <- fitSTARMAGARCH(fobj, Y_in, print = FALSE)

    # Use fallback if loglik not returned
    loglik <- if (!is.null(fit$loglik)) fit$loglik else -fit$optimization$objective

    # If loglik is invalid, skip
    if (is.nan(loglik) || is.infinite(loglik)) stop("LogLik is invalid")

    p_total <- length(fit$coefficients)
    n_obs   <- nrow(Y_in) * ncol(Y_in)  # total data points

    bic <- -2 * loglik + p_total * log(n_obs)

    data.frame(Weight = wname, LogLik = loglik, BIC = bic)
  }, error = function(e) {
    cat("⚠️ Failed for:", wname, "| Reason:", conditionMessage(e), "\n")
    NULL
  })

  # Store if successful
  if (!is.null(result)) {
    loglik_bic_list[[wname]] <- result
  }
}

# === Combine and save results ===
loglik_bic_combined <- bind_rows(loglik_bic_list)
print(loglik_bic_combined)


# Save results

save.image("loglik_bic_spMAGARCH_trainset.RData")


