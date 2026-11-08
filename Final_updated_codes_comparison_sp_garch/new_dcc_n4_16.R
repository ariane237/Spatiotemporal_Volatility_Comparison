# ==============================================================================
# DIMENSIONALITY ANALYSIS: DCC-GARCH
# n = 4, 5, ..., 16
# Fixed parameters, recursive one-step-ahead forecasts
# ==============================================================================
setwd("C:/Users/Ariane/OneDrive/Bureau/Analysis_geodata")


# ==============================================================================
# DCC-GARCH: DIMENSIONALITY ANALYSIS
# Dimensions: n = 4, 5, ..., 16
# Out-of-sample period: 252 observations
# One-step-ahead rolling forecasts without parameter re-estimation
# ==============================================================================

suppressPackageStartupMessages({
  library(rugarch)
  library(rmgarch)
})

# ==============================================================================
# 1. DATA AND EXPERIMENTAL DESIGN
# ==============================================================================

residuals <- as.matrix(readRDS("residuals_VAR1.rds"))
storage.mode(residuals) <- "double"

out.l <- 252L
train.l <- nrow(residuals) - out.l
dimensions <- 4:16

stopifnot(
  is.numeric(residuals),
  ncol(residuals) >= 16L,
  train.l >= 10L,
  all(is.finite(residuals))
)

# Use the same asset ordering for every model.
# If the ordering has already been saved, reuse it.

order_file <- "asset_order_dimensionality.rds"

if (file.exists(order_file)) {
  
  asset_order <- readRDS(order_file)
  
} else {
  
  set.seed(123)
  
  asset_order <- sample(
    seq_len(ncol(residuals)),
    size = 16L,
    replace = FALSE
  )
  
  saveRDS(asset_order, order_file)
}

stopifnot(
  length(asset_order) == 16L,
  all(!is.na(asset_order)),
  all(asset_order == as.integer(asset_order)),
  all(asset_order %in% seq_len(ncol(residuals))),
  !anyDuplicated(asset_order)
)

# Keep the selected assets in a fixed order.
# The first four assets form the core subset.

cat("Selected asset indices:\n")
print(asset_order)


# ==============================================================================
# 2. FORECAST ERROR METRICS
# ==============================================================================

# Mean of asset-specific RMSFEs.
# This preserves the aggregation used in the original analysis.

rmsfe_fn <- function(error_mat) {
  
  error_mat <- as.matrix(error_mat)
  
  mean(
    apply(
      error_mat,
      2,
      function(x) sqrt(mean(x^2))
    )
  )
}

mafe_fn <- function(error_mat) {
  
  error_mat <- as.matrix(error_mat)
  
  mean(
    apply(
      error_mat,
      2,
      function(x) mean(abs(x))
    )
  )
}


# ==============================================================================
# 3. VOLATILITY PROXIES
# ==============================================================================

compute_proxies <- function(Y, train.l, out.l, lambda = 0.94) {
  
  Y <- as.matrix(Y)
  
  n <- ncol(Y)
  
  idx_test <- seq.int(train.l + 1L, train.l + out.l)
  
  Y_test <- Y[idx_test, , drop = FALSE]
  
  # ---------------------------------------------------------------------------
  # Squared residuals
  # ---------------------------------------------------------------------------
  
  RV_proxy <- Y_test^2
  
  # ---------------------------------------------------------------------------
  # Five-observation rolling proxies
  # ---------------------------------------------------------------------------
  
  RV5_abs <- matrix(
    NA_real_,
    nrow = out.l,
    ncol = n
  )
  
  RV5_sq <- matrix(
    NA_real_,
    nrow = out.l,
    ncol = n
  )
  
  for (j in seq_len(out.l)) {
    
    t <- train.l + j
    
    window <- Y[
      (t - 4L):t,
      ,
      drop = FALSE
    ]
    
    RV5_abs[j, ] <- colMeans(abs(window))
    
    RV5_sq[j, ] <- sqrt(colMeans(window^2))
  }
  
  # ---------------------------------------------------------------------------
  # EWMA conditional variance proxy
  #
  # Initialize from the last training observation.
  # Update using observations through the current date.
  # ---------------------------------------------------------------------------
  
  EWMA_vol <- matrix(
    NA_real_,
    nrow = out.l,
    ncol = n
  )
  
  ewma_previous <- Y[train.l, ]^2
  
  for (j in seq_len(out.l)) {
    
    ewma_previous <-
      lambda * ewma_previous +
      (1 - lambda) * Y_test[j, ]^2
    
    EWMA_vol[j, ] <- ewma_previous
  }
  
  # ---------------------------------------------------------------------------
  # Log transformations
  # ---------------------------------------------------------------------------
  
  list(
    RV = log(pmax(RV_proxy, 1e-12)),
    RV5_abs = log(pmax(RV5_abs^2, 1e-12)),
    RV5_sq = log(pmax(RV5_sq^2, 1e-12)),
    EWMA = log(pmax(EWMA_vol, 1e-12))
  )
}


# ==============================================================================
# 4. DCC-GARCH ESTIMATION AND FORECASTING
# ==============================================================================

run_DCC <- function(Y, out.l = 252L) {
  
  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"
  
  n <- ncol(Y)
  T_obs <- nrow(Y)
  
  train.l <- T_obs - out.l
  
  stopifnot(
    n >= 2L,
    train.l >= 10L,
    all(is.finite(Y))
  )
  
  # ---------------------------------------------------------------------------
  # 4.1 Univariate GARCH(1,1) specification
  # ---------------------------------------------------------------------------
  
  spec_garch <- ugarchspec(
    variance.model = list(
      model = "sGARCH",
      garchOrder = c(1, 1)
    ),
    mean.model = list(
      armaOrder = c(0, 0),
      include.mean = FALSE
    ),
    distribution.model = "norm"
  )
  
  # ---------------------------------------------------------------------------
  # 4.2 DCC(1,1) specification
  # ---------------------------------------------------------------------------
  
  dcc_spec <- dccspec(
    uspec = multispec(
      replicate(n, spec_garch, simplify = FALSE)
    ),
    dccOrder = c(1, 1),
    distribution = "mvnorm"
  )
  
  # ---------------------------------------------------------------------------
  # 4.3 Estimation
  #
  # out.sample = 252 reserves the final 252 observations
  # for rolling out-of-sample forecasting.
  # ---------------------------------------------------------------------------
  
  start_est <- proc.time()[["elapsed"]]
  
  fit <- dccfit(
    spec = dcc_spec,
    data = Y,
    out.sample = out.l,
    solver = "solnp"
  )
  
  estimation_time <- proc.time()[["elapsed"]] - start_est
  
  # ---------------------------------------------------------------------------
  # 4.4 Convergence diagnostics
  # ---------------------------------------------------------------------------
  
  convergence_code <- tryCatch(
    as.integer(fit@mfit$convergence),
    error = function(e) NA_integer_
  )
  
  # ---------------------------------------------------------------------------
  # 4.5 Log-likelihood and BIC
  #
  # Confirmed from the user's installed rmgarch version:
  #
  # 1. Akaike
  # 2. Bayes
  # 3. Shibata
  # 4. Hannan-Quinn
  # ---------------------------------------------------------------------------
  
  loglik <- as.numeric(likelihood(fit))
  
  IC <- infocriteria(fit)
  
  if (length(IC) < 2L) {
    stop("Information criteria could not be extracted.")
  }
  
  BIC <- as.numeric(IC[2])
  
  # ---------------------------------------------------------------------------
  # 4.6 One-step-ahead rolling forecasts
  #
  # The model parameters remain fixed.
  # Forecast origins advance as new observations become available.
  # ---------------------------------------------------------------------------
  
  start_fc <- proc.time()[["elapsed"]]
  
  fc <- dccforecast(
    fit,
    n.ahead = 1L,
    n.roll = out.l - 1L
  )
  
  forecasting_time <- proc.time()[["elapsed"]] - start_fc
  
  # ---------------------------------------------------------------------------
  # 4.7 Extract marginal variance forecasts
  # ---------------------------------------------------------------------------
  
  sigma_array <- sigma(fc)
  
  expected_dim <- c(1L, n, out.l)
  
  actual_dim <- dim(sigma_array)
  
  if (
    length(actual_dim) != 3L ||
    !all(actual_dim == expected_dim)
  ) {
    stop(
      "Unexpected forecast dimensions. Expected: ",
      paste(expected_dim, collapse = " x "),
      "; obtained: ",
      paste(actual_dim, collapse = " x ")
    )
  }
  
  # sigma_array has dimensions:
  # [forecast horizon, asset, rolling forecast origin]
  #
  # Convert standard deviations to conditional variances.
  
  forecasts <- t(
    matrix(
      sigma_array[1, , ],
      nrow = n,
      ncol = out.l
    )
  )^2
  
  if (
    !all(dim(forecasts) == c(out.l, n)) ||
    any(!is.finite(forecasts)) ||
    any(forecasts <= 0)
  ) {
    stop("Invalid conditional variance forecasts.")
  }
  
  colnames(forecasts) <- colnames(Y)
  
  # ---------------------------------------------------------------------------
  # 4.8 Compute forecast errors
  # ---------------------------------------------------------------------------
  
  proxies <- compute_proxies(
    Y = Y,
    train.l = train.l,
    out.l = out.l
  )
  
  log_forecasts <- log(pmax(forecasts, 1e-12))
  
  errors <- lapply(proxies, function(proxy) {
    proxy - log_forecasts
  })
  
  # ---------------------------------------------------------------------------
  # 4.9 Forecast performance
  # ---------------------------------------------------------------------------
  
  metrics <- lapply(errors, function(error_mat) {
    
    c(
      RMSFE = rmsfe_fn(error_mat),
      MAFE = mafe_fn(error_mat)
    )
  })
  
  # ---------------------------------------------------------------------------
  # 4.10 Number of estimated parameters
  # ---------------------------------------------------------------------------
  
  n_parameters <- length(coef(fit))
  
  # ---------------------------------------------------------------------------
  # 4.11 Return results
  # ---------------------------------------------------------------------------
  
  list(
    n = n,
    fit = fit,
    forecasts = forecasts,
    metrics = metrics,
    loglik = loglik,
    BIC = BIC,
    n_parameters = n_parameters,
    convergence = convergence_code,
    estimation_time = estimation_time,
    forecasting_time = forecasting_time,
    total_time = estimation_time + forecasting_time
  )
}


# ==============================================================================
# 5. TEST THE FUNCTION AT n = 4
# ==============================================================================

cat("\nTesting DCC implementation at n = 4...\n")

Y_4 <- residuals[
  ,
  asset_order[1:4],
  drop = FALSE
]

test_DCC <- run_DCC(
  Y = Y_4,
  out.l = out.l
)

cat("\nBIC:\n")
print(test_DCC$BIC)

cat("\nConvergence code:\n")
print(test_DCC$convergence)

cat("\nForecast dimensions:\n")
print(dim(test_DCC$forecasts))

cat("\nForecast metrics:\n")
print(test_DCC$metrics)

stopifnot(
  all(dim(test_DCC$forecasts) == c(out.l, 4L)),
  all(is.finite(test_DCC$forecasts))
)


# ==============================================================================
# 6. DIMENSIONALITY EXPERIMENT
# ==============================================================================

results_DCC <- vector(
  "list",
  length(dimensions)
)

names(results_DCC) <- as.character(dimensions)

summary_rows <- list()

for (n in dimensions) {
  
  cat("\n==================================================\n")
  cat("DCC-GARCH | Dimension n =", n, "\n")
  cat("==================================================\n")
  
  selected_assets <- asset_order[seq_len(n)]
  
  Y_n <- residuals[
    ,
    selected_assets,
    drop = FALSE
  ]
  
  # Reuse the validated n = 4 result.
  
  if (n == 4L) {
    
    result <- test_DCC
    
  } else {
    
    result <- tryCatch(
      run_DCC(
        Y = Y_n,
        out.l = out.l
      ),
      error = function(e) {
        message(
          "DCC failed at n = ",
          n,
          ": ",
          conditionMessage(e)
        )
        NULL
      }
    )
  }
  
  if (is.null(result)) {
    next
  }
  
  results_DCC[[as.character(n)]] <- result
  
  m <- result$metrics
  
  summary_rows[[as.character(n)]] <- data.frame(
    Model = "DCC",
    n = n,
    
    LogLikelihood = result$loglik,
    BIC = result$BIC,
    
    Parameters = result$n_parameters,
    Convergence = result$convergence,
    
    EstimationTime = result$estimation_time,
    ForecastTime = result$forecasting_time,
    TotalTime = result$total_time,
    
    RMSFE_RV = unname(m$RV["RMSFE"]),
    MAFE_RV = unname(m$RV["MAFE"]),
    
    RMSFE_RV5_abs = unname(m$RV5_abs["RMSFE"]),
    MAFE_RV5_abs = unname(m$RV5_abs["MAFE"]),
    
    RMSFE_RV5_sq = unname(m$RV5_sq["RMSFE"]),
    MAFE_RV5_sq = unname(m$RV5_sq["MAFE"]),
    
    RMSFE_EWMA = unname(m$EWMA["RMSFE"]),
    MAFE_EWMA = unname(m$EWMA["MAFE"])
  )
  
  cat("\nResults for n =", n, "\n")
  print(summary_rows[[as.character(n)]])
  
  # Save after each completed dimension.
  
  saveRDS(
    results_DCC,
    "DCC_dimensionality_results.rds"
  )
  
  saveRDS(
    do.call(rbind, summary_rows),
    "DCC_dimensionality_summary.rds"
  )
}


# ==============================================================================
# 7. FINAL RESULTS
# ==============================================================================

summary_DCC <- if (length(summary_rows) > 0L) {
  do.call(rbind, summary_rows)
} else {
  data.frame()
}

rownames(summary_DCC) <- NULL

cat("\n==================================================\n")
cat("FINAL DCC DIMENSIONALITY RESULTS\n")
cat("==================================================\n")

print(summary_DCC, digits = 6)

write.csv(
  summary_DCC,
  "DCC_dimensionality_summary.csv",
  row.names = FALSE
)

saveRDS(
  results_DCC,
  "DCC_dimensionality_results.rds"
)

cat("\nDCC dimensionality experiment completed.\n")
