
# ==============================================================================
# DYNAMIC LOG-ARCH: DIMENSIONALITY EXPERIMENT
# n = 4, 5, ..., 16
# Network: Granger-filtered correlation
# Estimation: original GMM-2SLS function
# Forecasts: rolling one-step-ahead, without re-estimation
# ==============================================================================

suppressPackageStartupMessages({
  library(MASS)
  library(proxy)
  library(lmtest)
  library(tsdistributions)
})

#source("functions_dysparch.R")

# ==============================================================================
# 1. DATA AND ASSET ORDERING
# ==============================================================================

residuals <- as.matrix(readRDS("residuals_VAR1.rds"))
storage.mode(residuals) <- "double"

out.l <- 252L
train.l <- nrow(residuals) - out.l
dimensions <- 4:16

stopifnot(
  ncol(residuals) >= 16L,
  train.l > 20L,
  all(is.finite(residuals))
)

# IMPORTANT: Use exactly the same ordering as DCC.

asset_order <- readRDS("asset_order_dimensionality.rds")

stopifnot(
  length(asset_order) == 16L,
  all(!is.na(asset_order)),
  all(asset_order == as.integer(asset_order)),
  all(asset_order %in% seq_len(ncol(residuals))),
  !anyDuplicated(asset_order)
)

# ==============================================================================
# 2. LOG-SQUARED RESIDUALS
# ==============================================================================

# Preserve the original zero-handling convention:
# replace zero observations by the smallest positive
# squared observation for that asset.

compute_log_residuals <- function(Y) {
  
  Y <- as.matrix(Y)
  result <- matrix(
    NA_real_,
    nrow = nrow(Y),
    ncol = ncol(Y)
  )
  
  for (i in seq_len(ncol(Y))) {
    
    squared <- Y[, i]^2
    positive <- squared[squared > 0]
    
    if (length(positive) == 0L) {
      stop("An asset has only zero residuals.")
    }
    
    squared[squared == 0] <- min(positive)
    
    result[, i] <- log(squared)
  }
  
  colnames(result) <- colnames(Y)
  result
}

# ==============================================================================
# 3. GRANGER-FILTERED CORRELATION MATRIX
# ==============================================================================

construct_W_GF <- function(Y_train,
                           max_lag = 12L,
                           p_threshold = 0.05) {
  
  Y_train <- as.matrix(Y_train)
  n <- ncol(Y_train)
  
  # Original correlation-based distance definition.
  # Inverse distances are used as edge weights.
  
  # Correlation-based distance
  R <- cor(Y_train, method = "pearson")
  
  D <- sqrt(2 * pmax(1 - R, 0))
  diag(D) <- Inf
  
  if (any(!is.finite(D[!diag(n)])) ||
      any(D[!diag(n)] <= 0)) {
    stop("Invalid correlation distances.")
  }
  
  C <- 1 / D
  diag(C) <- 0
  
  # Granger filtering:
  # G[i,j] = 1 if series j Granger-causes series i.
  
  G <- matrix(0, nrow = n, ncol = n)
  
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      
      if (i == j) next
      
      test <- tryCatch(
        lmtest::grangertest(
          Y_train[, i] ~ Y_train[, j],
          order = max_lag
        ),
        error = function(e) NULL
      )
      
      if (is.null(test)) next
      
      pval <- test$`Pr(>F)`[2]
      
      if (is.finite(pval) && pval < p_threshold) {
        G[i, j] <- 1
      }
    }
  }
  
  # Filter and row-standardise.
  # Zero-degree rows remain zero.
  
  W <- C * G
  
  row_totals <- rowSums(W)
  nonzero <- row_totals > 0
  
  W[nonzero, ] <- sweep(
    W[nonzero, , drop = FALSE],
    1,
    row_totals[nonzero],
    "/"
  )
  
  W[!nonzero, ] <- 0
  diag(W) <- 0
  
  stopifnot(all(is.finite(W)))
  
  list(
    W = W,
    G = G,
    density = sum(G) / (n * (n - 1)),
    isolated = sum(rowSums(G) == 0)
  )
}

# ==============================================================================
# 4. FORECAST METRICS AND PROXIES
# ==============================================================================

rmsfe_fn <- function(E) {
  mean(apply(as.matrix(E), 2, function(x) sqrt(mean(x^2))))
}

mafe_fn <- function(E) {
  mean(apply(as.matrix(E), 2, function(x) mean(abs(x))))
}

compute_proxies <- function(Y, train.l, out.l, lambda = 0.74) {
  
  n <- ncol(Y)
  
  Y_test <- Y[
    (train.l + 1L):(train.l + out.l),
    ,
    drop = FALSE
  ]
  
  RV5_abs <- matrix(NA_real_, out.l, n)
  RV5_sq <- matrix(NA_real_, out.l, n)
  
  for (j in seq_len(out.l)) {
    
    t <- train.l + j
    
    window <- Y[(t - 4L):t, , drop = FALSE]
    
    RV5_abs[j, ] <- colMeans(abs(window))
    RV5_sq[j, ] <- sqrt(colMeans(window^2))
  }
  
  # EWMA variance, initialised from the final training date.
  
  EWMA <- matrix(NA_real_, out.l, n)
  previous <- Y[train.l, ]^2
  
  for (j in seq_len(out.l)) {
    
    previous <- lambda * previous +
      (1 - lambda) * Y_test[j, ]^2
    
    EWMA[j, ] <- previous
  }
  
  list(
    RV = log(pmax(Y_test^2, 1e-12)),
    RV5_abs = log(pmax(RV5_abs^2, 1e-12)),
    RV5_sq = log(pmax(RV5_sq^2, 1e-12)),
    EWMA = log(pmax(EWMA, 1e-12))
  )
}

# ==============================================================================
# 5. ESTIMATE AND FORECAST DSTARCH
# ==============================================================================

run_DSTARCH <- function(Y, out.l = 252L) {
  
  Y <- as.matrix(Y)
  
  n <- ncol(Y)
  train.l <- nrow(Y) - out.l
  
  # ---------------------------------------------------------------------------
  # 5.1 Construct network from training observations only
  # ---------------------------------------------------------------------------
  
  start_network <- proc.time()[["elapsed"]]
  
  network <- construct_W_GF(
    Y_train = Y[1:train.l, , drop = FALSE],
    max_lag = 12L,
    p_threshold = 0.05
  )
  
  network_time <- proc.time()[["elapsed"]] - start_network
  
  Wmat <- network$W
  
  # ---------------------------------------------------------------------------
  # 5.2 Prepare training data
  # ---------------------------------------------------------------------------
  
  log_Y <- compute_log_residuals(Y)
  
  y_train <- t(log_Y[1:train.l, , drop = FALSE])
  
  # ---------------------------------------------------------------------------
  # 5.3 Original GMM-2SLS estimator
  # ---------------------------------------------------------------------------
  
  start_est <- proc.time()[["elapsed"]]
  
  model <- GMM_SDPD_2SLS_ARCH_ind_timelags(
    Y = y_train,
    X = NULL,
    W = Wmat,
    info = list(
      ksy = 10,
      ksx = 10,
      stl = 0,
      tl = 1,
      ted = 0
    )
  )
  
  estimation_time <- proc.time()[["elapsed"]] - start_est
  
  theta <- as.numeric(model$theta)
  
  if (length(theta) != n + 1L ||
      any(!is.finite(theta))) {
    stop("Invalid DSTARCH coefficient vector.")
  }
  
  rho_hat <- theta[1]
  gamma_hat <- theta[2:(n + 1L)]
  
  # ---------------------------------------------------------------------------
  # 5.4 Correctly aligned individual intercepts
  # ---------------------------------------------------------------------------
  
  Y_current <- y_train[
    ,
    2:train.l,
    drop = FALSE
  ]
  
  Y_lagged <- y_train[
    ,
    1:(train.l - 1L),
    drop = FALSE
  ]
  
  mu_hat <- rowMeans(
    Y_current -
      rho_hat * (Wmat %*% Y_current) -
      sweep(Y_lagged, 1, gamma_hat, "*")
  )
  
  # ---------------------------------------------------------------------------
  # 5.5 Spatial multiplier
  # ---------------------------------------------------------------------------
  
  A <- diag(n) - rho_hat * Wmat
  
  reciprocal_condition <- rcond(A)
  
  if (!is.finite(reciprocal_condition) ||
      reciprocal_condition < 1e-10) {
    stop("Spatial multiplier is singular or ill-conditioned.")
  }
  
  # ---------------------------------------------------------------------------
  # 5.6 Rolling one-step-ahead forecasts
  # ---------------------------------------------------------------------------
  
  start_fc <- proc.time()[["elapsed"]]
  
  forecasts <- matrix(
    NA_real_,
    nrow = out.l,
    ncol = n
  )
  
  for (j in seq_len(out.l)) {
    
    # At forecast origin t-1, only the preceding
    # log-squared residual vector is used.
    
    lagged_y <- log_Y[train.l + j - 1L, ]
    
    rhs <- mu_hat + gamma_hat * lagged_y
    
    forecasts[j, ] <- as.numeric(
      solve(A, rhs)
    )
  }
  
  forecasting_time <- proc.time()[["elapsed"]] - start_fc
  
  if (any(!is.finite(forecasts))) {
    stop("Non-finite DSTARCH forecasts.")
  }
  
  # ---------------------------------------------------------------------------
  # 5.7 Forecast evaluation
  #
  # DSTARCH forecasts are already on the log scale.
  # ---------------------------------------------------------------------------
  
  proxies <- compute_proxies(Y, train.l, out.l)
  
  errors <- lapply(proxies, function(proxy) {
    proxy - forecasts
  })
  
  metrics <- lapply(errors, function(E) {
    c(
      RMSFE = rmsfe_fn(E),
      MAFE = mafe_fn(E)
    )
  })
  
  # Original log-squared residual target as well.
  
  error_logresid <- log_Y[
    (train.l + 1L):(train.l + out.l),
    ,
    drop = FALSE
  ] - forecasts
  
  metrics$LogResidual <- c(
    RMSFE = rmsfe_fn(error_logresid),
    MAFE = mafe_fn(error_logresid)
  )
  
  list(
    n = n,
    model = model,
    W = Wmat,
    G = network$G,
    theta = theta,
    mu = mu_hat,
    forecasts = forecasts,
    metrics = metrics,
    parameters = length(theta) + n,
    network_density = network$density,
    isolated_nodes = network$isolated,
    reciprocal_condition = reciprocal_condition,
    network_time = network_time,
    estimation_time = estimation_time,
    forecasting_time = forecasting_time,
    total_time = network_time +
      estimation_time + forecasting_time
  )
}

# ==============================================================================
# 6. PILOT TEST AT n = 4
# ==============================================================================

Y_4 <- residuals[
  ,
  asset_order[1:4],
  drop = FALSE
]

test_DSTARCH <- run_DSTARCH(Y_4, out.l = out.l)

cat("\nPilot test at n = 4:\n")
print(test_DSTARCH$theta)
print(test_DSTARCH$metrics)
print(dim(test_DSTARCH$forecasts))
print(test_DSTARCH$network_density)

stopifnot(
  all(dim(test_DSTARCH$forecasts) == c(out.l, 4L))
)

# ==============================================================================
# 7. FULL DIMENSIONALITY EXPERIMENT
# ==============================================================================

results_DSTARCH <- vector("list", length(dimensions))
names(results_DSTARCH) <- as.character(dimensions)

summary_rows <- list()

for (n in dimensions) {
  
  cat("\n============================================\n")
  cat("DSTARCH | Number of assets:", n, "\n")
  cat("============================================\n")
  
  if (n == 4L) {
    
    result <- test_DSTARCH
    
  } else {
    
    Y_n <- residuals[
      ,
      asset_order[seq_len(n)],
      drop = FALSE
    ]
    
    result <- tryCatch(
      run_DSTARCH(Y_n, out.l = out.l),
      error = function(e) {
        message(
          "DSTARCH failed at n = ",
          n, ": ", conditionMessage(e)
        )
        NULL
      }
    )
  }
  
  if (is.null(result)) next
  
  results_DSTARCH[[as.character(n)]] <- result
  
  m <- result$metrics
  
  summary_rows[[as.character(n)]] <- data.frame(
    Model = "DSTARCH",
    n = n,
    
    Parameters = result$parameters,
    NetworkDensity = result$network_density,
    IsolatedNodes = result$isolated_nodes,
    ReciprocalCondition = result$reciprocal_condition,
    
    NetworkTime = result$network_time,
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
    MAFE_EWMA = unname(m$EWMA["MAFE"]),
    
    RMSFE_LogResidual = unname(m$LogResidual["RMSFE"]),
    MAFE_LogResidual = unname(m$LogResidual["MAFE"])
  )
  
  print(summary_rows[[as.character(n)]])
  
  saveRDS(
    results_DSTARCH,
    "DSTARCH_dimensionality_results.rds"
  )
  
  saveRDS(
    do.call(rbind, summary_rows),
    "DSTARCH_dimensionality_summary.rds"
  )
}

# ==============================================================================
# 8. SAVE FINAL RESULTS
# ==============================================================================

summary_DSTARCH <- if (length(summary_rows) > 0L) {
  do.call(rbind, summary_rows)
} else {
  data.frame()
}

rownames(summary_DSTARCH) <- NULL

print(summary_DSTARCH, digits = 6)

write.csv(
  summary_DSTARCH,
  "DSTARCH_dimensionality_summary.csv",
  row.names = FALSE
)

saveRDS(
  results_DSTARCH,
  "DSTARCH_dimensionality_results.rds"
)

cat("\nDSTARCH dimensionality experiment completed.\n")
