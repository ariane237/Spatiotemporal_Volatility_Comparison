# ============================================================
# HOLLELAND STARMAGARCH: DIMENSIONALITY EXPERIMENT
# n = 4,...,16; 252 one-step-ahead forecasts
# ============================================================

library(TMB)

source("functions_spgarch_hol.R")

if (!("STARMAGARCH" %in% names(getLoadedDLLs()))) {
  dyn.load(TMB::dynlib("STARMAGARCH"))
}

residuals <- as.matrix(readRDS("residuals_VAR1.rds"))
asset_order <- readRDS("asset_order_dimensionality.rds")

out.l <- 252L
train.l <- nrow(residuals) - out.l
dimensions <- 4:16

stopifnot(
  train.l > 1L,
  length(asset_order) >= max(dimensions),
  !anyNA(residuals),
  all(is.finite(residuals))
)

# Same aggregation as DSTARCH and DCC
rmsfe_fn <- function(E) {
  mean(sqrt(colMeans(E^2)))
}

mafe_fn <- function(E) {
  mean(colMeans(abs(E)))
}


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

# ------------------------------------------------------------
# Fixed-parameter recursive forecasting
# ------------------------------------------------------------

forecast_holleland <- function(Y, W, par, init, train.l) {
  
  Y <- as.matrix(Y)
  W <- as.matrix(W)
  
  TT <- nrow(Y)
  n <- ncol(Y)
  out.l <- TT - train.l
  
  stopifnot(
    all(dim(W) == c(n, n)),
    length(init) == n,
    out.l > 0L
  )
  
  required <- c("mu", "phi", "theta",
                "omega", "alpha", "beta")
  
  if (!all(required %in% names(par))) {
    stop("Missing Holleland parameters.")
  }
  
  mu    <- unname(par["mu"])
  phi   <- unname(par["phi"])
  theta <- unname(par["theta"])
  omega <- unname(par["omega"])
  alpha <- unname(par["alpha"])
  beta  <- unname(par["beta"])
  
  # Rows = time; columns = assets
  x <- matrix(NA_real_, TT, n)
  h <- matrix(NA_real_, TT, n)
  
  # Exactly the C++ initialization
  x[1, ] <- 0
  h[1, ] <- init
  
  # Training recursion: observations 2,...,train.l
  for (t in 2:train.l) {
    
    x[t, ] <- (
      Y[t, ] - mu -
        phi * as.vector(W %*% (Y[t - 1, ] - mu)) -
        theta * as.vector(W %*% x[t - 1, ])
    )
    
    h[t, ] <- (
      omega +
        alpha * as.vector(W %*% (x[t - 1, ]^2)) +
        beta * as.vector(W %*% h[t - 1, ])
    )
  }
  
  # One-step-ahead variance forecasts
  forecasts <- matrix(NA_real_, out.l, n)
  
  for (k in seq_len(out.l)) {
    
    t <- train.l + k
    
    # Forecast h_t using information through t-1
    h[t, ] <- (
      omega +
        alpha * as.vector(W %*% (x[t - 1, ]^2)) +
        beta * as.vector(W %*% h[t - 1, ])
    )
    
    forecasts[k, ] <- h[t, ]
    
    # After forecasting, observe y_t and update x_t
    # This is used ONLY for the next forecast
    x[t, ] <- (
      Y[t, ] - mu -
        phi * as.vector(W %*% (Y[t - 1, ] - mu)) -
        theta * as.vector(W %*% x[t - 1, ])
    )
  }
  
  if (any(!is.finite(forecasts)) ||
      any(forecasts <= 0)) {
    stop("Non-positive or non-finite variance forecasts.")
  }
  
  list(
    forecasts = forecasts,
    innovations = x,
    conditional_variances = h
  )
}

# ------------------------------------------------------------
# Fit and evaluate one dimension
# ------------------------------------------------------------

run_Holleland <- function(Y, out.l = 252L) {
  
  Y <- as.matrix(Y)
  
  n <- ncol(Y)
  TT <- nrow(Y)
  train.l <- TT - out.l
  
  # Construct network from training data only
  t0 <- proc.time()[["elapsed"]]
  
  net <- construct_W_GF(
    Y_train = Y[1:train.l, , drop = FALSE]
  )
  
  W <- net$W
  
  network_time <- proc.time()[["elapsed"]] - t0
  
  W_arr <- array(W, dim = c(n, n, 1L))
  
  Y_train <- t(Y[1:train.l, , drop = FALSE])
  
  init_vec <- pmax(
    apply(Y_train, 1, var),
    1e-6
  )
  
  init_par <- list(
    mu = mean(Y_train),
    phi = matrix(0.7, 1, 1),
    theta = matrix(0.01, 1, 1),
    omega = 1,
    alpha = matrix(0.01, 1, 1),
    beta = matrix(0.01, 1, 1)
  )
  
  map <- parameterlist2maptemplate(init_par)
  
  # Estimate once
  t0 <- proc.time()[["elapsed"]]
  
  fobj <- CreateLikelihood(
    data = Y_train,
    W = W_arr,
    init = init_vec,
    parameters = init_par,
    map = map
  )
  
  fit <- fitSTARMAGARCH(
    fobj,
    Y_train,
    print = FALSE
  )
  
  estimation_time <- proc.time()[["elapsed"]] - t0
  
  if (fit$optimization$convergence != 0) {
    warning("Optimizer did not report convergence at n = ", n)
  }
  
  par <- fit$coefficients
  
  # Explicit report at the optimum
  report_opt <- fobj$report(fit$optimization$par)
  
  # Check that our recursion matches the TMB implementation
  train_check <- forecast_holleland(
    Y = Y[1:train.l, , drop = FALSE],
    W = W,
    par = par,
    init = init_vec,
    train.l = train.l - 1L
  )
  
  tmb_h <- t(as.matrix(report_opt$sigma))
  
  recursion_error <- max(abs(
    train_check$conditional_variances - tmb_h
  ))
  
  if (!is.finite(recursion_error) ||
      recursion_error > 1e-6) {
    stop(
      "R recursion does not match TMB at n = ",
      n, ". Maximum error: ", recursion_error
    )
  }
  
  # Genuine rolling OOS forecasts
  t0 <- proc.time()[["elapsed"]]
  
  fc <- forecast_holleland(
    Y = Y,
    W = W,
    par = par,
    init = init_vec,
    train.l = train.l
  )
  
  forecasting_time <- proc.time()[["elapsed"]] - t0
  
  log_fc <- log(pmax(fc$forecasts, 1e-12))
  
  Y_test <- Y[(train.l + 1L):TT, , drop = FALSE]
  
  # Primary proxy: log squared residuals
  log_RV <- log(pmax(Y_test^2, 1e-12))
  
  E <- log_fc - log_RV
  
  metrics <- list(
    RV = c(
      RMSFE = rmsfe_fn(E),
      MAFE = mafe_fn(E)
    )
  )
  
  list(
    n = n,
    parameters = length(par),
    coefficients = par,
    forecasts = fc$forecasts,
    metrics = metrics,
    network_density = net$density,
    isolated_nodes = net$isolated,
    network_time = network_time,
    estimation_time = estimation_time,
    forecasting_time = forecasting_time,
    total_time = network_time +
      estimation_time +
      forecasting_time,
    recursion_error = recursion_error,
    convergence = fit$optimization$convergence
  )
}


Y_4 <- residuals[
  , asset_order[1:4], drop = FALSE
]

test_Holleland <- run_Holleland(
  Y_4,
  out.l = out.l
)

print(test_Holleland$coefficients)
print(test_Holleland$metrics)
print(dim(test_Holleland$forecasts))
print(test_Holleland$recursion_error)
print(test_Holleland$convergence)


# ============================================================
# HOLLELAND DIMENSIONALITY EXPERIMENT
# ============================================================

results_Holleland <- vector("list", length(dimensions))
names(results_Holleland) <- as.character(dimensions)

summary_rows <- list()

for (n in dimensions) {
  
  cat("\nRunning Holleland for n =", n, "\n")
  
  Y_n <- residuals[
    , asset_order[seq_len(n)], drop = FALSE
  ]
  
  result <- if (n == 4L) {
    test_Holleland
  } else {
    tryCatch(
      run_Holleland(Y_n, out.l = out.l),
      error = function(e) {
        message(
          "Failed at n = ", n, ": ",
          conditionMessage(e)
        )
        NULL
      }
    )
  }
  
  if (is.null(result)) next
  
  results_Holleland[[as.character(n)]] <- result
  
  summary_rows[[as.character(n)]] <- data.frame(
    Model = "Holleland",
    n = n,
    Parameters = result$parameters,
    NetworkDensity = result$network_density,
    IsolatedNodes = result$isolated_nodes,
    Convergence = result$convergence,
    RecursionError = result$recursion_error,
    NetworkTime = result$network_time,
    EstimationTime = result$estimation_time,
    ForecastTime = result$forecasting_time,
    TotalTime = result$total_time,
    RMSFE_RV = unname(result$metrics$RV["RMSFE"]),
    MAFE_RV = unname(result$metrics$RV["MAFE"])
  )
  
  # Save intermediate results
  saveRDS(
    results_Holleland,
    "Holleland_dimensionality_results.rds"
  )
  
  saveRDS(
    do.call(rbind, summary_rows),
    "Holleland_dimensionality_summary.rds"
  )
}

summary_Holleland <- do.call(rbind, summary_rows)

write.csv(
  summary_Holleland,
  "Holleland_dimensionality_summary.csv",
  row.names = FALSE
)

print(summary_Holleland)



stopifnot(
  all(summary_Holleland$Convergence == 0),
  all(summary_Holleland$RecursionError < 1e-6),
  all(is.finite(summary_Holleland$RMSFE_RV)),
  all(is.finite(summary_Holleland$MAFE_RV))
)

print(summary_Holleland[
  , c("n", "Parameters", "EstimationTime",
      "ForecastTime", "RMSFE_RV", "MAFE_RV")
])
