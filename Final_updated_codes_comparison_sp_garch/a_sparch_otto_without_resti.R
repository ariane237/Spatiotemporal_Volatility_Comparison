# ================================================================
# DSTARCH: fixed-parameter, rolling one-step-ahead forecasts
# ================================================================
# Inputs: functions_dysparch.R, residuals_VAR1.rds, weight matrices
# Forecasts: log conditional variance / log squared residual scale
# No out-of-sample parameter re-estimation.

suppressPackageStartupMessages({
  library(MASS)    # ginv() used by functions_dysparch.R
  library(dplyr)
})
source("functions_dysparch.R")

Y <- as.matrix(readRDS("residuals_VAR1.rds"))
stopifnot(is.numeric(Y), all(is.finite(Y)), nrow(Y) > 260L)
out.l <- 252L
train.l <- nrow(Y) - out.l
n <- ncol(Y)
stopifnot(n >= 2L, train.l >= 10L)

weight_matrices <- list(
  W1 = readRDS("W1_eucl.rds"),
  W2 = readRDS("W1_corr.rds"),
  W3 = readRDS("W1_picollo.rds"),
  Wg1 = readRDS("W1_eucl_GF.rds"),
  Wg2 = readRDS("W2_corr_GF.rds"),
  Wg3 = readRDS("W3_picollo_GF.rds"),
  W1_5nn = readRDS("W1_eucl5nn.rds"),
  W2_5nn = readRDS("W2_corr5nn.rds"),
  W3_5nn = readRDS("W1_picollo5nn.rds"),
  Wg = readRDS("Wmat_granger.rds")
)

# ---------------------------------------------------------------
# 1. Log squared residuals: flooring chosen using training only
# ---------------------------------------------------------------
floor_by_asset <- vapply(seq_len(n), function(i) {
  sq <- Y[seq_len(train.l), i]^2
  positive <- sq[is.finite(sq) & sq > 0]
  if (!length(positive)) {
    stop("No nonzero squared training residual for asset ", i)
  }
  min(positive)
}, numeric(1))

log_sq <- function(x, floor_vector) {
  stopifnot(is.matrix(x), ncol(x) == length(floor_vector))
  log(sweep(x^2, 2L, floor_vector, FUN = pmax))
}
logY <- log_sq(Y, floor_by_asset)

# ---------------------------------------------------------------
# 2. Ex-post realised volatility targets
# ---------------------------------------------------------------
Ytest <- Y[(train.l + 1L):(train.l + out.l), , drop = FALSE]
RV5abs <- matrix(NA_real_, nrow = out.l, ncol = n)
RV5sq <- matrix(NA_real_, nrow = out.l, ncol = n)

for (k in seq_len(out.l)) {
  t <- train.l + k
  window <- Y[(t - 4L):t, , drop = FALSE]
  RV5abs[k, ] <- colMeans(abs(window))^2
  RV5sq[k, ] <- colMeans(window^2)
}

# EWMA variance (not standard deviation): initialise and propagate
# through the training period before extracting test-period values.
lambda <- 0.74
EWMA <- matrix(NA_real_, nrow = out.l, ncol = n)
prev <- Y[1L, ]^2
for (t in 2:nrow(Y)) {
  prev <- lambda * prev + (1 - lambda) * Y[t, ]^2
  if (t > train.l) EWMA[t - train.l, ] <- prev
}

# Identical first two targets are intentionally retained for
# compatibility with earlier output column names.
log_RV <- logY[(train.l + 1L):(train.l + out.l), , drop = FALSE]
proxies <- list(
  log_squared_residual = log_RV,
  RV = log_RV,
  RV5_abs = log(pmax(RV5abs, 1e-12)),
  RV5_sq = log(pmax(RV5sq, 1e-12)),
  EWMA = log(pmax(EWMA, 1e-12))
)

metric <- function(E) {
  if (!is.matrix(E) || any(!is.finite(E))) {
    stop("Non-finite or incorrectly shaped forecast errors")
  }
  c(RMSFE = mean(sqrt(colMeans(E^2))),
    MAFE = mean(colMeans(abs(E))))
}

# ---------------------------------------------------------------
# 3. GMM estimation and fixed-parameter forecasts by network
# ---------------------------------------------------------------
train <- t(logY[seq_len(train.l), , drop = FALSE])  # n x T
results <- list()
summaries <- list()
diagnostics <- list()

for (w_name in names(weight_matrices)) {
  cat("\nWeight matrix:", w_name, "\n")
  W <- as.matrix(weight_matrices[[w_name]])
  if (!identical(dim(W), c(n, n)) || any(!is.finite(W))) {
    stop("Invalid weight matrix: ", w_name)
  }

  model <- GMM_SDPD_2SLS_ARCH_ind_timelags(
    Y = train, X = NULL, W = W,
    info = list(ksy = 10, ksx = 10, stl = 0, tl = 1, ted = 0)
  )
  th <- as.numeric(model$theta)
  if (length(th) != n + 1L || any(!is.finite(th))) {
    stop("Invalid GMM coefficients for ", w_name)
  }
  rho <- th[1L]
  gamma <- th[2:(n + 1L)]
  A <- diag(n) - rho * W
  rc <- rcond(A)

  # Preserve a diagnostic record even when forecasts are skipped.
  diagnostics[[w_name]] <- data.frame(
    Weight = w_name, rho = rho, rcond = rc,
    temporal_radius = NA_real_, status = "not_evaluated"
  )

  if (!is.finite(rc) || rc < 1e-4) {
    diagnostics[[w_name]]$status <- "ill_conditioned_spatial_multiplier"
    warning(w_name, ": rcond(I-rho W)=", signif(rc, 5),
            "; forecasts omitted")
    results[[w_name]] <- list(model = model, rho = rho,
                              gamma = gamma, rcond = rc,
                              radius = NA_real_, forecasts = NULL)
    next
  }

  # y_t = mu + rho W y_t + diag(gamma) y_(t-1) + u_t
  # Recover asset-specific intercepts from training observations.
  current <- train[, 2:train.l, drop = FALSE]
  lagged <- train[, 1:(train.l - 1L), drop = FALSE]
  mu <- rowMeans(A %*% current - sweep(lagged, 1L, gamma, FUN = "*"))
  if (any(!is.finite(mu))) {
    diagnostics[[w_name]]$status <- "invalid_intercept"
    warning(w_name, ": non-finite intercept; forecasts omitted")
    results[[w_name]] <- list(model = model, rho = rho,
                              gamma = gamma, mu = mu,
                              rcond = rc, radius = NA_real_,
                              forecasts = NULL)
    next
  }

  B <- solve(A, diag(gamma, nrow = n))
  radius <- max(Mod(eigen(B, only.values = TRUE)$values))
  diagnostics[[w_name]]$temporal_radius <- radius

  if (!is.finite(radius) || radius >= 1) {
    diagnostics[[w_name]]$status <- "temporally_unstable"
    warning(w_name, ": temporal spectral radius=", signif(radius, 6),
            "; forecasts omitted; constrained re-estimation needed")
    results[[w_name]] <- list(model = model, rho = rho,
                              gamma = gamma, mu = mu,
                              rcond = rc, radius = radius,
                              forecasts = NULL)
    next
  }

  # Compute the spatial multiplier once; at each k, use only
  # the observation at t-1 for the forecast of time t.
  A_inv <- solve(A)
  forecasts <- matrix(NA_real_, nrow = out.l, ncol = n)
  for (k in seq_len(out.l)) {
    previous_log_y <- logY[train.l + k - 1L, ]
    forecasts[k, ] <- as.vector(A_inv %*% (mu + gamma * previous_log_y))
  }
  if (any(!is.finite(forecasts))) {
    stop("Non-finite forecasts for ", w_name)
  }

  scores <- lapply(proxies, function(target) metric(target - forecasts))
  summaries[[w_name]] <- data.frame(
    Weight = w_name,
    rho = rho,
    temporal_radius = radius,
    rcond = rc,
    RMSFE = unname(scores$log_squared_residual["RMSFE"]),
    MAFE = unname(scores$log_squared_residual["MAFE"]),
    RMSFE_RV = unname(scores$RV["RMSFE"]),
    MAFE_RV = unname(scores$RV["MAFE"]),
    RMSFE_RV5_abs = unname(scores$RV5_abs["RMSFE"]),
    MAFE_RV5_abs = unname(scores$RV5_abs["MAFE"]),
    RMSFE_RV5_sq = unname(scores$RV5_sq["RMSFE"]),
    MAFE_RV5_sq = unname(scores$RV5_sq["MAFE"]),
    RMSFE_EWMA = unname(scores$EWMA["RMSFE"]),
    MAFE_EWMA = unname(scores$EWMA["MAFE"])
  )
  diagnostics[[w_name]]$status <- "forecasted"
  results[[w_name]] <- list(
    model = model, mu = mu, rho = rho, gamma = gamma,
    rcond = rc, radius = radius, forecasts = forecasts,
    scores = scores
  )
  cat("  radius:", signif(radius, 6),
      " rcond:", signif(rc, 6),
      " RMSFE:", signif(scores$RV["RMSFE"], 6), "\n")
}

summary_table <- dplyr::bind_rows(summaries)
diagnostic_table <- dplyr::bind_rows(diagnostics)

cat("\n=== Stability diagnostics (all networks) ===\n")
print(diagnostic_table, row.names = FALSE, digits = 8)
cat("\n=== Forecast accuracy (stable networks only) ===\n")
print(summary_table, row.names = FALSE, digits = 8)

saveRDS(
  list(summary = summary_table, diagnostics = diagnostic_table,
       results = results, proxies = proxies,
       lambda = lambda, floor_by_asset = floor_by_asset,
       training_length = train.l, test_length = out.l),
  "DSTARCH_corrected_forecasting_results.rds"
)
write.csv(summary_table,
          "DSTARCH_corrected_forecasting_summary.csv",
          row.names = FALSE)
write.csv(diagnostic_table,
          "DSTARCH_corrected_forecasting_diagnostics.csv",
          row.names = FALSE)

for (name in names(weight_matrices)) {
  W <- as.matrix(weight_matrices[[name]])

  cat(
    name,
    "Row sums:", range(rowSums(W)),
    "Spectral radius:", max(Mod(eigen(W, only.values = TRUE)$values)),
    "\n"
  )
}