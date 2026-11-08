# =============================================================================
# CAPORIN-PARUOLO PROXIMITY-STRUCTURED BEKK
# Fixed-parameter, recursive one-step-ahead variance forecasts across networks
# =============================================================================
# Model:
# A = a0 I + a1 W; B = b0 I + b1 W
# C = (I - s1 W)^(-1) diag(exp(v)) (I - s1 W)^(-T)
# S_t = C + (A x_{t-1})(A x_{t-1})' + B S_{t-1} B'
# x_t = U_t - mean(U_train) (training-only mean)
#
# Primary evaluation uses UN-CENTRED U_t^2, as in the other model scripts.
# The original nonnegative dynamic-coefficient bounds are retained.
# Time units: seconds. Network construction/loading is not timed as estimation.
# =============================================================================

set.seed(20261016)

# ---------------------------- SETTINGS --------------------------------------
H <- 252L
nstarts <- 5L
maxit <- 400L
floor_value <- 1e-12
lambda_ewma <- 0.94
stability_margin <- 0.999
out_dir <- 'Caporin_BEKK_corrected_results'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

weight_files <- c(
  W1 = 'W1_eucl.rds',
  W2 = 'W1_corr.rds',
  W3 = 'W1_picollo.rds',
  Wg1 = 'W1_eucl_GF.rds',
  Wg2 = 'W2_corr_GF.rds',
  Wg3 = 'W3_picollo_GF.rds',
  W1_5nn = 'W1_eucl5nn.rds',
  W2_5nn = 'W2_corr5nn.rds',
  W3_5nn = 'W1_picollo5nn.rds',
  Wg = 'Wmat_granger.rds'
)

# ---------------------------- DATA ------------------------------------------
if (!file.exists('residuals_VAR1.rds')) stop('Missing residuals_VAR1.rds')
missing_files <- weight_files[!file.exists(weight_files)]
if (length(missing_files)) stop('Missing W files: ', paste(missing_files, collapse = ', '))

U <- as.matrix(readRDS('residuals_VAR1.rds'))
storage.mode(U) <- 'double'
if (!is.numeric(U) || any(!is.finite(U)) || ncol(U) < 2L || nrow(U) <= H + 101L)
  stop('Invalid residual matrix or insufficient training data.')

TT <- nrow(U)
n <- ncol(U)
Ttrain <- TT - H
idx <- seq.int(Ttrain + 1L, TT)
cat('n =', n, '; Ttrain =', Ttrain, '; H =', H, '\n')
cat('Median daily variance:', signif(median(diag(cov(U))), 4), '\n')

centre <- colMeans(U[seq_len(Ttrain), , drop = FALSE])
Z <- sweep(U, 2L, centre, '-')
Ztrain <- Z[seq_len(Ttrain), , drop = FALSE]

Wlist <- lapply(names(weight_files), function(nm) {
  W <- as.matrix(readRDS(weight_files[[nm]]))
  storage.mode(W) <- 'double'
  if (!identical(dim(W), c(n, n)) || any(!is.finite(W)))
    stop('Invalid W matrix: ', nm)
  if (!is.null(colnames(U))) {
    if (!is.null(rownames(W)) && !identical(rownames(W), colnames(U)))
      stop('W row ordering differs from U: ', nm)
    if (!is.null(colnames(W)) && !identical(colnames(W), colnames(U)))
      stop('W column ordering differs from U: ', nm)
  }
  rs <- rowSums(W)
  if (any(abs(rs - 1) > 1e-8 & abs(rs) > 1e-12))
    warning('Some rows of ', nm, ' sum to neither zero nor one.')
  W
})
names(Wlist) <- names(weight_files)

# ---------------------------- COMMON EVALUATION TARGETS ---------------------
# Use U, not Z, to match the uncentred residual-based proxies used elsewhere.
RV <- U[idx, , drop = FALSE]^2
RV5_abs <- RV5_sq <- matrix(NA_real_, H, n)
for (j in seq_len(H)) {
  tt <- Ttrain + j
  block <- U[seq.int(tt - 4L, tt), , drop = FALSE]
  RV5_abs[j, ] <- colMeans(abs(block))^2
  RV5_sq[j, ] <- colMeans(block^2)
}

# EWMA: training-only variance initialization, as in the Sp-EGARCH script.
EWMA <- matrix(NA_real_, H, n)
prev <- apply(U[seq_len(Ttrain), , drop = FALSE], 2L, var)
for (j in seq_len(H)) {
  prev <- lambda_ewma * prev + (1 - lambda_ewma) * U[Ttrain + j, ]^2
  EWMA[j, ] <- prev
}

targets <- list(RV = RV, RV5_abs = RV5_abs, RV5_sq = RV5_sq, EWMA = EWMA)
for (nm in names(targets)) colnames(targets[[nm]]) <- colnames(U)

# Existing comparison uses uncorrected log errors and asset-mean RMSFE.
# Extra pooled and Gaussian-offset metrics are diagnostic only.
OFFSET <- digamma(0.5) + log(2) # E[log(chi-square_1)]
score <- function(pred, target) {
  if (!identical(dim(pred), dim(target))) stop('Forecast/target dimension mismatch.')
  if (any(!is.finite(pred)) || any(pred <= 0) || any(!is.finite(target)) || any(target < 0))
    stop('Invalid forecast or target.')
  log_pred <- log(pmax(pred, floor_value))
  log_target <- log(pmax(target, floor_value))
  e <- log_pred - log_target
  ec <- e + OFFSET # log(h)-[log(RV)-OFFSET]
  p <- pmax(pred, floor_value)
  c(
    RMSFE = mean(sqrt(colMeans(e^2))),
    MAFE = mean(abs(e)),
    RMSFE_pooled = sqrt(mean(e^2)),
    RMSFE_corr = mean(sqrt(colMeans(ec^2))),
    MAFE_corr = mean(abs(ec)),
    QLIKE = mean(target / p + log(p)) # equivalent up to target-only constant
  )
}

# ---------------------------- MODEL -----------------------------------------
# n is fixed throughout this script; the parameter vector has length n+5.
matrices <- function(theta, W) {
  if (length(theta) != n + 5L || any(!is.finite(theta))) return(NULL)
  M <- diag(n) - theta[5L] * W
  rc <- tryCatch(rcond(M), error = function(e) NA_real_)
  if (!is.finite(rc) || rc < 1e-8) return(NULL)
  Minv <- tryCatch(solve(M), error = function(e) NULL)
  if (is.null(Minv)) return(NULL)
  v <- exp(theta[5L + seq_len(n)])
  if (any(!is.finite(v)) || any(v <= 0)) return(NULL)
  C <- tcrossprod(sweep(Minv, 2L, sqrt(v), '*'))
  list(
    A = theta[1L] * diag(n) + theta[2L] * W,
    B = theta[3L] * diag(n) + theta[4L] * W,
    C = (C + t(C)) / 2
  )
}

# Since A and B are polynomials in the same W, the eigenvalues of
# A (x) A + B (x) B are lambdaA_i*lambdaA_j + lambdaB_i*lambdaB_j.
rad <- function(theta, mu) {
  la <- theta[1L] + theta[2L] * mu
  lb <- theta[3L] + theta[4L] * mu
  max(Mod(outer(la, la) + outer(lb, lb)))
}

chol_eval <- function(S) {
  S <- (S + t(S)) / 2
  R <- tryCatch(chol(S), error = function(e) NULL)
  if (is.null(R)) return(NULL)
  list(R = R, logdet = 2 * sum(log(diag(R))))
}

initial_cov <- function(Utrain) {
  S <- cov(Utrain[seq_len(min(100L, nrow(Utrain))), , drop = FALSE])
  if (!identical(dim(S), c(n, n)) || any(!is.finite(S)))
    stop('Invalid initial covariance.')
  S <- (S + t(S)) / 2
  if (is.null(chol_eval(S))) {
    # Numerical ridge only if needed, with a documented scale.
    ridge <- max(1e-10, 1e-8 * mean(diag(cov(Utrain))))
    S <- S + diag(ridge, n)
    if (is.null(chol_eval(S))) stop('Initial covariance not positive definite.')
  }
  S
}

nll <- function(theta, Utrain, W, S0, mu) {
  radius <- rad(theta, mu)
  if (!is.finite(radius) || radius >= stability_margin) return(1e12)
  m <- matrices(theta, W)
  if (is.null(m)) return(1e12)
  S <- S0
  ans <- 0
  for (k in seq.int(2L, nrow(Utrain))) {
    xlag <- as.vector(m$A %*% Utrain[k - 1L, ])
    S <- m$C + tcrossprod(xlag) + m$B %*% S %*% t(m$B)
    d <- chol_eval(S)
    if (is.null(d)) return(1e12)
    z <- forwardsolve(t(d$R), Utrain[k, ])
    term <- 0.5 * (d$logdet + sum(z^2) + n * log(2 * pi))
    if (!is.finite(term)) return(1e12)
    ans <- ans + term
  }
  as.numeric(ans)
}

fit_network <- function(Utrain, W, nstarts = 5L, maxit = 400L) {
  S0 <- initial_cov(Utrain)
  mu <- eigen(W, only.values = TRUE)$values
  variances <- diag(cov(Utrain))
  if (any(!is.finite(variances)) || any(variances <= 0))
    stop('Nonpositive sample variance.')
  vbar <- log(variances)

  # Retain original positive dynamic-coefficient restrictions.
  lower <- c(rep(0.001, 4L), -0.5, vbar - 8)
  upper <- c(0.9, 0.9, 0.98, 0.9, 0.5, vbar + 1)

  starts <- vector('list', nstarts)
  starts[[1L]] <- c(0.12, 0.04, 0.65, 0.04, 0, vbar + log(0.2))
  if (nstarts > 1L) for (i in seq.int(2L, nstarts)) {
    starts[[i]] <- c(
      runif(2L, 0.02, 0.25),
      runif(2L, 0.02, 0.5),
      runif(1L, -0.2, 0.2),
      vbar + runif(n, -3, -1)
    )
  }

  fits <- lapply(starts, function(x) {
    tryCatch(
      optim(
        x, nll, Utrain = Utrain, W = W, S0 = S0, mu = mu,
        method = 'L-BFGS-B', lower = lower, upper = upper,
        control = list(maxit = maxit)
      ),
      error = function(e) NULL
    )
  })
  valid <- Filter(function(x) {
    !is.null(x) && is.finite(x$value) && x$value < 1e11 &&
      is.finite(rad(x$par, mu)) && rad(x$par, mu) < stability_margin &&
      !is.null(matrices(x$par, W))
  }, fits)
  if (!length(valid)) return(NULL)
  converged <- Filter(function(x) identical(x$convergence, 0L), valid)
  candidates <- if (length(converged)) converged else valid
  best <- candidates[[which.min(vapply(candidates, `[[`, numeric(1L), 'value'))]]
  best$lower <- lower
  best$upper <- upper
  best$mu <- mu
  best$n_valid <- length(valid)
  best$n_converged <- length(converged)
  best$selected_converged <- best$convergence == 0L
  best
}

# ---------------------------- FORECAST --------------------------------------
# At day k, S_k uses observations only through day k-1.
forecast_fixed <- function(theta, Z, W, Ttrain, H) {
  m <- matrices(theta, W)
  if (is.null(m)) stop('Invalid fitted matrices.')
  S <- initial_cov(Z[seq_len(Ttrain), , drop = FALSE])
  out <- matrix(NA_real_, H, n)
  for (k in seq.int(2L, Ttrain + H)) {
    xlag <- as.vector(m$A %*% Z[k - 1L, ])
    S <- m$C + tcrossprod(xlag) + m$B %*% S %*% t(m$B)
    if (is.null(chol_eval(S))) stop('Nonpositive forecast covariance at time ', k)
    if (k > Ttrain) out[k - Ttrain, ] <- diag(S)
  }
  if (any(!is.finite(out)) || any(out <= 0)) stop('Invalid forecast variances.')
  colnames(out) <- colnames(Z)
  out
}

# Diagnostic only: not included in forecast timing.
lookahead_check <- function(theta, W, reference_forecast) {
  if (H < 2L) return(invisible(TRUE))
  j <- min(100L, H - 1L)
  Z2 <- Z
  Z2[Ttrain + j, ] <- Z2[Ttrain + j, ] + 10 * pmax(abs(Z2[Ttrain + j, ]), 1e-4)
  altered <- forecast_fixed(theta, Z2, W, Ttrain, H)
  if (!isTRUE(all.equal(reference_forecast[j, ], altered[j, ], tolerance = 1e-10)))
    stop('Look-ahead detected: same-day forecast changed.')
  if (isTRUE(all.equal(reference_forecast[j + 1L, ], altered[j + 1L, ], tolerance = 1e-10)))
    warning('Next-day forecast unchanged after perturbation; inspect A and data.')
  invisible(TRUE)
}

# ---------------------------- RUN NETWORK -----------------------------------
run_one <- function(wname) {
  W <- Wlist[[wname]]
  tryCatch({
    t0 <- proc.time()[['elapsed']]
    fit <- fit_network(Ztrain, W, nstarts = nstarts, maxit = maxit)
    estimation_time <- proc.time()[['elapsed']] - t0
    if (is.null(fit)) stop('No valid optimisation.')

    theta <- fit$par
    radius <- rad(theta, fit$mu)
    kpar <- length(theta)
    Teff <- Ttrain - 1L
    at_bound <- sum(abs(theta - fit$lower) < 1e-6 |
                      abs(theta - fit$upper) < 1e-6)

    t0 <- proc.time()[['elapsed']]
    pred <- forecast_fixed(theta, Z, W, Ttrain, H)
    forecasting_time <- proc.time()[['elapsed']] - t0

    t0 <- proc.time()[['elapsed']]
    lookahead_check(theta, W, pred)
    diagnostic_time <- proc.time()[['elapsed']] - t0

    sc <- lapply(targets, function(target) score(pred, target))
    summ <- do.call(rbind, lapply(names(sc), function(target_name) {
      data.frame(Weight = wname, Target = target_name,
                 as.list(sc[[target_name]]), check.names = FALSE)
    }))
    rownames(summ) <- NULL

    params <- data.frame(
      Weight = wname, a0 = theta[1L], a1 = theta[2L],
      b0 = theta[3L], b1 = theta[4L], s1 = theta[5L],
      LogLik = -fit$value,
      AIC = 2 * fit$value + 2 * kpar,
      BIC = 2 * fit$value + log(Teff) * kpar,
      convergence = fit$convergence,
      n_valid_starts = fit$n_valid,
      n_converged_starts = fit$n_converged,
      BEKK_radius = radius,
      n_at_bound = at_bound,
      EstimationTime = estimation_time,
      ForecastTime = forecasting_time,
      DiagnosticTime = diagnostic_time,
      TotalModelTime = estimation_time + forecasting_time
    )
    if (!fit$selected_converged)
      warning(wname, ': no optimiser start converged; selected a non-converged fit.')
    if (at_bound > 0L)
      warning(wname, ': ', at_bound, ' parameter(s) at a bound.')

    result <- list(
      parameters = params,
      coefficients = setNames(theta, c('a0','a1','b0','b1','s1',
                                     paste0('v', seq_len(n)))),
      metrics = summ,
      forecasts = pred,
      timing = c(EstimationTime = estimation_time,
                 ForecastTime = forecasting_time,
                 DiagnosticTime = diagnostic_time,
                 TotalModelTime = estimation_time + forecasting_time)
    )
    saveRDS(result, file.path(out_dir, paste0('Caporin_', wname, '.rds')))
    message('DONE ', wname, ' | RMSFE=', signif(summ$RMSFE[summ$Target == 'RV'], 6),
            ' | MAFE=', signif(summ$MAFE[summ$Target == 'RV'], 6),
            ' | Estimation=', round(estimation_time, 2), 's')
    result
  }, error = function(e) {
    message('FAILED ', wname, ': ', conditionMessage(e))
    list(error = conditionMessage(e), Weight = wname)
  })
}

# ---------------------------- PARALLEL EXECUTION ----------------------------
wall_start <- Sys.time()

# Fixed seeds per network, independent of parallel scheduling.
network_seeds <- setNames(20261016L + seq_along(Wlist) - 1L, names(Wlist))
run_seeded <- function(wname) {
  set.seed(network_seeds[[wname]])
  run_one(wname)
}

if (.Platform$OS.type == 'windows') {
  cores <- 1L
} else {
  detected <- parallel::detectCores(logical = TRUE)
  if (is.na(detected)) detected <- 2L
  cores <- max(1L, min(length(Wlist), detected - 1L))
}
cat('Running with ', cores, ' worker(s).\n', sep = '')
if (cores == 1L) {
  results <- lapply(names(Wlist), run_seeded)
} else {
  results <- parallel::mclapply(names(Wlist), run_seeded,
                                mc.cores = cores, mc.preschedule = FALSE)
}
names(results) <- names(Wlist)
parallel_wall_seconds <- as.numeric(difftime(Sys.time(), wall_start, units = 'secs'))

ok <- Filter(function(x) is.list(x) && !is.null(x$metrics), results)
if (!length(ok)) stop('All network fits failed. Check diagnostic messages.')

summary_table <- do.call(rbind, lapply(ok, `[[`, 'metrics'))
parameter_table <- do.call(rbind, lapply(ok, `[[`, 'parameters'))
rownames(summary_table) <- NULL
rownames(parameter_table) <- NULL

cat('\nPrimary RV comparison:\n')
print(summary_table[summary_table$Target == 'RV',
                    c('Weight','RMSFE','MAFE')], row.names = FALSE, digits = 7)
cat('\nEstimation and forecasting times (seconds):\n')
print(parameter_table[, c('Weight','EstimationTime','ForecastTime',
                          'DiagnosticTime','TotalModelTime')],
      row.names = FALSE, digits = 7)
cat('\nParallel wall-clock time (minutes):', round(parallel_wall_seconds / 60, 3), '\n')

saveRDS(list(
  summary = summary_table,
  parameters = parameter_table,
  forecasts = lapply(ok, `[[`, 'forecasts'),
  coefficients = lapply(ok, `[[`, 'coefficients'),
  targets = targets,
  train_length = Ttrain,
  test_length = H,
  centre = centre,
  parallel_wall_seconds = parallel_wall_seconds,
  cores = cores,
  errors = lapply(Filter(function(x) !is.null(x$error), results), `[[`, 'error'),
  settings = list(nstarts = nstarts, maxit = maxit,
                  stability_margin = stability_margin,
                  primary_RMSFE = 'asset_mean',
                  proxy = 'uncentred residual squared')
), file.path(out_dir, 'Caporin_BEKK_corrected_results.rds'))

write.csv(summary_table, file.path(out_dir, 'Caporin_BEKK_corrected_summary.csv'),
          row.names = FALSE)
write.csv(parameter_table, file.path(out_dir, 'Caporin_BEKK_corrected_parameters.csv'),
          row.names = FALSE)
cat('Saved results to: ', normalizePath(out_dir), '\n', sep = '')
