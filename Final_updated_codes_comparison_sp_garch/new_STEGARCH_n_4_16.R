# ================================================================
# SpEGARCH dimensionality experiment: n = 4,...,16
# Estimate on training observations, produce 252 one-step-ahead
# log-variance forecasts, and compute asset-averaged RMSFE / MAFE.
# Run: source('run_full_spEGARCH_forecast.R')
# Requires functions_spegarch_cluster.R (matching rebuilt API).
# ================================================================

suppressPackageStartupMessages({
  library(Rsolnp)
  library(nleqslv)
  library(lmtest)
})
source('functions_spegarch_cluster.R')
required_functions <- c('construct_W_GF', 'estimate_spEGARCH', 'sp_invert',
                        'sp_evaluate', 'sp_stability')
stopifnot(all(vapply(required_functions, exists, logical(1), mode='function')))

# ---------- Settings ----------
dimensions <- 4:16
horizon <- 252L
n_random <- 8L
n_optim <- 2L
root_tol <- 1e-12
out_dir <- 'results_spegarch_full'
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
# Set TRUE only if you intentionally wish to overwrite completed files.
overwrite <- FALSE

# ---------- Data ----------
residuals <- as.matrix(readRDS('residuals_VAR1.rds'))
asset_order <- readRDS('asset_order_dimensionality.rds')
stopifnot(is.numeric(residuals), all(is.finite(residuals)),
          nrow(residuals) > horizon+5L, length(asset_order) >= max(dimensions))
if (is.character(asset_order)) {
  stopifnot(!is.null(colnames(residuals)),
            all(asset_order[seq_len(max(dimensions))] %in% colnames(residuals)))
} else {
  stopifnot(is.numeric(asset_order),
            all(asset_order[seq_len(max(dimensions))] %in% seq_len(ncol(residuals))))
}
train_length <- nrow(residuals)-horizon

# Previous n=4 candidate, included only for n=4.
p_previous_4 <- c(rho_0=0.0022975166, rho_1=0.0095046657,
                  theta_0=-0.2156919455, theta_1=-0.5941941106,
                  alpha=-7.8664769790, lambda_1=0.0011275458,
                  lambda_0=0.0003729631)

# ---------- Forecasting ----------
# Forecast target: log conditional variance (n x H).
# At forecast origin t-1, known states are eps_(t-1), log h_(t-1).
# E[g(eps_t) | F_(t-1)] is taken as zero, as in the model's
# centered-innovation one-step-ahead forecasting convention.
# Hence (I - lambda_0 W2) E[log h_t | F_(t-1)] =
#   alpha*1 + rho_1*g(eps_(t-1)) + lambda_1*log h_(t-1).
# After making each forecast, filter the observed test return y_t
# using sp_invert() to update the state for the NEXT forecast.
# This is a rolling one-step-ahead forecast with FIXED parameters.
forecast_spegarch <- function(fit, Y_test_n_by_H, W1, W2, tol=1e-12) {
  p <- fit$params
  needed <- c('rho_0','rho_1','theta_0','theta_1','alpha','lambda_1','lambda_0')
  stopifnot(all(needed %in% names(p)))
  p <- p[needed]
  Y_test_n_by_H <- as.matrix(Y_test_n_by_H)
  n <- nrow(Y_test_n_by_H); H <- ncol(Y_test_n_by_H)
  stopifnot(identical(dim(W1), c(n,n)), identical(dim(W2), c(n,n)),
            nrow(fit$eps)==n, nrow(fit$log_cond_var)==n,
            all(is.finite(Y_test_n_by_H)), all(Y_test_n_by_H != 0))
  A <- diag(n)-unname(p['lambda_0'])*W2
  if (!is.finite(rcond(A)) || rcond(A)<1e-12)
    stop('Singular spatial multiplier during forecasting')
  eps_prev <- as.numeric(fit$eps[,ncol(fit$eps)])
  logh_prev <- as.numeric(fit$log_cond_var[,ncol(fit$log_cond_var)])
  y_prev <- NULL
  pred <- matrix(NA_real_, n, H)
  updated_eps <- matrix(NA_real_, n, H)
  updated_logh <- matrix(NA_real_, n, H)
  residuals <- rep(NA_real_, H)
  for (j in seq_len(H)) {
    lagged_shock <- unname(p['theta_1'])*eps_prev +
      (abs(eps_prev)-sqrt(2/pi))
    rhs <- rep(unname(p['alpha']), n) +
      unname(p['rho_1'])*lagged_shock +
      unname(p['lambda_1'])*logh_prev
    pred[,j] <- as.vector(solve(A,rhs))
    # Update AFTER forecast, using current observed test returns.
    y_now <- Y_test_n_by_H[,j]
    if (j==1L) {
      # For sp_invert, Y_lag is needed in addition to eps_lag.
      # The training endpoint must be supplied through fit$last_train_y.
      stopifnot(!is.null(fit$last_train_y), length(fit$last_train_y)==n)
      y_lag <- fit$last_train_y
    } else {
      y_lag <- y_prev
    }
    out <- sp_invert(y=y_now, p=p, W1=W1, W2=W2,
                     eps_lag=eps_prev, Y_lag=y_lag, root_tol=tol)
    if (!isTRUE(out$ok))
      stop(sprintf('Test-period inversion failed at step %d: %s (residual %s)',
                   j, out$reason, format(out$residual)))
    updated_eps[,j] <- out$eps
    updated_logh[,j] <- out$lnh
    residuals[j] <- out$residual
    eps_prev <- out$eps
    logh_prev <- out$lnh
    y_prev <- y_now
  }
  list(log_forecast=pred, eps_test=updated_eps,
       logh_test=updated_logh, max_test_root_residual=max(residuals))
}

# ---------- Evaluation: average per-asset RMSFE and MAFE ----------
score_forecast <- function(log_forecast, Y_test_n_by_H) {
  stopifnot(identical(dim(log_forecast), dim(Y_test_n_by_H)))
  # Same realized proxy convention used previously.
  log_RV <- log(pmax(Y_test_n_by_H^2, 1e-12))
  errors <- log_forecast-log_RV
  per_asset_rmsfe <- sqrt(rowMeans(errors^2))
  per_asset_mafe <- rowMeans(abs(errors))
  list(RMSFE=mean(per_asset_rmsfe), MAFE=mean(per_asset_mafe),
       per_asset_RMSFE=per_asset_rmsfe, per_asset_MAFE=per_asset_mafe,
       log_RV=log_RV, errors=errors)
}

# ---------- One dimension ----------
run_one_n <- function(n) {
  start_clock <- Sys.time()
  cat('\n========== n =', n, 'started', format(start_clock), '==========\n')
  # T x n for Granger/correlation network construction.
  Y_all <- residuals[,asset_order[seq_len(n)],drop=FALSE]
  Y_train_Tn <- Y_all[seq_len(train_length),,drop=FALSE]
  Y_test_Tn <- Y_all[(train_length+1L):nrow(Y_all),,drop=FALSE]
  network <- construct_W_GF(Y_train_Tn, max_lag=12L, p_threshold=0.05)
  W <- network$W
  stopifnot(identical(dim(W),c(n,n)), all(is.finite(W)))
  # n x T for SpEGARCH estimation; n x 252 for forecast evaluation.
  Y_train <- t(Y_train_Tn)
  Y_test <- t(Y_test_Tn)
  previous <- if (n==4L) p_previous_4 else NULL
  fit <- estimate_spEGARCH(
    Y=Y_train, W_1=W, W_2=W,
    n_random=n_random, n_optim=n_optim,
    previous_params=previous, seed=101L+n,
    root_tol=root_tol, verbose=TRUE
  )
  fit$last_train_y <- as.numeric(Y_train[,ncol(Y_train)])
  # A feasible estimate is not necessarily an optimizer-converged solution.
  if (!identical(as.integer(fit$convergence), 0L))
    warning('n=',n,': optimizer convergence code ',fit$convergence,
            '; interpret forecasting scores as provisional.')
  fc <- forecast_spegarch(fit,Y_test,W,W,tol=root_tol)
  scores <- score_forecast(fc$log_forecast,Y_test)
  elapsed <- as.numeric(difftime(Sys.time(),start_clock,units='mins'))
  result <- list(n=n, params=fit$params, fit=fit, network=network,
                 forecasts=fc, scores=scores, Y_test=Y_test,
                 train_length=train_length, test_length=horizon,
                 elapsed_minutes=elapsed, completed_at=Sys.time(),
                 forecast_protocol='fixed-parameter rolling 1-step; update after observing y_t',
                 target='log(pmax(y_t^2, 1e-12))')
  path <- file.path(out_dir,sprintf('SpEGARCH_n%02d.rds',n))
  saveRDS(result,path)
  cat(sprintf('n=%d | NLL=%.8f | RMSFE=%.8f | MAFE=%.8f | code=%s | %.1f min\n',
              n,fit$objective,scores$RMSFE,scores$MAFE,
              as.character(fit$convergence),elapsed))
  invisible(result)
}

# ---------- Resume-friendly sequential run ----------
for (n in dimensions) {
  result_path <- file.path(out_dir,sprintf('SpEGARCH_n%02d.rds',n))
  if (file.exists(result_path) && !overwrite) {
    cat('Skipping n=',n,' (already saved).\n')
    next
  }
  tryCatch(run_one_n(n), error=function(e) {
    msg <- sprintf('n=%d failed: %s',n,conditionMessage(e))
    message(msg)
    writeLines(msg,file.path(out_dir,sprintf('error_n%02d.txt',n)))
  })
}

# ---------- Summary, including forecast scores ----------
rows <- list()
for (n in dimensions) {
  path <- file.path(out_dir,sprintf('SpEGARCH_n%02d.rds',n))
  if (!file.exists(path)) next
  r <- tryCatch(readRDS(path),error=function(e) NULL)
  if (is.null(r) || is.null(r$scores)) next
  p <- r$params
  rows[[length(rows)+1L]] <- data.frame(
    n=n, model='SpEGARCH', RMSFE=r$scores$RMSFE, MAFE=r$scores$MAFE,
    NLL=r$fit$objective, convergence=r$fit$convergence,
    spatial_radius=unname(r$fit$stability['spatial']),
    temporal_radius=unname(r$fit$stability['temporal']),
    max_train_root_residual=r$fit$max_root_residual,
    max_test_root_residual=r$forecasts$max_test_root_residual,
    density=r$network$density, isolated=r$network$isolated,
    rho_0=unname(p['rho_0']), rho_1=unname(p['rho_1']),
    theta_0=unname(p['theta_0']), theta_1=unname(p['theta_1']),
    alpha=unname(p['alpha']), lambda_1=unname(p['lambda_1']),
    lambda_0=unname(p['lambda_0']), elapsed_minutes=r$elapsed_minutes
  )
}
if (length(rows)) {
  summary_results <- do.call(rbind,rows)
  write.csv(summary_results,file.path(out_dir,'SpEGARCH_comparison_n4_n16.csv'),
            row.names=FALSE)
  print(summary_results[,c('n','RMSFE','MAFE','NLL','convergence',
                          'spatial_radius','temporal_radius')],row.names=FALSE)
} else {
  message('No complete results to summarize.')
}
