# =============================================================================
# SpEGARCH: fixed-parameter recursive one-step-ahead forecasting, all networks
# Requires functions_spegarch.R in the working directory.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(Rsolnp)
  library(nleqslv)
})
source('functions_spegarch.R')

# ------------------------------ SETTINGS ------------------------------------
residual_file <- 'residuals_VAR1.rds'
out_l <- 252L
n_init <- 30L                 # random candidate points evaluated
n_starts <- 5L                # candidates actually optimised
burn_in <- 5L
base_seed <- 123L
out_dir <- 'SpEGARCH_corrected_results'
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

weight_files <- c(
  W1='W1_eucl.rds',
  W2='W1_corr.rds',
  W3='W1_picollo.rds',
  Wg1='W1_eucl_GF.rds',
  Wg2='W2_corr_GF.rds',
  Wg3='W3_picollo_GF.rds',
  W1_5nn='W1_eucl5nn.rds',
  W2_5nn='W2_corr5nn.rds',
  W3_5nn='W1_picollo5nn.rds',
  Wg='Wmat_granger.rds'
)

# --------------------------- LOAD AND VALIDATE -------------------------------
if (!file.exists(residual_file)) stop('Missing: ', residual_file)
missing_files <- weight_files[!file.exists(weight_files)]
if (length(missing_files)) stop('Missing weight files: ', paste(missing_files,collapse=', '))
residuals <- as.matrix(readRDS(residual_file))  # rows=time, columns=assets
if (!is.numeric(residuals) || any(!is.finite(residuals)))
  stop('Residuals must be a finite numeric matrix.')
if (any(residuals == 0))
  stop('Exact zero residuals detected. Current inversion requires nonzero observations.')
TT <- nrow(residuals)
n <- ncol(residuals)
train_l <- TT-out_l
if (n < 2L || train_l <= burn_in+1L) stop('Insufficient data.')
Y_train <- t(residuals[seq_len(train_l),,drop=FALSE])
Y_test <- residuals[seq.int(train_l+1L,TT),,drop=FALSE]

# IMPORTANT: weight files must have been constructed using training data only.
# For meaningful model comparisons, the residuals, ordering, splits, and proxies
# must match those used for DSTARCH, DCC and STGARCH.

# ----------------------------- LOG PROXIES ----------------------------------
# Primary proxy: log of squared observed residual (variance proxy).
log_RV <- log(pmax(Y_test^2,1e-12))

# Five-day trailing absolute-return proxy, squared to yield a variance proxy.
# For the first four test dates, the window includes training observations.
RV5_abs <- matrix(NA_real_,out_l,n)
RV5_sq <- matrix(NA_real_,out_l,n)
for (j in seq_len(out_l)) {
  idx <- seq.int(train_l+j-4L,train_l+j)
  block <- residuals[idx,,drop=FALSE]
  RV5_abs[j,] <- colMeans(abs(block))^2
  RV5_sq[j,] <- colMeans(block^2)
}
log_RV5_abs <- log(pmax(RV5_abs,1e-12))
log_RV5_sq <- log(pmax(RV5_sq,1e-12))

# EWMA is ALREADY a variance; do not square it again.
# Initialize using training-only sample variance; update with realized return.
EWMA_var <- matrix(NA_real_,out_l,n)
v_prev <- apply(residuals[seq_len(train_l),,drop=FALSE],2,var)
for (j in seq_len(out_l)) {
  v_prev <- .94*v_prev+.06*Y_test[j,]^2
  EWMA_var[j,] <- v_prev
}
log_EWMA <- log(pmax(EWMA_var,1e-12))

# ------------------------------ EVALUATION ----------------------------------
eval_metric <- function(log_h,log_proxy) {
  if (!identical(dim(log_h),dim(log_proxy))) stop('Forecast/proxy dimension mismatch.')
  if (any(!is.finite(log_h)) || any(!is.finite(log_proxy)))
    stop('Nonfinite forecast or proxy values.')
  e <- log_h-log_proxy
  c(RMSFE=sqrt(mean(e^2)),MAFE=mean(abs(e)))
}

# ---------------------------- FORECAST ROUTINE -------------------------------
forecast_one_network <- function(fit,W_1,W_2) {
  p <- fit$params
  A <- diag(n)-as.numeric(p['lambda_0'])*W_2
  if (!is.finite(rcond(A)) || rcond(A)<1e-10)
    stop('Near-singular spatial multiplier.')
  # One factorisation per network, not a matrix inversion each period.
  eps_current <- as.numeric(fit$eps[,ncol(fit$eps)])
  logh_current <- as.numeric(fit$log_cond_var[,ncol(fit$log_cond_var)])
  log_forecasts <- matrix(NA_real_,out_l,n)
  filtered_eps <- matrix(NA_real_,out_l,n)
  filtered_logh <- matrix(NA_real_,out_l,n)
  for (j in seq_len(out_l)) {
    # Predict t+1 using only filtered state through t.
    lag_innov <- as.numeric(p['theta_1'])*eps_current +
      (abs(eps_current)-sqrt(2/pi))
    rhs <- as.numeric(p['alpha']) +
      as.numeric(p['rho_1'])*lag_innov +
      as.numeric(p['lambda_1'])*logh_current
    logh_pred <- as.numeric(solve(A,rhs))
    if (any(!is.finite(logh_pred))) stop('Nonfinite forecast at j=',j)
    log_forecasts[j,] <- logh_pred

    # Observe Y_{t+1} ONLY after forecast is stored, then filter next state.
    filtered <- g(
      y=as.numeric(Y_test[j,]),
      alpha=as.numeric(p['alpha']),
      rhoW_1=as.numeric(p['rho_0'])*W_1,
      lambdaW_2=as.numeric(p['lambda_0'])*W_2,
      rho_1=as.numeric(p['rho_1']),
      lambda_1=as.numeric(p['lambda_1']),
      theta_0=as.numeric(p['theta_0']),
      theta_1=as.numeric(p['theta_1']),
      xi=1,
      eps_lag=eps_current,
      logh_lag=logh_current
    )
    if (is.null(filtered)) stop('OOS filtering failed at j=',j)
    eps_current <- as.numeric(filtered$eps)
    logh_current <- as.numeric(filtered$lnh)
    filtered_eps[j,] <- eps_current
    filtered_logh[j,] <- logh_current
  }
  colnames(log_forecasts) <- colnames(residuals)
  colnames(filtered_eps) <- colnames(residuals)
  colnames(filtered_logh) <- colnames(residuals)
  list(log_forecasts=log_forecasts,
       forecasts=exp(log_forecasts),
       filtered_eps=filtered_eps,
       filtered_logh=filtered_logh)
}

# ------------------------------ ALL NETWORKS --------------------------------
results_list <- list()
summary_rows <- list()
for (k in seq_along(weight_files)) {
  w_name <- names(weight_files)[k]
  message('\n======= Processing ',w_name,' (',k,'/',length(weight_files),') =======')
  W_1 <- as.matrix(readRDS(weight_files[[k]]))
  W_2 <- W_1 # Original restricted specification: same network for both terms.
  if (!is.numeric(W_1) || !identical(dim(W_1),c(n,n)) ||
      any(!is.finite(W_1))) stop('Invalid weight matrix: ',w_name)
  # Optional: if W has dimnames, check ordering against residual columns.
  if (!is.null(colnames(residuals)) && !is.null(rownames(W_1)) &&
      !identical(rownames(W_1),colnames(residuals)))
    stop('Asset ordering mismatch in rows of ',w_name)
  if (!is.null(colnames(residuals)) && !is.null(colnames(W_1)) &&
      !identical(colnames(W_1),colnames(residuals)))
    stop('Asset ordering mismatch in columns of ',w_name)

  # Avoid losing results if a later network fails. Each network has own checkpoint.
  result <- tryCatch({
    fit <- estimate_spEGARCH(
      Y=Y_train,W_1=W_1,W_2=W_2,
      n_init=n_init,n_starts=n_starts,
      seed=base_seed+k-1L,burn_in=burn_in
    )
    pred <- forecast_one_network(fit,W_1,W_2)
    met <- list(
      RV=eval_metric(pred$log_forecasts,log_RV),
      RV5_abs=eval_metric(pred$log_forecasts,log_RV5_abs),
      RV5_sq=eval_metric(pred$log_forecasts,log_RV5_sq),
      EWMA=eval_metric(pred$log_forecasts,log_EWMA)
    )
    list(fit=fit,pred=pred,metrics=met)
  },error=function(e) {
    message('FAILED ',w_name,': ',conditionMessage(e))
    NULL
  })
  if (is.null(result)) next
  results_list[[w_name]] <- result
  saveRDS(result,file.path(out_dir,paste0('SpEGARCH_',w_name,'.rds')))
  m <- result$metrics
  summary_rows[[w_name]] <- data.frame(
    Weight=w_name,
    RMSFE_RV=unname(m$RV['RMSFE']),MAFE_RV=unname(m$RV['MAFE']),
    RMSFE_RV5_abs=unname(m$RV5_abs['RMSFE']),
    MAFE_RV5_abs=unname(m$RV5_abs['MAFE']),
    RMSFE_RV5_sq=unname(m$RV5_sq['RMSFE']),
    MAFE_RV5_sq=unname(m$RV5_sq['MAFE']),
    RMSFE_EWMA=unname(m$EWMA['RMSFE']),
    MAFE_EWMA=unname(m$EWMA['MAFE']),
    NLL=result$fit$nll,
    convergence=result$fit$convergence,
    n_success=result$fit$n_success,
    spatial_radius=result$fit$spatial_radius,
    temporal_radius=result$fit$temporal_radius
  )
  summary_table <- dplyr::bind_rows(summary_rows)
  write.csv(summary_table,file.path(out_dir,'SpEGARCH_all_weights_summary.csv'),row.names=FALSE)
  message('DONE ',w_name,' RMSFE=',round(m$RV['RMSFE'],5),
          ' MAFE=',round(m$RV['MAFE'],5))
}
if (!length(summary_rows)) stop('All network estimations failed; see messages above.')
summary_table <- dplyr::bind_rows(summary_rows)
print(summary_table,row.names=FALSE,digits=6)
saveRDS(results_list,file.path(out_dir,'SpEGARCH_all_weights_full_results.rds'))
cat('\nCompleted ',length(results_list),'/',length(weight_files),' networks.\n',sep='')
cat('Results saved to: ',normalizePath(out_dir),'\n',sep='')
