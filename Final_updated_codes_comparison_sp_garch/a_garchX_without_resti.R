# SpGARCH-X: fixed-parameter one-step-ahead conditional variance forecasts
# 252 out-of-sample observations, same ten original networks.
# The spatial lag of lagged squared innovations enters the VARIANCE equation.
# Requires rugarch and the RDS files listed below.

suppressPackageStartupMessages(library(rugarch))

Y <- as.matrix(readRDS('residuals_VAR1.rds'))
storage.mode(Y) <- 'double'
stopifnot(all(is.finite(Y)), ncol(Y) >= 2L, nrow(Y) > 260L)
N <- ncol(Y); H <- 252L; T_train <- nrow(Y) - H
if (is.null(colnames(Y))) colnames(Y) <- paste0('Asset_', seq_len(N))

weight_matrices <- list(
  W1=readRDS('W1_eucl.rds'), W2=readRDS('W1_corr.rds'),
  W3=readRDS('W1_picollo.rds'), Wg1=readRDS('W1_eucl_GF.rds'),
  Wg2=readRDS('W2_corr_GF.rds'), Wg3=readRDS('W3_picollo_GF.rds'),
  W1_5nn=readRDS('W1_eucl5nn.rds'), W2_5nn=readRDS('W2_corr5nn.rds'),
  W3_5nn=readRDS('W1_picollo5nn.rds'), Wg=readRDS('Wmat_granger.rds')
)
for (name in names(weight_matrices)) {
  W <- as.matrix(weight_matrices[[name]])
  stopifnot(identical(dim(W), c(N,N)), all(is.finite(W)))
  weight_matrices[[name]] <- W  # no network alteration
}

# Common ex-post targets: compare identical definitions across all models.
eps <- 1e-12
Y_test <- Y[(T_train+1L):nrow(Y), , drop=FALSE]
RV5_abs <- RV5_sq <- matrix(NA_real_, H, N)
for (j in seq_len(H)) {
  window <- Y[(T_train+j-4L):(T_train+j), , drop=FALSE]
  RV5_abs[j, ] <- colMeans(abs(window))^2
  RV5_sq[j, ] <- colMeans(window^2)
}
lambda <- 0.94  # must match DCC, DSTARCH and other models
EWMA <- matrix(NA_real_, H, N)
prev <- Y[1L,]^2
for (t in 2:nrow(Y)) {
  prev <- lambda*prev + (1-lambda)*Y[t,]^2
  if (t > T_train) EWMA[t-T_train,] <- prev
}
targets <- list(
  RV=log(pmax(Y_test^2, eps)),
  RV5_abs=log(pmax(RV5_abs, eps)),
  RV5_sq=log(pmax(RV5_sq, eps)),
  EWMA=log(pmax(EWMA, eps))
)
metric <- function(E) {
  stopifnot(all(is.finite(E)))
  c(RMSFE=mean(sqrt(colMeans(E^2))), MAFE=mean(colMeans(abs(E))))
}

# X[t,i] = sum_j W[i,j] * Y[t-1,j]^2, for t=2,...,T.
# The first observation is excluded consistently in estimation and filtering.
results <- list(); summary_rows <- list(); diagnostic_rows <- list()
for (w_name in names(weight_matrices)) {
  cat('\nNetwork:', w_name, '\n')
  W <- weight_matrices[[w_name]]
  X <- (Y[1:(nrow(Y)-1L),,drop=FALSE]^2) %*% t(W)
  OOS_var <- matrix(NA_real_, H, N, dimnames=list(NULL,colnames(Y)))
  fit_list <- vector('list', N)
  statuses <- character(N)

  for (i in seq_len(N)) {
    x_all <- matrix(X[,i], ncol=1L)
    x_train <- matrix(X[1:(T_train-1L),i], ncol=1L)
    spec_train <- ugarchspec(
      variance.model=list(model='sGARCH',garchOrder=c(1L,1L),
                          external.regressors=x_train),
      mean.model=list(armaOrder=c(0L,0L),include.mean=FALSE),
      distribution.model='norm')
    fit <- tryCatch(ugarchfit(spec=spec_train,
                             data=Y[2:T_train,i], solver='hybrid',
                             fit.control=list(scale=0)),
                    error=function(e) e)
    if (inherits(fit,'error')) {
      statuses[i] <- paste('fit error:',conditionMessage(fit)); next
    }
    if (!is.null(fit@fit$convergence) && fit@fit$convergence != 0L) {
      statuses[i] <- paste('nonconverged:',fit@fit$convergence); next
    }
    fit_list[[i]] <- fit
    # Re-filter full series using training-estimated parameters only.
    # The filtered sigma[t] is a prediction based on information through t-1.
    # n.old preserves training-based GARCH recursion initialisation.
    spec_full <- ugarchspec(
      variance.model=list(model='sGARCH',garchOrder=c(1L,1L),
                          external.regressors=x_all),
      mean.model=list(armaOrder=c(0L,0L),include.mean=FALSE),
      distribution.model='norm',
      fixed.pars=as.list(coef(fit)))
    filtered <- tryCatch(ugarchfilter(spec_full,
                                     data=Y[2:nrow(Y),i],
                                     n.old=T_train-1L),
                         error=function(e) e)
    if (inherits(filtered,'error')) {
      statuses[i] <- paste('filter error:',conditionMessage(filtered)); next
    }
    sig <- as.numeric(sigma(filtered))
    if (length(sig) != nrow(Y)-1L || any(!is.finite(sig)) || any(sig<=0)) {
      statuses[i] <- 'invalid filtered sigma'; next
    }
    # Y[2:T] corresponds to filtered positions 1:(T-1).
    OOS_var[,i] <- sig[T_train:(nrow(Y)-1L)]^2
    statuses[i] <- 'ok'
  }
  diagnostic_rows[[w_name]] <- data.frame(Weight=w_name,
    Asset=colnames(Y), Status=statuses, stringsAsFactors=FALSE)
  if (any(statuses!='ok')) {
    warning(w_name,': incomplete forecasts; see diagnostic CSV. Metrics omitted.')
    results[[w_name]] <- list(variance_forecasts=OOS_var,
                              fits=fit_list, statuses=statuses)
    next
  }
  log_fc <- log(pmax(OOS_var, eps))
  scores <- lapply(targets, function(target) metric(target-log_fc))
  summary_rows[[w_name]] <- data.frame(
    Weight=w_name,
    RMSFE_RV=unname(scores$RV['RMSFE']),MAFE_RV=unname(scores$RV['MAFE']),
    RMSFE_RV5_abs=unname(scores$RV5_abs['RMSFE']),
    MAFE_RV5_abs=unname(scores$RV5_abs['MAFE']),
    RMSFE_RV5_sq=unname(scores$RV5_sq['RMSFE']),
    MAFE_RV5_sq=unname(scores$RV5_sq['MAFE']),
    RMSFE_EWMA=unname(scores$EWMA['RMSFE']),
    MAFE_EWMA=unname(scores$EWMA['MAFE']))
  results[[w_name]] <- list(variance_forecasts=OOS_var,
    log_variance_forecasts=log_fc, fits=fit_list, scores=scores,
    statuses=statuses)
  cat('RMSFE RV:',signif(scores$RV['RMSFE'],6),'\n')
}

summary_table <- do.call(rbind,summary_rows)
if (is.null(summary_table)) summary_table <- data.frame()
diagnostic_table <- do.call(rbind,diagnostic_rows)
row.names(summary_table) <- NULL
row.names(diagnostic_table) <- NULL
print(summary_table, row.names=FALSE, digits=8)
write.csv(summary_table,'SpGARCH_X_corrected_summary.csv',row.names=FALSE)
write.csv(diagnostic_table,'SpGARCH_X_corrected_diagnostics.csv',row.names=FALSE)
saveRDS(list(summary=summary_table, diagnostics=diagnostic_table,
             results=results, targets=targets, lambda=lambda,
             training_length=T_train, test_length=H),
        'SpGARCH_X_corrected_results.rds')

print(diagnostic_table, row.names = FALSE, digits = 6)
