# Holleland STARMAGARCH: fixed-parameter genuine one-step-ahead forecasts
# Requires functions_spgarch_hol.R, STARMAGARCH compiled TMB library and RDS inputs.
suppressPackageStartupMessages({library(TMB); library(dplyr)})
source('functions_spgarch_hol.R')
if (!is.loaded('objective_function')) dyn.load(TMB::dynlib('STARMAGARCH'))
Y <- as.matrix(readRDS('residuals_VAR1.rds'))
stopifnot(is.numeric(Y), all(is.finite(Y)), nrow(Y) > 257L)
N <- ncol(Y); Tn <- nrow(Y); H <- 252L; train_n <- Tn-H
weights <- list(
  W1=readRDS('W1_eucl.rds'), W2=readRDS('W1_corr.rds'),
  W3=readRDS('W1_picollo.rds'), Wg1=readRDS('W1_eucl_GF.rds'),
  Wg2=readRDS('W2_corr_GF.rds'), Wg3=readRDS('W3_picollo_GF.rds'),
  W1_5nn=readRDS('W1_eucl5nn.rds'), W2_5nn=readRDS('W2_corr5nn.rds'),
  W3_5nn=readRDS('W1_picollo5nn.rds'), Wg=readRDS('Wmat_granger.rds')
)
# One matrix / one temporal lag in each coefficient block, as in original code.
get_W <- function(w) {
  if (is.matrix(w)) w <- array(w, c(N,N,1L))
  stopifnot(length(dim(w))==3L, identical(as.integer(dim(w)),c(N,N,1L)), all(is.finite(w)))
  w
}
# Exact C++ recursion: x[,1]=0, h[,1]=init; for t >= 2,
# x_t=(y_t-mu)-phi*W*(y_(t-1)-mu)-theta*W*x_(t-1)
# h_t=omega+alpha*W*(x_(t-1)^2)+beta*W*h_(t-1).
# Forecast for date t is h_t computed using information through t-1.
filter_model <- function(y, W, par, init) {
  m <- nrow(y); tt <- ncol(y); W <- W[,,1]
  x <- h <- matrix(NA_real_, m, tt)
  x[,1] <- 0; h[,1] <- init
  mu <- as.numeric(par$mu); phi <- as.numeric(par$phi)
  theta <- as.numeric(par$theta); omega <- as.numeric(par$omega)
  alpha <- as.numeric(par$alpha); beta <- as.numeric(par$beta)
  for (t in 2:tt) {
    x[,t] <- (y[,t]-mu) - phi*drop(W %*% (y[,t-1]-mu)) -
      theta*drop(W %*% x[,t-1])
    h[,t] <- omega + alpha*drop(W %*% (x[,t-1]^2)) +
      beta*drop(W %*% h[,t-1])
    if (any(!is.finite(h[,t])) || any(h[,t] <= 0))
      stop('Nonpositive or nonfinite conditional variance at time ',t)
  }
  list(x=x,h=h)
}
# Proxy definitions consistent with original script; five-day windows use
# observations ending on the forecast target date (evaluation only).
test_dates <- (train_n+1L):Tn
rv <- Y[test_dates,,drop=FALSE]^2
rv5a <- rv5s <- matrix(NA_real_,H,N)
for (k in seq_len(H)) {
  window <- Y[(train_n+k-4L):(train_n+k),,drop=FALSE]
  rv5a[k,] <- colMeans(abs(window))^2
  rv5s[k,] <- colMeans(window^2)
}
# EWMA is a VARIANCE proxy, not squared variance. Initialise on training sample.
lambda <- 0.94
EW <- matrix(NA_real_,H,N)
prev <- apply(Y[seq_len(train_n),,drop=FALSE],2,var)
for(k in seq_len(H)) {
  prev <- lambda*prev+(1-lambda)*Y[train_n+k,]^2
  EW[k,] <- prev
}
proxies <- lapply(list(RV=rv,RV5_abs=rv5a,RV5_sq=rv5s,EWMA=EW),
                  function(z) log(pmax(z,1e-12)))
# Choose the aggregation used by ALL models in the comparison.
# 'pooled' replicates the original Holleland evaluation script.
# 'asset_mean' matches the previously discussed DSTARCH dimensionality script.
aggregation <- 'asset_mean'
score <- function(err) {
  if (aggregation=='asset_mean') {
    c(RMSFE=mean(sqrt(colMeans(err^2))), MAFE=mean(colMeans(abs(err))))
  } else {
    c(RMSFE=sqrt(mean(err^2)), MAFE=mean(abs(err)))
  }
}
fits <- list(); summaries <- list()
for (name in names(weights)) {
  cat('Fitting',name,'\n')
  W <- get_W(weights[[name]])
  train <- t(Y[seq_len(train_n),,drop=FALSE])
  init <- pmax(apply(train,1,var),1e-6)
  parameters <- list(mu=mean(train),phi=matrix(0.7,1,1),
                     theta=matrix(0.01,1,1),omega=1,
                     alpha=matrix(0.01,1,1),beta=matrix(0.01,1,1))
  obj <- CreateLikelihood(train,W,init=init,parameters=parameters,
                          map=parameterlist2maptemplate(parameters))
  # Parameter ordering: mu, phi, theta, omega, alpha, beta.
  # Bounds from original fitSTARMAGARCH(), but constrain omega > 0
  # because it enters the conditional variance directly.
  stopifnot(length(obj$par)==6L)
  lower <- c(-1e4,-5,-5,1e-8,1e-8,1e-8)
  upper <- c( 1e4, 5, 5,1e4,1,1)
  opt <- nlminb(obj$par,obj$fn,obj$gr,lower=lower,upper=upper)
  if(opt$convergence != 0) warning(name,': nlminb convergence code ',opt$convergence)
  par <- obj$env$parList(opt$par)
  # Crucial verification: same recursion as TMB, evaluated at OPTIMISED parameters.
  rep <- obj$report(opt$par)
  check <- filter_model(train,W,par,init)
  delta <- max(abs(check$h-rep$sigma))
  delta_x <- max(abs(check$x-rep$x))
  if (!is.finite(delta) || !is.finite(delta_x) ||
      delta>1e-6 || delta_x>1e-6)
    stop(name,': R recursion does not reproduce TMB (h=',delta,', x=',delta_x,')')
  # Full-sample filtering with FIXED training estimates: h[,t] depends only
  # on information through t-1, even though x[,t] uses observed y[,t].
  filt <- filter_model(t(Y),W,par,init)
  forecast <- t(filt$h[,test_dates,drop=FALSE])
  log_forecast <- log(pmax(forecast,1e-12))
  metrics <- lapply(proxies,function(target) score(log_forecast-target))
  summaries[[name]] <- data.frame(Weight=name,
    RMSFE_1d=metrics$RV['RMSFE'],MAFE_1d=metrics$RV['MAFE'],
    RMSFE_R5_a=metrics$RV5_abs['RMSFE'],MAFE_R5_a=metrics$RV5_abs['MAFE'],
    RMSFE_R5_s=metrics$RV5_sq['RMSFE'],MAFE_R5_s=metrics$RV5_sq['MAFE'],
    RMSFE_EW=metrics$EWMA['RMSFE'],MAFE_EW=metrics$EWMA['MAFE'])
  fits[[name]] <- list(par=par,optimizer=opt,variance_forecast=forecast,
                       log_variance_forecast=log_forecast,metrics=metrics,
                       max_training_recursion_error=c(h=delta,x=delta_x))
}
summary_table <- dplyr::bind_rows(summaries)
print(summary_table,digits=8)
saveRDS(list(summary=summary_table,fits=fits,proxies=proxies,
             aggregation=aggregation,lambda=lambda),
        'STGARCH_Holleland_corrected_results.rds')
write.csv(summary_table,'STGARCH_Holleland_corrected_summary.csv',row.names=FALSE)