# === Libraries ===
# This will install the package only if it's not already available.
if(!require("BEKKs")) install.packages("BEKKs")
if(!require("doParallel")) install.packages("doParallel")
if(!require("foreach")) install.packages("foreach")
if(!require("parallel")) install.packages("parallel")
library(BEKKs)
library(doParallel)
library(foreach)
library(parallel)


# --- Data ---
residuals <- as.matrix(readRDS("residuals_VAR1.rds"))
Ttot <- nrow(residuals); N <- ncol(residuals)
train.l <- Ttot - 252; out.l <- 252

# --- Realized-volatility proxies ---
RV_proxy <- residuals[(train.l + 1):(train.l + out.l), ]^2
RV5_mean_sq <- t(sapply(1:out.l, function(i) colMeans(residuals[(train.l + i - 4):(train.l + i), ]^2)))
RV5_mean_abs <- t(sapply(1:out.l, function(i) colMeans(abs(residuals[(train.l + i - 4):(train.l + i), ]))))
ewma_var <- function(x, lambda = 0.94) {
  out <- numeric(length(x)); out[1] <- x[1]^2
  for (t in 2:length(x)) out[t] <- lambda * out[t - 1] + (1 - lambda) * x[t]^2
  out
}
EWMA_var <- apply(residuals[(train.l + 1):(train.l + out.l), ], 2, ewma_var)

# --- Log-transformed versions ---
eps <- 1e-12
log_RV <- log(pmax(RV_proxy, eps))
log_RV5_mean_sq <- log(pmax(RV5_mean_sq, eps))
log_RV5_mean_abs <- log(pmax(RV5_mean_abs, eps))
log_EWMA_var <- log(pmax(EWMA_var, eps))

# --- Metrics functions ---
rmsfe_fn <- function(err) mean(apply(err, 2, function(x) sqrt(mean(x^2, na.rm=TRUE))), na.rm=TRUE)
mafe_fn  <- function(err) mean(apply(err, 2, function(x) mean(abs(x), na.rm=TRUE)),    na.rm=TRUE)

# --- BEKK model spec ---
bekk_spec_obj <- bekk_spec(model = list(type = "dbekk", asymmetric = FALSE))
extract_var1 <- function(pred) {
  sds <- as.numeric(pred$volatility_forecast[grepl("Conditional standard deviation", names(pred$volatility_forecast))])
  sds^2
}

# --- Single estimation using the training window ---
cat("Estimating BEKK model once using the training window...\n")
fit <- bekk_fit(bekk_spec_obj, data = residuals[1:train.l, ],
                QML_t_ratios = FALSE, max_iter = 200, crit = 1e-8)

# --- Forecasts with fixed model over rolling test set ---
cat("Performing parallelized one-step-ahead rolling forecast...\n")
n.cores <- min(124, detectCores())
cl <- makeCluster(n.cores, type = "FORK")
registerDoParallel(cl)

OOSf_bekk <- foreach(j = 1:out.l, .combine = "rbind", .packages = "BEKKs",
                     .export = c("fit", "extract_var1")) %dopar% {
                       # Prepare data: rolling residuals up to time t-1
                       test_data <- residuals[j:(j + train.l - 1), ]
                       test_data <- scale(test_data, center = TRUE, scale = FALSE)  # de-mean
                       
                       # Use fitted model to forecast (re-centering may affect)
                       pred <- tryCatch(predict(fit, n.ahead = 1, newdata = test_data), error = function(e) NULL)
                       if (is.null(pred)) rep(NA_real_, N) else extract_var1(pred)
                     }

stopCluster(cl)

# --- In-sample Fit ---
loglik_bekk <- fit$log_likelihood
bic_bekk <- fit$BIC

# --- Forecast errors ---
log_OOS <- log(pmax(OOSf_bekk, eps))
metrics_list <- list(
  LogLikelihood = loglik_bekk,
  BIC           = bic_bekk,
  RMSFE_RV          = rmsfe_fn(log_RV - log_OOS),          MAFE_RV          = mafe_fn(log_RV - log_OOS),
  RMSFE_RV5_mean_sq = rmsfe_fn(log_RV5_mean_sq - log_OOS), MAFE_RV5_mean_sq = mafe_fn(log_RV5_mean_sq - log_OOS),
  RMSFE_RV5_mean_abs= rmsfe_fn(log_RV5_mean_abs - log_OOS),MAFE_RV5_mean_abs= mafe_fn(log_RV5_mean_abs - log_OOS),
  RMSFE_EWMA        = rmsfe_fn(log_EWMA_var - log_OOS),    MAFE_EWMA        = mafe_fn(log_EWMA_var - log_OOS)
)
print(metrics_list)

save.image("symmetric_bekk_model_rolling_summary_noreestimate.RData")
