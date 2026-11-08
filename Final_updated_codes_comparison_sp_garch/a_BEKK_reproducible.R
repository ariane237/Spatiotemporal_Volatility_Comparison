# ==============================================================================
# SYMMETRIC DIAGONAL BEKK(1,1)
# REPRODUCIBLE FIXED-PARAMETER OUT-OF-SAMPLE EVALUATION
# ==============================================================================
#
# Required input:
#   residuals_VAR1.rds
#
# Required package:
#   BEKKs
#
# Forecast horizon:
#   H = 252 trading days
#
# Forecasting:
#   Fixed parameters
#   Sequential one-step-ahead forecasts
#
# Evaluation proxies:
#   RV, RV5_abs, RV5_sq, EWMA
#
# Primary evaluation:
#   Pooled RMSFE
#   MAFE
#
# Additional diagnostic:
#   Asset-averaged RMSFE
#
# Outputs:
#   Fitted model
#   Conditional variance forecasts
#   Forecast errors
#   AIC and BIC
#   Computation times
#   R session information
#
# ==============================================================================


# ==============================================================================
# 0. PACKAGES AND SETTINGS
# ==============================================================================

suppressPackageStartupMessages(library(BEKKs))

options(digits = 10)

H <- 252L
eps <- 1e-12
lambda <- 0.94

outdir <- "BEKK_symmetric_reproducible_results"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

start_all <- Sys.time()
total_clock <- proc.time()[["elapsed"]]

cat("\n============================================================\n")
cat("       SYMMETRIC DIAGONAL BEKK(1,1)\n")
cat("============================================================\n\n")

cat("Start:", format(start_all), "\n")
cat("BEKKs version:", as.character(packageVersion("BEKKs")), "\n\n")


# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

U <- as.matrix(readRDS("residuals_VAR1.rds"))

storage.mode(U) <- "double"

stopifnot(
  is.numeric(U),
  all(is.finite(U)),
  nrow(U) > H + 4L,
  ncol(U) >= 2L
)

if (is.null(colnames(U))) {
  colnames(U) <- paste0("Asset_", seq_len(ncol(U)))
}

N <- ncol(U)
Ttot <- nrow(U)
Ttrain <- Ttot - H

Utrain <- U[seq_len(Ttrain), , drop = FALSE]

Utest <- U[
  (Ttrain + 1L):Ttot,
  ,
  drop = FALSE
]

cat("DATA INFORMATION\n")
cat("-----------------------------\n")
cat("Assets:               ", N, "\n")
cat("Total observations:   ", Ttot, "\n")
cat("Training observations:", Ttrain, "\n")
cat("Forecast observations:", H, "\n")
cat("EWMA lambda:          ", lambda, "\n\n")

# Fingerprint of the input file
input_md5 <- unname(tools::md5sum("residuals_VAR1.rds"))

cat("Input RDS MD5:", input_md5, "\n\n")


# ==============================================================================
# 2. MODEL SPECIFICATION
# ==============================================================================

spec <- bekk_spec(
  model = list(
    type = "dbekk",
    asymmetric = FALSE
  )
)


# ==============================================================================
# 3. MODEL ESTIMATION
# ==============================================================================

cat("============================================================\n")
cat("                 MODEL ESTIMATION\n")
cat("============================================================\n")

cat("Fitting symmetric diagonal BEKK...\n")

estimation_start <- proc.time()[["elapsed"]]

fit <- bekk_fit(
  spec,
  data = Utrain,
  QML_t_ratios = FALSE,
  max_iter = 200,
  crit = 1e-8
)

estimation_sec <- proc.time()[["elapsed"]] - estimation_start

# Save immediately to preserve expensive estimation results.
saveRDS(
  fit,
  file.path(outdir, "BEKK_symmetric_fitted_model.rds")
)

cat("\nEstimation completed.\n")
cat("Estimation time:", round(estimation_sec, 3), "seconds\n")
cat("Estimation time:", round(estimation_sec / 60, 3), "minutes\n\n")

if (!isTRUE(fit$BEKK_valid)) {
  warning(
    "BEKK_valid is not TRUE. The fitted model requires investigation."
  )
}


# ==============================================================================
# 4. EXTRACT ESTIMATED MATRICES
# ==============================================================================

C0 <- as.matrix(fit$C0)
A <- as.matrix(fit$A)
G <- as.matrix(fit$G)

stopifnot(
  identical(dim(C0), c(N, N)),
  identical(dim(A), c(N, N)),
  identical(dim(G), c(N, N)),
  all(is.finite(C0)),
  all(is.finite(A)),
  all(is.finite(G))
)

# Verify diagonal restrictions on A and G.

offdiag <- function(M) {
  diag(M) <- 0
  max(abs(M))
}

stopifnot(
  offdiag(A) < 1e-8,
  offdiag(G) < 1e-8
)


# ==============================================================================
# 5. EXTRACT LAST IN-SAMPLE CONDITIONAL COVARIANCE
# ==============================================================================

Hhist <- fit$H_t

if (
  is.matrix(Hhist) &&
  ncol(Hhist) == N * N &&
  nrow(Hhist) == Ttrain
) {
  
  Hprev <- matrix(
    as.numeric(Hhist[Ttrain, ]),
    nrow = N,
    ncol = N
  )
  
} else if (
  is.list(Hhist) &&
  length(Hhist) == Ttrain
) {
  
  Hprev <- as.matrix(Hhist[[Ttrain]])
  
} else if (
  length(dim(Hhist)) == 3L &&
  identical(as.integer(dim(Hhist)), c(N, N, Ttrain))
) {
  
  Hprev <- Hhist[, , Ttrain]
  
} else {
  
  stop(
    "Unexpected fit$H_t structure. ",
    "Inspect the saved fitted model."
  )
}

stopifnot(
  identical(dim(Hprev), c(N, N)),
  all(is.finite(Hprev))
)

Hprev <- (Hprev + t(Hprev)) / 2


# ==============================================================================
# 6. BEKK COVARIANCE RECURSION
# ==============================================================================

# H_(t+1) = C0 C0' + A' u_t u_t' A + G' H_t G

step_cov <- function(u_last, H_last) {
  
  u <- matrix(
    as.numeric(u_last),
    ncol = 1L
  )
  
  result <-
    tcrossprod(C0) +
    crossprod(A, tcrossprod(u) %*% A) +
    crossprod(G, H_last %*% G)
  
  (result + t(result)) / 2
}


# ==============================================================================
# 7. VERIFY FIRST FORECAST AGAINST BEKKs
# ==============================================================================

cat("\n============================================================\n")
cat("              FIRST FORECAST VERIFICATION\n")
cat("============================================================\n")

first_manual <- step_cov(
  U[Ttrain, ],
  Hprev
)

pf <- predict(
  fit,
  n.ahead = 1L
)

package_h <- pf$H_t_forecast

if (is.null(package_h)) {
  stop("BEKKs predict() returned no H_t_forecast.")
}

if (
  is.matrix(package_h) &&
  ncol(package_h) == N * N
) {
  
  first_package <- matrix(
    as.numeric(package_h[nrow(package_h), ]),
    N,
    N
  )
  
} else if (
  is.list(package_h) &&
  length(package_h) >= 1L
) {
  
  first_package <- as.matrix(
    package_h[[length(package_h)]]
  )
  
} else if (
  length(dim(package_h)) == 3L &&
  dim(package_h)[1] == N &&
  dim(package_h)[2] == N
) {
  
  first_package <- package_h[
    ,
    ,
    dim(package_h)[3]
  ]
  
} else {
  
  stop(
    "Unrecognised BEKKs forecast structure. ",
    "Inspect predict(fit)."
  )
}

first_package <- (first_package + t(first_package)) / 2

check_abs <- max(
  abs(first_manual - first_package)
)

check_rel <- check_abs /
  max(1, max(abs(first_package)))

cat("Maximum absolute difference:", check_abs, "\n")
cat("Relative difference:", check_rel, "\n")

if (
  !isTRUE(
    all.equal(
      first_manual,
      first_package,
      tolerance = 1e-7,
      check.attributes = FALSE
    )
  )
) {
  
  stop(
    "Manual covariance recursion does not match ",
    "the BEKKs first-step forecast."
  )
}

cat("First forecast verification passed.\n")


# ==============================================================================
# 8. SEQUENTIAL ONE-STEP-AHEAD FORECASTS
# ==============================================================================

cat("\n============================================================\n")
cat("                    FORECASTING\n")
cat("============================================================\n")

forecast_start <- proc.time()[["elapsed"]]

forecast_var <- matrix(
  NA_real_,
  nrow = H,
  ncol = N,
  dimnames = dimnames(Utest)
)

for (j in seq_len(H)) {
  
  # First forecast uses the last training observation.
  # Subsequent forecasts use the most recent observed residual.
  
  last_idx <- Ttrain + j - 1L
  
  Hnext <- step_cov(
    U[last_idx, ],
    Hprev
  )
  
  d <- diag(Hnext)
  
  if (
    any(!is.finite(d)) ||
    any(d <= 0)
  ) {
    stop(
      "Invalid forecast variance at OOS step ",
      j
    )
  }
  
  forecast_var[j, ] <- d
  
  # Update filtered covariance.
  # The next forecast will incorporate the next observed residual.
  Hprev <- Hnext
}

forecast_sec <- proc.time()[["elapsed"]] - forecast_start

cat("Forecasting completed.\n")
cat("Forecasting time:", round(forecast_sec, 3), "seconds\n\n")


# ==============================================================================
# 9. CONSTRUCT EVALUATION PROXIES
# ==============================================================================

cat("Constructing evaluation proxies...\n")

evaluation_start <- proc.time()[["elapsed"]]

RV <- Utest^2

RV5_abs <- matrix(
  NA_real_,
  H,
  N,
  dimnames = dimnames(Utest)
)

RV5_sq <- matrix(
  NA_real_,
  H,
  N,
  dimnames = dimnames(Utest)
)

for (j in seq_len(H)) {
  
  t <- Ttrain + j
  
  w <- U[
    (t - 4L):t,
    ,
    drop = FALSE
  ]
  
  RV5_abs[j, ] <- colMeans(abs(w))^2
  
  RV5_sq[j, ] <- colMeans(w^2)
}


# EWMA proxy

EWMA <- matrix(
  NA_real_,
  H,
  N,
  dimnames = dimnames(Utest)
)

EWMA_prev <- U[1L, ]^2

for (t in 2:Ttot) {
  
  EWMA_prev <-
    lambda * EWMA_prev +
    (1 - lambda) * U[t, ]^2
  
  if (t > Ttrain) {
    EWMA[t - Ttrain, ] <- EWMA_prev
  }
}


# Log evaluation targets

targets <- list(
  RV = log(pmax(RV, eps)),
  RV5_abs = log(pmax(RV5_abs, eps)),
  RV5_sq = log(pmax(RV5_sq, eps)),
  EWMA = log(pmax(EWMA, eps))
)

log_forecast <- log(pmax(forecast_var, eps))


# ==============================================================================
# 10. FORECAST EVALUATION
# ==============================================================================

errors <- lapply(
  targets,
  function(target) {
    target - log_forecast
  }
)


# Primary pooled RMSFE

pooled_rmsfe <- function(E) {
  sqrt(mean(E^2))
}


# Additional diagnostic: asset-averaged RMSFE

asset_mean_rmsfe <- function(E) {
  mean(sqrt(colMeans(E^2)))
}


# MAFE

mafe <- function(E) {
  mean(abs(E))
}


error_table <- data.frame(
  Proxy = names(errors),
  
  RMSFE_pooled = vapply(
    errors,
    pooled_rmsfe,
    numeric(1)
  ),
  
  MAFE = vapply(
    errors,
    mafe,
    numeric(1)
  ),
  
  RMSFE_asset_mean_diagnostic = vapply(
    errors,
    asset_mean_rmsfe,
    numeric(1)
  ),
  
  row.names = NULL
)

evaluation_sec <- proc.time()[["elapsed"]] - evaluation_start


# ==============================================================================
# 11. INFORMATION CRITERIA
# ==============================================================================

llh <- as.numeric(fit$log_likelihood)
aic_pkg <- as.numeric(fit$AIC)
bic_pkg <- as.numeric(fit$BIC)

if (
  length(llh) != 1L ||
  length(aic_pkg) != 1L ||
  length(bic_pkg) != 1L ||
  any(!is.finite(c(llh, aic_pkg, bic_pkg)))
) {
  
  stop(
    "Missing or invalid BEKK log-likelihood, AIC or BIC."
  )
}


# Theoretical number of parameters for symmetric diagonal BEKK:
# N(N+1)/2 free elements in C0, plus N in A and N in G.

k_theory <- N * (N + 1) / 2 + 2 * N


cat("\n============================================================\n")
cat("                 INFORMATION CRITERIA\n")
cat("============================================================\n")

cat("Log-likelihood:", llh, "\n")
cat("AIC (BEKKs):", aic_pkg, "\n")
cat("BIC (BEKKs):", bic_pkg, "\n")
cat("Theoretical parameter count:", k_theory, "\n")
cat("BEKK_valid:", fit$BEKK_valid, "\n")
cat("Iterations:", fit$iter, "\n")


# ==============================================================================
# 12. PRINT FORECAST ERRORS
# ==============================================================================

cat("\n============================================================\n")
cat("                  FORECAST ERRORS\n")
cat("============================================================\n")

print(
  error_table,
  digits = 9,
  row.names = FALSE
)


# ==============================================================================
# 13. COMPUTATION TIMES
# ==============================================================================

total_sec <- proc.time()[["elapsed"]] - total_clock

timing <- data.frame(
  Stage = c(
    "Estimation",
    "Forecasting",
    "Evaluation and proxies",
    "Total"
  ),
  
  Seconds = c(
    estimation_sec,
    forecast_sec,
    evaluation_sec,
    total_sec
  )
)

timing$Minutes <- timing$Seconds / 60

cat("\n============================================================\n")
cat("                 COMPUTATION TIMES\n")
cat("============================================================\n")

print(
  timing,
  digits = 7,
  row.names = FALSE
)


# ==============================================================================
# 14. FINAL SUMMARY TABLE
# ==============================================================================

summary_table <- data.frame(
  Model = "BEKK_diagonal_symmetric",
  
  N = N,
  T_train = Ttrain,
  H = H,
  
  LogLikelihood = llh,
  AIC = aic_pkg,
  BIC = bic_pkg,
  
  N_parameters_theoretical = k_theory,
  
  BEKK_valid = isTRUE(fit$BEKK_valid),
  Iterations = as.integer(fit$iter),
  
  RMSFE_RV = error_table$RMSFE_pooled[1],
  MAFE_RV = error_table$MAFE[1],
  
  RMSFE_RV5_abs = error_table$RMSFE_pooled[2],
  MAFE_RV5_abs = error_table$MAFE[2],
  
  RMSFE_RV5_sq = error_table$RMSFE_pooled[3],
  MAFE_RV5_sq = error_table$MAFE[3],
  
  RMSFE_EWMA = error_table$RMSFE_pooled[4],
  MAFE_EWMA = error_table$MAFE[4],
  
  Estimation_seconds = estimation_sec,
  Forecast_seconds = forecast_sec,
  Evaluation_seconds = evaluation_sec,
  Total_seconds = total_sec
)


cat("\n============================================================\n")
cat("                    FINAL RESULTS\n")
cat("============================================================\n")

print(
  summary_table,
  digits = 9,
  row.names = FALSE
)


# ==============================================================================
# 15. SAVE RESULTS
# ==============================================================================

saveRDS(
  list(
    summary = summary_table,
    errors = errors,
    error_table = error_table,
    timing = timing,
    
    forecast_variance = forecast_var,
    log_variance_forecast = log_forecast,
    
    log_targets = targets,
    
    first_forecast_check = list(
      max_abs = check_abs,
      relative = check_rel
    ),
    
    input_md5 = input_md5,
    
    training_length = Ttrain,
    test_length = H,
    
    lambda = lambda,
    eps = eps,
    
    fit = fit
  ),
  
  file.path(
    outdir,
    "full_results.rds"
  )
)


write.csv(
  summary_table,
  file.path(outdir, "summary.csv"),
  row.names = FALSE
)

write.csv(
  error_table,
  file.path(outdir, "forecast_errors.csv"),
  row.names = FALSE
)

write.csv(
  timing,
  file.path(outdir, "computation_times.csv"),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(outdir, "session_info.txt")
)


cat("\n============================================================\n")
cat("                 COMPUTATION COMPLETED\n")
cat("============================================================\n")

cat("Results saved in:", outdir, "\n")
cat("Total computation time:", round(total_sec / 60, 3), "minutes\n")
cat("End time:", format(Sys.time()), "\n")