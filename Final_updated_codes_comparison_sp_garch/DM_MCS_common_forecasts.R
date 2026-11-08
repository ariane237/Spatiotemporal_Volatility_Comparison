# ============================================================================
# DIEBOLD-MARIANO AND MODEL CONFIDENCE SET: COMMON FORECAST EVALUATION
# Requires: forecast_audit_output/audit_all_results.rds
# No model re-estimation. No changes to saved forecasts.
# ============================================================================
options(stringsAsFactors = FALSE)

ROOT <- getwd()
AUDIT_FILE <- file.path(ROOT, "forecast_audit_output", "audit_all_results.rds")
OUT <- file.path(ROOT, "forecast_audit_output", "significance_tests")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(AUDIT_FILE)) stop("Run audit_all_volatility_forecasts.R first; missing: ", AUDIT_FILE)

# Analysis choices
PROXY <- "RV"                     # "RV", "RV5_abs", "RV5_sq", "EWMA"
LOSS <- "squared"                 # "squared" or "absolute" log forecast error
NETWORK_PROTOCOL <- "best_ex_post" # "best_ex_post" or "all_networks"
ALPHA <- 0.10                     # 90% Model Confidence Set
B <- 5000L                        # bootstrap replications for MCS
BLOCK_LENGTH <- 10L               # block length for MCS (sensitivity: 5, 20)
DM_HAC_LAG <- 0L                  # 1-step forecasts: lag 0; sensitivity: 5, 10
SEED <- 20261020L

# The MCS package is the published MCS implementation (Hansen, Lunde, Nason).
if (!requireNamespace("MCS", quietly = TRUE)) {
  stop("Please install the MCS package first: install.packages('MCS')")
}

obj <- readRDS(AUDIT_FILE)
if (is.null(obj$errors) || is.null(obj$scores)) stop("Unexpected audit file structure")
score_table <- obj$scores
keycol <- paste0(PROXY, "_RMSFE_asset_mean")
if (!(keycol %in% names(score_table))) stop("Proxy not found in score table: ", PROXY)

# IMPORTANT: network selection on OOS errors is retrospective.
# Treat the best_ex_post results as exploratory, not an unbiased selection test.
if (NETWORK_PROTOCOL == "best_ex_post") {
  pieces <- split(score_table, score_table$Model)
  chosen <- do.call(rbind, lapply(pieces, function(d) {
    d[which.min(d[[keycol]]), , drop = FALSE]
  }))
} else if (NETWORK_PROTOCOL == "all_networks") {
  chosen <- score_table
} else stop("Unknown NETWORK_PROTOCOL")
rownames(chosen) <- NULL
chosen <- chosen[order(chosen$Model, chosen$Network), , drop = FALSE]

# Build one daily cross-sectional loss per forecast specification.
# Every model uses the SAME date, assets, proxy, and loss transformation.
series <- list()
for (r in seq_len(nrow(chosen))) {
  model <- chosen$Model[r]
  network <- chosen$Network[r]
  key <- paste(model, network, PROXY, sep = "::")
  e <- obj$errors[[key]]
  if (is.null(e)) stop("Missing stored forecast errors: ", key)
  e <- as.matrix(e)
  if (!is.numeric(e) || any(!is.finite(e))) stop("Invalid errors: ", key)
  if (LOSS == "squared") {
    daily <- rowMeans(e^2)
  } else if (LOSS == "absolute") {
    daily <- rowMeans(abs(e))
  } else stop("LOSS must be squared or absolute")
  nm <- if (NETWORK_PROTOCOL == "all_networks") paste(model, network, sep = "__") else model
  if (nm %in% names(series)) stop("Duplicate specification label: ", nm)
  series[[nm]] <- daily
}
Tn <- unique(lengths(series))
if (length(Tn) != 1L || Tn < 30L) stop("Loss series have unequal or insufficient lengths")
L <- do.call(cbind, series)
if (any(!is.finite(L))) stop("Loss matrix contains non-finite values")
if (ncol(L) < 2L) stop("At least two models required")

cat("Observations:", nrow(L), " | Models/specifications:", ncol(L), "\n")
cat("Proxy:", PROXY, " | Loss:", LOSS, " | Network protocol:", NETWORK_PROTOCOL, "\n")
cat("CAUTION: Best-network selection based on OOS loss is exploratory.\n")
write.csv(chosen, file.path(OUT, "selected_specifications.csv"), row.names = FALSE)
write.csv(data.frame(DateIndex = seq_len(nrow(L)), L, check.names = FALSE),
          file.path(OUT, "daily_loss_matrix.csv"), row.names = FALSE)

# ----------------------------------------------------------------------------
# 1. Diebold-Mariano test (two-sided, asymptotic normal, Bartlett HAC)
# d_t = loss(Model_1,t) - loss(Model_2,t)
# Negative mean d => Model_1 has LOWER mean loss.
# The standard DM statistic is valid for 1-step forecasts under appropriate
# weak dependence assumptions. No Harvey small-sample adjustment applied.
# ----------------------------------------------------------------------------
hac_long_run_var <- function(d, lag) {
  n <- length(d)
  z <- d - mean(d)
  v <- mean(z^2)
  if (lag > 0L) for (k in seq_len(lag)) {
    gamma <- sum(z[(k + 1L):n] * z[seq_len(n-k)]) / n
    v <- v + 2 * (1 - k / (lag + 1)) * gamma
  }
  v
}

names_m <- colnames(L)
pairs <- combn(seq_along(names_m), 2L)
dm_rows <- vector("list", ncol(pairs))
pvals <- matrix(NA_real_, ncol(L), ncol(L), dimnames = list(names_m, names_m))
statmat <- pvals
mean_diff_mat <- pvals
diag(pvals) <- 1
diag(statmat) <- 0
diag(mean_diff_mat) <- 0
for (j in seq_len(ncol(pairs))) {
  a <- pairs[1L,j]; b <- pairs[2L,j]
  d <- L[,a] - L[,b]
  md <- mean(d)
  lrv <- hac_long_run_var(d, DM_HAC_LAG)
  if (is.finite(lrv) && lrv > 0) {
    se <- sqrt(lrv / length(d))
    stat <- md / se
    pv <- 2 * pnorm(-abs(stat))
  } else {
    stat <- NA_real_; pv <- NA_real_
    warning("Zero/invalid long-run variance for ", names_m[a], " vs ", names_m[b])
  }
  dm_rows[[j]] <- data.frame(
    Model_1 = names_m[a], Model_2 = names_m[b],
    Mean_loss_1 = mean(L[,a]), Mean_loss_2 = mean(L[,b]),
    Mean_loss_difference_1_minus_2 = md,
    DM_statistic = stat, p_value_two_sided = pv,
    Lower_mean_loss = if (md < 0) names_m[a] else if (md > 0) names_m[b] else "tie",
    Significant_at_5pct = is.finite(pv) && pv < 0.05,
    Significant_at_10pct = is.finite(pv) && pv < 0.10
  )
  pvals[a,b] <- pvals[b,a] <- pv
  statmat[a,b] <- stat; statmat[b,a] <- -stat
  mean_diff_mat[a,b] <- md; mean_diff_mat[b,a] <- -md
}
dm <- do.call(rbind, dm_rows)
write.csv(dm, file.path(OUT, "DM_pairwise_results.csv"), row.names = FALSE)
write.csv(pvals, file.path(OUT, "DM_pvalues_matrix.csv"))
write.csv(statmat, file.path(OUT, "DM_statistics_matrix.csv"))
write.csv(mean_diff_mat, file.path(OUT, "DM_mean_loss_differences_matrix.csv"))

# For interpretation: mean squared loss determines expected loss, NOT the
# average of asset-specific RMSE values used in the descriptive table.
ranked <- data.frame(Model = names_m, Mean_daily_loss = colMeans(L))
ranked <- ranked[order(ranked$Mean_daily_loss), , drop = FALSE]
ranked$Rank <- seq_len(nrow(ranked))
ranked <- ranked[,c("Rank", "Model", "Mean_daily_loss")]
write.csv(ranked, file.path(OUT, "mean_loss_ranking.csv"), row.names = FALSE)

cat("\n=== MEAN DAILY LOSS RANKING ===\n")
print(ranked, row.names = FALSE, digits = 8)
cat("\n=== DIEBOLD-MARIANO PAIRWISE RESULTS ===\n")
print(dm[,c("Model_1", "Model_2", "Mean_loss_difference_1_minus_2", "DM_statistic", "p_value_two_sided")],
      row.names = FALSE, digits = 5)

# ----------------------------------------------------------------------------
# 2. Model Confidence Set (Hansen, Lunde & Nason)
# MCS::MCSprocedure uses block bootstrap and sequential elimination.
# Statistic = Tmax (maximum standardized relative loss).
# Set seed to make results reproducible.
# ----------------------------------------------------------------------------
set.seed(SEED)
cat("\nRunning MCS with B=", B, ", k=", BLOCK_LENGTH, ", alpha=", ALPHA, "\n", sep = "")
mcs <- tryCatch(
  MCS::MCSprocedure(Loss = L, alpha = ALPHA, B = B,
                    statistic = "Tmax", k = BLOCK_LENGTH),
  error = function(e) stop("MCSprocedure failed: ", conditionMessage(e))
)
saveRDS(mcs, file.path(OUT, "MCS_result_object.rds"))
capture.output(show(mcs), file = file.path(OUT, "MCS_printed_results.txt"))
cat("\n=== MODEL CONFIDENCE SET ===\n")
show(mcs)

# Save full reproducibility configuration.
config <- list(proxy = PROXY, loss = LOSS, network_protocol = NETWORK_PROTOCOL,
               alpha = ALPHA, bootstrap_replications = B,
               block_length = BLOCK_LENGTH, dm_hac_lag = DM_HAC_LAG,
               seed = SEED, n_forecast_dates = nrow(L),
               n_assets = ncol(obj$errors[[paste(chosen$Model[1], chosen$Network[1], PROXY, sep = "::")]]),
               session_info = capture.output(sessionInfo()))
saveRDS(config, file.path(OUT, "analysis_configuration.rds"))
cat("\nSaved all results in: ", OUT, "\n", sep = "")
cat("IMPORTANT: Ex-post best network selection can inflate significance.\n")
cat("For a robustness analysis, repeat with NETWORK_PROTOCOL='all_networks'.\n")
cat("For inference, avoid claiming statistically validated network selection.\n")

