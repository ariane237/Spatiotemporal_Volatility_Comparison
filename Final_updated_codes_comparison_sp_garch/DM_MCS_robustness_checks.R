# ============================================================================
# ROBUSTNESS CHECKS: DM HAC LAGS AND MCS BOOTSTRAP BLOCK LENGTHS
# Reads forecasts from forecast_audit_output/audit_all_results.rds
# Does NOT re-estimate models or modify saved forecasts.
# ============================================================================
options(stringsAsFactors = FALSE)

ROOT <- getwd()
INPUT <- file.path(ROOT, 'forecast_audit_output', 'audit_all_results.rds')
OUT <- file.path(ROOT, 'forecast_audit_output', 'significance_robustness')
if (!file.exists(INPUT)) stop('Missing file: ', INPUT)
if (!requireNamespace('MCS', quietly = TRUE)) {
  stop("Install package first: install.packages('MCS')")
}
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ----- Settings -------------------------------------------------------------
PROXY <- 'RV'                  # RV, RV5_abs, RV5_sq, EWMA
LOSS <- 'squared'              # squared or absolute
PROTOCOLS <- c('best_ex_post', 'all_networks')
DM_LAGS <- c(0L, 5L, 10L)
MCS_BLOCKS <- c(5L, 10L, 20L)
MCS_B <- 5000L
MCS_ALPHA <- 0.10
MCS_STATISTIC <- 'Tmax'
SEED <- 20261020L

obj <- readRDS(INPUT)
stopifnot(is.list(obj), !is.null(obj$errors), !is.null(obj$scores))
score_table <- obj$scores
score_column <- paste0(PROXY, '_RMSFE_asset_mean')
if (!score_column %in% names(score_table)) stop('Missing score: ', score_column)

make_loss_matrix <- function(protocol) {
  if (protocol == 'best_ex_post') {
    parts <- split(score_table, score_table$Model)
    chosen <- do.call(rbind, lapply(parts, function(z) {
      z[which.min(z[[score_column]]), , drop = FALSE]
    }))
  } else if (protocol == 'all_networks') {
    chosen <- score_table
  } else stop('Unknown protocol: ', protocol)
  rownames(chosen) <- NULL
  chosen <- chosen[order(chosen$Model, chosen$Network), , drop = FALSE]

  series <- vector('list', nrow(chosen))
  labels <- character(nrow(chosen))
  for (j in seq_len(nrow(chosen))) {
    model <- chosen$Model[j]
    network <- chosen$Network[j]
    key <- paste(model, network, PROXY, sep = '::')
    e <- obj$errors[[key]]
    if (is.null(e)) stop('Missing error matrix: ', key)
    e <- as.matrix(e)
    if (!is.numeric(e) || any(!is.finite(e))) stop('Invalid errors: ', key)
    if (LOSS == 'squared') {
      series[[j]] <- rowMeans(e^2)
    } else if (LOSS == 'absolute') {
      series[[j]] <- rowMeans(abs(e))
    } else stop('Unsupported LOSS')
    labels[j] <- if (protocol == 'all_networks') {
      paste(model, network, sep = '__')
    } else model
  }
  if (anyDuplicated(labels)) stop('Duplicate model labels')
  if (length(unique(lengths(series))) != 1L) stop('Different forecast horizons')
  L <- do.call(cbind, series)
  colnames(L) <- labels
  if (nrow(L) < 30L || ncol(L) < 2L || any(!is.finite(L))) {
    stop('Invalid loss matrix')
  }
  list(L = L, chosen = chosen)
}

# Bartlett-kernel HAC long-run variance of the loss differential.
hac_lrv <- function(d, lag) {
  n <- length(d)
  if (lag < 0L || lag >= n) stop('Invalid HAC lag')
  z <- d - mean(d)
  v <- sum(z^2) / n
  if (lag > 0L) {
    for (k in seq_len(lag)) {
      gamma_k <- sum(z[(k + 1L):n] * z[1L:(n-k)]) / n
      v <- v + 2 * (1 - k / (lag + 1)) * gamma_k
    }
  }
  v
}

run_dm <- function(L, lag) {
  pairs <- combn(seq_len(ncol(L)), 2L)
  rows <- lapply(seq_len(ncol(pairs)), function(j) {
    a <- pairs[1L, j]; b <- pairs[2L, j]
    d <- L[, a] - L[, b]
    md <- mean(d)
    lrv <- hac_lrv(d, lag)
    stat <- if (is.finite(lrv) && lrv > 0) md / sqrt(lrv / length(d)) else NA_real_
    pv <- if (is.finite(stat)) 2 * pnorm(-abs(stat)) else NA_real_
    data.frame(
      Model_1 = colnames(L)[a], Model_2 = colnames(L)[b],
      HAC_lag = lag, Mean_loss_1 = mean(L[, a]),
      Mean_loss_2 = mean(L[, b]), Mean_loss_difference_1_minus_2 = md,
      DM_statistic = stat, p_value_two_sided = pv,
      Lower_mean_loss = if (md < 0) colnames(L)[a] else if (md > 0) colnames(L)[b] else 'tie',
      Significant_at_5pct = !is.na(pv) && pv < 0.05,
      Significant_at_10pct = !is.na(pv) && pv < 0.10
    )
  })
  do.call(rbind, rows)
}

all_dm <- list()
run_status <- list()
for (protocol in PROTOCOLS) {
  cat('\n================ ', protocol, ' ================\n', sep = '')
  folder <- file.path(OUT, protocol)
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  dat <- make_loss_matrix(protocol)
  L <- dat$L
  write.csv(dat$chosen, file.path(folder, 'selected_specifications.csv'), row.names = FALSE)
  write.csv(data.frame(DateIndex = seq_len(nrow(L)), L, check.names = FALSE),
            file.path(folder, 'daily_loss_matrix.csv'), row.names = FALSE)
  rank <- data.frame(Specification = colnames(L), Mean_daily_loss = colMeans(L))
  rank <- rank[order(rank$Mean_daily_loss), , drop = FALSE]
  write.csv(rank, file.path(folder, 'mean_loss_ranking.csv'), row.names = FALSE)
  cat('Forecast dates:', nrow(L), '  Specifications:', ncol(L), '\n')

  dm <- do.call(rbind, lapply(DM_LAGS, function(lag) run_dm(L, lag)))
  dm$Protocol <- protocol
  write.csv(dm, file.path(folder, 'DM_all_HAC_lags.csv'), row.names = FALSE)
  all_dm[[protocol]] <- dm

  # Additional focused table: DSTARCH vs every other candidate.
  dst <- dm[dm$Model_1 == 'DSTARCH' | dm$Model_2 == 'DSTARCH' |
              grepl('^DSTARCH__', dm$Model_1) | grepl('^DSTARCH__', dm$Model_2), , drop = FALSE]
  write.csv(dst, file.path(folder, 'DM_comparisons_involving_DSTARCH.csv'), row.names = FALSE)

  for (k in MCS_BLOCKS) {
    cat('Running MCS: protocol=', protocol, ', block=', k,
        ', B=', MCS_B, '\n', sep = '')
    # Reproducible but different random streams for each setting.
    set.seed(SEED + k + if (protocol == 'all_networks') 1000L else 0L)
    logfile <- file.path(folder, paste0('MCS_block_', k, '.txt'))
    result <- tryCatch({
      transcript <- capture.output({
        mcs <- MCS::MCSprocedure(
          Loss = L, alpha = MCS_ALPHA, B = MCS_B,
          statistic = MCS_STATISTIC, k = k
        )
        show(mcs)
      })
      writeLines(transcript, logfile)
      saveRDS(mcs, file.path(folder, paste0('MCS_block_', k, '.rds')))
      data.frame(Protocol = protocol, Block_length = k, Status = 'ok',
                 Details = '', stringsAsFactors = FALSE)
    }, error = function(e) {
      writeLines(paste('MCS FAILED:', conditionMessage(e)), logfile)
      data.frame(Protocol = protocol, Block_length = k, Status = 'failed',
                 Details = conditionMessage(e), stringsAsFactors = FALSE)
    })
    run_status[[paste(protocol, k)]] <- result
  }
}

write.csv(do.call(rbind, all_dm), file.path(OUT, 'DM_all_protocols_all_lags.csv'), row.names = FALSE)
write.csv(do.call(rbind, run_status), file.path(OUT, 'MCS_run_status.csv'), row.names = FALSE)
saveRDS(list(PROXY = PROXY, LOSS = LOSS, PROTOCOLS = PROTOCOLS,
             DM_LAGS = DM_LAGS, MCS_BLOCKS = MCS_BLOCKS,
             MCS_B = MCS_B, MCS_ALPHA = MCS_ALPHA,
             MCS_STATISTIC = MCS_STATISTIC, SEED = SEED,
             session_info = capture.output(sessionInfo())),
        file.path(OUT, 'analysis_configuration.rds'))
cat('\nDone. Outputs saved to: ', OUT, '\n', sep = '')
cat('IMPORTANT: best_ex_post uses OOS network selection; inference is exploratory.\n')
cat('All-networks MCS tests the set of network-specific specifications;\n')
cat('it does not correct for every aspect of model/network selection.\n')


res <- readRDS("Caporin_BEKK_corrected_results.rds")

# Check the stored objects
names(res)
names(res$parameters)

# Number of estimated parameters
k <- 21

# Set this to the actual number of training observations
T_train <- 6244

# Calculate BIC assuming LogLik contains maximised log-likelihoods
bic_results <- res$parameters[, c("Weight", "LogLik")]

bic_results$BIC <- -2 * bic_results$LogLik +
  k * log(T_train)

# Keep only the RV forecast evaluation
rv_results <- subset(
  res$summary,
  Target == "RV",
  select = c("Weight", "RMSFE", "MAFE")
)

# Combine estimation and forecasting results
final_results <- merge(
  bic_results,
  rv_results,
  by = "Weight",
  sort = FALSE
)

# Preserve the original network order
final_results <- final_results[
  match(res$parameters$Weight, final_results$Weight),
]

print(final_results, row.names = FALSE)
