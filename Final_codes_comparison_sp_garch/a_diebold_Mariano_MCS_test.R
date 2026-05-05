######### Diebold–Mariano Test + Model Ranking + Model Confidence Set #########
# Hansen, Lunde & Nason (2011, Econometrica 79:453-497)
# MCS package version 0.2.0

library(zoo)
library(sandwich)
library(dplyr)
library(MCS)    
library(xtable)

# === Set working directory ===


# ============================================================
# SECTION 0: Load data
# ============================================================

e_sq_bekk         <- readRDS("e_sq_bekk.rds")
e_sq_bekk_asy     <- readRDS("e_sq_bekk_asy.rds")
e_sq_dcc          <- readRDS("e_sq_dcc.rds")
e_sq_dyarchotto   <- readRDS("e_sq_dyarchotto.rds")
e_sq_garchx_borov <- readRDS("e_sq_garchx_borov.rds")
e_sq_sp_egarch    <- readRDS("e_sq_sp_egarch.rds")
e_sq_proximi_capo <- readRDS("e_sq_proximi_capo.rds")
e_sq_garch_hol    <- readRDS("e_sq_garch_hol.rds")

# Named list — e_sq_* objects are T x n matrices of squared forecast errors
models_losses <- list(
  BEKK         = e_sq_bekk,
  BEKK_Asy     = e_sq_bekk_asy,
  DCC          = e_sq_dcc,
  DYARCH_Otto  = e_sq_dyarchotto,
  GARCHX_Borov = e_sq_garchx_borov,
  SpEGARCH     = e_sq_sp_egarch,
  Proximity    = e_sq_proximi_capo,
  ST_GARCH_Hol = e_sq_garch_hol
)

# ============================================================
# SECTION 1: Diebold–Mariano Test
# ============================================================

perform_dm_test <- function(err1, err2, loss = "SE", h = 1) {
  # Average across assets per time period
  agg_err1 <- rowMeans(err1, na.rm = TRUE)
  agg_err2 <- rowMeans(err2, na.rm = TRUE)
  
  # Loss differential
  # NOTE: e_sq already contains squared errors — use directly for SE
  #       For AE, take sqrt first to recover absolute errors
  if (loss == "SE") {
    d_t <- agg_err1 - agg_err2
  } else if (loss == "AE") {
    d_t <- sqrt(agg_err1) - sqrt(agg_err2)
  } else {
    stop("Loss must be 'SE' or 'AE'")
  }
  
  d_t   <- d_t[is.finite(d_t)]
  T     <- length(d_t)
  d_bar <- mean(d_t)
  lag   <- h - 1
  gamma0 <- var(d_t)
  
  if (lag > 0) {
    gamma_sum <- 0
    for (l in 1:lag) {
      w_l     <- 1 - l / (lag + 1)
      gamma_l <- sum((d_t[(l+1):T] - mean(d_t)) *
                       (d_t[1:(T-l)]  - mean(d_t))) / T
      gamma_sum <- gamma_sum + 2 * w_l * gamma_l
    }
    d_var <- gamma0 + gamma_sum
  } else {
    d_var <- gamma0
  }
  
  DM_stat <- d_bar / sqrt(d_var / T)
  p_val   <- 2 * (1 - pnorm(abs(DM_stat)))   # two-sided
  return(list(statistic = DM_stat, p_value = p_val))
}

# --- Pairwise DM matrix ---
model_names <- names(models_losses)
n_models    <- length(model_names)
dm_matrix   <- matrix(NA, nrow = n_models, ncol = n_models,
                      dimnames = list(model_names, model_names))

for (i in 1:n_models) {
  for (j in 1:n_models) {
    if (i != j) {
      res            <- perform_dm_test(models_losses[[i]], models_losses[[j]])
      dm_matrix[i, j] <- res$p_value
    }
  }
}

cat("\nPairwise Diebold–Mariano p-values (two-sided):\n")
print(round(dm_matrix, 4))

# --- Average loss ranking ---
# e_sq already squared — mean gives MSFE per model
avg_losses    <- sapply(models_losses, function(x) mean(x, na.rm = TRUE))
ranking_table <- data.frame(
  Model     = names(avg_losses),
  Mean_Loss = avg_losses
) %>%
  arrange(Mean_Loss) %>%
  mutate(Rank = rank(Mean_Loss, ties.method = "min"))

cat("\nAverage Loss Ranking (lower = better):\n")
print(ranking_table)

# --- Significance map ---
sig_matrix <- ifelse(!is.na(dm_matrix) & dm_matrix < 0.05, "*", "")
cat("\nSignificance Map (* = p < 0.05):\n")
print(sig_matrix)

# ============================================================
# SECTION 2: Build T x M loss matrices for MCS
# ============================================================
# MCSprocedure() requires a T x M matrix where each column is the
# per-period loss series for one model (T = forecast periods, M = models).
#
# e_sq_* are T x n (n = 16 assets).
# We average across assets to get a single loss per period per model.
#
# SE loss : mean squared error across assets at each t
#           e_sq already stores squared errors → rowMeans directly
#
# AE loss : mean absolute error across assets at each t
#           absolute error = sqrt(squared error) → rowMeans(sqrt(x))

loss_matrix_SE <- sapply(models_losses, function(x) {
  rowMeans(x, na.rm = TRUE)          # already squared — no double-squaring
})

loss_matrix_AE <- sapply(models_losses, function(x) {
  rowMeans(sqrt(x), na.rm = TRUE)    # recover absolute errors first
})

cat("\nLoss matrix dimensions (T x M):", dim(loss_matrix_SE), "\n")
cat("Expected: 252 x", n_models, "\n")

# Sanity check — column means should match avg_losses above
cat("\nColumn means of SE loss matrix (should match avg_losses):\n")
print(round(colMeans(loss_matrix_SE), 6))

# ============================================================
# SECTION 3: Model Confidence Set (MCS package v0.2.0)
# ============================================================
# MCSprocedure() arguments in v0.2.0:
#   Loss      : T x M numeric matrix of losses
#   alpha     : significance level (0.10 → 90% MCS)
#   B         : bootstrap replications (5000 recommended)
#   statistic : "Tmax" (range, recommended) or "TR" (semi-quadratic)
#   k         : NULL = test all models
#   min.k     : stop when k models remain (default 3)
#   seed      : integer seed for reproducibility (replaces set.seed())
#   verbose   : TRUE to print elimination steps

# --- MCS under SE loss (RMSFE-equivalent) ---
cat("\n===================================================\n")
cat(" Model Confidence Set — Squared Error (SE) Loss\n")
cat("===================================================\n")

mcs_SE <- MCSprocedure(
  Loss      = loss_matrix_SE,
  alpha     = 0.10,
  B         = 5000,
  statistic = "Tmax",
  seed      = 42,
  verbose   = TRUE
)

cat("\nMCS results (SE loss):\n")
print(mcs_SE)

# --- MCS under AE loss (MAFE-equivalent) ---
cat("\n===================================================\n")
cat(" Model Confidence Set — Absolute Error (AE) Loss\n")
cat("===================================================\n")

mcs_AE <- MCSprocedure(
  Loss      = loss_matrix_AE,
  alpha     = 0.10,
  B         = 5000,
  statistic = "Tmax",
  seed      = 42,
  verbose   = TRUE
)

cat("\nMCS results (AE loss):\n")
print(mcs_AE)

# ============================================================
# SECTION 4: Extract and display MCS results
# ============================================================

extract_mcs_summary <- function(mcs_obj, loss_label, ranking_table) {
  
  # Results are in @show: rows = models, cols = Avg.Loss, p-Value, MCS p-Value
  show_df <- as.data.frame(mcs_obj@show)
  
  # Models in the superior set come directly from @Info$included
  in_mcs <- mcs_obj@Info$included
  
  summary_df <- data.frame(
    Model      = rownames(show_df),
    Avg_Loss   = round(show_df[["Avg.Loss"]], 4),
    MCS_pvalue = round(show_df[["MCS p-Value"]], 4),
    In_MCS     = ifelse(rownames(show_df) %in% in_mcs, "Yes", "No"),
    stringsAsFactors = FALSE
  ) %>%
    left_join(ranking_table[, c("Model", "Rank")], by = "Model") %>%
    arrange(Rank)
  
  cat(sprintf("\n--- MCS Summary: %s Loss ---\n", loss_label))
  cat(sprintf("Models in the 90%% MCS: %s\n\n",
              paste(in_mcs, collapse = ", ")))
  print(summary_df)
  
  return(summary_df)
}

mcs_summary_SE <- extract_mcs_summary(mcs_SE, "Squared Error (SE)", ranking_table)
mcs_summary_AE <- extract_mcs_summary(mcs_AE, "Absolute Error (AE)", ranking_table)
# ============================================================
# SECTION 5: Combined summary table (DM p-values + MCS)
# ============================================================

dm_rounded <- round(dm_matrix, 4)
dm_display <- apply(dm_rounded, 2, function(x)
  ifelse(is.na(x), "--",
         ifelse(x < 0.001, "<0.001",
                formatC(x, digits = 4, format = "f"))))

# Align MCS membership to dm_display row order
mcs_se_col <- mcs_summary_SE$In_MCS[
  match(rownames(dm_display), mcs_summary_SE$Model)]
mcs_ae_col <- mcs_summary_AE$In_MCS[
  match(rownames(dm_display), mcs_summary_AE$Model)]

final_summary <- cbind(
  ranking_table[match(rownames(dm_display), ranking_table$Model),
                c("Mean_Loss", "Rank")],
  In_MCS_SE = mcs_se_col,
  In_MCS_AE = mcs_ae_col,
  dm_display
)

cat("\n=== Combined Table: Mean Loss + MCS Membership + DM p-values ===\n")
print(final_summary, right = FALSE, row.names = TRUE)

#write.csv(final_summary, "DMtest_MCS_full_summary.csv", row.names = TRUE)
#cat("\nSaved: DMtest_MCS_full_summary.csv\n")
