######### Diebold-Mariono test for equal predictive perfomance#######


######### Diebold–Mariano Test + Model Ranking #########

library(zoo)
library(sandwich)
library(dplyr)

# === Set working directory ===

# === Import squared forecast errors (L = e^2) ===
e_sq_bekk         <- readRDS("e_sq_bekk.rds")
e_sq_bekk_asy     <- readRDS("e_sq_bekk_asy.rds")
e_sq_dcc          <- readRDS("e_sq_dcc.rds")
e_sq_dyarchotto   <- readRDS("e_sq_dyarchotto.rds")
e_sq_garchx_borov <- readRDS("e_sq_garchx_borov.rds")
e_sq_sp_egarch    <- readRDS("e_sq_sp_egarch.rds")
e_sq_proximi_capo <- readRDS("e_sq_proximi_capo.rds")
e_sq_garch_hol    <- readRDS("e_sq_garch_hol.rds")

# === Store in a named list ===
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

# === Diebold–Mariano Test Function ===
perform_dm_test <- function(err1, err2, loss = "SE", h = 1) {
  # Average across assets per time
  agg_err1 <- rowMeans(err1, na.rm = TRUE)
  agg_err2 <- rowMeans(err2, na.rm = TRUE)
  
  # Loss differential
  if (loss == "SE") {
    d_t <- agg_err1^2 - agg_err2^2
  } else if (loss == "AE") {
    d_t <- abs(agg_err1) - abs(agg_err2)
  } else {
    stop("Loss must be 'SE' or 'AE'")
  }
  
  # Drop NA values
  d_t <- d_t[is.finite(d_t)]
  T <- length(d_t)
  d_bar <- mean(d_t)
  
  # === Manual Newey–West variance ===
  # lag truncation = h - 1 (usually 0 for one-step ahead)
  lag <- h - 1
  gamma0 <- var(d_t)
  
  if (lag > 0) {
    gamma_sum <- 0
    for (l in 1:lag) {
      w_l <- 1 - l / (lag + 1)  # Bartlett weights
      gamma_l <- sum((d_t[(l + 1):T] - mean(d_t)) * (d_t[1:(T - l)] - mean(d_t))) / T
      gamma_sum <- gamma_sum + 2 * w_l * gamma_l
    }
    d_var <- gamma0 + gamma_sum
  } else {
    d_var <- gamma0
  }
  
  # === DM statistic and p-value ===
  DM_stat <- d_bar / sqrt(d_var / T)
  p_val <- 2 * (1 - pnorm(abs(DM_stat)))  # two-sided
  
  return(list(statistic = DM_stat, p_value = p_val))
}


# === 1. Compute pairwise DM test p-values ===
model_names <- names(models_losses)
n_models <- length(model_names)
dm_matrix <- matrix(NA, nrow = n_models, ncol = n_models,
                    dimnames = list(model_names, model_names))

for (i in 1:n_models) {
  for (j in 1:n_models) {
    if (i != j) {
      res <- perform_dm_test(models_losses[[i]], models_losses[[j]])
      dm_matrix[i, j] <- res$p_value
    } else {
      dm_matrix[i, j] <- NA
    }
  }
}

cat("\nPairwise Diebold–Mariano p-values (two-sided):\n")
print(round(dm_matrix, 4))

# === 2. Compute average loss per model (lower = better) ===
avg_losses <- sapply(models_losses, function(x) mean(x, na.rm = TRUE))
ranking_table <- data.frame(
  Model = names(avg_losses),
  Mean_Loss = avg_losses
) %>%
  arrange(Mean_Loss) %>%
  mutate(Rank = rank(Mean_Loss, ties.method = "min"))

cat("\nAverage Loss Ranking (lower = better):\n")
print(ranking_table)

# === 3. Identify significant differences (p < 0.05) ===
sig_matrix <- ifelse(dm_matrix < 0.05, "*", "")
sig_matrix[is.na(sig_matrix)] <- ""

cat("\nSignificance Map (* means p < 0.05):\n")
print(sig_matrix)



# === Combine ranking and DM matrix ===
# Assumes dm_matrix and ranking_table are already in memory

# Round DM p-values for readability
dm_rounded <- round(dm_matrix, 4)

# Convert p-values < 0.001 to "<0.001" for nice display
dm_display <- apply(dm_rounded, 2, function(x) ifelse(x < 0.001, "<0.001", formatC(x, digits=4, format="f")))

# Merge mean loss and rank as the leftmost columns
final_summary <- cbind(
  ranking_table[match(rownames(dm_display), ranking_table$Model), c("Mean_Loss", "Rank")],
  dm_display
)

# Clean up names
colnames(final_summary)[1:2] <- c("Mean_Loss", "Rank")

# Print nicely
cat("\n=== Combined DM Test (p-values) + Mean Loss Ranking ===\n")
print(final_summary, right = FALSE, row.names = TRUE)

# Optional: save to CSV
write.csv(final_summary, "DMtest_full_summary_with_pvalues.csv", row.names = TRUE)

# === LaTeX-ready DM Test Summary Table ===
library(xtable)

# Round and label significance
dm_disp <- round(dm_matrix, 4)
dm_disp_fmt <- apply(dm_disp, 2, function(x) {
  sapply(x, function(val) {
    if (is.na(val)) return("--")
    if (val < 0.001) return("\\textbf{<0.001}")
    if (val < 0.01)  return(sprintf("\\textbf{%.3f}", val))
    if (val < 0.05)  return(sprintf("\\textit{%.3f}", val))
    sprintf("%.3f", val)
  })
})

# Merge mean losses and ranks
dm_final <- cbind(
  Model = ranking_table$Model,
  `Mean Loss` = sprintf("%.3f", ranking_table$Mean_Loss),
  Rank = ranking_table$Rank,
  dm_disp_fmt[match(ranking_table$Model, rownames(dm_disp_fmt)), ]
)

# Create LaTeX table
tab <- xtable(dm_final, align = rep("c", ncol(dm_final) + 1),
              caption = "Pairwise Diebold--Mariano (DM) test p-values for equality of predictive accuracy. 
              Smaller p-values indicate significant differences. 
              The two leftmost columns report each model’s average loss and rank (lower = better). 
              Bold entries denote strong significance (p < 0.01); italics denote moderate significance (p < 0.05).")

print(tab, include.rownames = FALSE, sanitize.text.function = identity,
      caption.placement = "top", size = "\\footnotesize", booktabs = TRUE)





