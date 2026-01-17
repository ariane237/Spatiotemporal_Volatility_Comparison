library("rio")
library("BatchGetSymbols")
library("cluster")
library("TSclust")
library("lgarch")
library("MASS")
library("doParallel")
library(ggplot2)
library(reshape2)
library(Rsolnp)
library(classInt)
library(nleqslv)
library(Matrix)
library("igraph")
library("wesanderson")
library(RColorBrewer)
library(scales)  # for rescale()
library(spdep) 
library(MTS)
library(gplots)
library(rugarch)
library(gridExtra)
library(dplyr)
library(vars)
library(car)  # For grangertest function
library(lmtest)


# === Load Data ===
residuals <- readRDS("residuals_VAR1.rds")
train.l <- nrow(residuals) - 252
out.l   <- 252
test.l <- out.l



# Prepare initial training dataset
Datar <- residuals
logDatar2 <- Datar[1:train.l, ]

cat("Standard deviation of returns over the considered period:\n")
print(apply(Datar, 2, sd))

# Replace zero values with small random noise
total_zero_counts <- numeric(ncol(Datar))
for (i in seq_len(ncol(Datar))) {
  total_zero_counts[i] <- sum(Datar[, i] == 0)
  Datar[, i] <- ifelse(Datar[, i] == 0, rnorm(1, 0, 0.001), Datar[, i])
}

columns_with_zeros <- total_zero_counts[total_zero_counts > 0]
names(columns_with_zeros) <- colnames(Datar)[total_zero_counts > 0]

cat("\nStocks with returns equal to zero on some days:\n")
print(columns_with_zeros)

cat("\nPercentage of zero returns in the dataset: ",
    round((sum(total_zero_counts) / (nrow(Datar) * ncol(Datar))) * 100, 3), "%\n")

cat("\nDimensions of Datar: ", paste(dim(Datar), collapse = " x "), "\n")

# Output training and test lengths
cat("Training window: ", train.l, " days\n")
cat("Test window: ", test.l, " days\n")


# Create a data frame with an additional 'Time' column to represent the sequence of dates or time points
df <- data.frame(Time = 1:nrow(Datar), Datar)

# Melt the data frame to long format for ggplot
df_melted <- melt(df, id.vars = "Time")


# Plot using ggplot2
ggplot(df_melted, aes(x = Time, y = value)) +
  geom_line(color = "blue", alpha = 0.8) +
  labs(title = "Log-Return residuals from VAR model: series 16 stocks from S&P 500",
       x = "Time",
       y = "Residuals") +
  theme_minimal() +
  facet_wrap(~ variable, scales = "free_y", ncol = 5)  # Adjust ncol to control the number of columns




# === Load Weight Matrices ===
weight_matrices <- list(
  W1      = readRDS("W1_eucl.rds"),
  W2      = readRDS("W1_corr.rds"),
  W3      = readRDS("W1_picollo.rds"),
  Wg1     = readRDS("W1_eucl_GF.rds"),
  Wg2     = readRDS("W2_corr_GF.rds"),
  Wg3     = readRDS("W3_picollo_GF.rds"),
  W1_5nn  = readRDS("W1_eucl5nn.rds"),
  W2_5nn  = readRDS("W2_corr5nn.rds"),
  W3_5nn  = readRDS("W1_picollo5nn.rds"),
  Wg      = readRDS("Wmat_granger.rds")
)

# === Node names and sector mapping ===
tickers <- names(logDatar2)  # or: tickers <- colnames(weight_matrices$W1)
names <- tickers

sectors <- c(
  AMZN = "Consumer Discretionary", BA = "Industrials", BAC = "Financials",
  CAT = "Industrials", CMCSA = "Communication Services", CSCO = "Information Technology",
  DIS = "Communication Services", INTC = "Information Technology", JNJ = "Health Care",
  JPM = "Financials", KO = "Consumer Staples", MCD = "Consumer Discretionary",
  MSFT = "Information Technology", PFE = "Health Care", T = "Communication Services",
  WMT = "Consumer Staples"
)

sector_levels <- unique(sectors)
sector_palette <- brewer.pal(n = length(sector_levels), name = "Set2")
sector_colors <- setNames(sector_palette, sector_levels)

# === Helper Function for Reuse ===
plot_network <- function(W, title, directed = FALSE, scale_size = TRUE, is_5nn = FALSE) {
  edges <- which(W > 0, arr.ind = TRUE)
  g <- graph(edges = as.vector(t(edges)), n = nrow(W), directed = directed)
  
  V(g)$label <- names
  V(g)$name <- names
  V(g)$color <- sector_colors[sectors[V(g)$name]]
  V(g)$label.color <- "black"
  V(g)$frame.color <- "black"
  
  if (scale_size) {
    conn <- rowSums(W) + colSums(W)
    V(g)$size <- rescale(conn, to = c(8, 20))
  } else {
    V(g)$size <- 14
  }
  
  if (!is_5nn) {
    classes <- classIntervals(W[edges], 8)
    E(g)$color <- findColours(classes, pal = gray(seq(0.9, 0.1, length = 8)))
  } else {
    E(g)$color <- "gray"
  }
  
  plot(
    g, layout = layout_with_kk(g),
    main = title,
    edge.arrow.size = 0.4,
    edge.width = 0.9,
    vertex.label.cex = 0.8
  )
}

# === Section A: Fully Connected Networks ===
plot_network(weight_matrices$W1, "S&P 500 stocks Euclidean-Based Network")
plot_network(weight_matrices$W2, "S&P 500 Correlation-Based Network")
plot_network(weight_matrices$W3, "S&P 500 Piccolo-Based Network")

# === Section B: 5-NN Networks ===
plot_network(weight_matrices$W1_5nn, "S&P 500 5-NN Euclidean Network", directed = TRUE, is_5nn = TRUE)
plot_network(weight_matrices$W2_5nn, "S&P 500 5-NN Correlation Network", directed = TRUE, is_5nn = TRUE)
plot_network(weight_matrices$W3_5nn, "S&P 500 5-NN Piccolo Network", directed = TRUE, is_5nn = TRUE)

# === Section C: Granger-Filtered Networks ===
plot_network(weight_matrices$Wg1, "S&P 500 Granger-Filtered Euclidean Network", directed = TRUE, is_5nn = TRUE)
plot_network(weight_matrices$Wg2, "S&P 500 Granger-Filtered Correlation Network", directed = TRUE, is_5nn = TRUE)
plot_network(weight_matrices$Wg3, "S&P 500 Granger-Filtered Piccolo Network", directed = TRUE, is_5nn = TRUE)

# === Section D: EGARCH-Based Spillover Network ===
plot_network(weight_matrices$Wg, "S&P 500 Granger-Based Volatility Spillover Network", directed = TRUE)


# === Final: Shared Sector Legend ===
# Create an empty plot with no axes or points
plot.new()
legend("center",
       legend = names(sector_colors),
       col = sector_colors,
       pch = 19,
       pt.cex = 1.5,
       bty = "n",
       title = "Sectors")

