##Required packages
# Packages

library("rio")
library("BatchGetSymbols")
library("cluster")
library("TSclust")
library("lgarch")
library("MASS")
library(ggplot2)
library(reshape2)
library("lgarch")
library(dplyr)
library(classInt)
library(Matrix)
library("igraph")
library("wesanderson")
library("classInt")
library(spdep)
library(MTS)
library(stats)
library(ggcorrplot)
library(tidyverse)




# Import the Excel file containing selected S&P 500 tickers
SP500.dat <- import("sp500_selected_tickers.xlsx")

# Extract the tickers
tickers <- SP500.dat$Symbol

# Define the date range
start.date <- as.Date("1998-12-22")
end.date <- as.Date("2024-10-20")

# Retrieve stock data

l.out <- BatchGetSymbols(tickers = tickers,
                         first.date = start.date, 
                         last.date =  end.date, 
                         thresh.bad.data = 1, 
                         bench.ticker = "^GSPC", # S&P 500 index
                         type.return = "log",
                         freq.data = "daily") 



#export(l.out,"l.out.RDS")
Datar <- unstack(l.out$df.tickers, l.out$df.tickers$ret.adjusted.prices~l.out$df.tickers$ticker)
Datar <- Datar[-1,]
head(Datar)
dim(Datar)

# Remove NA values
Datar <- na.omit(Datar)

# Convert to matrix format for modeling
Datar_matrix <- as.matrix(Datar)

# Check summary statistics
summary(Datar_matrix)



# Summary statistics (mean, median, std deviation, min, max, skewness, kurtosis)
summary_stats <- Datar %>%
  summarise_all(list(
    Mean = mean, 
    Median = median, 
    SD = sd, 
    Min = min, 
    Max = max, 
    Skewness = ~e1071::skewness(.), 
    Kurtosis = ~e1071::kurtosis(.)
  ))

# Display summary statistics
print(summary_stats)



library(e1071)  # Load library for kurtosis computation

# Compute kurtosis for each stock (excluding Date column)
kurtosis_values <- apply(Datar, 2, kurtosis) 

# Count stocks with kurtosis > 3
leptokurtic_count <- sum(kurtosis_values > 3)

# Compute percentage
percentage_leptokurtic <- (leptokurtic_count / length(kurtosis_values)) * 100

# Print results
cat("Leptokurtic Series Count:", leptokurtic_count, "\n")
cat("Percentage of Leptokurtic Series:", round(percentage_leptokurtic, 2), "%\n")


# Create a dataframe for kurtosis values
kurtosis_df <- data.frame(
  variable = names(kurtosis_values),
  kurtosis = round(kurtosis_values, 2)  # Round for better readability
)

# Define colors manually: Red for kurtosis > 3, Black otherwise
kurtosis_df$color <- ifelse(kurtosis_df$kurtosis > 3, "red", "black")


# Ensure "Date" column exists (if not, create one)
if (!"Date" %in% colnames(Datar)) {
  Datar$Date <- seq(as.Date("1998-12-21"), by = "days", length.out = nrow(Datar))
}

# Melt data into long format for boxplot
Datar_melted <- melt(Datar, id.vars = "Date")

# Boxplot with Kurtosis Annotations
p1 <- ggplot(Datar_melted, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot(outlier.shape = TRUE) +
  geom_text(data = kurtosis_df, 
            aes(x = variable, y = max(Datar_melted$value, na.rm = TRUE) + 0.01, 
                label = kurtosis), 
            color = kurtosis_df$color,
            angle = 90, size = 2, hjust = 0) +
  theme_minimal() +
  labs(title = "Boxplot with Kurtosis", x = "Stock", y = "Log returns") +
  guides(fill = FALSE) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  )

# Time Series Plot of Log Returns
p2 <- ggplot(Datar_melted, aes(x = Date, y = value, color = variable)) +
  geom_line(alpha = 0.7, size = 0.7) +
  theme_minimal() +
  labs(title = "Log returns series", x = "Date", y = "Log returns", color = "S&P 500 stocks") +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),        # size of the boxes in legend
    legend.spacing.y = unit(0.2, "cm"),       # vertical spacing between items
    legend.box.margin = margin(5, 5, 5, 5),   # margin around the whole legend
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  )




# --- ACF Plots ---
# Convert to long format
Datar_long <- Datar %>%
  mutate(Obs = row_number()) %>%
  pivot_longer(-Obs, names_to = "Stock", values_to = "LogReturn")

# Compute ACFs up to lag 20 for each stock
acf_df <- Datar_long %>%
  group_by(Stock) %>%
  summarise(acf_values = list(acf(LogReturn, lag.max = 20, plot = FALSE)$acf)) %>%
  unnest(acf_values) %>%
  group_by(Stock) %>%
  mutate(Lag = 0:(n() - 1))

T <- nrow(Datar)
conf_level <- 1.96 / sqrt(T)

p3 <- ggplot(acf_df %>% filter(Lag > 0), aes(x = Lag, y = acf_values)) +
  geom_bar(stat = "identity", fill = "blue") +
  geom_hline(yintercept = c(conf_level, -conf_level), linetype = "dashed", color = "red") +
  facet_wrap(~Stock, scales = "free_y") +
  theme_minimal() +
  labs(title = "ACF of log returns", x = "Lag", y = "ACF") +
  theme(
    strip.text = element_text(size = 9),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  )



# --- Correlation matrix ---
# Compute correlation matrix
cor_mat <- cor_mat <- cor(Datar[, sapply(Datar, is.numeric)])


# Plot
p4 <- ggcorrplot(cor_mat, hc.order = TRUE, type = "lower",
                 lab = TRUE, lab_size = 3, method = "square",
                 colors = c("blue", "white", "red"),
                 title = "Correlation Matrix of log returns",
                 ggtheme = theme_minimal()) +
  theme(
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )




# --- VAR model and residuals ---
library(vars)

# Clean Datar
Datar_clean <- Datar %>%
  dplyr::select(where(is.numeric)) %>%
  na.omit()

# VAR lag selection and estimation

p <- 1
var_model <- VAR(Datar_clean, p = p, type = "const")

# Extract residuals
residuals_df <- as.data.frame(residuals(var_model))
residuals_df$Date <- Datar$Date[(p + 1):nrow(Datar)]

# Melt for plotting
res_melted <- melt(residuals_df, id.vars = "Date")

# 1. Boxplot with kurtosis
kurtosis_vals <- apply(dplyr::select(residuals_df, -Date), 2, kurtosis)
kurtosis_df_res <- data.frame(
  variable = names(kurtosis_vals),
  kurtosis = round(kurtosis_vals, 2),
  color = ifelse(kurtosis_vals > 3, "red", "black")
)

p1_res <- ggplot(res_melted, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot(outlier.shape = TRUE) +
  geom_text(data = kurtosis_df_res,
            aes(x = variable, y = max(res_melted$value, na.rm = TRUE) + 0.01,
                label = kurtosis, color = color),
            angle = 90, size = 2, hjust = 0) +
  scale_color_identity() +
  labs(title = "Boxplot with Kurtosis (VAR residuals)", x = "Stock", y = "Residuals") +
  theme_minimal() +
  guides(fill = FALSE) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  )
# 2. Time Series
p2_res <- ggplot(res_melted, aes(x = Date, y = value, color = variable)) +
  geom_line(alpha = 0.7, size = 0.7) +
  theme_minimal() +
  labs(title = "Time series of VAR residuals", x = "Date", y = "Residuals", color = "S&P 500 Stocks"  ) + # This sets the legend title
  theme(
    legend.position = "right",
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),        # size of the boxes in legend
    legend.spacing.y = unit(0.2, "cm"),       # vertical spacing between items
    legend.box.margin = margin(5, 5, 5, 5),   # margin around the whole legend
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  )


# 3. ACF Plot
acf_df_res <- residuals_df %>%
  dplyr::select(-Date) %>%  # force use of dplyr's select
  tidyr::pivot_longer(everything(), names_to = "Stock", values_to = "Residual") %>%
  dplyr::group_by(Stock) %>%
  dplyr::summarise(acf_values = list(acf(Residual, lag.max = 20, plot = FALSE)$acf)) %>%
  tidyr::unnest(acf_values) %>%
  dplyr::group_by(Stock) %>%
  dplyr::mutate(Lag = 0:(n() - 1))


T <- nrow(residuals_df)
conf_level <- 1.96 / sqrt(T)

p3_res <- ggplot(acf_df_res %>% filter(Lag > 0), aes(x = Lag, y = acf_values)) +
  geom_bar(stat = "identity", fill = "blue") +
  geom_hline(yintercept = c(conf_level, -conf_level), linetype = "dashed", color = "red") +
  facet_wrap(~Stock, scales = "free_y") +
  theme_minimal() +
  labs(title = "ACF of VAR residuals", x = "Lag", y = "ACF") +
  theme(
    strip.text = element_text(size = 9),            # Stock labels on facets
    axis.text = element_text(size = 8),             # Tick labels
    axis.title = element_text(size = 10),           # Axis titles
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)  # Main title
  )

# 4. Correlation matrix
cor_mat_res <- cor(dplyr::select(residuals_df, -Date))
p4_res <- ggcorrplot(cor_mat_res, hc.order = TRUE, type = "lower",
                     lab = TRUE, lab_size = 3, method = "square",
                     colors = c("blue", "white", "red"),
                     title = "Correlation matrix of VAR residuals",
                     ggtheme = theme_minimal()) +
  theme(
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

  

## evaluating asymmetry in data :

library(rugarch)

# Define EGARCH(1,1) model specification
spec_egarch <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1,1)),
                          mean.model = list(armaOrder = c(0,0)),
                          distribution.model = "std")  # Using Student-t for heavy tails

# Function to fit EGARCH model and extract asymmetry parameter
fit_egarch <- function(series) {
  fit <- tryCatch(ugarchfit(spec = spec_egarch, data = series), error = function(e) NULL)
  if (!is.null(fit)) {
    coef_fit <- coef(fit)
    gamma <- coef_fit["gamma1"]  # Extract asymmetry parameter
    p_value <- fit@fit$matcoef["gamma1", 4]  # Extract p-value
    return(c(gamma, p_value))
  } else {
    return(c(NA, NA))
  }
}


# Select only numeric columns dynamically
Datar_numeric <- Datar[sapply(Datar, is.numeric)]

# Apply the EGARCH fitting function only on numeric columns to all stocks
egarch_results <- as.data.frame(t(apply(Datar_numeric, 2, fit_egarch)))

colnames(egarch_results) <- c("Gamma", "P_Value")

# Count number of stocks with significant asymmetry (p < 0.05)
asymmetric_count <- sum(egarch_results$P_Value < 0.05, na.rm = TRUE)

# Compute percentage
percentage_asymmetric <- (asymmetric_count / nrow(egarch_results)) * 100

# Print results
cat("Number of stocks with asymmetric volatility:", asymmetric_count, "\n")
cat("Percentage of series with significant asymmetry:", round(percentage_asymmetric, 2), "%\n")


# Identify stocks with non-significant asymmetry (p >= 0.05)
non_asymmetric_stocks <- rownames(egarch_results[egarch_results$P_Value >= 0.05 | is.na(egarch_results$P_Value), ])

# Print them
cat("Stocks without significant EGARCH asymmetry:\n")
print(non_asymmetric_stocks)



# Select top 5 stocks with the highest absolute gamma (asymmetry)
top_asym_stocks <- rownames(egarch_results[order(abs(egarch_results$Gamma), decreasing = TRUE), ])[1:5]

# Initialize an empty list to store fitted volatilities
vol_data <- list()

# Loop through selected stocks and extract conditional volatility
for (stock in top_asym_stocks) {
  fit <- ugarchfit(spec = spec_egarch, data = Datar[[stock]])
  vol_data[[stock]] <- sigma(fit)  # Extract conditional volatility
}

# Convert list to dataframe
vol_df <- as.data.frame(vol_data)
vol_df$Date <- Datar$Date  # Add time column
vol_df_melted <- melt(vol_df, id.vars = "Date")  # Melt for ggplot

# Plot EGARCH-fitted volatilities
ggplot(vol_df_melted, aes(x = Date, y = value, color = variable)) +
  geom_line(size = 1) +
  theme_minimal() +
  labs(title = "EGARCH-Fitted Volatility for Selected Stocks",
       x = "Date", y = "Conditional Volatility",
       color = "Stock") +

  theme(legend.position = "top")
