##Required packages
# Packages

library("rio")
library("BatchGetSymbols")
library("cluster")
library("TSclust")
library("lgarch")
library("MASS")
library(ggplot2)
library(gridExtra)
library(grid)
library(reshape2)
library("lgarch")
library(classInt)
library(nleqslv)
library(Matrix)
library("igraph")
library("wesanderson")
library("classInt")
library(spdep)
library(MTS)
library(rugarch)



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



library(vars)
library(lmtest)

# Ensure data is numeric and complete
Datar_clean <- na.omit(Datar)

# Fit the VAR(1) model
var.model <- VAR(Datar_clean, p = 1, type = "const")

# View summary
summary(var.model)


# Get residuals (εt)
residuals.var <- resid(var.model)

# Convert to data frame if needed
residuals.df <- as.data.frame(residuals.var)
head(residuals.df)

# Save residuals as RDS
saveRDS(residuals.df, file = "residuals_VAR1.rds")

# Load the residuals
residuals <- readRDS("residuals_VAR1.rds")

# Check it loaded correctly
head(residuals)


train_dat <- residuals[1:(nrow(residuals) - 252),]


test_dat <- residuals[(nrow(residuals) - 251):nrow(residuals),]

train.l <- (nrow(residuals) - 252)

out.l <- nrow(test_dat) 


### construction of weight matrices 


# Construct the W matrix:
logDatar2 <- train_dat


# Alternative A: standard Euclidean distance
Wmat1 <- as.matrix(diss(t(logDatar2), "EUCL"))
Wmat1 <- 1/Wmat1
diag(Wmat1) <- 0
Wmat1 <- Wmat1 / max(eigen(Wmat1, only.values = TRUE)$values)

saveRDS(Wmat1 , "W1_eucl.rds")

# Alternative B: correlation-based distance
Wmat2 <- as.matrix(diss(t(logDatar2), "COR"))
Wmat2 <- 1/Wmat2
diag(Wmat2) <- 0
Wmat2 <- Wmat2 / max(eigen(Wmat2, only.values = TRUE)$values)

saveRDS(Wmat2 , "W1_corr.rds")

# Alternative C: log-ARCH approach

Wmat3 <- as.matrix(diss(t(logDatar2), "AR.PIC"))
Wmat3 <- 1/Wmat3
diag(Wmat3) <- 0
Wmat3 <- Wmat3 / max(eigen(Wmat3, only.values = TRUE)$values)

saveRDS(Wmat3 , "W1_picollo.rds")


#### Granger filtering matrix

# Function to compute a single Granger causality filtering matrix
granger_causality_matrix <- function(residuals, max_lag = 12, p_threshold = 0.05) {
  n <- ncol(residuals)
  causality_matrix <- matrix(0, n, n)  # Initialize matrix
  
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        test_result <- tryCatch({
          grangertest(residuals[,i] ~ residuals[,j], order = max_lag)
        }, error = function(e) return(NULL))
        
        if (!is.null(test_result)) {
          p_value <- test_result$`Pr(>F)`[2]  # Extract p-value
          if (p_value < p_threshold) {
            causality_matrix[i, j] <- 1  # Mark as significant
          }
        }
      }
    }
  }
  return(causality_matrix)
}

# Compute the single Granger causality matrix
GrangerMat <- granger_causality_matrix(logDatar2)
# Assume G is your binary Granger matrix (1 = j Granger causes i)
# For filtering, element-wise multiply:


W1_f <- Wmat1 * GrangerMat
W1_f <- W1_f / rowSums(W1_f)
diag(W1_f) <- 0
# Replace NaNs (from 0/0) with 0
W1_f[!is.finite(W1_f)] <- 0

saveRDS(W1_f, "W1_eucl_GF.rds")

W2_f <- Wmat2 * GrangerMat
W2_f <- W2_f / rowSums(W2_f)
diag(W2_f) <- 0
W2_f[!is.finite(W2_f)] <- 0

saveRDS(W2_f, "W2_corr_GF.rds")

W3_f <- Wmat3 * GrangerMat
W3_f <- W3_f / rowSums(W3_f)
diag(W3_f) <- 0
W3_f[!is.finite(W3_f)] <- 0

saveRDS(W3_f, "W3_picollo_GF.rds")



# Approach : k-nearest neighbours

# Alternative A: standard Euclidean distance

dist_mat1 <- as.matrix(diss(t(logDatar2), "EUCL"))

# 5NN
k <- 5 # Choose 3, 5 or 10
Wmat1_10    <- t(sapply(1:ncol(dist_mat1), function(i) ifelse(dist_mat1[i,] < sort(dist_mat1[i,])[k+2] & dist_mat1[i,] > 0, 1/k, 0))) # k+2, because of the diagonal zero entry and the strict inequality
diag(Wmat1_10) <- 0

saveRDS(Wmat1_10 , "W1_eucl5nn.rds")

# Alternative B: correlation-based distance

dist_mat2 <-  as.matrix(diss(t(logDatar2), "COR")) 
#5NN
k <- 5 # Choose 3, 5 or 10
Wmat2_10    <- t(sapply(1:ncol(dist_mat2), function(i) ifelse(dist_mat2[i,] < sort(dist_mat2[i,])[k+2] & dist_mat2[i,] > 0, 1/k, 0))) # k+2, because of the diagonal zero entry and the strict inequality
diag(Wmat2_10) <- 0

saveRDS(Wmat2_10 , "W2_corr5nn.rds")

# Alternative C: log-ARCH approach

dist_mat3 <- as.matrix(diss(t(logDatar2), "AR.PIC"))

#5NN
k <- 5 # Choose 3, 5 or 10
Wmat3_10    <- t(sapply(1:ncol(dist_mat3), function(i) ifelse(dist_mat3[i,] < sort(dist_mat3[i,])[k+2] & dist_mat3[i,] > 0, 1/k, 0))) # k+2, because of the diagonal zero entry and the strict inequality
diag(Wmat3_10) <- 0

saveRDS(Wmat3_10 , "W1_picollo5nn.rds")


## construction of the weight matrix based on volatility spillovers#########
construct_granger_weight_matrix <- function(logDatar2) {
  
  # Step 1: Fit EGARCH(1,1) models and extract volatilities
  spec_egarch <- ugarchspec(
    variance.model = list(model = "eGARCH", garchOrder = c(1,1)),
    mean.model = list(armaOrder = c(0,0)),
    distribution.model = "std"
  )
  
  volatility_data <- list()  # Store EGARCH volatility series
  
  for (stock in colnames(logDatar2)) {
    fit <- tryCatch(
      ugarchfit(spec = spec_egarch, data = logDatar2[[stock]], solver.control = list(trace = 0)),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      volatility_data[[stock]] <- rugarch::sigma(fit)
    } else {
      volatility_data[[stock]] <- rep(NA, nrow(logDatar2))
    }
    
  }
  
  # Convert volatility data to a dataframe
  vol_df <- as.data.frame(volatility_data)
  
  # Step 2: Compute spillovers using Granger causality
  spillover_matrix <- matrix(0, ncol = ncol(vol_df), nrow = ncol(vol_df))
  colnames(spillover_matrix) <- rownames(spillover_matrix) <- colnames(vol_df)
  
  for (i in 1:ncol(vol_df)) {
    for (j in 1:ncol(vol_df)) {
      if (i != j && all(!is.na(vol_df[[i]])) && all(!is.na(vol_df[[j]]))) {
        test <- tryCatch(grangertest(vol_df[[i]], vol_df[[j]], order = 2), error = function(e) NULL)
        if (!is.null(test)) {
          spillover_matrix[i, j] <- test$`Pr(>F)`[2]  # Store p-value
        }
      }
    }
  }
  
  # Convert p-values to weights (lower p-value = stronger spillover)
  adj_matrix <- ifelse(spillover_matrix < 0.05, 1 - spillover_matrix, 0)
  
  # Step 3: Create a financial network graph
  graph <- graph.adjacency(adj_matrix, mode = "directed", weighted = TRUE, diag = FALSE)
  
  # Step 4: Extract the weight matrix for further modeling
  weight_matrix <- as.matrix(as_adjacency_matrix(graph, attr = "weight"))
  
  # Normalize rows (each row sums to 1)
  row_sums <- rowSums(weight_matrix)
  weight_matrix <- sweep(weight_matrix, 1, row_sums, FUN = "/")
  
  # Ensure diagonal is zero
  diag(weight_matrix) <- 0
  
  return(weight_matrix)
}


# Compute the Granger-based weight matrix
Wmat_granger <- construct_granger_weight_matrix(logDatar2)


saveRDS(Wmat_granger , "Wmat_granger.rds")



