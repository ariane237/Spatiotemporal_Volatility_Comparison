# Function to estimate Spatio-Temporal EGARCH model and compute conditional variances
estimate_spEGARCH <- function(Y, W_1, W_2, n_init = 20, seed = 123) {
  set.seed(seed)
  
  # Define likelihood wrapper
  extra_args <- list(Y = Y, W_1 = W_1, W_2 = W_2)
  
  # Generate initial parameter guesses
  params_list <- lapply(1:n_init, function(i) generate_random_params(Y))
  neg_ll_values <- sapply(params_list, function(p) likelihood(p, extra_args))
  
  # Choose best starting point
  best_index <- which.min(neg_ll_values)
  best_params <- params_list[[best_index]]
  
  # Constraint: |lambda_0| + |lambda_1| < 1
  constraint_function <- function(params) abs(params["lambda_0"]) + abs(params["lambda_1"]) - 1
  
  # Optimize
  result <- solnp(
    pars = best_params,
    fun = function(params) likelihood(params, extra_args),
    ineqfun = constraint_function,
    ineqLB = -Inf,
    ineqUB = 0,
    LB = c(0.0001, 0.0001, -1, -1, -10, 0.0001, 0.0001),
    UB = c(0.99, 0.999999999999, 1, 1, 10, 1, 0.95),
    control = list(trace = 0)
  )
  
  estimate <- result$pars
  
  # Compute conditional variances h_t
  eps_hat <- array(0.01, dim = dim(Y))
  h_hat <- array(0.01, dim = dim(Y))
  
  for (t in 1:ncol(Y)) {
    if (t == 1) {
      out <- g(Y[, t], alpha = estimate[5], rhoW_1 = estimate[1] * W_1, lambdaW_2 = estimate[7] * W_2,
               theta_0 = estimate[3], theta_1 = estimate[4], xi = 1, rho_1 = estimate[2],
               lambda_1 = estimate[6], eps_lag = rep(0.001, nrow(Y)), Y_lag = rep(0.001, nrow(Y)))
    } else {
      out <- g(Y[, t], alpha = estimate[5], rhoW_1 = estimate[1] * W_1, lambdaW_2 = estimate[7] * W_2,
               theta_0 = estimate[3], theta_1 = estimate[4], xi = 1, rho_1 = estimate[2],
               lambda_1 = estimate[6], eps_lag = eps_hat[, t-1], Y_lag = Y[, t-1])
    }
    eps_hat[, t] <- out$eps
    h_hat[, t]  <- exp(out$lnh)
  }
  
  return(list(params = estimate, cond_var = h_hat, eps = eps_hat))
}


mean_abs_eps <- sqrt(2/pi)
#gamma <- 0.57721


# Numerical version to solve for epsilon
#define the function mapping y to epsilon
g <- function(y, alpha, rhoW_1, lambdaW_2, rho_1, lambda_1, theta_0, theta_1, xi, eps_lag, Y_lag){ # mapping from y to eps
  n     <- length(y)
  function_y_eps <- function(x, y, alpha, rhoW_1, lambdaW_2, rho_1, lambda_1, theta_0, theta_1, xi){
    n <- length(x)
    eps <- x
    return(    (diag(n) - lambdaW_2) %*% log(abs(y)^2) -
                 (rep(alpha,n) +
                    rhoW_1 %*% (theta_0 * eps + xi*(abs(eps) - mean_abs_eps)) +
                    rho_1 * (theta_1 * eps_lag + xi*(abs(eps_lag) - mean_abs_eps)) +
                    lambda_1 * (log(abs(Y_lag)^2) - log(abs(eps_lag)^2)) +
                    (diag(n) - lambdaW_2) %*% log(abs(eps)^2))  )
  }
  out <- nleqslv(x = y, fn = function_y_eps, alpha = alpha, rhoW_1 = rhoW_1, lambdaW_2 = lambdaW_2,
                 rho_1 = rho_1, lambda_1 = lambda_1, theta_0 = theta_0, theta_1=theta_1, xi = xi, y = y, control = list(ftol = 1e-10))
  eps <- sign(y) * abs(out$x)   #out$x
  lnh <- log(abs(y)^2) - log(abs(eps)^2)
  return(list(eps = eps, lnh =  lnh))
}


# Define the derivative of h wrt epsilon function

d_h_d_eps <- function(eps, h, alpha, theta_0, xi, rhoW_1, lambdaW_2){
  tau_eps_prime <- function(eps, theta_0, xi){
    return(theta_0 + xi*sign(eps))
  }
  n     <- length(eps)
  aux   <- solve(diag(n) - Matrix(lambdaW_2)) %*% Matrix(rhoW_1)
  eps_j <- t(array(rep(eps, n), dim = c(n, n)))
  h_i   <- array(rep(h, n), dim = c(n, n))
  return( aux %*% tau_eps_prime(eps_j, theta_0, xi) * (h_i))
}

#defining the likelihood

likelihood <- function(params,  extra_args) {
  
  Y <- extra_args$Y
  W_2 <- extra_args$W_2
  W_1 <- extra_args$W_1
  
  rho_0 <- params[1]
  rho_1 <- params[2]
  theta_0 <- params[3]
  theta_1 <- params[4]
  alpha <- params[5]
  xi <- 1
  lambda_1 <- params[6]
  lambda_0 <- params[7]
  n <- nrow(Y)
  Total <- ncol(Y)
  eps_hat <- array(0.01, dim = c(n, Total))
  
  # Compute log of the determinant of the jacobian matrix
  log_det_J <- numeric(Total)
  for (t in 1:Total) {
    if(t == 1){
      out <- g(Y[, t], alpha = alpha, rhoW_1 = rho_0 * W_1, lambdaW_2 = lambda_0 * W_2, theta_0 = theta_0, theta_1 = theta_1, xi = xi, rho_1 = rho_1, lambda_1 = lambda_1,
               eps_lag = rep(0.001, length(Y[, t])), Y_lag = rep(0.001, length(Y[, t])))
    } else {
      out <- g(Y[, t], alpha = alpha, rhoW_1 = rho_0 * W_1, lambdaW_2 = lambda_0 * W_2, theta_0 = theta_0, theta_1 = theta_1, xi = xi, rho_1 = rho_1, lambda_1 = lambda_1,
               eps_lag = eps_hat[,t-1], Y_lag = Y[,t-1])
    }
    eps_hat[,t] <- out$eps
    
    h_hat       <- exp(out$lnh)
    J <- 0.5 * (out$eps / sqrt(h_hat)) * d_h_d_eps(out$eps, h_hat, alpha, theta_0, xi, rho_0 * W_1, lambda_0 * W_2) + diag(sqrt(h_hat))
    log_det_J[t] <- determinant(J, logarithm = TRUE)$modulus 
  }
  
  likelihood <-  sum(log(dnorm(eps_hat[,-(1:5)]))) - sum(log_det_J[-(1:5)]) #sum(log(dnorm(t(Y)%*%diag((sqrt(1/h_hat))))))
  
  return(-likelihood)  # Return negative likelihood for maximization
}


# Function to generate random parameters
generate_random_params <- function(Y) {
  c(
    rho_0 = runif(1, 0, 0.3),
    rho_1 = runif(1, 0, 0.5),
    theta_0 = runif(1, -0.5, 0.51),
    theta_1 = runif(1, -0.5, 0.51),
    alpha = sd(Y), #runif(1, -0.2, 0.2),
    lambda_1 = runif(1, 0, 0.5),
    lambda_0 = runif(1, 0, 0.5)
  )
}
