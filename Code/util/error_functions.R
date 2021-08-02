# Log odds Analytical MSE
compute_log_odds_bias2 <- function(n, p, lambda) {
  # common values
  x1 <- (p[1, 2] + p[2, 1]) * (1 / p[1, 1] + 1 / p[2, 2]) - 
    (p[1, 1] + p[2, 2]) * (1 / p[1, 2] + 1 / p[2, 1])
  x2 <- 1 / p[1, 1] - 1 / p[1, 2] - 1 / p[2, 1] + 1 / p[2, 2]
  
  # coefficients
  a1 <- x1
  a2 <- p[2, 2] / p[1, 1] - p[2, 1] / p[1, 2] - p[1, 2] / p[2, 1] + p[1, 1] / p[2, 2] - x1
  a3 <- -x2 / 2
  a4 <-  x1 + 2 * x2
  a5 <- 1 / 2 * (
    (p[1, 1] + p[2, 2])^2 * (1 / p[1, 2]^2 + 1 / p[2, 1]^2) - (p[1, 2] + p[2, 1])^2 * (1 / p[1, 1]^2 + 1 / p[2, 2]^2)
    ) + 2 * x1
  
  # squared bias
  output <- a1^2 * lambda^2 + 
    2 * a1 * a3 * lambda / n + 
    2 * (a1 * a4 + a2 * a3 + a3 * a5) * lambda^2 / n
  return(output)
}

compute_log_odds_var <- function(n, p, lambda) {
  # common values
  x1 <- sum(1 / p)
  x2 <- (p[1, 2] + p[2, 1]) * (1 / p[1, 1] + 1 / p[2, 2]) - (p[1, 1] + p[2, 2]) * (1 / p[1, 2] + 1 / p[2, 1])
  
  # coefficients
  b0 <- x1
  b1 <- -8 * x1
  b2 <- (p[1, 1] + p[2, 2]) * (1 / p[1, 2] + 1 / p[2, 1])^2 + (p[1, 2] + p[2, 1]) * (1 / p[1, 1] + 1 / p[2, 2])^2 - 
    x2^2 + 20 * x1
  
  # variance
  output <- b0 / n + b1 * lambda / n + b2 * lambda^2 / n
  
  return(output)
}

compute_log_odds_mse <- function(n, p, lambda) {
  return(compute_log_odds_bias2(n, p, lambda) + compute_log_odds_var(n, p, lambda))
}

# Log odds Empirical MSE
compute_empirical_log_odds_error <- function(omega, samples, lambda) {
  if(length(lambda) > 1) return(sapply(lambda, compute_empirical_log_odds_error, omega=omega, samples=samples))
  
  estimates <- apply(samples, 3, compute_smoothed_log_odds, lambda=lambda)
  bias2_ <- (mean(estimates) - log(omega))^2
  var_ <- (n - 1) / n * var(estimates)
  mse_ <- bias2_ + var_
  
  return(matrix(c(bias2_, var_, mse_), nrow=3))
}

# Prob MISE
compute_prob_bias_ij <- function(n, p, lambda, i, j, use_exact=FALSE) {
  i2 <- 3 - i
  j2 <- 3 - j
  output <- lambda * (p[i, j2] + p[i2, j] - 2 * p[i, j])
  if (use_exact) output <- output + lambda^2 * (p[i, j] - p[i, j2] - p[i2, j] + p[i2, j2])  
  return(output)
}

compute_prob_bias_matrix <- function(n, p, lambda, use_exact=FALSE) {
  bias_matrix <- matrix(0, 2, 2)
  for(i in 1:2){
    for(j in 1:2){
      bias_matrix[i, j] <- compute_prob_bias_ij(n, p, lambda, i, j, use_exact=use_exact)
    }
  }
  return(bias_matrix)
}

compute_prob_var_ij <- function(n, p, lambda, i, j, use_exact=FALSE) {
  i2 <- 3 - i
  j2 <- 3 - j
  
  if(use_exact) {
    # variance terms
    a1 <- p[i, j] * (1 - p[i, j]) * (1 - lambda)^4 / n
    a2 <- p[i, j2] * (1 - p[i, j2]) * lambda^2 * (1 - lambda)^2 / n
    a3 <- p[i2, j] * (1 - p[i2, j]) * lambda^2 * (1 - lambda)^2 / n
    a4 <- p[i2, j2] * (1 - p[i2, j2]) * lambda^4 / n
    
    # covariance terms
    b12 <- -2 * p[i, j] * p[i, j2] * lambda * (1 - lambda)^3 / n
    b13 <- -2 * p[i, j] * p[i2, j] * lambda * (1 - lambda)^3 / n
    b14 <- -2 * p[i, j] * p[i2, j2] * lambda^2 * (1 - lambda)^2 / n
    b23 <- -2 * p[i, j2] * p[i2, j] * lambda^2 * (1 - lambda)^2 / n
    b24 <- -2 * p[i, j2] * p[i2, j2] * lambda^3 * (1 - lambda) / n
    b34 <- -2 * p[i2, j] * p[i2, j2] * lambda^3 * (1 - lambda) / n
    
    output <- a1 + a2 + a3 + a4 + b12 + b13 + b14 + b23 + b24 + b34
    return(output)
  }
  
  b0 <- p[i, j] * (1 - p[i, j])
  b1 <- -2 * p[i, j] * (2 * (1 - p[i, j]) + p[i, j2] + p[i2, j])
  b2 <- 6 * p[i, j] * (1 - p[i, j] + p[i, j2] + p[i2, j]) + 
    p[i, j2] * (1 - p[i, j2]) + p[i2, j] * (1 - p[i2, j]) -
    2 * p[i, j] * p[i2, j2] - 2 * p[i, j2] * p[i2, j]
  
  output <- b0 / n + b1 * lambda / n + b2 * lambda^2 / n
  return(output)
}

compute_prob_var_matrix <- function(n, p, lambda, use_exact=FALSE) {
  var_matrix <- matrix(0, 2, 2)
  for(i in 1:2){
    for(j in 1:2){
      var_matrix[i, j] <- compute_prob_var_ij(n, p, lambda, i, j, use_exact=use_exact)
    }
  }
  return(var_matrix)
}

compute_prob_mse_matrix <- function(n, p, lambda, use_exact=FALSE) {
  output <- compute_prob_bias_matrix(n, p, lambda, use_exact)^2 + compute_prob_var_matrix(n, p, lambda, use_exact)
  return(output)
}

compute_prob_mise <- function(n, p, lambda, use_exact=FALSE) {
  if(use_exact == FALSE) {
    x1 <- sum(p^2)
    x2 <- (p[1, 1] + p[2, 2]) * (p[1, 2] + p[2, 1])
    x3 <- p[1, 1] * p[2, 2] + p[1, 2] * p[2, 1]
    
    a0 <- 1 - x1
    a1 <- -4 + 4 * x1 - 4 * x2
    a2 <- 8 - 8 * x1 + 12 * x2 - 8 * x3
    a3 <- 6 * x1 - 8 * x2 + 4 * x3
    
    output <- a0 / n + a1 * lambda / n + a2 * lambda^2 / n + a3 * lambda^2
    return(output)
  }
  else {
    output <- sapply(lambda, function(x) sum(compute_prob_mse_matrix(n, p, x, TRUE)))
    return(output)
  }
}

# Empirical Prob MISE - for verification purposes only
compute_empirical_prob_mise <- function(p, samples, lambda, n=NULL) {
  if(length(lambda) > 1) return(sapply(lambda, compute_empirical_prob_error, p=p, samples=samples, n=n))
  
  estimates <- apply(samples, 3, compute_smoothed_p, lambda=lambda, n=n)
  output <- mean(apply((estimates - c(p))^2, 2, sum))
  return(output)
}

# Prob ISE, CV
compute_cv_ise <- function(counts, lambda, use_exact=FALSE) {
  n <- sum(counts)
  p <- counts / n
  
  x1 <- sum(p^2)
  x2 <- (p[1, 1] + p[2, 2]) * (p[1, 2] + p[2, 1])
  x3 <- (p[1, 1] * p[2, 2] + p[1, 2] * p[2, 1])
  
  if(use_exact) {
    term1 <- sum(compute_smoothed_p(counts, lambda)^2)
  } else {
    a0 <- x1
    a1 <- -4 * x1 + 4 * x2
    a2 <- 8 * x1 - 12 * x2 + 8 * x3
    term1 <- a0 + a1 * lambda + a2 * lambda^2
  }
  
  b0 <- -1 / (n - 1) + n / (n - 1) * x1
  b1 <- 2 / (n - 1) + 2 * n / (n- 1) * (-x1 + x2)
  b2 <- -1 / (n - 1) + n / (n - 1) * (x1 - 2 * x2 + 2 * x3)
  term2 <- b0 + b1 * lambda + b2 * lambda^2
  
  output <- term1 - 2 * term2
  return(output)
}

# Prob ISE, True - for verification purposes only
compute_true_ise <- function(counts, p, lambda) {
  if(length(lambda) > 1) return(sapply(lambda, compute_true_ise, counts=counts, p=p))
  p_smoothed <- compute_smoothed_p(counts, lambda)
  
  return(sum(p_smoothed^2 - 2 * p * p_smoothed))
}