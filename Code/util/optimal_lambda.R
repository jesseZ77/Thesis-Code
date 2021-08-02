# ANA
compute_optimal_lambda_analytical <- function(n, p) {
  # common values
  x1 <- (p[1, 2] + p[2, 1]) * (1 / p[1, 1] + 1 / p[2, 2]) - (p[1, 1] + p[2, 2]) * (1 / p[1, 2] + 1 / p[2, 1])
  x2 <- 1 / p[1, 1] - 1 / p[1, 2] - 1 / p[2, 1] + 1 / p[2, 2]
  x3 <- sum(1 / p)
  
  # coefficients
  a1 <- x1
  a2 <- p[2, 2] / p[1, 1] - p[2, 1] / p[1, 2] - p[1, 2] / p[2, 1] + p[1, 1] / p[2, 2] - x1
  a3 <- -x2 / 2
  a4 <-  x1 + 2 * x2
  a5 <- 1 / 2 * (
    (p[1, 1] + p[2, 2])^2 * (1 / p[1, 2]^2 + 1 / p[2, 1]^2) - (p[1, 2] + p[2, 1])^2 * (1 / p[1, 1]^2 + 1 / p[2, 2]^2)
  ) + 2 * x1
  
  b1 <- -8 * x3
  b2 <- (p[1, 1] + p[2, 2]) * (1 / p[1, 2] + 1 / p[2, 1])^2 + (p[1, 2] + p[2, 1]) * (1 / p[1, 1] + 1 / p[2, 2])^2 - 
    x1^2 + 20 * x3
  
  num <- -2 * a1 * a3 - b1
  den <- 2 * n * a1^2 + 4 * (a1 * a4 + a2 * a3 + a3 * a5) + 2 * b2
  output <- num / den
  output <- min(ifelse(output <= 0, 1 / n, output), 0.25)
  return(output)
}

# BOOT
compute_optimal_lambda_bootstrap <- function(n, p, m=1000, samples=NULL) {
  log_omega <- compute_log_odds(p)
  
  if(is.null(samples)) samples <- generate_samples(m, n, p)
  
  f <- function(lambda){
    if(length(lambda) > 1) return(sapply(lambda, f))
    estimates <- apply(samples, 3, compute_smoothed_log_odds, lambda=lambda)
    output <- mean((estimates - log_omega)^2)
    return(output)
  }

  output <- optim(0.01, f, lower=0, upper=0.25, method='Brent')$par
  output <- min(ifelse(output <= 0, 1 / n, output), 0.25)
  return(output)
}

# MISE
compute_optimal_lambda_prob <- function(n, p, use_exact=FALSE) {
  if(use_exact) {
    f <- function(x) compute_prob_mise(n, p, x, use_exact=TRUE)
    output <- optim(0.01, f, method="Brent", lower=0, upper=0.25)$par
    output <- min(ifelse(output <= 0, 1 / n, output), 0.25)
    return(output)
  } 
  
  x1 <- sum(p^2)
  x2 <- (p[1, 1] + p[2, 2]) * (p[1, 2] + p[2, 1])
  x3 <- p[1, 1] * p[2, 2] + p[1, 2] * p[2, 1]
  
  a1 <- -4 + 4 * x1 - 4 * x2
  a2 <- 8 - 8 * x1 + 12 * x2 - 8 * x3
  a3 <- 6 * x1 - 8 * x2 + 4 * x3
  
  output <- -a1 / (2 * (a2 + n * a3))
  output <- min(ifelse(output <= 0, 1 / n, output), 0.25)
  return(output)
}

# CV
compute_optimal_lambda_cv <- function(counts, use_exact=FALSE) {
  if (use_exact) {
    f <- function(x) compute_cv_ise(counts, x, use_exact=TRUE)
    output <- optim(0.01, f, method="Brent", lower=0, upper=0.25)$par
    output <- min(ifelse(output <= 0, 1 / n, output), 0.25)
    return(output)
  }
  
  n <- sum(counts)
  p <- counts / n
  
  x1 <- sum(p^2)
  x2 <- (p[1, 1] + p[2, 2]) * (p[1, 2] + p[2, 1])
  x3 <- (p[1, 1] * p[2, 2] + p[1, 2] * p[2, 1])
  
  a1 <- -4 * x1 + 4 * x2
  a2 <- 8 * x1 - 12 * x2 + 8 * x3
  b1 <- 2 / (n - 1) + 2 * n / (n- 1) * (-x1 + x2)
  b2 <- -1 / (n - 1) + n / (n - 1) * (x1 - 2 * x2 + 2 * x3)
  
  num <- 2 * b1 - a1
  den <- 2 * a2 - 4 * b2
  output <- num / den
  output <- min(ifelse(output <= 0, 1 / n, output), 0.25)
  return(output)
}