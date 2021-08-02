generate_samples <- function(m, n, p) {
  p_vec <- c(p)
  output <- array(rmultinom(m, n, p_vec), c(2, 2, m))
  return(output)
}

generate_p_matrix <- function(omega, p0x, px0) {
  if (omega == 1) {
    p00 <- p0x * px0
    p01 <- p0x * (1 - px0)
    p10 <- (1 - p0x) * px0
    p11 <- (1 - p0x) * (1 - px0)
  } 
  else {
    a <- (omega - 1)
    b <- (1 - omega) * p0x + (1 - omega) * px0 - 1
    c <- omega * p0x * px0
    p00 <- (-b - sqrt(b^2 - 4 * a * c)) / (2 * a)
    p01 <- p0x - p00
    p10 <- px0 - p00
    p11 <- 1 - p00 - p01 - p10
    p <- c(p00, p01, p10, p11)
  }
  
  p <- c(p00, p01, p10, p11)
  p_mat <- matrix(p, 2, 2, byrow = T)
  return(p_mat)
}

compute_log_odds <- function(p) {
  return(log(p[1, 1] * p[2, 2]) - log(p[1, 2] * p[2, 1]))
}

compute_smoothed_p <- function(counts, lambda=0, n=NULL){
  if(is.null(n)) n <- sum(counts)
  smoothed_p_matrix <- matrix(0, 2, 2)
  for(i1 in 1:2) {
    for(j1 in 1:2) {
      i2 <- 3 - i1
      j2 <- 3 - j1
      smoothed_p_matrix[i1, j1] <- (1 - lambda)^2 * counts[i1, j1] / n + 
        lambda * (1 - lambda) * (counts[i1, j2] / n + counts[i2, j1] / n) + 
        lambda^2 * counts[i2, j2] / n
    }
  }
  return(smoothed_p_matrix)
}

compute_smoothed_log_odds <- function(counts, lambda=0) { 
  smoothed_p <- compute_smoothed_p(counts, lambda)
  output <- compute_log_odds(smoothed_p)
  return(output)
}