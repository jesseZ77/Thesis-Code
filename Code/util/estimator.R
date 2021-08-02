# Wrapper for computing optimal lambda
compute_optimal_lambda <- function(n=NULL, p=NULL, counts=NULL, method=NULL, 
                                   m=1000, samples=NULL) {
  if(method == 1) return(compute_optimal_lambda_analytical(n, p))
  if(method == 2) return(compute_optimal_lambda_bootstrap(n, p, m, samples))
  if(method == 3) return(compute_optimal_lambda_prob(n, p))
  if(method == 4) return(compute_optimal_lambda_cv(counts))
}

# Compute "basic" estimates, ie. MLE, HA
estimator_basic <- function(counts, method=1, c_=0) {
  if (!is.na(dim(counts)[3])) return(apply(counts, 3, estimator_basic, method=method, c_=c_))
  
  if (method == 1) {
    # MLE / Woolf
    return(compute_log_odds(counts))
  } 
  else if (method == 2) {
    # Haldane-Anscombe
    counts_smoothed <- counts + c_
    return(compute_log_odds(counts_smoothed))
  } 
  else if(method == 3) {
    counts_smoothed <- pmax(counts, c_)
    return(compute_log_odds(counts_smoothed))
  }
  else if(method == 4) {
    # Jewell
    log_omega_hat <- compute_log_odds(counts)
    if (is.nan(log_omega_hat)) return(NaN)
    
    if (log_omega_hat >= 0) {
      counts_smoothed <- counts + (1 - diag(2))
      return(compute_log_odds(counts_smoothed))
    }
    
    counts_smoothed <- counts + diag(2)
    return(compute_log_odds(counts_smoothed))
  }
  else if(method == 5) {
    # Agresti
    n <- sum(counts)
    p_hat <- counts / n
    p_smoothed <- matrix(0, 2, 2)
    
    for(i1 in 1:2) {
      for(j1 in 1:2) {
        ni <- sum(counts[i1, ])
        nj <- sum(counts[ , j1])
        p_smoothed[i1, j1] <- n / (n + 4 * c_) * p_hat[i1, j1] + 4 * c_ / (n + 4 * c_) * ni * nj / n^2    
      }
    }
    return(compute_log_odds(p_smoothed))
  }
}

# Compute new estimates, ie. ANA, BOOT etc.
estimator_complex <- function(counts, plugin_method=-1, c_=0, lambda=NULL, lambda_method=-1, m=1000, pre_smooth=FALSE) {
  
  if(!is.na(dim(counts)[3])) return(apply(counts, 3, estimator_complex, plugin_method, c_, lambda, lambda_method, m, pre_smooth))
  
  counts2 <- counts
  if(pre_smooth) counts2 <- pmax(counts, 0.5)
  
  if (lambda_method == 4) {
    lambda <- compute_optimal_lambda(counts=counts, method=4)
    return(c(compute_smoothed_log_odds(counts2, lambda), lambda))
  }
  
  if (plugin_method == 0) {
    return(c(estimate=compute_smoothed_log_odds(counts2, lambda), lambda))
  }
  
  if (plugin_method == 1) {
    n_plugin <- sum(counts)
    p_plugin <- counts / n_plugin
    lambda <- compute_optimal_lambda(n=n_plugin, p=p_plugin, method=lambda_method, m=m)
    return(c(compute_smoothed_log_odds(counts2, lambda), lambda))
  } 
  
  if (plugin_method == 2) {
    counts_smoothed <- counts + c_
    n_plugin <- sum(counts_smoothed)
    p_plugin <- counts_smoothed / n_plugin
    lambda <- compute_optimal_lambda(n=n_plugin, p=p_plugin, method=lambda_method, m=m)
    return(c(compute_smoothed_log_odds(counts2, lambda), lambda))
  } 
  if (plugin_method == 3) {
    counts_smoothed <- pmax(counts, c_)
    n_plugin <- sum(counts_smoothed)
    p_plugin <- counts_smoothed / n_plugin
    lambda <- compute_optimal_lambda(n=n_plugin, p=p_plugin, method=lambda_method, m=m)
    return(c(compute_smoothed_log_odds(counts2, lambda), lambda))
  }
}

# Performs simulations, except BOOT (since this takes a long time)
simulation_runner <- function(omega, p0x, px0, n, m=1000, samples=NULL) {
  # non-bootstrap simulations

  p <- generate_p_matrix(omega, p0x, px0)
  
  estimate_df <- data.frame(omega=rep(omega, m), p0x, px0, n)
  lambda_df <- data.frame(omega=rep(omega, m), p0x, px0, n)
  
  if(is.null(samples)) {    
    samples <- generate_samples(m, n, p)
  }
  
  # 0. Oracle lambda
  for (method in c(1, 3)) {
    lambda <- compute_optimal_lambda(n=n, p=p, method=method)
    est <- estimator_complex(samples, plugin_method=0, lambda=lambda, pre_smooth=TRUE)
    
    name <- paste0("oracle_", method)
    estimate_df[name] <- est[1, ]
    lambda_df[name] <- est[2, ]
  }
  
  # 1. Simple lambda
  est <- estimator_complex(samples, plugin_method=0, lambda=1/n, pre_smooth=TRUE)
  estimate_df['simple_n'] <- est[1, ]
  lambda_df['simple_n'] <- est[2, ]
  
  est <- estimator_complex(samples, plugin_method=0, lambda=1/sqrt(n), pre_smooth=TRUE)
  estimate_df['simple_sqrt_n'] <- est[1, ]
  lambda_df['simple_sqrt_n'] <- est[2, ]
  
  # 2. Basic
  estimate_df['basic_1'] <- estimator_basic(samples)
  
  for (method in 2:3) {
    for (c_ in c(0.5, 1)) {
      name <- paste0("basic_", method, "_", c_)
      estimate_df[name] <- estimator_basic(samples, method=method, c_=c_)
    }
  }
  
  estimate_df['basic_4'] <- estimator_basic(samples, method=4)
  estimate_df['basic_5'] <- estimator_basic(samples, method=5, c_=0.5)
  
  # 3. Plug-in
  ## directly plug-in counts, only for lambda_method=3
  est <- estimator_complex(samples, plugin_method=1, lambda_method=3, pre_smooth=TRUE)
  estimate_df["lambda_3_plugin_1"] <- est[1, ]
  lambda_df["lambda_3_plugin_1"] <- est[2, ]
  
  ## use 2 to include bootstrap
  for (lambda_method in c(1, 3)) {
    for (plugin_method in 2:3) {
      for (c_ in c(0.5, 1)) {
        est <- estimator_complex(samples, plugin_method=plugin_method, c_=c_, 
                                 lambda_method=lambda_method, m=100, pre_smooth=TRUE)
        name <- paste0("lambda_", lambda_method, "_plugin_", plugin_method, "_", c_)
        estimate_df[name] <- est[1, ]
        lambda_df[name] <- est[2, ]
      }
    }
  }
  
  # 3. Cross validation
  est <- estimator_complex(samples, lambda_method=4, pre_smooth=TRUE)
  estimate_df["lambda_4"] <- est[1, ]
  lambda_df["lambda_4"] <- est[2, ]
  
  return(list(estimate=estimate_df, lambda=lambda_df, samples=samples))
}

# Performs simulations for BOOT only
simulation_runner_bootstrap <- function(omega, p0x, px0, n, samples=NULL, m=1000) {
  p <- generate_p_matrix(omega, p0x, px0)

  if(is.null(samples)) {    
    samples <- generate_samples(m, n, p)
  }

  estimate_df <- data.frame(omega=rep(omega, m), p0x, px0, n)
  lambda_df <- data.frame(omega=rep(omega, m), p0x, px0, n)

  # oracle
  method <- 2
  lambda <- compute_optimal_lambda(n=n, p=p, method=method, m=1000)
  est <- estimator_complex(samples, plugin_method=0, lambda=lambda, pre_smooth=TRUE)
  
  name <- paste0("oracle_", method)
  estimate_df[name] <- est[1, ]
  lambda_df[name] <- est[2, ]

  # plug-in
  lambda_method <- 2
  for (plugin_method in 2:3) {
    for (c_ in c(0.5, 1)) {
      est <- estimator_complex(samples[, , 1:1000], plugin_method=plugin_method, c_=c_, 
                                 lambda_method=lambda_method, m=100, pre_smooth=TRUE)
      name <- paste0("lambda_", lambda_method, "_plugin_", plugin_method, "_", c_)
      estimate_df[name] <- est[1, ]
      lambda_df[name] <- est[2, ]
    }
  }

  return(list(estimate=estimate_df, lambda=lambda_df, samples=samples))
}