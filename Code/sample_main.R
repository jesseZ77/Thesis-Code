source("Code/util/basic_functions.R")
source("Code/util/error_functions.R")
source("Code/util/optimal_lambda.R")
source("Code/util/estimator.R")

m <- 1000
n <- 1000

p0x <- 0.3
px0 <- 0.7
omega <- 5

p <- generate_p_matrix(omega, p0x, px0)
samples <- generate_samples(m, n, p)

# Log-odds MSE Comparison

lambda_analytical <- compute_optimal_lambda_analytical(n, p)
lambda_empirical <- compute_optimal_lambda_bootstrap(n, p, samples=samples)

bias2_analytical <- compute_log_odds_bias2(n, p, lambda_analytical)
var_analytical <- compute_log_odds_var(n, p, lambda_analytical)
mse_analytical <- bias2_analytical + var_analytical

error_empirical <- compute_empirical_log_odds_error(omega, samples, lambda_analytical)
bias2_empirical <- c(error_empirical[1]^2)
var_empirical <- c(error_empirical[2])
mse_empirical <- c(error_empirical[3])

row_log_odds <- data.frame(omega, p0x, px0, n, lambda_analytical, lambda_empirical,
                       bias2_analytical, var_analytical, mse_analytical,
                       bias2_empirical, var_empirical, mse_empirical)

row_log_odds


# Estimator Comparison, all computation is handled by simulation_runner
output <- simulation_runner(omega, p0x, px0, n, m=m, samples=samples)
head(output$estimate[, 1:8])
head(output$lambda[, 1:8])
output$samples[, , 1:3]
