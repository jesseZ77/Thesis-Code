# Thesis-Code
Simulation code and full results. See description below for contents of each folder.

## Code
- sample_main.R: sample code for simulations (ie. not entire code), imports various functions from the files in Code/util

### Code/util
- contains utility files for simulation and estimation
- basic_functions.R: compute p matrix from given odds and margins, simulates multinomial data etc.
- error_functions.R: calculates log-odds MSE, prob MISE etc.
- optimal_lambda.R: compute optimal bandwidths using different methods (ANA, BOOT etc.)
- estimator.R: computes log-odds estimates using different methods; simulation runner (ie. wrapper to perform estimator comparison)

## Results

### Results/MSE approximation
- Plots and error values (sq. bias, var, MSE) for different distributions and sample sizes

### Results/Estimator comparison
- full MSE results and bandwidth distribution (mean, sd)
- includes all plug-in methods and some other comparison methods not included in the report
- Note that estimator names are in the form lambda_{lambda_method}_plugin{plugin_method}_{c}
  - lambda_method: 1 = ANA, 2 = BOOT, 3 = MISE, 4 = CV
  - plugin_method: 1 = ML, 2 = add c, 3 = add c to 0 counts only
- Basic estimators numbered as follows: 1 = ML, 2 = HA, 3 = add c to 0 counts only, 4 = Jewell, 5 = Agresti
