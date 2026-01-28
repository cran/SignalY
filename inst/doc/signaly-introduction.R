## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----install, eval=FALSE------------------------------------------------------
# # Install from CRAN (when available)
# install.packages("SignalY")
# 
# # Install development version from GitHub
# # remotes::install_github("username/SignalY")
# 
# # For Bayesian methods (Horseshoe, HP-GC), install cmdstanr:
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# cmdstanr::install_cmdstan()

## ----quickstart, eval=FALSE---------------------------------------------------
# library(SignalY)
# 
# # Generate example data
# set.seed(42)
# n <- 100
# p <- 20
# 
# # Create correlated predictors with factor structure
# factors <- matrix(rnorm(n * 2), n, 2)
# loadings <- matrix(runif(p * 2, -1, 1), p, 2)
# X <- factors %*% t(loadings) + matrix(rnorm(n * p, 0, 0.5), n, p)
# colnames(X) <- paste0("X", 1:p)
# 
# # True signal depends on 5 predictors
# true_beta <- c(rep(1, 5), rep(0, 15))
# Y <- X %*% true_beta + rnorm(n, 0, 0.5)
# 
# # Combine into data frame
# data <- data.frame(Y = Y, X)
# 
# # Run comprehensive analysis
# result <- signal_analysis(
#   data = data,
#   y_formula = "Y",
#   methods = "all",
#   verbose = TRUE
# )
# 
# # View results
# print(result)
# summary(result)
# plot(result)

## ----wavelet, eval=FALSE------------------------------------------------------
# # Apply wavelet decomposition
# Y <- sin(seq(0, 4*pi, length.out = 128)) + rnorm(128, 0, 0.3)
# 
# wav_result <- filter_wavelet(
#   x = Y,
#   filter = "la8",      # Daubechies least asymmetric, 8 vanishing moments
#   levels = c(3, 4),    # Combine D3 + D4 (business cycle frequencies)
#   first_diff = FALSE
# )
# 
# # Examine results
# names(wav_result)
# # [1] "combined" "smooth" "details" "level_variance" "filter" "levels"
# 
# # Combined signal captures 8-32 period cycles (for annual data)
# plot(wav_result$combined, type = "l", main = "Wavelet-Filtered Signal")

## ----emd, eval=FALSE----------------------------------------------------------
# emd_result <- filter_emd(
#   x = Y,
#   max_imf = 10,
#   boundary = "wave"
# )
# 
# # IMFs are ordered from highest to lowest frequency
# matplot(emd_result$imf[, 1:3], type = "l",
#         main = "First 3 Intrinsic Mode Functions")
# 
# # Residue captures the trend
# plot(emd_result$residue, type = "l", main = "EMD Residue (Trend)")

## ----hpgc, eval=FALSE---------------------------------------------------------
# # Requires cmdstanr
# hpgc_result <- filter_hpgc(
#   x = Y,
#   lambda = NULL,  # Auto-select via DIC
#   n_chains = 4,
#   n_iter = 2000
# )
# 
# # Compare to standard HP filter
# plot(Y, type = "l", col = "gray")
# lines(hpgc_result$trend, col = "blue", lwd = 2)
# legend("topright", c("Original", "Bayesian Trend"),
#        col = c("gray", "blue"), lwd = c(1, 2))

## ----horseshoe, eval=FALSE----------------------------------------------------
# # Fit Horseshoe regression
# hs_fit <- fit_horseshoe(
#   y = Y,
#   X = X,
#   p0 = 5,               # Expected non-zeros (can be NULL for auto)
#   n_chains = 4,
#   n_iter = 2000,
#   n_warmup = 1000,
#   adapt_delta = 0.95,   # Target acceptance (increase if divergences)
#   use_qr = TRUE         # QR decomposition for multicollinearity
# )
# 
# # Examine shrinkage
# print(hs_fit)
# 
# # Key outputs:
# # - beta: Coefficient estimates
# # - kappa: Shrinkage factors (0 = no shrinkage, 1 = full shrinkage)
# # - m_eff: Effective number of non-zeros
# 
# # Select variables by shrinkage threshold
# selection <- select_by_shrinkage(hs_fit, threshold = 0.5)
# which(selection$selected)  # Indices of selected predictors

## ----pca, eval=FALSE----------------------------------------------------------
# pca_result <- pca_bootstrap(
#   X = X,
#   n_components = NULL,      # Auto-select
#   rotation = "varimax",     # Or "none", "oblimin"
#   n_boot = 1000,
#   block_length = NULL,      # Auto: sqrt(T)
#   significance_level = 0.05
# )
# 
# # Key outputs
# pca_result$variance_explained  # Proportion by component
# pca_result$loadings            # Variable loadings
# pca_result$entropy             # Entropy of loadings (concentration measure)
# pca_result$loadings_significant  # Bootstrap significance
# 
# # Interpretation
# # Low entropy = concentrated loadings (few variables dominate)
# # High entropy = diffuse loadings (many variables contribute)

## ----dfm, eval=FALSE----------------------------------------------------------
# dfm_result <- estimate_dfm(
#   X = X,
#   n_factors = NULL,         # Auto-select via Bai-Ng IC
#   max_factors = 10,
#   var_lags = 1,             # VAR lags for factor dynamics
#   ic_criterion = "bai_ng_2"
# )
# 
# # Key outputs
# dfm_result$n_factors         # Optimal number
# dfm_result$factors           # Estimated factors (T x r)
# dfm_result$loadings          # Factor loadings (p x r)
# dfm_result$variance_explained
# dfm_result$var_coefficients  # VAR transition matrix

## ----unitroot, eval=FALSE-----------------------------------------------------
# ur_result <- test_unit_root(
#   x = Y,
#   tests = c("adf", "ers", "kpss", "pp")
# )
# 
# # Synthesized conclusion
# ur_result$synthesis
# # $conclusion: "stationary", "trend_stationary", "difference_stationary", or "inconclusive"
# # $confidence: "high", "medium", "low"
# # $evidence: Detailed reasoning
# 
# # Individual test results
# ur_result$tests$adf_none
# ur_result$tests$ers_dfgls

## ----master, eval=FALSE-------------------------------------------------------
# result <- signal_analysis(
#   data = data,
#   y_formula = "Y",              # Or formula: Y ~ X1 + X2 + X3
#   time_var = NULL,              # Time variable name
#   group_var = NULL,             # Panel group variable
#   methods = "all",              # Or c("wavelet", "horseshoe", "pca")
# 
#   # Method-specific configuration
#   filter_config = list(
#     wavelet_filter = "la8",
#     wavelet_levels = c(3, 4),
#     hpgc_n_chains = 4
#   ),
# 
#   horseshoe_config = list(
#     p0 = NULL,                  # Auto-calibrate
#     n_chains = 4,
#     adapt_delta = 0.95,
#     kappa_threshold = 0.5
#   ),
# 
#   pca_config = list(
#     rotation = "none",
#     n_boot = 1000
#   ),
# 
#   dfm_config = list(
#     ic_criterion = "bai_ng_2"
#   ),
# 
#   # General options
#   na_action = "interpolate",    # Or "omit", "fail"
#   standardize = TRUE,
#   first_difference = FALSE,
#   verbose = TRUE,
#   seed = 42
# )
# 
# # Access components
# result$filters$wavelet
# result$horseshoe
# result$pca
# result$dfm
# result$unitroot
# result$interpretation

## ----interp, eval=FALSE-------------------------------------------------------
# result$interpretation$signal_characteristics
# # $smoothness: Variance of second differences
# # $smoothness_interpretation: "Very smooth", "Moderately volatile", etc.
# # $trend_share: Proportion of variance from trend
# 
# result$interpretation$variable_selection
# # $sparsity_ratio: Proportion shrunk to zero
# # $n_selected: Number of selected predictors
# # $top_predictors: Data frame of top variables
# 
# result$interpretation$factor_structure
# # $pc1_entropy: Shannon entropy of PC1 loadings
# # $topology_interpretation: "Concentrated", "Diffuse", etc.
# 
# result$interpretation$persistence
# # $conclusion: Stationarity type
# # $interpretation: Plain-language description
# 
# result$interpretation$overall_summary
# # Combined narrative synthesis

## ----explore, eval=FALSE------------------------------------------------------
# # Visualize the data first
# plot(Y, type = "l")
# acf(Y)
# pacf(Y)
# 
# # Check for obvious non-stationarity
# ur_quick <- test_unit_root(Y, tests = "adf")

## ----selective, eval=FALSE----------------------------------------------------
# # For trend extraction only
# result <- signal_analysis(data, "Y", methods = c("wavelet", "hpgc"))
# 
# # For variable selection only
# result <- signal_analysis(data, "Y", methods = "horseshoe")
# 
# # For factor analysis only
# result <- signal_analysis(data, "Y", methods = c("pca", "dfm"))

## ----diag, eval=FALSE---------------------------------------------------------
# # Check Horseshoe convergence
# result$horseshoe$diagnostics
# # Look for: n_divergent = 0, max_rhat < 1.01, ESS > 400
# 
# # Check LOO-CV
# result$horseshoe$loo
# # Pareto k > 0.7 indicates problematic observations

## ----memory, eval=FALSE-------------------------------------------------------
# # Use fewer bootstrap replications
# pca_config = list(n_boot = 500)
# 
# # Use fewer MCMC iterations
# horseshoe_config = list(n_iter = 1000, n_warmup = 500)
# 
# # Analyze subsets
# result <- signal_analysis(data[1:500, ], ...)

