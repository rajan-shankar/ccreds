# ccreds <a href="https://rajan-shankar.github.io/ccreds/"><img src="man/figures/logo.png" align="right" height="138" alt="ccreds website" /></a>

<!-- badges: start -->
<!-- badges: end -->

**ccreds** is an R package for performing Censored-Covariate Regression with Dual Selection. ccreds allows users to fit penalised censored-covariate regression models with two separate L1 penalties --- one for selecting features associated with the response variable and one for selecting features associated with the censored covariate. Tuning parameters are selected via cross-validation. Functions for producing visual diagnostics such as cross-validation heat maps, coefficient paths, and plotting estimated densities, are provided.

## Installation

Install from GitHub with:

```r
# install.packages("pak")
pak::pak("rajan-shankar/ccreds")
```

OR if using Windows:

```r
# install.packages("devtools")
devtools::install_github("rajan-shankar/ccreds")
```

## Usage

The main function is `ccreds()`. Input data must be a data frame with columns in this order:

| Column | Description |
|--------|-------------|
| `case` | Integer observation ID |
| `y` | Response variable |
| `w` | Observed value of the censored covariate, i.e. $x$ if completely observed or $c$ if right-censored
| `d` | Censoring indicator (1 = observed, 0 = censored) |
| `z1`, `z2`, ..., `zp` | Features |

```r
library(ccreds)

results <- ccreds(
  data = my_data,
  lambda_1s = c(0.01, 0.1, 1),
  lambda_2s = c(0.01, 0.1, 1),
  k = 5
)
```

The output is a tibble with one row per (lambda_1, lambda_2) pair, containing cross-validation metrics and the full model fit. Use `extract_fit()` to retrieve the best model:

```r
# Minimum CV MAE rule (default)
fit <- extract_fit(results)

# One-standard-error rule (most parsimonious model within 1 SE of the minimum)
fit <- extract_fit(results, rule = "1se")

fit$beta_0    # Intercept (beta_0)
fit$alpha     # Coefficient of censored covariate x (alpha)
fit$beta      # Feature coefficients for y (beta vector)
fit$sigma_sq  # Error variance (sigma^2)
fit$gamma_0   # Intercept for x model (gamma_0)
fit$gamma     # Feature coefficients for x (gamma vector)
```

For parallel computation over the lambda grid:

```r
results <- ccreds(
  data = my_data,
  lambda_1s = c(0.01, 0.1, 1),
  lambda_2s = c(0.01, 0.1, 1),
  cores = 4
)
```

This requires the `future` and `furrr` packages.
