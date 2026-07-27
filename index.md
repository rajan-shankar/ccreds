# ccreds

Feature selection in censored-covariate regression models, as described
in:

> Shankar R, Garcia T, Ormerod J, Tarr G (2026). “Feature Selection in
> Censored-Covariate Regression Models.” *Statistics in Medicine*.

The package fits a penalised censored-covariate regression model using
an EM algorithm with two separate L1 penalties — one for features
associated with the response variable and one for features associated
with the censored covariate. Tuning parameters are selected via
cross-validation.

## Installation

Install from GitHub with:

``` r

# install.packages("pak")
pak::pak("rajan-shankar/ccreds")
```

OR if using Windows:

``` r

# install.packages("devtools")
devtools::install_github("rajan-shankar/ccreds")
```

## Usage

The main function is
[`ccreds()`](https://rajan-shankar.github.io/ccreds/reference/ccreds.md).
Input data must be a data frame with columns in this order:

| Column              | Description                                         |
|---------------------|-----------------------------------------------------|
| `case`              | Integer observation ID                              |
| `y`                 | Response variable                                   |
| `w`                 | Observed value of x (or censoring time if censored) |
| `d`                 | Censoring indicator (1 = observed, 0 = censored)    |
| `z1`, `z2`, …, `zp` | Features                                            |

``` r

library(ccreds)

results <- ccreds(
  data = my_data,
  lambda_1s = c(0.01, 0.1, 1),
  lambda_2s = c(0.01, 0.1, 1),
  k = 5
)
```

The output is a tibble with one row per (lambda_1, lambda_2) pair,
containing cross-validation metrics and the full model fit. To extract
the best model:

``` r

best <- results |> dplyr::slice_min(CV_MAE_y)
fit <- best$full_fit[[1]]

fit$beta_0    # Intercept (beta_0)
fit$alpha     # Coefficient of censored covariate x (alpha)
fit$beta      # Feature coefficients for y (beta vector)
fit$sigma_sq  # Error variance (sigma^2)
fit$gamma_0   # Intercept for x model (gamma_0)
fit$gamma     # Feature coefficients for x (gamma vector)
```

For parallel computation over the lambda grid:

``` r

results <- ccreds(
  data = my_data,
  lambda_1s = c(0.01, 0.1, 1),
  lambda_2s = c(0.01, 0.1, 1),
  cores = 4
)
```

This requires the `future` and `furrr` packages.
