# Censored-Covariate Regression with Dual Selection

Fits penalised censored-covariate regression models over a grid of
tuning parameters \\(\lambda_1, \lambda_2)\\ using cross-validation to
evaluate prediction performance. For each pair, an EM algorithm
estimates the model parameters while applying \\L_1\\ penalties to
perform feature selection for both the response-variable model component
and the censored-covariate model component.

## Usage

``` r
ccreds(
  data,
  lambda_1s,
  lambda_2s,
  force_active_beta = integer(0),
  force_active_gamma = integer(0),
  k = 5,
  cores = 1
)
```

## Arguments

- data:

  A data frame with columns in this order: `case` (integer ID), `y`
  (response), `w` (observed value of x or censoring time), `d`
  (censoring indicator: 1 = observed, 0 = censored), `z1`, `z2`, ...,
  `zp` (features).

- lambda_1s:

  Numeric vector of candidate values for \\\lambda_1\\, the penalty on
  \\\pmb{\beta}\\.

- lambda_2s:

  Numeric vector of candidate values for \\\lambda_2\\, the penalty on
  \\\pmb{\gamma}\\.

- force_active_beta:

  Integer vector of feature indices (referring to z column numbers) that
  should not be penalised in \\\pmb{\beta}\\. Default is `integer(0)`.

- force_active_gamma:

  Integer vector of feature indices that should not be penalised in
  \\\pmb{\gamma}\\. Default is `integer(0)`.

- k:

  Number of cross-validation folds. Default is 5.

- cores:

  Number of parallel workers. If greater than 1, uses future and furrr
  for parallel computation. Default is 1.

## Value

A tibble with one row per \\(\lambda_1, \lambda_2)\\ pair and columns:

- lambda_1, lambda_2:

  The tuning parameter values.

- full_fit:

  A list containing the fitted model on the full data, with elements:

  beta_0

  :   Intercept \\\beta_0\\.

  alpha

  :   Coefficient \\\alpha\\ of the censored covariate \\x\\.

  beta

  :   Feature coefficient vector \\\pmb{\beta}\\ of length \\p\\.

  sigma_sq

  :   Error variance \\\sigma^2\\.

  gamma_0

  :   Intercept \\\gamma_0\\.

  gamma

  :   Feature coefficient vector \\\pmb{\gamma}\\ of length \\p\\.

  beta_selected

  :   Number of non-zero entries in \\\hat{\pmb{\beta}}\\ (excluding
      forced-active features).

  gamma_selected

  :   Number of non-zero entries in \\\hat{\pmb{\gamma}}\\ (excluding
      forced-active features).

  iterations

  :   Number of EM iterations until convergence.

- CV_MAE_y:

  Cross-validated mean absolute error for \\y\\.

- CV_MAE_y_SE:

  Standard error of CV_MAE_y across folds.

- CV_MAE_x:

  Cross-validated mean absolute error for \\x\\ (computed on complete
  observations only).

- CV_MAE_x_SE:

  Standard error of CV_MAE_x across folds.

- time:

  Elapsed time for fitting that lambda pair.

## Model

The censored-covariate regression model consists of two components:

**Response-variable model:** \$\$y_i = \beta_0 + \alpha x_i +
\pmb{\beta}^\top \mathbf{z}\_i + \varepsilon_i, \quad \varepsilon_i \sim
N(0, \sigma^2)\$\$

**Censored-covariate model:** \$\$X_i \mid \mathbf{z}\_i \sim
\chi^2(\mu_i), \quad \mu_i = \exp(\gamma_0 + \pmb{\gamma}^\top
\mathbf{z}\_i)\$\$

Here \\\beta_0\\ and \\\gamma_0\\ are intercepts, \\\alpha\\ is the
coefficient of the censored covariate \\x\\, and \\\pmb{\beta}\\ and
\\\pmb{\gamma}\\ are feature coefficient vectors subject to \\L_1\\
penalisation.

## References

Shankar R, Garcia T, Ormerod J, Tarr G (2026). "Feature Selection in
Censored-Covariate Regression Models." *Statistics in Medicine*.

## Examples

``` r
if (FALSE) { # \dontrun{
results <- ccreds(
  data = my_data,
  lambda_1s = c(0.01, 0.1, 1),
  lambda_2s = c(0.01, 0.1, 1),
  k = 5
)
} # }
```
