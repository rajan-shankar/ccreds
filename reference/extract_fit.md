# Extract a fitted model from ccreds results

Extracts the fitted model from a
[`ccreds()`](https://rajan-shankar.github.io/ccreds/reference/ccreds.md)
results tibble using either the minimum CV MAE rule or the
one-standard-error (1-SE) rule.

## Usage

``` r
extract_fit(results, rule = c("min", "1se"))
```

## Arguments

- results:

  A tibble returned by
  [`ccreds()`](https://rajan-shankar.github.io/ccreds/reference/ccreds.md).

- rule:

  Character string specifying the selection rule. `"min"` (default)
  selects the \\(\lambda_1, \lambda_2)\\ pair with the smallest
  `CV_MAE_y`. `"1se"` selects the most parsimonious model whose
  `CV_MAE_y` is within one standard error of the minimum — that is,
  among all models with `CV_MAE_y <= min(CV_MAE_y) + CV_MAE_y_SE` at the
  minimum, the model with the fewest total selected features is chosen.

## Value

A named list containing the fitted model parameters:

- beta_0:

  Intercept \\\beta_0\\.

- alpha:

  Coefficient \\\alpha\\ of the censored covariate \\x\\.

- beta:

  Feature coefficient vector \\\pmb{\beta}\\ of length \\p\\.

- sigma_sq:

  Error variance \\\sigma^2\\.

- gamma_0:

  Intercept \\\gamma_0\\.

- gamma:

  Feature coefficient vector \\\pmb{\gamma}\\ of length \\p\\.

- beta_selected:

  Number of non-zero entries in \\\hat{\pmb{\beta}}\\.

- gamma_selected:

  Number of non-zero entries in \\\hat{\pmb{\gamma}}\\.

- iterations:

  Number of EM iterations until convergence.

## Examples

``` r
if (FALSE) { # \dontrun{
results <- ccreds(data, lambda_1s, lambda_2s)

# Best model by minimum CV MAE
fit <- extract_fit(results)

# Most parsimonious model within 1 SE of the minimum
fit_1se <- extract_fit(results, rule = "1se")

fit$beta
fit$gamma
} # }
```
