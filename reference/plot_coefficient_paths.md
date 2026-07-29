# Plot coefficient paths

Shows how the estimated feature coefficients \\\pmb{\beta}\\ or
\\\pmb{\gamma}\\ change as the corresponding penalty parameter varies.
For `component = "beta"`, the plot varies \\\lambda_1\\ while holding
\\\lambda_2\\ fixed, and vice versa for `component = "gamma"`.

## Usage

``` r
plot_coefficient_paths(
  results,
  component = c("beta", "gamma"),
  lambda_fixed = NULL,
  feature_names = NULL,
  highlight = TRUE
)
```

## Arguments

- results:

  A tibble returned by
  [`ccreds()`](https://rajan-shankar.github.io/ccreds/reference/ccreds.md).

- component:

  Either `"beta"` or `"gamma"`.

- lambda_fixed:

  The value of the other penalty to hold constant while the plotted
  penalty varies. When `component = "beta"`, this fixes \\\lambda_2\\;
  when `component = "gamma"`, this fixes \\\lambda_1\\. If `NULL`
  (default), uses the value from the minimum `CV_MAE_y` model.

- feature_names:

  Optional character vector of length \\p\\ giving display names for the
  features. If `NULL`, uses `z1`, `z2`, etc.

- highlight:

  Logical; if `TRUE` (default), highlights features that are non-zero at
  the minimum-CV model and adds a dashed vertical line at the chosen
  penalty value.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) { # \dontrun{
results <- ccreds(data, lambda_1s, lambda_2s)
plot_coefficient_paths(results, "beta")
plot_coefficient_paths(results, "gamma",
  feature_names = c("age", "weight", "height")
)
} # }
```
