# Plot cross-validation metric heatmap

Displays a heatmap of a cross-validation metric over the \\(\lambda_1,
\lambda_2)\\ tuning parameter grid. The pair that minimises the chosen
metric is marked.

## Usage

``` r
plot_cv(results, metric = "CV_MAE_y")
```

## Arguments

- results:

  A tibble returned by
  [`ccreds()`](https://rajan-shankar.github.io/ccreds/reference/ccreds.md).

- metric:

  Character string naming the column to plot. Default is `"CV_MAE_y"`.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) { # \dontrun{
results <- ccreds(data, lambda_1s, lambda_2s)
plot_cv(results)
plot_cv(results, metric = "CV_MAE_x")
} # }
```
