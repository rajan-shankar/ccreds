# Plot predictive densities

Plots the joint density \\f(y, x \mid \mathbf{z})\\ as a filled contour,
or, when `y_given` is supplied, overlays the marginal density \\f(x \mid
\mathbf{z})\\ and the conditional density \\f(x \mid \mathbf{z}, y)\\
for a specific value of \\y\\.

## Usage

``` r
plot_density(
  fit,
  z,
  x_range = NULL,
  y_range = NULL,
  y_given = NULL,
  n_grid = 200
)
```

## Arguments

- fit:

  A single fitted model list (e.g. one element of the `full_fit` column
  from
  [`ccreds()`](https://rajan-shankar.github.io/ccreds/reference/ccreds.md)
  output), containing `beta_0`, `alpha`, `beta`, `sigma_sq`, `gamma_0`,
  `gamma`.

- z:

  Numeric vector of feature values (length \\p\\).

- x_range:

  Numeric vector of length 2 giving the range of \\x\\ values to plot.
  Default is computed from the chi-squared distribution implied by the
  model.

- y_range:

  Numeric vector of length 2 giving the range of \\y\\ values (used for
  the joint contour plot). Default is computed from the model.

- y_given:

  Optional scalar. If supplied, produces a density overlay of \\f(x \mid
  \mathbf{z})\\ and \\f(x \mid \mathbf{z}, y = y\_{given})\\ instead of
  the joint contour.

- n_grid:

  Number of grid points per axis. Default is 200.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) { # \dontrun{
results <- ccreds(data, lambda_1s, lambda_2s)
fit <- results$full_fit[[1]]
z_vals <- c(0.5, -0.3, 1.2)

# Joint density contour
plot_density(fit, z_vals)

# Marginal vs conditional density of x given y = 60
plot_density(fit, z_vals, y_given = 60)
} # }
```
