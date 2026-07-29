# Getting Started with ccreds

This vignette walks through a complete `ccreds` workflow: simulating
censored-covariate regression data, fitting the model with
cross-validated penalty selection, and visualising the results.

## Setup

``` r

library(ccreds)
library(MASS)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following object is masked from 'package:MASS':
#> 
#>     select
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
```

## Simulate data

We generate data from the censored-covariate regression model with \\n =
200\\ observations and \\p = 30\\ features.

**Response-variable model:** \\y_i = \beta_0 + \alpha x_i +
\pmb{\beta}^\top \mathbf{z}\_i + \varepsilon_i, \quad \varepsilon_i \sim
N(0, \sigma^2)\\

**Censored-covariate model:** \\X_i \mid \mathbf{z}\_i \sim
\chi^2(\mu_i), \quad \mu_i = \exp(\gamma_0 + \pmb{\gamma}^\top
\mathbf{z}\_i)\\

We set \\\pmb{\beta}\\ and \\\pmb{\gamma}\\ to be sparse with partially
overlapping support, so that some features affect only the response,
some affect only the censored covariate, and some affect both. We set
the censoring rate to be approximately 50%.

``` r

set.seed(8371)

n <- 200
p <- 30

# Correlated features (AR(0.5) structure)
rho <- 0.5
Sigma <- 9*toeplitz(rho^(0:(p - 1)))
z <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)

# True parameters
beta_0 <- 5
alpha <- 10
beta <- c(4, 4, 2, rep(0, p-5), -2, -4)
sigma_sq <- 100

gamma_0 <- 2
gamma <- c(0.04, 0.04, 0.02, 0.02, rep(0, p-4))

# Generate the censored covariate x
x_mean <- exp(cbind(1, z) %*% c(gamma_0, gamma)) |> drop()
x <- rchisq(n, df = x_mean)

# Generate censoring times to achieve ~50% censoring
compute_c_mean <- function(censoring_rate, x_mean) {
  f <- function(c_mean) {
    1 - pbeta(0.5, x_mean / 2, c_mean / 2) - censoring_rate
  }
  uniroot(f, lower = 1e-12, upper = 1000)$root
}
c_mean <- vapply(x_mean, \(m) compute_c_mean(0.5, m), numeric(1))
cens <- rchisq(n, df = c_mean)

# Generate response
e <- rnorm(n, mean = 0, sd = sqrt(sigma_sq))
y <- cbind(1, x, z) %*% c(beta_0, alpha, beta) + e

# Observed covariate and censoring indicator
w <- pmin(x, cens)
d <- as.numeric(x <= cens)

# Assemble the data frame in the format ccreds() expects
z_df <- as.data.frame(z)
colnames(z_df) <- paste0("z", 1:p)
dat <- tibble(
  case = 1:n,
  y = drop(y),
  w = w,
  d = d
) |>
  bind_cols(z_df)

dat
#> # A tibble: 200 × 34
#>     case     y     w     d        z1     z2      z3        z4      z5    z6
#>    <int> <dbl> <dbl> <dbl>     <dbl>  <dbl>   <dbl>     <dbl>   <dbl> <dbl>
#>  1     1  70.9  4.58     0 -3.65     -1.37   5.27   10.2       4.54    7.54
#>  2     2 157.  12.6      0  2.78      4.35   3.78    5.21      7.00    1.56
#>  3     3 156.  12.0      1  1.80     -0.182 -1.09   -3.97     -2.81    1.44
#>  4     4  26.4  4.13     1  0.227    -1.15  -1.44   -6.43     -5.97   -1.22
#>  5     5  45.3  5.27     1 -1.56     -2.30   0.606   0.970     6.37    5.09
#>  6     6 122.   2.02     0 -0.000930 -2.10   0.476  -2.46      5.04    4.09
#>  7     7  24.9  1.37     0 -3.73     -4.34  -1.58   -0.000185  0.0260  3.69
#>  8     8 128.   7.74     1 -0.0911    3.58   2.09    2.92      5.13    1.11
#>  9     9  96.6  7.02     1 -0.101     2.45   2.64    1.69      1.21    2.16
#> 10    10 113.   5.54     0 -4.11     -1.80   0.0936 -4.47     -1.48   -1.56
#> # ℹ 190 more rows
#> # ℹ 24 more variables: z7 <dbl>, z8 <dbl>, z9 <dbl>, z10 <dbl>, z11 <dbl>,
#> #   z12 <dbl>, z13 <dbl>, z14 <dbl>, z15 <dbl>, z16 <dbl>, z17 <dbl>,
#> #   z18 <dbl>, z19 <dbl>, z20 <dbl>, z21 <dbl>, z22 <dbl>, z23 <dbl>,
#> #   z24 <dbl>, z25 <dbl>, z26 <dbl>, z27 <dbl>, z28 <dbl>, z29 <dbl>, z30 <dbl>
```

``` r

# Censoring rate should be approximately 50%
mean(dat$d == 0)
#> [1] 0.475
```

## Fit the model

We search over a small grid of penalty parameters. In practice you would
use a finer grid and potentially run in parallel with `cores > 1`.

``` r

lambda_1s <- c(0.0001, 0.001, 0.01, 0.1, 1)
lambda_2s <- c(0.0001, 0.001, 0.01, 0.1, 1)

results <- ccreds(
  data = dat,
  lambda_1s = lambda_1s,
  lambda_2s = lambda_2s,
  k = 5
)

results
#> # A tibble: 25 × 8
#>    lambda_1 lambda_2 full_fit          CV_MAE_y CV_MAE_y_SE CV_MAE_x CV_MAE_x_SE
#>       <dbl>    <dbl> <list>               <dbl>       <dbl>    <dbl>       <dbl>
#>  1   0.0001   0.0001 <named list [13]>     32.7        1.77     2.92       0.267
#>  2   0.0001   0.001  <named list [13]>     32.6        1.76     2.92       0.265
#>  3   0.0001   0.01   <named list [13]>     32.1        1.58     2.88       0.243
#>  4   0.0001   0.1    <named list [13]>     29.8        1.52     2.70       0.160
#>  5   0.0001   1      <named list [13]>     33.4        1.23     3.00       0.204
#>  6   0.001    0.0001 <named list [13]>     32.7        1.79     2.92       0.268
#>  7   0.001    0.001  <named list [13]>     32.6        1.77     2.92       0.266
#>  8   0.001    0.01   <named list [13]>     32.0        1.59     2.88       0.244
#>  9   0.001    0.1    <named list [13]>     29.7        1.50     2.69       0.160
#> 10   0.001    1      <named list [13]>     33.4        1.22     3.00       0.203
#> # ℹ 15 more rows
#> # ℹ 1 more variable: time <drtn>

round(results$time, digits = 2)
#> Time differences in secs
#>  [1] 7.30 6.89 6.78 6.10 5.90 5.85 5.88 5.92 5.25 5.27 4.01 4.08 3.94 3.60 3.83
#> [16] 2.91 2.91 2.91 2.66 2.72 2.92 2.76 2.81 2.58 2.55
```

## Identify the best model

Use
[`extract_fit()`](https://rajan-shankar.github.io/ccreds/reference/extract_fit.md)
to retrieve the fitted model. The default `rule = "min"` selects the
\\(\lambda_1, \lambda_2)\\ pair that minimises `CV_MAE_y`. The
`rule = "1se"` alternative selects the most parsimonious model whose
`CV_MAE_y` is within one standard error of the minimum.

``` r

fit <- extract_fit(results)
```

We can inspect the first few estimated coefficients and compare them to
the truth:

``` r

# Response-model coefficients
tibble(
  feature = paste0("z", 1:p),
  true_beta = beta,
  estimated_beta = round(fit$beta, 3)
)
#> # A tibble: 30 × 3
#>    feature true_beta estimated_beta
#>    <chr>       <dbl>          <dbl>
#>  1 z1              4          3.23 
#>  2 z2              4          4.34 
#>  3 z3              2          1.54 
#>  4 z4              0          0    
#>  5 z5              0          0    
#>  6 z6              0          0.329
#>  7 z7              0          0    
#>  8 z8              0          0    
#>  9 z9              0          0    
#> 10 z10             0          0    
#> # ℹ 20 more rows

# Censored-covariate coefficients
tibble(
  feature = paste0("z", 1:p),
  true_gamma = gamma,
  estimated_gamma = round(fit$gamma, 5)
)
#> # A tibble: 30 × 3
#>    feature true_gamma estimated_gamma
#>    <chr>        <dbl>           <dbl>
#>  1 z1            0.04         0.0297 
#>  2 z2            0.04         0.0418 
#>  3 z3            0.02         0.0341 
#>  4 z4            0.02         0      
#>  5 z5            0            0      
#>  6 z6            0            0      
#>  7 z7            0           -0.00484
#>  8 z8            0            0      
#>  9 z9            0            0      
#> 10 z10           0            0      
#> # ℹ 20 more rows
```

## Visualise: Cross-validation heatmap

[`plot_cv()`](https://rajan-shankar.github.io/ccreds/reference/plot_cv.md)
shows how the CV metric varies over the penalty grid. The white cross
marks the optimal pair.

``` r

plot_cv(results)
```

![Heatmap of cross-validated MAE over the lambda
grid.](ccreds-workflow_files/figure-html/plot-cv-1.png)

## Visualise: Coefficient paths

[`plot_coefficient_paths()`](https://rajan-shankar.github.io/ccreds/reference/plot_coefficient_paths.md)
shows how the estimated coefficients change as the penalty parameter
varies, holding the other penalty fixed at the value from the best
model.

``` r

plot_coefficient_paths(results, component = "beta", highlight = TRUE)
```

![Coefficient paths for beta as lambda_1
varies.](ccreds-workflow_files/figure-html/coef-paths-beta-1.png)

``` r

plot_coefficient_paths(results, component = "gamma", highlight = TRUE)
```

![Coefficient paths for gamma as lambda_2
varies.](ccreds-workflow_files/figure-html/coef-paths-gamma-1.png)

## Visualise: Predictive densities

[`plot_density()`](https://rajan-shankar.github.io/ccreds/reference/plot_density.md)
displays the joint or conditional density of \\(y, x)\\ for a given
feature vector \\\mathbf{z}\\.

``` r

z_new <- rep(0, p)
plot_density(fit, z = z_new, x_range = c(0,20), y_range = c(0, 200))
```

![Joint density contour plot of y and
x.](ccreds-workflow_files/figure-html/density-joint-1.png)

Supplying `y_given` shows the marginal density \\f(x \mid \mathbf{z})\\
alongside the conditional density \\f(x \mid \mathbf{z}, y)\\:

``` r

plot_density(fit, z = z_new, y_given = 50)
```

![Marginal vs conditional density of x given z and
y.](ccreds-workflow_files/figure-html/density-conditional-1.png)

## Summary

The typical `ccreds` workflow is:

1.  **Prepare data** in the required format: `case`, `y`, `w`, `d`,
    `z1`, …, `zp`.
2.  **Fit** with
    [`ccreds()`](https://rajan-shankar.github.io/ccreds/reference/ccreds.md)
    over a grid of \\(\lambda_1, \lambda_2)\\ values.
3.  **Select** the best model with
    [`extract_fit()`](https://rajan-shankar.github.io/ccreds/reference/extract_fit.md).
4.  **Visualise** with
    [`plot_cv()`](https://rajan-shankar.github.io/ccreds/reference/plot_cv.md),
    [`plot_coefficient_paths()`](https://rajan-shankar.github.io/ccreds/reference/plot_coefficient_paths.md),
    and
    [`plot_density()`](https://rajan-shankar.github.io/ccreds/reference/plot_density.md).
