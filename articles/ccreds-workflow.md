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

set.seed(2026)

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
#>     case      y     w     d     z1      z2     z3       z4     z5      z6    z7
#>    <int>  <dbl> <dbl> <dbl>  <dbl>   <dbl>  <dbl>    <dbl>  <dbl>   <dbl> <dbl>
#>  1     1  -6.33  3.72     0 -3.80  -3.24   -1.85  -0.0231  -3.34  -1.43    4.59
#>  2     2  23.5   3.46     1 -4.53  -0.0690 -3.93   2.98     7.22   5.55    2.89
#>  3     3  -3.79  2.85     0  0.698 -7.17   -6.90  -1.07     2.17  -1.17   -1.55
#>  4     4  -2.67  3.63     1 -5.12  -0.847  -0.182 -4.16    -2.85  -4.43   -4.50
#>  5     5 159.    7.11     0  3.32  -0.365  -1.42  -4.99     0.946  0.200   2.99
#>  6     6  95.5   3.86     0  1.23  -1.26    3.49   0.00796  3.45   0.0594  2.18
#>  7     7   5.01  4.32     1 -3.31   2.51    3.72   4.71     3.38   4.61    3.32
#>  8     8 173.    6.34     0  1.74   0.959   0.681  1.35     4.43   5.00    2.85
#>  9     9 226.   11.2      0  5.11   4.75    2.92  -0.642    0.462 -1.39   -5.84
#> 10    10  54.8   8.13     1  1.18  -1.45    2.00  -0.840   -0.143  2.43    2.77
#> # ℹ 190 more rows
#> # ℹ 23 more variables: z8 <dbl>, z9 <dbl>, z10 <dbl>, z11 <dbl>, z12 <dbl>,
#> #   z13 <dbl>, z14 <dbl>, z15 <dbl>, z16 <dbl>, z17 <dbl>, z18 <dbl>,
#> #   z19 <dbl>, z20 <dbl>, z21 <dbl>, z22 <dbl>, z23 <dbl>, z24 <dbl>,
#> #   z25 <dbl>, z26 <dbl>, z27 <dbl>, z28 <dbl>, z29 <dbl>, z30 <dbl>
```

``` r

# Censoring rate should be approximately 50%
mean(dat$d == 0)
#> [1] 0.46
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
#>  1   0.0001   0.0001 <named list [13]>     38.9        2.64     3.19       0.401
#>  2   0.0001   0.001  <named list [13]>     38.9        2.63     3.19       0.399
#>  3   0.0001   0.01   <named list [13]>     38.7        2.48     3.15       0.380
#>  4   0.0001   0.1    <named list [13]>     37.3        1.85     2.93       0.300
#>  5   0.0001   1      <named list [13]>     40.1        1.28     2.97       0.178
#>  6   0.001    0.0001 <named list [13]>     38.9        2.68     3.18       0.400
#>  7   0.001    0.001  <named list [13]>     38.8        2.67     3.18       0.398
#>  8   0.001    0.01   <named list [13]>     38.7        2.52     3.14       0.378
#>  9   0.001    0.1    <named list [13]>     37.2        1.93     2.91       0.300
#> 10   0.001    1      <named list [13]>     40.1        1.21     2.96       0.178
#> # ℹ 15 more rows
#> # ℹ 1 more variable: time <drtn>

round(results$time, digits = 2)
#> Time differences in secs
#>  [1] 10.50  9.97  9.63  8.85  8.06  8.78  8.70  8.68  7.90  7.19  5.59  5.67
#> [13]  5.49  4.97  4.50  4.01  3.98  3.92  3.23  2.93  4.39  4.21  4.17  3.64
#> [25]  3.18
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
#>  1 z1              4          3.51 
#>  2 z2              4          4.34 
#>  3 z3              2          1.38 
#>  4 z4              0          0    
#>  5 z5              0          0.161
#>  6 z6              0          0    
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
#>  1 z1            0.04         0.0288 
#>  2 z2            0.04         0.0341 
#>  3 z3            0.02         0.0363 
#>  4 z4            0.02         0.00365
#>  5 z5            0            0      
#>  6 z6            0            0      
#>  7 z7            0            0      
#>  8 z8            0            0      
#>  9 z9            0            0      
#> 10 z10           0            0.0126 
#> # ℹ 20 more rows
```

## Visualise: Cross-validation heatmap

[`plot_cv()`](https://rajan-shankar.github.io/ccreds/reference/plot_cv.md)
shows how the CV metric varies over the penalty grid. The red cross
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
