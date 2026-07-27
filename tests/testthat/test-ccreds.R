test_that("ccreds runs and returns expected structure", {
  set.seed(4781)
  n <- 60
  p <- 5

  z <- matrix(rnorm(n * p), nrow = n)
  colnames(z) <- paste0("z", 1:p)

  gamma_0_true <- log(10)
  gamma_true <- c(0.3, -0.2, rep(0, p - 2))
  mu <- exp(gamma_0_true + z %*% gamma_true)
  x <- rchisq(n, df = mu)

  beta_0_true <- 1
  alpha_true <- 0.5
  beta_true <- c(0.4, 0, 0, -0.3, 0)
  y <- beta_0_true + alpha_true * x + z %*% beta_true + rnorm(n, sd = 0.5)

  # Introduce some censoring
  c_vals <- rexp(n, rate = 1 / mean(x))
  d <- as.integer(x <= c_vals)
  w <- ifelse(d == 1, x, c_vals)

  data <- data.frame(case = 1:n, y = as.numeric(y), w = w, d = d, z)

  results <- ccreds(
    data = data,
    lambda_1s = c(0.1),
    lambda_2s = c(0.1),
    k = 2
  )

  expect_s3_class(results, "tbl_df")
  expect_named(
    results,
    c(
      "lambda_1", "lambda_2", "full_fit",
      "CV_MAE_y", "CV_MAE_y_SE", "CV_MAE_x", "CV_MAE_x_SE", "time"
    )
  )

  fit <- results$full_fit[[1]]
  expect_type(fit, "list")
  expect_named(
    fit[c(
      "beta_0", "alpha", "beta", "sigma_sq",
      "gamma_0", "gamma", "beta_selected", "gamma_selected"
    )],
    c(
      "beta_0", "alpha", "beta", "sigma_sq",
      "gamma_0", "gamma", "beta_selected", "gamma_selected"
    )
  )
  expect_length(fit$beta, p)
  expect_length(fit$gamma, p)
})

test_that("ccreds rejects empty lambda grids", {
  data <- data.frame(
    case = 1:10, y = rnorm(10), w = rexp(10),
    d = c(rep(1, 5), rep(0, 5)),
    z1 = rnorm(10)
  )
  expect_error(
    ccreds(data = data, lambda_1s = numeric(0), lambda_2s = c(0.1)),
    "lambda_1 AND lambda_2"
  )
})
