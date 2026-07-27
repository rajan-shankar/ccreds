# --- helpers to build minimal test fixtures ---

make_test_fit <- function(p = 3) {
  list(
    beta_0 = 1,
    alpha = 0.5,
    beta = c(0.3, 0, -0.2),
    sigma_sq = 1,
    gamma_0 = log(10),
    gamma = c(0.2, 0, 0)
  )
}

make_test_results <- function() {
  fit <- make_test_fit()
  p <- length(fit$beta)

  # Small 2x2 grid
  tibble::tibble(
    lambda_1 = c(0.1, 0.1, 1, 1),
    lambda_2 = c(0.1, 1, 0.1, 1),
    full_fit = lapply(1:4, function(i) {
      f <- fit
      f$beta <- fit$beta * (1 / i)
      f$gamma <- fit$gamma * (1 / i)
      f
    }),
    CV_MAE_y = c(2.1, 2.5, 3.0, 3.5),
    CV_MAE_y_SE = c(0.1, 0.2, 0.3, 0.4),
    CV_MAE_x = c(1.5, 1.8, 2.0, 2.2),
    CV_MAE_x_SE = c(0.1, 0.1, 0.2, 0.2),
    time = c(1, 1, 1, 1)
  )
}

# --- tests ---

test_that("plot_cv returns a ggplot", {
  results <- make_test_results()
  p <- plot_cv(results)
  expect_s3_class(p, "ggplot")
})

test_that("plot_cv errors on invalid metric", {
  results <- make_test_results()
  expect_error(plot_cv(results, metric = "nonexistent"), "not found")
})

test_that("plot_coefficient_paths returns a ggplot for beta and gamma", {
  results <- make_test_results()
  p_beta <- plot_coefficient_paths(results, "beta")
  p_gamma <- plot_coefficient_paths(results, "gamma")
  expect_s3_class(p_beta, "ggplot")
  expect_s3_class(p_gamma, "ggplot")
})

test_that("plot_coefficient_paths accepts feature_names", {
  results <- make_test_results()
  p <- plot_coefficient_paths(results, "beta",
    feature_names = c("feat_a", "feat_b", "feat_c")
  )
  expect_s3_class(p, "ggplot")
})

test_that("plot_coefficient_paths errors on wrong feature_names length", {
  results <- make_test_results()
  expect_error(
    plot_coefficient_paths(results, "beta", feature_names = c("a", "b")),
    "feature_names"
  )
})

test_that("plot_density returns a ggplot for joint and conditional", {
  fit <- make_test_fit()
  z <- c(0.5, -0.3, 1.0)

  p_joint <- plot_density(fit, z)
  expect_s3_class(p_joint, "ggplot")

  p_cond <- plot_density(fit, z, y_given = 5)
  expect_s3_class(p_cond, "ggplot")
})

