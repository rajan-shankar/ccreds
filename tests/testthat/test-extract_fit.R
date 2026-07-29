make_fit <- function(beta_selected = 3, gamma_selected = 2) {
  list(
    beta_0 = 1, alpha = 0.5,
    beta = c(1, 2, 0), sigma_sq = 1,
    gamma_0 = 0.5, gamma = c(0.1, 0, 0),
    beta_selected = beta_selected,
    gamma_selected = gamma_selected,
    iterations = 10
  )
}

make_results <- function() {
  tibble::tibble(
    lambda_1 = c(0.01, 0.1, 1.0),
    lambda_2 = c(0.01, 0.1, 1.0),
    full_fit = list(
      make_fit(beta_selected = 5, gamma_selected = 3),
      make_fit(beta_selected = 3, gamma_selected = 2),
      make_fit(beta_selected = 1, gamma_selected = 0)
    ),
    CV_MAE_y = c(10.0, 9.0, 9.5),
    CV_MAE_y_SE = c(1.0, 0.8, 0.5),
    CV_MAE_x = c(5, 4, 6),
    CV_MAE_x_SE = c(0.5, 0.4, 0.6),
    time = c(1, 1, 1)
  )
}

test_that("extract_fit with rule='min' returns the minimum CV_MAE_y model", {
  results <- make_results()
  fit <- extract_fit(results, rule = "min")
  expect_type(fit, "list")
  expect_equal(fit$beta_selected, 3)
  expect_equal(fit$gamma_selected, 2)
})

test_that("extract_fit defaults to rule='min'", {
  results <- make_results()
  fit_default <- extract_fit(results)
  fit_min <- extract_fit(results, rule = "min")
  expect_identical(fit_default, fit_min)
})

test_that("extract_fit with rule='1se' picks the most parsimonious eligible model", {
  results <- make_results()
  # min CV_MAE_y = 9.0 at row 2, SE = 0.8, threshold = 9.8
  # Eligible: row 2 (9.0) and row 3 (9.5)
  # Row 3 has fewer selected features (1 + 0 = 1) vs row 2 (3 + 2 = 5)
  fit <- extract_fit(results, rule = "1se")
  expect_equal(fit$beta_selected, 1)
  expect_equal(fit$gamma_selected, 0)
})

test_that("extract_fit with rule='1se' breaks ties by highest penalty", {
  results <- make_results()
  # Make rows 2 and 3 equally parsimonious
  results$full_fit[[2]] <- make_fit(beta_selected = 1, gamma_selected = 0)
  fit <- extract_fit(results, rule = "1se")
  # Both eligible and equally parsimonious; row 3 has higher lambdas
  expect_equal(fit$beta_selected, 1)
})

test_that("extract_fit errors on missing columns", {
  bad <- tibble::tibble(x = 1)
  expect_error(extract_fit(bad), "missing required columns")
})

test_that("extract_fit errors when all CV_MAE_y are NA", {
  results <- make_results()
  results$CV_MAE_y <- NA_real_
  expect_error(extract_fit(results), "No valid models")
})

test_that("extract_fit skips rows with NA CV_MAE_y", {
  results <- make_results()
  results$CV_MAE_y[2] <- NA
  # Min should now be row 3 (CV_MAE_y = 9.5)
  fit <- extract_fit(results, rule = "min")
  expect_equal(fit$beta_selected, 1)
})
