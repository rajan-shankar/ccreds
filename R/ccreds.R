#' Censored-Covariate Regression with Dual Selection
#'
#' Fits penalised censored-covariate regression models over a grid of
#' tuning parameters \eqn{(\lambda_1, \lambda_2)} using cross-validation
#' to evaluate prediction performance. For each pair, an EM algorithm
#' estimates the model parameters while applying \eqn{L_1} penalties
#' to perform feature selection for both the response-variable model
#' component and the censored-covariate model component.
#'
#' @section Model:
#' The censored-covariate regression model consists of two components:
#'
#' \strong{Response-variable model:}
#' \deqn{y_i = \beta_0 + \alpha x_i + \pmb{\beta}^\top
#'   \mathbf{z}_i + \varepsilon_i, \quad \varepsilon_i \sim N(0,
#'   \sigma^2)}
#'
#' \strong{Censored-covariate model:}
#' \deqn{X_i \mid \mathbf{z}_i \sim \chi^2(\mu_i), \quad \mu_i =
#'   \exp(\gamma_0 + \pmb{\gamma}^\top \mathbf{z}_i)}
#'
#' Here \eqn{\beta_0} and \eqn{\gamma_0} are intercepts, \eqn{\alpha}
#' is the coefficient of the censored covariate \eqn{x}, and
#' \eqn{\pmb{\beta}} and \eqn{\pmb{\gamma}} are feature
#' coefficient vectors subject to \eqn{L_1} penalisation.
#'
#' @param data A data frame with columns in this order: `case` (integer
#'   ID), `y` (response), `w` (observed value of x or censoring time),
#'   `d` (censoring indicator: 1 = observed, 0 = censored), `z1`,
#'   `z2`, ..., `zp` (features).
#' @param lambda_1s Numeric vector of candidate values for
#'   \eqn{\lambda_1}, the penalty on \eqn{\pmb{\beta}}.
#' @param lambda_2s Numeric vector of candidate values for
#'   \eqn{\lambda_2}, the penalty on \eqn{\pmb{\gamma}}.
#' @param force_active_beta Integer vector of feature indices (referring
#'   to z column numbers) that should not be penalised in
#'   \eqn{\pmb{\beta}}. Default is `integer(0)`.
#' @param force_active_gamma Integer vector of feature indices that
#'   should not be penalised in \eqn{\pmb{\gamma}}. Default is
#'   `integer(0)`.
#' @param k Number of cross-validation folds. Default is 5.
#' @param cores Number of parallel workers. If greater than 1, uses
#'   \pkg{future} and \pkg{furrr} for parallel computation. Default
#'   is 1.
#'
#' @return A tibble with one row per \eqn{(\lambda_1, \lambda_2)} pair
#'   and columns:
#'   \describe{
#'     \item{lambda_1, lambda_2}{The tuning parameter values.}
#'     \item{full_fit}{A list containing the fitted model on the full
#'       data, with elements:
#'       \describe{
#'         \item{beta_0}{Intercept \eqn{\beta_0}.}
#'         \item{alpha}{Coefficient \eqn{\alpha} of the censored
#'           covariate \eqn{x}.}
#'         \item{beta}{Feature coefficient vector
#'           \eqn{\pmb{\beta}} of length \eqn{p}.}
#'         \item{sigma_sq}{Error variance \eqn{\sigma^2}.}
#'         \item{gamma_0}{Intercept \eqn{\gamma_0}.}
#'         \item{gamma}{Feature coefficient vector
#'           \eqn{\pmb{\gamma}} of length \eqn{p}.}
#'         \item{beta_selected}{Number of non-zero entries in
#'           \eqn{\hat{\pmb{\beta}}} (excluding forced-active
#'           features).}
#'         \item{gamma_selected}{Number of non-zero entries in
#'           \eqn{\hat{\pmb{\gamma}}} (excluding forced-active
#'           features).}
#'         \item{iterations}{Number of EM iterations until
#'           convergence.}
#'       }
#'     }
#'     \item{CV_MAE_y}{Cross-validated mean absolute error for
#'       \eqn{y}.}
#'     \item{CV_MAE_y_SE}{Standard error of CV_MAE_y across folds.}
#'     \item{CV_MAE_x}{Cross-validated mean absolute error for
#'       \eqn{x} (computed on complete observations only).}
#'     \item{CV_MAE_x_SE}{Standard error of CV_MAE_x across folds.}
#'     \item{time}{Elapsed time for fitting that lambda pair.}
#'   }
#'
#' @references
#' Shankar R, Garcia T, Ormerod J, Tarr G (2026). "Feature Selection
#' in Censored-Covariate Regression Models." \emph{Statistics in
#' Medicine}.
#'
#' @examples
#' \dontrun{
#' results <- ccreds(
#'   data = my_data,
#'   lambda_1s = c(0.01, 0.1, 1),
#'   lambda_2s = c(0.01, 0.1, 1),
#'   k = 5
#' )
#' }
#'
#' @export
ccreds <- function(data,
                   lambda_1s = 10^seq(-4, 0, length.out = 9),
                   lambda_2s = 10^seq(-4, 0, length.out = 9),
                   force_active_beta = integer(0),
                   force_active_gamma = integer(0),
                   k = 5,
                   cores = 1) {
  # Enforce minimum w for censored obs
  min_w <- data |>
    dplyr::filter(.data$d == 0) |>
    dplyr::pull(.data$w) |>
    min()
  if (min_w < 0.01) {
    data <- data |>
      dplyr::mutate(w = dplyr::case_when(
        .data$d == 0 & .data$w < 0.01 ~ 0.01,
        TRUE ~ .data$w
      ))
    message(
      "some censored w values are < 0.01; ",
      "these have been changed to 0.01"
    )
  }

  J1 <- length(lambda_1s)
  J2 <- length(lambda_2s)

  if (any(c(J1, J2) < 1)) {
    stop("Number of values for lambda_1 AND lambda_2 must be >= 1")
  }

  p <- data |> dplyr::select(dplyr::starts_with("z")) |> ncol()

  folds <- get_folds(data, k)

  order_fit1 <- rep(1:J1, each = J2)
  order_fit2 <- rep(1:J2, times = J1)

  results <- tibble::tibble(
    lambda_1 = lambda_1s[order_fit1],
    lambda_2 = lambda_2s[order_fit2]
  )

  lambda_pair_fn <- function(lambda_1, lambda_2) {
    t1 <- Sys.time()

    error_output <- list(
      full_fit = NA,
      CV_MAE_y = NA,
      CV_MAE_y_SE = NA,
      CV_MAE_x = NA,
      CV_MAE_x_SE = NA,
      time = NA
    )

    full_fit <- tryCatch(
      fit_model(
        init_alpha = 1,
        init_beta_0 = 0,
        init_beta = rep(0, p),
        init_sigma_sq = 1,
        init_gamma_0 = 0,
        init_gamma = rep(0, p),
        data = data,
        lambda_1 = lambda_1,
        lambda_2 = lambda_2,
        force_active_beta = force_active_beta,
        force_active_gamma = force_active_gamma
      ),
      error = function(e) e
    )

    if (inherits(full_fit, "error")) {
      error_output$full_fit <- list(
        occurred = "main full_fit",
        error = full_fit
      )
      return(error_output)
    }

    CV_output <- tibble::tibble(
      fold = numeric(),
      y_errors = list(),
      x_errors = list()
    )

    fit <- NULL
    for (i in 1:k) {
      train_data <- data |> dplyr::slice(purrr::list_c(folds[-i]))
      test_data <- data |>
        dplyr::slice(folds[[i]]) |>
        dplyr::arrange(dplyr::desc(.data$d))

      fit <- tryCatch(
        fit_model(
          init_alpha = 1,
          init_beta_0 = 0,
          init_beta = rep(0, p),
          init_sigma_sq = 1,
          init_gamma_0 = 0,
          init_gamma = rep(0, p),
          data = train_data,
          lambda_1 = lambda_1,
          lambda_2 = lambda_2,
          force_active_beta = force_active_beta,
          force_active_gamma = force_active_gamma
        ),
        error = function(e) e
      )

      if (inherits(fit, "error")) break

      # Predict x (chi-squared mean)
      test_z_complete <- test_data |>
        dplyr::filter(.data$d == 1) |>
        dplyr::select(dplyr::starts_with("z")) |>
        as.matrix()
      test_z_censored <- test_data |>
        dplyr::filter(.data$d == 0) |>
        dplyr::select(dplyr::starts_with("z")) |>
        as.matrix()

      x_pred_complete <- exp(drop(
        cbind(1, test_z_complete) %*% c(fit$gamma_0, fit$gamma)
      ))
      x_pred_censored <- exp(drop(
        cbind(1, test_z_censored) %*% c(fit$gamma_0, fit$gamma)
      ))

      # Predict y
      test_z_all <- test_data |>
        dplyr::select(dplyr::starts_with("z")) |>
        as.matrix()
      x_pred_all <- c(x_pred_complete, x_pred_censored)
      y_pred <- drop(
        cbind(1, x_pred_all, test_z_all) %*%
          c(fit$beta_0, fit$alpha, fit$beta)
      )

      y_errors <- y_pred - test_data$y
      x_errors <- x_pred_complete - test_data$w[test_data$d == 1]
      CV_output <- CV_output |>
        tibble::add_row(
          fold = i,
          y_errors = list(y_errors),
          x_errors = list(x_errors)
        )
    }

    if (inherits(fit, "error")) {
      error_output$full_fit <- list(occurred = "CV loop", error = fit)
      return(error_output)
    }

    CV_MAE_y <- mean(
      purrr::map_dbl(CV_output$y_errors, \(e) mean(abs(e)))
    )
    CV_MAE_y_SE <- sd(
      purrr::map_dbl(CV_output$y_errors, \(e) mean(abs(e)))
    ) / sqrt(k)
    CV_MAE_x <- mean(
      purrr::map_dbl(CV_output$x_errors, \(e) mean(abs(e)))
    )
    CV_MAE_x_SE <- sd(
      purrr::map_dbl(CV_output$x_errors, \(e) mean(abs(e)))
    ) / sqrt(k)

    t2 <- Sys.time()

    list(
      full_fit = full_fit,
      CV_MAE_y = CV_MAE_y,
      CV_MAE_y_SE = CV_MAE_y_SE,
      CV_MAE_x = CV_MAE_x,
      CV_MAE_x_SE = CV_MAE_x_SE,
      time = t2 - t1
    )
  }

  # Fit over lambda grid
  if (cores > 1) {
    if (!requireNamespace("future", quietly = TRUE) ||
      !requireNamespace("furrr", quietly = TRUE)) {
      stop(
        "Packages 'future' and 'furrr' are required for ",
        "parallel computation. Install them with:\n",
        "  install.packages(c(\"future\", \"furrr\"))"
      )
    }
    future::plan(future::multisession, workers = cores)
    on.exit(future::plan(future::sequential), add = TRUE)

    results <- results |>
      dplyr::mutate(result = furrr::future_map2(
        .data$lambda_1, .data$lambda_2,
        .f = lambda_pair_fn,
        .progress = TRUE,
        .options = furrr::furrr_options(seed = TRUE)
      )) |>
      tidyr::unnest_wider("result")
  } else {
    results <- results |>
      dplyr::mutate(result = purrr::map2(
        .data$lambda_1, .data$lambda_2,
        .f = lambda_pair_fn,
        .progress = TRUE
      )) |>
      tidyr::unnest_wider("result")
  }

  results
}
