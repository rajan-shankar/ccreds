#' Extract a fitted model from ccreds results
#'
#' Extracts the fitted model from a [ccreds()] results tibble using
#' either the minimum CV MAE rule or the one-standard-error (1-SE)
#' rule.
#'
#' @param results A tibble returned by [ccreds()].
#' @param rule Character string specifying the selection rule.
#'   `"min"` (default) selects the \eqn{(\lambda_1, \lambda_2)} pair
#'   with the smallest `CV_MAE_y`. `"1se"` selects the most
#'   parsimonious model whose `CV_MAE_y` is within one standard error
#'   of the minimum --- that is, among all models with
#'   `CV_MAE_y <= min(CV_MAE_y) + CV_MAE_y_SE` at the minimum, the
#'   model with the fewest total selected features is chosen.
#'
#' @return A named list containing the fitted model parameters:
#'   \describe{
#'     \item{beta_0}{Intercept \eqn{\beta_0}.}
#'     \item{alpha}{Coefficient \eqn{\alpha} of the censored
#'       covariate \eqn{x}.}
#'     \item{beta}{Feature coefficient vector
#'       \eqn{\pmb{\beta}} of length \eqn{p}.}
#'     \item{sigma_sq}{Error variance \eqn{\sigma^2}.}
#'     \item{gamma_0}{Intercept \eqn{\gamma_0}.}
#'     \item{gamma}{Feature coefficient vector
#'       \eqn{\pmb{\gamma}} of length \eqn{p}.}
#'     \item{beta_selected}{Number of non-zero entries in
#'       \eqn{\hat{\pmb{\beta}}}.}
#'     \item{gamma_selected}{Number of non-zero entries in
#'       \eqn{\hat{\pmb{\gamma}}}.}
#'     \item{iterations}{Number of EM iterations until convergence.}
#'   }
#'
#' @examples
#' \dontrun{
#' results <- ccreds(data, lambda_1s, lambda_2s)
#'
#' # Best model by minimum CV MAE
#' fit <- extract_fit(results)
#'
#' # Most parsimonious model within 1 SE of the minimum
#' fit_1se <- extract_fit(results, rule = "1se")
#'
#' fit$beta
#' fit$gamma
#' }
#'
#' @export
extract_fit <- function(results, rule = c("min", "1se")) {
  rule <- match.arg(rule)

  required_cols <- c(
    "lambda_1", "lambda_2", "full_fit", "CV_MAE_y", "CV_MAE_y_SE"
  )
  missing_cols <- setdiff(required_cols, names(results))
  if (length(missing_cols) > 0) {
    stop(
      "results is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      ". Did you pass a tibble returned by ccreds()?"
    )
  }

  valid <- results[!is.na(results$CV_MAE_y), ]
  if (nrow(valid) == 0) {
    stop("No valid models found in results (all CV_MAE_y values are NA).")
  }

  if (rule == "min") {
    min_idx <- which.min(valid$CV_MAE_y)
    return(valid$full_fit[[min_idx]])
  }

  # 1-SE rule
  min_idx <- which.min(valid$CV_MAE_y)
  threshold <- valid$CV_MAE_y[min_idx] + valid$CV_MAE_y_SE[min_idx]

  eligible <- valid[valid$CV_MAE_y <= threshold, ]

  total_selected <- vapply(
    eligible$full_fit,
    \(f) f$beta_selected + f$gamma_selected,
    numeric(1)
  )
  total_penalty <- eligible$lambda_1 + eligible$lambda_2

  # Sort by: fewest selected features, then lowest CV_MAE_y, then

  # highest total penalty (most regularisation)
  ord <- order(total_selected, eligible$CV_MAE_y, -total_penalty)

  eligible$full_fit[[ord[1]]]
}
