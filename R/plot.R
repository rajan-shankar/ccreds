#' Plot cross-validation metric heatmap
#'
#' Displays a heatmap of a cross-validation metric over the
#' \eqn{(\lambda_1, \lambda_2)} tuning parameter grid. The pair
#' that minimises the chosen metric is marked.
#'
#' @param results A tibble returned by [ccreds()].
#' @param metric Character string naming the column to plot.
#'   Default is `"CV_MAE_y"`.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' results <- ccreds(data, lambda_1s, lambda_2s)
#' plot_cv(results)
#' plot_cv(results, metric = "CV_MAE_x")
#' }
#'
#' @export
plot_cv <- function(results, metric = "CV_MAE_y") {
  if (!metric %in% names(results)) {
    stop("Column '", metric, "' not found in results.")
  }

  best <- results[which.min(results[[metric]]), ]

  lambda_1_label <- lambda_label("lambda_1")
  lambda_2_label <- lambda_label("lambda_2")

  ggplot2::ggplot(results, ggplot2::aes(
    x = .data$lambda_1,
    y = .data$lambda_2
  )) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data[[metric]])) +
    ggplot2::geom_point(
      data = best,
      shape = 4, size = 3, stroke = 1.2, colour = "red"
    ) +
    ggplot2::scale_x_log10(expand = c(0, 0)) +
    ggplot2::scale_y_log10(expand = c(0, 0)) +
    ggplot2::scale_fill_viridis_c(direction = -1) +
    ggplot2::labs(
      x = lambda_1_label,
      y = lambda_2_label,
      fill = metric
    ) +
    ggplot2::theme_bw()
}


#' Plot coefficient paths
#'
#' Shows how the estimated feature coefficients
#' \eqn{\pmb{\beta}} or \eqn{\pmb{\gamma}} change
#' as the corresponding penalty parameter varies. For
#' `component = "beta"`, the plot varies \eqn{\lambda_1} while
#' holding \eqn{\lambda_2} fixed, and vice versa for
#' `component = "gamma"`.
#'
#' @param results A tibble returned by [ccreds()].
#' @param component Either `"beta"` or `"gamma"`.
#' @param lambda_fixed The value of the other penalty to hold
#'   constant while the plotted penalty varies. When
#'   `component = "beta"`, this fixes \eqn{\lambda_2}; when
#'   `component = "gamma"`, this fixes \eqn{\lambda_1}. If `NULL`
#'   (default), uses the value from the minimum `CV_MAE_y` model.
#' @param feature_names Optional character vector of length
#'   \eqn{p} giving display names for the features. If `NULL`,
#'   uses `z1`, `z2`, etc.
#' @param highlight Logical; if `TRUE` (default), highlights
#'   features that are non-zero at the minimum-CV model and adds
#'   a dashed vertical line at the chosen penalty value.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' results <- ccreds(data, lambda_1s, lambda_2s)
#' plot_coefficient_paths(results, "beta")
#' plot_coefficient_paths(results, "gamma",
#'   feature_names = c("age", "weight", "height")
#' )
#' }
#'
#' @export
plot_coefficient_paths <- function(results,
                                   component = c("beta", "gamma"),
                                   lambda_fixed = NULL,
                                   feature_names = NULL,
                                   highlight = TRUE) {
  component <- match.arg(component)

  best_row <- results[which.min(results$CV_MAE_y), ]

  if (component == "beta") {
    fixed_col <- "lambda_2"
    varying_col <- "lambda_1"
    if (is.null(lambda_fixed)) lambda_fixed <- best_row$lambda_2
    best_lambda <- best_row$lambda_1
  } else {
    fixed_col <- "lambda_1"
    varying_col <- "lambda_2"
    if (is.null(lambda_fixed)) lambda_fixed <- best_row$lambda_1
    best_lambda <- best_row$lambda_2
  }

  # Filter to rows with the fixed lambda value
  subset <- results[results[[fixed_col]] == lambda_fixed, ]

  if (nrow(subset) == 0) {
    stop(
      "No rows found with ", fixed_col, " == ", lambda_fixed,
      ". Available values: ",
      paste(sort(unique(results[[fixed_col]])), collapse = ", ")
    )
  }

  # Extract coefficients from each fit
  coefs <- lapply(subset$full_fit, function(fit) {
    if (is.list(fit) && component %in% names(fit)) {
      fit[[component]]
    } else {
      NA
    }
  })

  # Drop rows where fit failed

  valid <- !is.na(coefs)
  subset <- subset[valid, ]
  coefs <- coefs[valid]

  if (nrow(subset) == 0) {
    stop("No valid fitted models found for the chosen lambda slice.")
  }

  p <- length(coefs[[1]])
  if (is.null(feature_names)) {
    feature_names <- paste0("z", seq_len(p))
  }
  if (length(feature_names) != p) {
    stop(
      "feature_names has length ", length(feature_names),
      " but there are ", p, " features."
    )
  }

  # Build long-format data
  coef_mat <- do.call(rbind, coefs)
  colnames(coef_mat) <- feature_names

  path_data <- data.frame(
    lambda = subset[[varying_col]],
    coef_mat,
    check.names = FALSE
  )
  path_long <- tidyr::pivot_longer(
    path_data,
    cols = -"lambda",
    names_to = "feature",
    values_to = "value"
  )
  # Identify which features are selected at the best lambda
  best_fit <- best_row$full_fit[[1]]
  best_coefs <- best_fit[[component]]
  selected <- feature_names[best_coefs != 0]
  path_long$selected <- path_long$feature %in% selected

  # Order legend by coefficient value (highest to lowest)
  feature_order <- feature_names[order(best_coefs, decreasing = TRUE)]
  path_long$feature <- factor(path_long$feature, levels = feature_order)

  lambda_label_text <- lambda_label(varying_col)

  plt <- ggplot2::ggplot(
    path_long,
    ggplot2::aes(
      x = .data$lambda,
      y = .data$value,
      colour = .data$feature,
      group = .data$feature
    )
  )

  component_label <- ifelse(
    component == "beta",
    expression(hat(beta)),
    expression(hat(gamma))
  )

  if (highlight) {
    # Unselected lines: grey and excluded from legend
    plt <- plt +
      ggplot2::geom_line(
        data = path_long[!path_long$selected, ],
        colour = "grey80", linewidth = 0.6,
        show.legend = FALSE
      ) +
      ggplot2::geom_line(
        data = path_long[path_long$selected, ],
        linewidth = 0.8
      )
  } else {
    plt <- plt +
      ggplot2::geom_line(linewidth = 0.7)
  }

  plt +
    ggplot2::geom_vline(
      xintercept = best_lambda,
      linetype = "dashed"
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      x = lambda_label_text,
      y = component_label,
      colour = "Feature"
    ) +
    ggplot2::theme_bw()
}


#' Plot predictive densities
#'
#' Plots the joint density \eqn{f(y, x \mid \mathbf{z})} as a
#' filled contour, or, when `y_given` is supplied, overlays the
#' marginal density \eqn{f(x \mid \mathbf{z})} and the
#' conditional density \eqn{f(x \mid \mathbf{z}, y)} for a
#' specific value of \eqn{y}.
#'
#' @param fit A single fitted model list (e.g. one element of the
#'   `full_fit` column from [ccreds()] output), containing
#'   `beta_0`, `alpha`, `beta`, `sigma_sq`, `gamma_0`, `gamma`.
#' @param z Numeric vector of feature values (length \eqn{p}).
#' @param x_range Numeric vector of length 2 giving the range
#'   of \eqn{x} values to plot. Default is computed from the
#'   chi-squared distribution implied by the model.
#' @param y_range Numeric vector of length 2 giving the range
#'   of \eqn{y} values (used for the joint contour plot).
#'   Default is computed from the model.
#' @param y_given Optional scalar. If supplied, produces a
#'   density overlay of \eqn{f(x \mid \mathbf{z})} and
#'   \eqn{f(x \mid \mathbf{z}, y = y_{given})} instead of the
#'   joint contour.
#' @param n_grid Number of grid points per axis. Default is 200.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' results <- ccreds(data, lambda_1s, lambda_2s)
#' fit <- results$full_fit[[1]]
#' z_vals <- c(0.5, -0.3, 1.2)
#'
#' # Joint density contour
#' plot_density(fit, z_vals)
#'
#' # Marginal vs conditional density of x given y = 60
#' plot_density(fit, z_vals, y_given = 60)
#' }
#'
#' @export
plot_density <- function(fit,
                         z,
                         x_range = NULL,
                         y_range = NULL,
                         y_given = NULL,
                         n_grid = 200) {
  # Compute distribution parameters from the model
  mu_x <- exp(fit$gamma_0 + sum(fit$gamma * z))
  mu_y <- fit$beta_0 + fit$alpha * mu_x + sum(fit$beta * z)
  sd_y <- sqrt(fit$sigma_sq)
  
  if (is.null(x_range)) {
    x_range <- c(
      max(0.01, qchisq(0.001, df = mu_x)),
      qchisq(0.999, df = mu_x)
    )
    }
  if (is.null(y_range)) {
    y_range <- c(mu_y - 4 * sd_y, mu_y + 4 * sd_y)
  }
  
  x_seq <- seq(x_range[1], x_range[2], length.out = n_grid)
  y_seq <- seq(y_range[1], y_range[2], length.out = n_grid)

  # Log-density helper functions
  log_f_y_given_x <- function(y, x) {
    mean_y <- fit$beta_0 + fit$alpha * x + sum(fit$beta * z)
    dnorm(y, mean = mean_y, sd = sd_y, log = TRUE)
  }

  log_f_x <- function(x) {
    dchisq(x, df = mu_x, log = TRUE)
  }

  log_f_joint <- function(y, x) {
    log_f_y_given_x(y, x) + log_f_x(x)
  }

  if (is.null(y_given)) {
    # Joint density contour
    grid <- expand.grid(x = x_seq, y = y_seq)
    grid$density <- exp(log_f_joint(grid$y, grid$x))

    ggplot2::ggplot(grid, ggplot2::aes(
      x = .data$x,
      y = .data$y,
      z = .data$density
    )) +
      ggplot2::geom_contour_filled(bins = 8) +
      ggplot2::scale_fill_grey(start = 1, end = 0.2) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::labs(x = "x", y = "y") +
      ggplot2::theme(legend.position = "none")
  } else {
    # Marginal f(x | z) and conditional f(x | z, y)
    diff_x <- x_seq[2] - x_seq[1]

    marginal <- exp(log_f_x(x_seq))
    joint_at_y <- exp(log_f_joint(y_given, x_seq))
    conditional <- joint_at_y / (sum(joint_at_y) * diff_x)

    dens_data <- data.frame(
      x = rep(x_seq, 2),
      density = c(conditional, marginal),
      type = rep(
        c("conditional", "marginal"),
        each = n_grid
      )
    )
    dens_data$type <- factor(
      dens_data$type,
      levels = c("conditional", "marginal")
    )

    cond_label <- paste0("f(x | z, y = ", y_given, ")")
    marg_label <- "f(x | z)"

    ggplot2::ggplot(dens_data, ggplot2::aes(
      x = .data$x, y = .data$density
    )) +
      ggplot2::geom_line(ggplot2::aes(colour = .data$type)) +
      ggplot2::geom_area(
        ggplot2::aes(fill = .data$type),
        alpha = 0.2, position = "identity"
      ) +
      ggplot2::scale_colour_manual(
        labels = c("conditional" = cond_label, "marginal" = marg_label),
        values = c("conditional" = "orange", "marginal" = "black")
      ) +
      ggplot2::scale_fill_manual(
        labels = c("conditional" = cond_label, "marginal" = marg_label),
        values = c("conditional" = "orange", "marginal" = "black")
      ) +
      ggplot2::labs(x = "x", y = "Density", colour = NULL, fill = NULL) +
      ggplot2::theme_bw()
  }
}


# Helper: produce a lambda label, using latex2exp if available
lambda_label <- function(name) {
  tex <- switch(name,
    lambda_1 = "$\\lambda_1$",
    lambda_2 = "$\\lambda_2$"
  )
  if (requireNamespace("latex2exp", quietly = TRUE)) {
    latex2exp::TeX(tex)
  } else {
    gsub("\\$", "", tex)
  }
}
