# Internal helper functions (not exported)

unstandardise_fit <- function(fit, scaled_features) {
  p <- length(fit$beta)
  center_z <- attr(scaled_features, "scaled:center")[-1]
  scale_z <- attr(scaled_features, "scaled:scale")[-1]
  scale_y <- attr(scaled_features, "scaled:scale")[1]
  center_y <- attr(scaled_features, "scaled:center")[1]


  # Back-transform beta (z coefficients) and beta_0, alpha

  beta_orig <- fit$beta / scale_z * scale_y
  alpha_orig <- fit$alpha * scale_y
  beta_0_orig <- (fit$beta_0 - sum(fit$beta * center_z / scale_z)) * scale_y +
    center_y

  # Back-transform gamma (z coefficients) and gamma_0
  gamma_orig <- fit$gamma / scale_z
  gamma_0_orig <- fit$gamma_0 - sum(fit$gamma * center_z / scale_z)

  # Back-transform sigma_sq
  sigma_sq_orig <- fit$sigma_sq * scale_y^2

  fit$beta_0 <- beta_0_orig
  fit$alpha <- alpha_orig
  fit$beta <- beta_orig
  fit$gamma_0 <- gamma_0_orig
  fit$gamma <- gamma_orig
  fit$sigma_sq <- sigma_sq_orig

  fit
}


compute_expectations <- function(a, b, c, w, tol = 1e-5) {
  m <- s <- ell <- NA

  log_u <- function(x, a, b, c) {
    -a * (x - b)^2 + (c - 1) * log(x)
  }

  discriminant <- 4 * (a * b)^2 + 8 * a * (c - 1)
  M <- ifelse(
    discriminant >= 0,
    yes = max(w, b / 2 + 1 / (4 * a) * sqrt(discriminant)),
    no = w
  )
  log_u_M <- log_u(M, a, b, c)

  Ik_f <- function(x, a, b, c, log_u_M) {
    exp(log_u(x, a, b, c) - log_u_M)
  }
  Ilog_f <- function(x, a, b, c, log_u_M) {
    log(x) * Ik_f(x, a, b, c, log_u_M)
  }

  I0 <- I1 <- I2 <- Ilog <- NA
  x_values_I0 <- NA

  for (k in 0:2) {
    w_fval <- Ik_f(w, a, b, c + k, log_u_M)

    discriminant <- 4 * a^2 * b^2 + 8 * a * (c + k - 1)
    argmax_fval <- ifelse(
      discriminant >= 0,
      yes = max(w, b / 2 + 1 / (4 * a) * sqrt(discriminant)),
      no = w
    )
    max_fval <- Ik_f(argmax_fval, a, b, c + k, log_u_M)

    integral_return_value <- NA
    lower_bound <- w

    if (w_fval < tol) {
      if (max_fval < tol) {
        integral_return_value <- 0
      } else if (w < argmax_fval) {
        lower_bound <- uniroot(
          f = function(x, a, b, c, log_u_M, threshold, max_fval) {
            Ik_f(x, a, b, c, log_u_M) / max_fval - threshold
          },
          interval = c(w, argmax_fval),
          a = a, b = b, c = c + k, log_u_M = log_u_M,
          threshold = tol, max_fval = max_fval,
          extendInt = "upX", tol = 1e-16
        )$root
      } else {
        stop("w > argmax_fval")
      }
    }

    if (is.na(integral_return_value)) {
      upper_bound <- uniroot(
        f = function(x, a, b, c, log_u_M, threshold, max_fval) {
          Ik_f(x, a, b, c, log_u_M) / max_fval - threshold
        },
        interval = c(argmax_fval, 10000),
        a = a, b = b, c = c + k, log_u_M = log_u_M,
        threshold = tol, max_fval = max_fval,
        extendInt = "downX", tol = 1e-16
      )$root

      if (lower_bound > upper_bound) {
        stop("Lower bound is greater than upper bound in integral")
      }

      x_values <- seq(lower_bound, upper_bound, length.out = 21)
      y_values <- Ik_f(x_values, a, b, c + k, log_u_M) / max_fval
      integral_return_value <- pracma::trapz(x_values, y_values) * max_fval
    }

    if (k == 0) {
      I0 <- integral_return_value
      x_values_I0 <- x_values
    } else if (k == 1) {
      I1 <- integral_return_value
    } else if (k == 2) {
      I2 <- integral_return_value
    }
  }

  y_values <- Ilog_f(x_values_I0, a, b, c, log_u_M)
  Ilog <- pracma::trapz(x_values_I0, y_values)

  if (any(is.na(c(I0, I1, I2, Ilog)))) {
    stop("One or more integrals evaluated to NA/NaN")
  }

  m <- I1 / I0
  s <- sqrt(I2 / I0 - m^2)
  ell <- Ilog / I0

  if (any(is.na(c(m, s, ell)))) {
    stop("One or more of m, s, ell evaluated to NA/NaN")
  }

  if (m < w) {
    stop("m is being estimated < w")
  }

  list(
    m = m,
    s = s,
    ell = ell,
    log_u_M_plus_log_I0 = log_u_M + log(I0)
  )
}


get_folds <- function(data, k) {
  complete_indices <- data |> dplyr::filter(.data$d == 1) |> dplyr::pull(.data$case)
  censored_indices <- data |> dplyr::filter(.data$d == 0) |> dplyr::pull(.data$case)

  complete_folds <- split(
    sample(complete_indices),
    rep(1:k, length = length(complete_indices))
  )
  censored_folds <- split(
    sample(censored_indices),
    rep(1:k, length = length(censored_indices))
  )

  folds <- list()
  for (i in 1:k) {
    folds[[i]] <- c(complete_folds[[i]], censored_folds[[i]])
  }

  folds
}
