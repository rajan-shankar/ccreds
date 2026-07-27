# Internal function: fit model for a single (lambda_1, lambda_2) pair

fit_model <- function(lambda_1,
                      lambda_2,
                      data,
                      init_alpha = 1,
                      init_beta_0 = 0,
                      init_beta = NA,
                      init_sigma_sq = 1,
                      init_gamma_0 = 0,
                      init_gamma = NA,
                      force_active_beta = integer(0),
                      force_active_gamma = integer(0),
                      maxit = 100,
                      eps = 1e-02,
                      verbose = FALSE) {
  # Data set variables must be in this order:
  #   case, y, w, d, z1, z2, ..., zp

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
    message("some censored w values are < 0.01; these have been changed to 0.01")
  }

  if (length(force_active_gamma) > 0) {
    data <- data |>
      dplyr::select(
        "case", "y", "w", "d",
        paste0("z", force_active_gamma),
        dplyr::everything()
      )
  }

  z_order <- data |>
    dplyr::select(dplyr::starts_with("z")) |>
    colnames() |>
    stringr::str_remove("z") |>
    as.numeric()

  if (length(force_active_beta) > 0) {
    force_active_beta <- purrr::map_dbl(
      force_active_beta,
      \(x) which(x == z_order)
    )
  }

  # Complete observations before censored observations
  data <- data |> dplyr::arrange(dplyr::desc(.data$d))

  # Standardise predictors
  scaled_features <- data |>
    dplyr::select("y", dplyr::starts_with("z")) |>
    as.matrix() |>
    scale()

  data <- data |>
    dplyr::select(-"y", -dplyr::starts_with("z")) |>
    dplyr::bind_cols(tibble::as_tibble(scaled_features))

  # Extract data components
  y <- data$y
  w <- data$w
  d <- data$d
  z <- data |>
    dplyr::select(dplyr::starts_with("z")) |>
    as.matrix()
  Z <- cbind(1, z)

  data_1 <- data |> dplyr::filter(.data$d == 1)
  y1 <- data_1$y
  w1 <- data_1$w
  z1 <- data_1 |>
    dplyr::select(dplyr::starts_with("z")) |>
    as.matrix()

  data_0 <- data |> dplyr::filter(.data$d == 0)
  y0 <- data_0$y
  w0 <- data_0$w
  z0 <- data_0 |>
    dplyr::select(dplyr::starts_with("z")) |>
    as.matrix()

  n <- nrow(data)
  p <- ncol(z)

  gelnet_d <- rep(1, p)
  if (length(force_active_beta) > 0) gelnet_d[force_active_beta] <- 0
  # Penalty weights: 0 for beta_0, 0 for alpha, then gelnet_d for beta
  gelnet_d <- c(0, 0, gelnet_d)

  # Control
  diff <- NA
  observed_data_loss_fn <- NA

  # Parameters (using paper notation internally)
  beta_0_est <- init_beta_0
  alpha_est <- init_alpha
  beta_est <- init_beta
  sigma_sq_est <- init_sigma_sq
  gamma_0_est <- init_gamma_0
  gamma_est <- init_gamma

  if (length(force_active_beta) > 0) {
    beta_est <- beta_est[z_order]
  }
  if (length(force_active_gamma) > 0) {
    gamma_est <- gamma_est[z_order]
  }

  # OWL-QN functions for gamma
  obj_gamma_smooth_fn <- function(x, a_vec, Z, n) {
    e_Zx <- exp(drop(Z %*% x))
    drop(
      1 / (2 * n) * t(a_vec) %*% e_Zx +
        1 / n * sum(lgamma(1 / 2 * e_Zx))
    )
  }

  obj_gamma_smooth_grad <- function(x, a_vec, Z, n) {
    e_Zx <- exp(drop(Z %*% x))
    drop(
      1 / (2 * n) * t(Z * e_Zx) %*%
        (a_vec + digamma(1 / 2 * e_Zx))
    )
  }

  # EM algorithm
  m <- w0
  for (i in 1:maxit) {
    if (abs(alpha_est) < 1e-5) {
      stop(
        "alpha being estimated or initialised too close to 0; ",
        "consider whether your x variable really affects your y variable."
      )
    }

    # --- E-step ---
    m <- ell <- s <- log_int_t_tilde <- NA
    for (j in seq_along(w0)) {
      # nu = y_i - beta_0 - beta^T z_i (residual without x contribution)
      nu <- y0[j] - beta_0_est - drop(t(beta_est) %*% z0[j, ])
      mu <- exp(gamma_0_est + drop(t(gamma_est) %*% z0[j, ]))

      a_var <- drop(alpha_est^2 / (2 * sigma_sq_est))
      b_var <- drop(nu / alpha_est - sigma_sq_est / (2 * alpha_est^2))
      c_var <- drop(mu / 2)

      expectations <- compute_expectations(a_var, b_var, c_var, w0[j])
      m[j] <- expectations$m
      s[j] <- expectations$s
      ell[j] <- expectations$ell

      log_A <- -1 / 2 * (nu / alpha_est - sigma_sq_est / (4 * alpha_est^2)) -
        1 / 2 * log(2 * pi * sigma_sq_est) -
        mu / 2 * log(2) - lgamma(mu / 2)
      log_int_t_tilde[j] <- log_A + expectations$log_u_M_plus_log_I0
    }

    # Penalised observed-data loss function (negative, divided by n)
    observed_data_loss_fn[i] <- -1 / n * sum(dnorm(
      y1,
      mean = beta_0_est + alpha_est * w1 + z1 %*% beta_est,
      sd = sqrt(sigma_sq_est),
      log = TRUE
    ) + dchisq(
      w1,
      df = exp(gamma_0_est + z1 %*% gamma_est),
      log = TRUE
    )) - 1 / n * sum(log_int_t_tilde)

    observed_data_loss_fn[i] <- observed_data_loss_fn[i] +
      lambda_1 / sigma_sq_est * sum(abs(beta_est)) +
      lambda_2 * sum(abs(gamma_est))

    # --- M-step (beta_0, alpha, beta, sigma_sq) ---
    # Design matrix: [1, w^(k), Z] where w^(k) = [x_complete, m_censored]
    A <- cbind(1, c(w1, m), z)

    # gelnet expects the coefficient vector in the order [beta_0, alpha, beta]
    gelnet_w_init <- c(beta_0_est, alpha_est, beta_est)

    res_beta <- gelnet::gelnet(
      A, y,
      l1 = lambda_1,
      l2 = 1 / n,
      d = gelnet_d,
      P = diag(c(0, sum(s^2), rep(0, p))),
      w.init = gelnet_w_init,
      b.init = 0,
      fix.bias = TRUE,
      silent = TRUE
    )

    obj_beta <- gelnet::gelnet.lin.obj(
      w = res_beta$w,
      b = res_beta$b,
      X = A,
      z = y,
      l1 = lambda_1,
      l2 = 1 / n,
      d = gelnet_d,
      P = diag(c(0, sum(s^2), rep(0, p)))
    )

    # --- M-step (gamma_0, gamma) ---
    a_vec <- c(log(2) - log(w1), log(2) - ell)

    gamma_vec_init <- c(gamma_0_est, gamma_est)
    res_gamma <- lbfgs::lbfgs(
      call_eval = obj_gamma_smooth_fn,
      call_grad = obj_gamma_smooth_grad,
      vars = gamma_vec_init,
      a_vec = a_vec,
      Z = Z,
      n = n,
      orthantwise_c = lambda_2,
      orthantwise_start = 1 + length(force_active_gamma),
      invisible = 1
    )

    # Update estimates
    beta_0_est <- res_beta$w[1]
    alpha_est <- res_beta$w[2]
    beta_est <- res_beta$w[3:(p + 2)]
    sigma_sq_est <- drop(2 * obj_beta)
    gamma_0_est <- res_gamma$par[1]
    gamma_est <- res_gamma$par[2:(p + 1)]

    if (i > 1) {
      diff[i] <- (observed_data_loss_fn[i] - observed_data_loss_fn[i - 1]) /
        observed_data_loss_fn[i - 1] * 100
      if (abs(diff[i]) < eps) break
    }

    if (verbose) {
      cat(
        i,
        "\tbeta_0:", round(beta_0_est, 3),
        "\talpha:", round(alpha_est, 3),
        "\tsigma_sq:", round(sigma_sq_est, 3),
        "\trel_diff:", round(diff[i], 3),
        "beta_sel:", sum(beta_est != 0),
        "gamma_sel:", sum(gamma_est != 0),
        "obj:", round(observed_data_loss_fn[i], 3),
        "\n"
      )
    }
  }

  res <- list(
    lambda_1 = lambda_1,
    lambda_2 = lambda_2,
    beta_0 = beta_0_est,
    alpha = alpha_est,
    beta = beta_est,
    sigma_sq = sigma_sq_est,
    gamma_0 = gamma_0_est,
    gamma = gamma_est,
    beta_selected = sum(beta_est != 0) - length(force_active_beta),
    gamma_selected = sum(gamma_est != 0) - length(force_active_gamma),
    iterations = i,
    observed_data_loss_fn = observed_data_loss_fn,
    diff = diff
  )

  # Back-transform from standardised to original scale

  res <- unstandardise_fit(res, scaled_features)

  # Restore original z ordering
  if (length(force_active_beta) > 0) {
    res$beta <- res$beta[order(z_order)]
  }
  if (length(force_active_gamma) > 0) {
    res$gamma <- res$gamma[order(z_order)]
  }

  res
}
