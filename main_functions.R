
unstandardise_fit = function(fit, scaled_features) {
  p = length(fit$beta) - 2
  fit$beta = c(
    fit$beta[1] - sum(fit$beta[3:(p+2)] * attr(scaled_features, "scaled:center")[-1] / attr(scaled_features, "scaled:scale")[-1]), 
    fit$beta[2],
    fit$beta[3:(p+2)] / attr(scaled_features, "scaled:scale")[-1]
  ) * attr(scaled_features, "scaled:scale")[1]
  
  fit$beta[1] = fit$beta[1] + attr(scaled_features, "scaled:center")[1]
  
  fit$gamma = c(
    fit$gamma[1] - sum(fit$gamma[2:(p+1)] * attr(scaled_features, "scaled:center")[-1] / attr(scaled_features, "scaled:scale")[-1]),
    fit$gamma[2:(p+1)] / attr(scaled_features, "scaled:scale")[-1]
  )
  
  fit$sigma_sq = fit$sigma_sq * attr(scaled_features, "scaled:scale")[1]^2
  
  return(fit)
}



compute_expectations = function(a, b, c, w, tol = 1e-5) {
  
  # Outputs
  m = s = ell = NA
  
  # Compute maximum value of log u(x) over the interval [w, infinity)
  log_u = function(x, a, b, c) {-a*(x-b)^2 + (c-1)*log(x)}
  discriminant = 4*(a*b)^2 + 8*a*(c-1)
  M = ifelse(
    discriminant >= 0,
    yes = max(w, b/2 + 1/(4*a) * sqrt(discriminant)),
    no = w
  )
  log_u_M = log_u(M, a, b, c)  # log u(M), where M is the mode
  
  # Integrand functions
  Ik_f = function(x, a, b, c, log_u_M) {exp(log_u(x, a, b, c) - log_u_M)}
  Ilog_f = function(x, a, b, c, log_u_M) {log(x) * Ik_f(x, a, b, c, log_u_M)}
  
  # Compute integrals for I0, I1, I2
  I0 = I1 = I2 = Ilog = NA
  x_values_I0 = NA
  for (k in 0:2) {
    # Value at w
    w_fval = Ik_f(w, a, b, c+k, log_u_M)
    
    # Maximum value over the interval [w, infinity)
    discriminant = 4*a^2*b^2 + 8*a*(c+k-1)
    argmax_fval = ifelse(
      discriminant >= 0,
      yes = max(w, b/2 + 1/(4*a) * sqrt(discriminant)),
      no = w
    )
    max_fval = Ik_f(argmax_fval, a, b, c+k, log_u_M)
    
    integral_return_value = NA
    lower_bound = w
    # Lower bound should be w unless value at w is not large enough
    if (w_fval < tol) {
      if (max_fval < tol) {
        # Function is too small; only happens to I1_f and I2_f, AND only when 'action' is between 0 and 1. I0_f is normalised so its max_fval = 1
        integral_return_value = 0
      } else if (w < argmax_fval) {
        # w is too far away from the 'action'; choose a more suitable lower bound closer to the action
        lower_bound = uniroot(
          f = function(x, a, b, c, log_u_M, threshold, max_fval) {
            Ik_f(x, a, b, c, log_u_M) / max_fval - threshold
          }, 
          interval = c(w, argmax_fval), 
          a = a,
          b = b,
          c = c+k,
          log_u_M = log_u_M,
          threshold = tol,
          max_fval = max_fval,
          extendInt = "upX", 
          tol = 1e-16
        )$root
      } else {
        stop("w > argmax_fval")
      }
    }
    
    if (is.na(integral_return_value)) {
      # Choose suitable upper bound
      upper_bound = uniroot(
        f = function(x, a, b, c, log_u_M, threshold, max_fval) {
          Ik_f(x, a, b, c, log_u_M) / max_fval - threshold
        }, 
        interval = c(argmax_fval, 10000), 
        a = a,
        b = b,
        c = c+k,
        log_u_M = log_u_M,
        threshold = tol,
        max_fval = max_fval,
        extendInt = "downX", 
        tol = 1e-16
      )$root
      
      if (lower_bound > upper_bound) {
        stop("Lower bound is greater than upper bound in integral")
      }
      
      # Integrate with trapezoidal rule
      x_values = seq(lower_bound, upper_bound, length.out = 21)
      y_values = Ik_f(x_values, a, b, c+k, log_u_M) / max_fval
      integral_return_value = pracma::trapz(x_values, y_values) * max_fval
    }
    
    if (k == 0) {
      I0 = integral_return_value
      x_values_I0 = x_values  # save for computing Ilog
    } else if (k == 1) {
      I1 = integral_return_value
    } else if (k == 2) {
      I2 = integral_return_value
    }
  }
  
  # Compute Ilog with bounds used in computing I0
  y_values = Ilog_f(x_values_I0, a, b, c, log_u_M)
  integral_return_value = pracma::trapz(x_values_I0, y_values)
  Ilog = integral_return_value
  
  if (any(is.na(c(I0, I1, I2, Ilog)))) {
    stop("One or more integrals evaluated to NA/NaN")
  }
  
  # Expectation quantities
  m = I1 / I0
  s = sqrt(I2 / I0 - m^2)
  ell = Ilog / I0
  
  if (any(is.na(c(m, s, ell)))) {
    stop("One or more of m, s, ell evaluated to NA/NaN")
  }
  
  if (m < w) {
    stop("m is being estimated < w")
  }
  
  return(list(
    # debug = list(I0, I1, I2, Ilog),
    m = m,
    s = s,
    ell = ell,
    log_u_M_plus_log_I0 = log_u_M + log(I0)
  ))
}


fit_model = function(lambda_1,
                     lambda_2,
                     data,
                     init_beta = NA,
                     init_sigma_sq = NA,
                     init_gamma = NA,
                     force_active_beta = integer(0),
                     force_active_gamma = integer(0),
                     maxit = 100,
                     eps = 1e-02,
                     verbose = FALSE) {
  
  # Data set variables must be in this order: 
  #   case, y, w, d, z1, z2, ..., zp
  
  # Make minimum value of w be 0.01:
  min_w = data |> filter(d == 0) |> pull(w) |> min()
  if (min_w < 0.01) {
    data = data |> 
      mutate(w = case_when(d == 0 & w < 0.01 ~ 0.01,
                           TRUE ~ w))
    message("some censored w values are < 0.01; these have been changed to 0.01")
  }
  
  
  if (length(force_active_gamma) > 0) {
    data = data |> 
      select(case, y, w, d, 
             paste0("z", force_active_gamma),
             everything())
  }
  
  z_order = data |> 
    select(starts_with("z")) |> 
    colnames() |> 
    str_remove("z") |> 
    as.numeric()
  
  if (length(force_active_beta) > 0) {
    # Update force_active_beta to new order
    force_active_beta = map_dbl(
      force_active_beta, 
      ~ which(. == z_order)
    )
  }
  
  # IMPORTANT! Complete observations before censored observations
  data = data |> arrange(desc(d))
  
  data_unstandardised = data
  
  # Standardise predictors
  scaled_features = data |>
    select(y, starts_with("z")) |>
    as.matrix() |>
    scale()
  
  data = data |>
    select(-y, -starts_with("z")) |>
    bind_cols(as_tibble(scaled_features))
  
  # Data
  y = data$y
  w = data$w
  d = data$d
  z = data |> 
    select(starts_with("z")) |> 
    as.matrix()
  Z = cbind(1, z)
  
  data_1 = data |> filter(d == 1)
  y1 = data_1$y
  w1 = data_1$w
  d1 = data_1$d
  z1 = data_1 |> 
    select(starts_with("z")) |> 
    as.matrix()
  
  data_0 = data |> filter(d == 0)
  y0 = data_0$y
  w0 = data_0$w
  d0 = data_0$d
  z0 = data_0 |> 
    select(starts_with("z")) |> 
    as.matrix()
  
  n = nrow(data)
  p = ncol(z)
  
  gelnet_d = rep(1, p)
  if (length(force_active_beta) > 0) {gelnet_d[force_active_beta] = 0}
  gelnet_d = c(0, 0, gelnet_d)
  
  # Control
  diff = NA
  observed_data_loss_fn = NA
  
  # Parameters
  beta_est = init_beta
  sigma_sq_est = init_sigma_sq
  gamma_est = init_gamma
  
  if(length(force_active_beta) > 0) {
    beta_est = c(beta_est[1:2], beta_est[3:(p+2)][z_order])
  }
  if(length(force_active_gamma) > 0) {
    gamma_est = c(gamma_est[1], gamma_est[2:(p+1)][z_order])
  }
  
  # OWL-QN functions
  obj_gamma_smooth_fn = function(x, a_vec, Z, n) {
    e_Zx = exp(drop(Z %*% x))
    drop(
      1/(2*n) * t(a_vec) %*% e_Zx + 1/n * sum(lgamma(1/2 * e_Zx))
    )
  }
  
  obj_gamma_smooth_grad = function(x, a_vec, Z, n) {
    e_Zx = exp(drop(Z %*% x))
    drop(
      1/(2*n) * t(Z * e_Zx) %*% (a_vec + digamma(1/2 * e_Zx))
    )
    # Note that doing Z * e_Zx does 'column-wise' element-wise multiplication of e_Zx to Z (order doesn't matter, it's always column-wise not row-wise)
  }
  
  # Fitting
  m = w0  # for computing observed_data_loss_fn in first iteration
  for (i in 1:maxit) {
    
    beta_0_est = beta_est[1]
    beta_1_est = beta_est[2]
    beta_2_est = beta_est[3:length(beta_est)]
    
    gamma_0_est = gamma_est[1]
    gamma_1_est = gamma_est[2:length(gamma_est)]
    
    if (abs(beta_1_est) < 1e-5) {
      stop("beta_1 being estimated or initialised too close to 0; consider whether your x variable really affects your y variable.")
    }
    
    # Expectation step
    m = NA
    ell = NA
    s = NA
    log_int_t_tilde = NA
    for (j in 1:length(w0)) {
      nu = y0[j] - beta_0_est - t(beta_2_est) %*% z0[j,]
      mu = exp(gamma_0_est + t(gamma_1_est) %*% z0[j,])
      
      a_var = drop(beta_1_est^2 / (2*sigma_sq_est))
      b_var = drop(nu/beta_1_est - sigma_sq_est/(2*beta_1_est^2))
      c_var = drop(mu/2)
      
      expectations = compute_expectations(a_var, b_var, c_var, w0[j])
      m[j] = expectations$m
      s[j] = expectations$s
      ell[j] = expectations$ell
      
      log_A = -1/2 * (nu / beta_1_est - sigma_sq_est / (4*beta_1_est^2)) - 1/2*log(2*pi*sigma_sq_est) - mu/2*log(2) - lgamma(mu/2)
      log_int_t_tilde[j] = log_A + expectations$log_u_M_plus_log_I0
    }
    
    # Observed data log-likelihood (penalised, negative, divided by n):
    observed_data_loss_fn[i] = -1/n * sum(dnorm(
      y1, 
      mean = beta_0_est + beta_1_est*w1 + z1 %*% beta_2_est,
      sd = sqrt(sigma_sq_est),
      log = TRUE
    ) + dchisq(
      w1,
      df = exp(gamma_0_est + z1 %*% gamma_1_est),
      log = TRUE
    )) - 1/n * sum(
      log_int_t_tilde
    )
    
    
    observed_data_loss_fn[i] = observed_data_loss_fn[i] + lambda_1/sigma_sq_est*sum(abs(beta_2_est)) +  lambda_2*sum(abs(gamma_1_est))
    
    
    # Maximisation step (beta, sigma_sq)
    A = cbind(1, c(w1,m), z)
    
    res_beta = gelnet(A, y,
                      l1 = lambda_1,
                      l2 = 1/n,
                      d = gelnet_d,
                      P = diag(c(0, sum(s^2), rep(0, length(beta_2_est)))),
                      w.init = beta_est,
                      b.init = 0,
                      fix.bias = TRUE,
                      silent = TRUE)
    
    obj_beta = gelnet.lin.obj(
      w = res_beta$w,
      b = res_beta$b,
      X = A,
      z = y,
      l1 = lambda_1,
      l2 = 1/n,
      d = gelnet_d,
      P = diag(c(0, sum(s^2), rep(0, length(beta_2_est)))),
    )
    
    # Maximisation step (gamma)
    a_vec = c(log(2) - log(w1), 
              log(2) - ell)
    
    res_gamma = lbfgs::lbfgs(
      call_eval = obj_gamma_smooth_fn,
      call_grad = obj_gamma_smooth_grad,
      vars = gamma_est,
      a_vec = a_vec,
      Z = Z,
      n = n,
      orthantwise_c = lambda_2,
      orthantwise_start = 1 + length(force_active_gamma),
      invisible = 1
    )
    
    # Get new estimates
    old_est = c(beta_est, sigma_sq_est, gamma_est)
    beta_est = res_beta$w
    sigma_sq_est = drop(2 * obj_beta)
    gamma_est = res_gamma$par
    new_est = c(beta_est, sigma_sq_est, gamma_est)
    
    if (i > 1) {
      diff[i] = (observed_data_loss_fn[i] - observed_data_loss_fn[i-1]) / observed_data_loss_fn[i-1] * 100
      if (abs(diff[i]) < eps) {break}
    }
    
    if (verbose == TRUE) {
      cat(i,
          "\tbeta0:", beta_est[1] |> round(3),
          "\tbeta1:", beta_est[2] |> round(3),
          "\tsigma^2:", sigma_sq_est |> round(3),
          "\trel_diff:", diff[i] |> round(3),
          "beta_sel:", sum(beta_est[3:length(beta_est)] != 0),
          "gamma_sel:", sum(gamma_est[2:length(gamma_est)] != 0),
          "obj:", observed_data_loss_fn[i] |> round(3),
          "\n")
    }
  }
  
  res = list(
    lambda_1 = lambda_1,
    lambda_2 = lambda_2,
    beta = beta_est,
    sigma_sq = sigma_sq_est,
    gamma = gamma_est,
    beta_sel = sum(beta_est[3:length(beta_est)] != 0) - length(force_active_beta),
    gamma_sel = sum(gamma_est[2:length(gamma_est)] != 0) - length(force_active_gamma),
    iterations = i,
    observed_data_loss_fn = observed_data_loss_fn,
    diff = diff
  )
  
  # Remember that observed_data_loss_fn is computed and remains on the *standardised* scale.
  res = unstandardise_fit(res, scaled_features)
  
  if(length(force_active_beta) > 0) {
    res$beta = c(res$beta[1:2], res$beta[3:(p+2)][order(z_order)])
  }
  if(length(force_active_gamma) > 0) {
    res$gamma = c(res$gamma[1], res$gamma[2:(p+1)][order(z_order)])
  }
  
  return(res)
}


get_folds = function(data, k) {
  
  complete_indices = data |> filter(d == 1) |> pull(case)
  censored_indices = data |> filter(d == 0) |> pull(case)
  
  complete_folds = split(sample(complete_indices), 
                         rep(1:k, length = length(complete_indices)))
  censored_folds = split(sample(censored_indices), 
                         rep(1:k, length = length(censored_indices)))
  
  folds = list()
  for (i in 1:k) {
    folds[[i]] = c(complete_folds[[i]], censored_folds[[i]])
  }
  
  return(folds)
}

fit_tuned_model = function(
    lambda_1s, 
    lambda_2s,
    force_active_beta = integer(0),
    force_active_gamma = integer(0),
    data,
    k = 5,
    cores = 1
) {
  
  # Make minimum value of w be 0.01:
  min_w = data |> filter(d == 0) |> pull(w) |> min()
  if (min_w < 0.01) {
    data = data |> 
      mutate(w = case_when(d == 0 & w < 0.01 ~ 0.01,
                           TRUE ~ w))
    message("some censored w values are < 0.01; these have been changed to 0.01")
  }
  
  # Number of predictors
  p = data |> select(starts_with("z")) |> ncol()
  
  # Get strata for stratified cross-validation
  folds = get_folds(data, k)
  
  # Order in which to traverse the lambda grid
  J1 = length(lambda_1s)
  J2 = length(lambda_2s)
  
  if (any(c(J1, J2) < 1)) {
    stop("Number of values for lambda_1 AND lambda_2 need to be >= 1")
  }
  
  order_fit1 = rep(1:J1, each = J2)
  order_fit2 = rep(1:J2, times = J1)
  
  # Initialise results output
  results = tibble(
    lambda_1 = lambda_1s[order_fit1],
    lambda_2 = lambda_2s[order_fit2]
  )
  
  # Function to fit a model and compute its cross-validation metrics
  lambda_pair_fn = function(lambda_1, lambda_2) {
    
    t1 = Sys.time()
    
    # In case any errors occur
    error_output = list("full_fit" = NA,
                        "CV_MAE_y" = NA,
                        "CV_MAE_y_SE" = NA,
                        "CV_MAE_x" = NA,
                        "CV_MAE_x_SE" = NA,
                        "time" = NA)
    
    # Fit model
    full_fit = tryCatch({ 
      fit_model(
        init_beta = c(0,1,rep(0,p)),
        init_sigma_sq = 1,
        init_gamma = c(0,rep(0,p)),
        data = data,
        lambda_1 = lambda_1,
        lambda_2 = lambda_2,
        force_active_beta = force_active_beta,
        force_active_gamma = force_active_gamma
      )},
      error = function(e) e
    )
    
    if (inherits(full_fit, "error")) {
      error_output$full_fit = list(
        "occurred" = "main full_fit",
        "error" = full_fit
      )
      return(error_output)
    }
    
    # Initialise output for cross-validation
    CV_output = tibble(
      fold = numeric(),
      y_errors = list(),
      x_errors = list()
    )
    
    # Cross-validation
    for (i in 1:k) {
      train_data = data |> slice(list_c(folds[-i]))
      test_data = data |> slice(folds[[i]]) |> arrange(desc(d))
      
      # Fit model excluding i^th fold
      fit = tryCatch({ 
        fit_model(
          init_beta = c(0,1,rep(0,p)),
          init_sigma_sq = 1,
          init_gamma = c(0,rep(0,p)),
          data = train_data,
          lambda_1 = lambda_1,
          lambda_2 = lambda_2,
          force_active_beta = force_active_beta,
          force_active_gamma = force_active_gamma
        )},
        error = function(e) e
      )
      
      if (inherits(fit, "error")) {
        break
      }
      
      # Make predictions for x
      x_pred_complete = exp(drop(cbind(
        1, test_data |> 
          filter(d == 1) |> 
          select(starts_with("z")) |> 
          as.matrix()
      ) %*% 
        fit$gamma))
      
      x_pred_censored = exp(drop(cbind(
        1, test_data |> 
          filter(d == 0) |> 
          select(starts_with("z")) |> 
          as.matrix()
      ) %*% 
        fit$gamma))
      
      # Make predictions for y
      y_pred = cbind(
        1, c(x_pred_complete, x_pred_censored),
        test_data |> 
          select(starts_with("z")) |> 
          as.matrix()
      ) %*% 
        fit$beta
      
      # Gather prediction errors on test fold
      y_errors = y_pred - test_data$y
      x_errors = x_pred_complete - test_data$w[test_data$d == 1]
      CV_output = CV_output |> 
        add_row(
          fold = i,
          y_errors = list(y_errors),
          x_errors = list(x_errors)
        )
    }
    
    if (inherits(fit, "error")) {
      error_output$full_fit = list(
        "occurred" = "CV loop",
        "error" = fit
      )
      return(error_output)
    }
    
    # Compute cross-validation metrics
    CV_MAE_y = mean(CV_output$y_errors |> map_dbl(~ mean(abs(.))))
    CV_MAE_y_SE = sd(CV_output$y_errors |> map_dbl(~ mean(abs(.)))) / sqrt(k)
    CV_MAE_x = mean(CV_output$x_errors |> map_dbl(~ mean(abs(.))))
    CV_MAE_x_SE = sd(CV_output$x_errors |> map_dbl(~ mean(abs(.)))) / sqrt(k)
    
    t2 = Sys.time()
    
    # Output
    result = list("full_fit" = full_fit,
                  "CV_MAE_y" = CV_MAE_y,
                  "CV_MAE_y_SE" = CV_MAE_y_SE,
                  "CV_MAE_x" = CV_MAE_x,
                  "CV_MAE_x_SE" = CV_MAE_x_SE,
                  "time" = t2 - t1)
    
    return(result)
  }
  
  # Fit models over lambda grid
  if (cores > 1) {
    
    future::plan(future::multisession, workers = cores)
    
    results = results |> 
      mutate(result = furrr::future_map2(
        lambda_1, lambda_2, .f = lambda_pair_fn,
        .progress = TRUE,
        .options = furrr::furrr_options(seed = TRUE)
      )) |> 
      unnest_wider(result)
    
    future::plan(future::sequential)
    
  } else {
    
    results = results |> 
      mutate(result = purrr::map2(
        lambda_1, lambda_2, .f = lambda_pair_fn,
        .progress = TRUE
      )) |> 
      unnest_wider(result)
  }
  
  # Return results output
  return(results)
}

