SETTING_NAME = "some similarity"
GRID_SIZE = 30
NCPUS = 34
NSAMPLES = 100
# Guideline: Takes 10 hours for 1 CPU to fit over a 30x30 grid

library(pracma)
library(gelnet)
library(lbfgs)
library(emdbook)
library(future)
library(furrr)
library(purrr)
library(dplyr)
library(tibble)
library(utf8)
library(withr)
library(stringr)

settings = tibble(
  setting = character(),
  beta_0 = numeric(),
  beta_1 = numeric(),
  beta_2 = list(),
  sigma_sq = numeric(),
  gamma_0 = numeric(),
  gamma_1 = list(),
) |> 
  add_row(
    setting = "exactly same",
    beta_0 = 5,
    beta_1 = 10,
    beta_2 = list(c(2,2,1,1, rep(0,22), -2,-2,-1,-1) * 2),
    sigma_sq = 100,
    gamma_0 = 2,
    gamma_1 = list(c(2,2,1,1, rep(0,22), -2,-2,-1,-1) / 100 * 2)
  ) |> 
  add_row(
    setting = "some similarity",
    beta_0 = 5,
    beta_1 = 10,
    beta_2 = list(c(2,2,1,1, rep(0,22), -2,-2,-1,-1) * 2),
    sigma_sq = 100,
    gamma_0 = 2,
    gamma_1 = list(c(rep(0,2), 2,2,1,1, rep(0,18), -2,-2,-1,-1, rep(0,2)) / 100 * 2)
  ) |>  
  add_row(
    setting = "mutually exclusive",
    beta_0 = 5,
    beta_1 = 10,
    beta_2 = list(c(2,2,1,1, rep(0,22), -2,-2,-1,-1) * 2),
    sigma_sq = 100,
    gamma_0 = 2,
    gamma_1 = list(c(rep(0,4), 2,2,1,1, rep(0,14), -2,-2,-1,-1, rep(0,4)) / 100 * 2)
  ) |> 
  add_row(
    setting = "all gamma 0",
    beta_0 = 5,
    beta_1 = 10,
    beta_2 = list(c(2,2,1,1, rep(0,22), -2,-2,-1,-1) * 2),
    sigma_sq = 100,
    gamma_0 = 2,
    gamma_1 = list(c(rep(0,30)) / 100 * 2)
  )


# Set seed for reproducibility
set.seed(128)

data_sets = tibble(setting = character(),
                   sample = numeric(),
                   data_set = list())

for (setting_name in settings$setting) {
  # for (setting_name in "some similarity") {
  for (sample_number in 1:NSAMPLES) {
    
    # Dimensions
    n = 1000
    p = 30
    
    # Predictors
    rho = 0.5
    Sigma = 9*toeplitz(rho^(0:(p - 1)))
    z = MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
    
    # Parameters
    current_setting = settings |> 
      filter(setting == setting_name)
    
    beta = c(current_setting$beta_0, 
             current_setting$beta_1, 
             current_setting$beta_2[[1]])
    
    sigma_sq = current_setting$sigma_sq * 1
    
    gamma = c(current_setting$gamma_0, 
              current_setting$gamma_1[[1]] / 1)
    
    # Generate data
    x_mean = exp(cbind(1, z) %*% gamma) |> drop()  # mean parameter
    
    x = rchisq(n, df = x_mean)  # to-be-censored covariate
    
    compute_c_mean = function(censoring_rate, x_mean) {
      f = function(c_mean) {
        1 - pbeta(0.5, x_mean/2, c_mean/2) - censoring_rate
      }
      
      c_mean = uniroot(f, lower = 1e-12, upper = 1000)$root
      return(c_mean)
    }
    
    censoring_rate = 0.5
    c_mean = map_dbl(x_mean, ~ compute_c_mean(censoring_rate, .))
    c = rchisq(n, df = c_mean)
    
    
    #--- Start: Strong signal gamma
    
    # if (setting_name == "strong signal") {
    #   x = NA; c = NA
    #   x = pmax(rnorm(n, mean = x_mean, sd = 0.1), 0.1)
    #   c_start = 0.5*x
    #   c_end = 1.5*x
    #   c = runif(n, c_start, c_end)
    # }
    
    #--- End
    
    
    e = rnorm(n, mean = 0, sd = sqrt(sigma_sq))  # random errors
    y = cbind(1, x, z) %*% beta + e  # continuous outcome
    w = pmin(x, c)  # observed covariate value
    d = as.numeric(x <= c)  # "event" indicator
    
    # Construct data set
    z = as.data.frame(z)
    colnames(z) = paste0("z", 1:p)
    random_right_dat = data.frame(y, w, d) |> cbind(z)
    # random_right_dat = data.frame(y, w, x, d) |> cbind(z)
    
    random_right_dat = as_tibble(random_right_dat) |> 
      arrange(desc(d)) |>
      mutate(case = 1:n) |> 
      select(case, everything())
    
    
    data_sets = data_sets |> 
      add_row(setting = setting_name,
              sample = sample_number,
              data_set = list(random_right_dat))
  }
}


set.seed(129)
plan(multisession, workers = NCPUS)
  
res_tune = data_sets |> 
  filter(setting == SETTING_NAME) |> 
  mutate(
    model = furrr::future_map(
      data_set, 
      .f = function(data_set) {
        fit_tuned_model(
          lambda_1s = c(emdbook::lseq(1e-3 * 1, 1, length.out = GRID_SIZE)),
          lambda_2s = c(emdbook::lseq(1e-3 * 1, 1, length.out = GRID_SIZE)),
          data = data_set)
      },
      .progress = FALSE, 
      .options = furrr::furrr_options(seed = TRUE)))

saveRDS(res_tune, 
        paste0("sim_results/res_tune_B",
               NSAMPLES,
               "_",
               gsub(" ", "_", SETTING_NAME), 
               ".rds"))

plan(sequential)

