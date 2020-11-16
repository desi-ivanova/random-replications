if (!require("tidyverse")) install.packages("tidyverse")
if (!require("leaps")) install.packages("leaps")
if (!require("lars")) install.packages("lars")

library(tidyverse)
source("utilities.R")

run_variable_selection <- function(rho, n_vars = 4, N = 300, d = 100, seed=-1, 
                                   nreps = 50, all_fixed = FALSE) {
  # if only y is to be sampled set `all_fixed` to TRUE
  
  if (seed > 0) {
    set.seed(seed)
  }
  
  backward_mses <- list()
  forward_mses <- list()
  lars_mses <- list()
  # TODO
  lasso_mses <- list()
  ridge_mses <- list()
  
  sim_data <- generate_data(rho, n_vars, N, d)
  data <- sim_data$data %>% 
    mutate_all(~. - mean(.))
  true_betas <- sim_data$true_betas
  y_mean <- sim_data$y_mean
  
  full_mod <- lm(y_observed ~ . + 0, data=data)
  # could it be the case that they've used fitted betas instead of true betas?
  #true_betas$value <- coef(full_mod)
  
  for (i in 1:nreps) {
    if (all_fixed) {
      # only change y observed by adding noise to the mean
      data$y_observed <- y_mean + rnorm(N, 0, sqrt(6.25))
    } else {
      # substitute all
      sim_data <- generate_data(rho, n_vars, N, d)
      data <- sim_data$data
      true_betas <- sim_data$true_betas
    }
    

    # run variable selections
    backward_fits <- leaps::regsubsets(
      y_observed ~ ., data = data, nbest = 1, nvmax = d, method = "backward", intercept = FALSE
    ) %>% extract_coefs_leaps()
    forward_fits <- leaps::regsubsets(
      y_observed ~ ., data = data, nbest = 1, nvmax = d, method = "forward", intercept = FALSE
    ) %>% extract_coefs_leaps()
    lars_fits <- lars::lars(
      x = data %>% dplyr::select(-y_observed) %>% as.matrix(), 
      y = data %>% pull(y_observed), 
      type = "lar", intercept = FALSE
    ) %>% coef() %>% as_tibble() %>% 
      filter(row_number() > 1)
    
    # calculate MSEs
    bb <- sum(true_betas$value^2)
    backward_mses[[i]] <- calculate_mse_coefs(backward_fits, true_betas, i) 
    forward_mses[[i]] <- calculate_mse_coefs(forward_fits, true_betas, i) 
    lars_mses[[i]] <- calculate_mse_coefs(lars_fits, true_betas, i)
  }
  
  print(bb)
  bind_rows(list(
    "backward_selection" = reduce(backward_mses, inner_join, by = "subset_size"), 
    "forward_selection" = reduce(forward_mses, inner_join, by = "subset_size"),
    "lars_selection" = reduce(lars_mses, inner_join, by = "subset_size")
  ), .id = "method") 
}



# testing
if (FALSE) {
  res <-  run_variable_selection(rho = 0, d = 31, N = 300, n_vars = 10, 
                                 seed = -1, nreps = 50, all_fixed = TRUE) 
  
  res %>% 
    gather(var, val, -c(method, subset_size)) %>% 
    group_by(method, subset_size) %>% 
    summarise(mean_mse = mean(val), .groups = 'drop') %>% 
    ggplot(aes(subset_size, mean_mse, col = method)) + 
    geom_line() +
    geom_point()
  
  
  av <- c()
  for (i in 1:1000) {
    betas <- rep(0, 31)
    betas[1:10] <- rnorm(10, 0, sqrt(0.4))
    av <- append(av, sum(betas^2))
  }
  
  
}



