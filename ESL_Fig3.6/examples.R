if (!require("tidyverse")) install.packages("tidyverse")
if (!require("leaps")) install.packages("leaps")
if (!require("lars")) install.packages("lars")

library(tidyverse)
source("utilities.R")

run_variable_selection <- function(rho, n_vars = 4, N = 200, d = 100, seed=-1, nreps = 50) {
  if (seed > 1) {
    set.seed(seed)
  }
  
  backward_mses <- list()
  forward_mses <- list()
  lars_mses <- list()
  # TODO
  lasso_mses <- list()
  ridge_mses <- list()
  
  for (i in 1:nreps) {
    sim_data <- generate_data(rho, n_vars, N, d)
    data <- sim_data$data
    true_betas <- sim_data$true_betas
    
    # run variable selections
    backward_fits <- leaps::regsubsets(
      y_observed ~ ., data = data, nbest = 1, nvmax = d, method = "backward"
    ) %>% extract_coefs_leaps()
    forward_fits <- leaps::regsubsets(
      y_observed ~ ., data = data, nbest = 1, nvmax = d, method = "forward"
    ) %>% extract_coefs_leaps()
    lars_fits <- lars::lars(
      x = data %>% select(-y_observed) %>% as.matrix, 
      y = data %>% pull(y_observed), type = "lar"
    ) %>% coef() %>% as_tibble() %>% 
      filter(row_number() > 1)
    
    # calculate MSEs
    backward_mses[[i]] <- calculate_mse_coefs(backward_fits, true_betas, i) 
    forward_mses[[i]] <- calculate_mse_coefs(forward_fits, true_betas, i) 
    lars_mses[[i]] <- calculate_mse_coefs(lars_fits, true_betas, i)
  }
  
  bind_rows(list(
    "backward_selection" = reduce(backward_mses, inner_join, by = "subset_size"), 
    "forward_selection" = reduce(forward_mses, inner_join, by = "subset_size"),
    "lars_selection" = reduce(lars_mses, inner_join, by = "subset_size")
  ), .id = "method")
}

# testing
if (FALSE) {
  res <- run_variable_selection(rho = 0, d = 15, N = 300, n_vars = 3,  seed = 420, nreps = 50) 
  
  res %>% 
    gather(var, val, -c(method, subset_size)) %>% 
    group_by(method, subset_size) %>% 
    summarise(mean_mse = mean(val), .groups = 'drop') %>% 
    ggplot(aes(subset_size, mean_mse, col = method)) + 
    geom_line() +
    geom_point()
  
}



