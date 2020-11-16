if (!require("mvtnorm")) install.packages("mvtnorm")
library(tidyverse)


generate_data <- function(rho, n_vars, N, d) {
  ### Following the example in ESL p 59 figure 3.6
  # ----------------------------------------------
  # generate betas
  true_betas <- rep(0, d)
  non_zero_betas <- sample(1:d, n_vars)
  true_betas[non_zero_betas] <- rnorm(n_vars, 0, sqrt(0.4))
  names(true_betas) <- paste0("feature_", seq(1:d))
  
  # generate features
  sigma <- matrix(rho, nrow = d, ncol = d)
  diag(sigma) <- 1
  X_features <- mvtnorm::rmvnorm(N, mean = rep(0, d), sigma = sigma) 
  colnames(X_features) <- paste0("feature_", seq(1:d))
  X_df <- as_tibble(X_features)
  
  # calculate y and add noise
  y_mean <- X_features %*% true_betas 
  y_observed <- y_mean + rnorm(N, mean = 0, sd = sqrt(6.25))
  colnames(y_observed) <- "y_observed"
  y_df <- as_tibble(y_observed)
  
  # bind in df
  data <- bind_cols(y_df, X_df)
  true_betas_df <- true_betas %>% t() %>% as_tibble() %>% gather(variable_name, value)
  return(list("data" = data, "true_betas" = true_betas_df, "y_mean" = as.numeric(y_mean)))
}

extract_coefs_leaps <- function(leaps_output) {
  d <- leaps_output$np
  estimated_betas <- coef(leaps_output, 1:d) %>% 
    tibble() %>% 
    set_names("res") %>%
    mutate(res = map(res, ~as_tibble(t(.x)))) %>% 
    unnest(res) 
  
  estimated_betas
} 

calculate_mse_coefs <- function(coefs_leaps, beta_true, i) {
  d <- nrow(coefs_leaps)
  
  coefs_leaps %>% 
    mutate(subset_size = 1:d) %>% 
    gather(variable_name, estimate, -subset_size) %>% 
    mutate(estimate = replace_na(estimate, 0)) %>% 
    left_join(beta_true, by = "variable_name") %>% 
    mutate(ss = (estimate - value)^2) %>% 
    group_by(subset_size) %>% 
    summarise(mse = sum(ss), .groups = 'drop') %>% 
    dplyr::select(subset_size, !!paste0("mse_", i) := mse)
}
 
