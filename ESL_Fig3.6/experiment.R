library(tidyverse)
source("utilities.R")
source("examples.R")

theme_set(theme_bw(base_size = 14))

nreps <- 100
N <- 200
rhos <- c(0, 0.3, 0.5, 0.7, 0.9)
dim_features<- c(30, 100, 200)
n_vars <- c(5, 10, 15)

results <- list()
set.seed(420)

for (d in dim_features) {
  for (rho in rhos) {
    for (n_var in n_vars) {
      print(paste ("d =", d, "rho =", rho, "n_var =", n_var))
      res <- run_variable_selection(rho = rho, d = d, N = N, n_vars = n_var, nreps = nreps) 
      readr::write_csv(res, paste0("res_rho", rho, "_n_var", n_var, "_dim_features", d, ".csv"))
      results[[as.character(rho)]][[as.character(n_var)]][[as.character(d)]] <- res
    }
  }
}

results <- tibble()
for (d in dim_features) {
  for (rho in rhos) {
    for (n_var in n_vars) {
      suppressMessages(
        res <- read_csv(paste0("res_rho", rho, "_n_var", n_var, "_dim_features", d, ".csv")) 
      )
      res$rho <- rho
      res$n_var <- n_var
      res$d <- d
      results <- bind_rows(results, res)
    }
  }
}

results_long <- results  %>% 
  gather(var, val, -c(method, subset_size, rho, d, n_var)) %>% 
  group_by(method, rho, n_var, d) %>% 
  # step-wise procedure can differ in number of steps.
  mutate(subset_size = d * (subset_size - min(subset_size)) / (max(subset_size) - min(subset_size))) %>% 
  ungroup()


plot_d <- function(dd) {
  results_long %>% 
    filter(d == dd, rho != 0.9, rho != 0.3) %>% 
    group_by(method, subset_size, rho, n_var) %>% 
    summarise(val = mean(val)) %>% 
    ungroup() %>%
    mutate(
      method = case_when(
        method == "backward_stepwise_selection"  ~ "Backward step-wise",
        method == "forward_stepwise_selection"   ~ "Forward step-wise",
        method == "lars_selection"  ~ "LARS",
        method == "lasso_selection"  ~ "LASSO",
        method == "stagewise_selection"  ~ "Stage-wise",
      )
    ) %>% 
    mutate(rho = paste0("Corr = ", rho*100, "%")) %>% 
    mutate(n_var = paste0("Non-zero betas = ", n_var)) %>% 
    mutate(n_var = forcats::fct_relevel(n_var, "Non-zero betas = 5", after=0)) %>% 
    ggplot(aes(subset_size, val, col = method, group)) +
    facet_grid(n_var~rho, scales = "free_y") +
    geom_line(size = 1) +
    # geom_point() +
    scale_color_brewer(type = "qual") +
    scale_x_continuous(breaks = seq(0, dd, by = as.integer(dd / 10))) +
    theme(legend.position = "bottom") +
    labs(x = "Subset size", y = paste0("MSE"), col ="")
}

plot_d(30) 
