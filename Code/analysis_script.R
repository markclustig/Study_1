# STUDY 1  #
# Analysis #

# Setup ####
here::i_am("Code/analysis_script.R")

library(tidyverse)
library(DescTools)
library(here)

# make output dir for plots
if (!dir.exists(here("Plots"))) {
  dir.create(here("Plots"))
}

# search for data files and read newest one
data_files <- list.files(path = here("Data"), pattern = "simulation_results_")
data_sim <- read_csv(here("Data", last(data_files))) %>% 
  mutate(condition = as.factor(condition),
         n_cases = as.factor(n_cases))

design_matrix <- read_csv(here("Data", "design_matrix.csv")) %>% 
  mutate(across(where(is.numeric), as.factor))

# elongate data (just a tad)
data_sim_longer <- data_sim %>%
  pivot_longer(
    cols = starts_with("mod_"), 
    names_pattern = "(mod_[a-zA-Z]+_)(.*)", 
    names_to = c("model", ".value")
  ) %>% 
  mutate(model = str_sub(model, end = -2))


# RESULTS #### 
## Bias ####
bias <- function(est, true) {
  mean(est - true)
}

results_bias <- data_sim_longer %>%
  group_by(condition, method, model, n_cases) %>% 
  summarize(
    bias.int = bias(`est_int~1`, int),
    bias.int_var = bias(`est_int~~int`, int_var),
    bias.slope = bias(`est_slope~1`, slope),
    bias.slope_var = bias(`est_slope~~slope`, slope_var),
    bias.lat_cov = bias(`est_int~~slope`, lat_cov),
    bias.res_var_1 = bias(`est_t1~~t1`, res_var),
    bias.res_var_2 = bias(`est_t2~~t2`, res_var),
    bias.res_var_3 = bias(`est_t3~~t3`, res_var),
    bias.res_var_4 = bias(`est_t4~~t4`, res_var),
    bias.res_var_5 = bias(`est_t5~~t5`, res_var),
    bias.res_var_6 = bias(`est_t6~~t6`, res_var),
    bias.res_var_7 = bias(`est_t7~~t7`, res_var),
    bias.res_var_8 = bias(`est_t8~~t8`, res_var),
    bias.res_cor_12 = bias(`est_t1~~t2`, res_cor),
    bias.res_cor_23 = bias(`est_t2~~t3`, res_cor),
    bias.res_cor_34 = bias(`est_t3~~t4`, res_cor),
    bias.res_cor_45 = bias(`est_t4~~t5`, res_cor),
    bias.res_cor_56 = bias(`est_t5~~t6`, res_cor),
    bias.res_cor_67 = bias(`est_t6~~t7`, res_cor),
    bias.res_cor_78 = bias(`est_t7~~t8`, res_cor)
  ) %>% 
  left_join(design_matrix)


## Coverage ####
coverage <- function(est, true, se) {
  ci_lo <- est - qnorm(.975) * se
  ci_up <- est + qnorm(.975) * se
  true %[]% cbind(ci_lo, ci_up) %>% 
    mean()
}

results_coverage <- data_sim_longer %>% 
  group_by(condition, method, model, n_cases) %>% 
  summarize(
    cover.int = coverage(`est_int~1`, int, `se_int~1`),
    cover.int_var = coverage(`est_int~~int`, int_var, `se_int~~int`),
    cover.slope = coverage(`est_slope~1`, slope, `se_slope~1`),
    cover.slope_var = coverage(`est_slope~~slope`, slope_var, `se_slope~~slope`),
    cover.lat_cov = coverage(`est_int~~slope`, lat_cov, `se_int~~slope`),
    cover.res_var_1 = coverage(`est_t1~~t1`, res_var, `se_t1~~t1`),
    cover.res_var_2 = coverage(`est_t2~~t2`, res_var, `se_t2~~t2`),
    cover.res_var_3 = coverage(`est_t3~~t3`, res_var, `se_t3~~t3`),
    cover.res_var_4 = coverage(`est_t4~~t4`, res_var, `se_t4~~t4`),
    cover.res_var_5 = coverage(`est_t5~~t5`, res_var, `se_t5~~t5`),
    cover.res_var_6 = coverage(`est_t6~~t6`, res_var, `se_t6~~t6`),
    cover.res_var_7 = coverage(`est_t7~~t7`, res_var, `se_t7~~t7`),
    cover.res_var_8 = coverage(`est_t8~~t8`, res_var, `se_t8~~t8`),
    cover.res_cor_12 = coverage(`est_t1~~t2`, res_cor, `se_t1~~t2`),
    cover.res_cor_23 = coverage(`est_t2~~t3`, res_cor, `se_t2~~t3`),
    cover.res_cor_34 = coverage(`est_t3~~t4`, res_cor, `se_t3~~t4`),
    cover.res_cor_45 = coverage(`est_t4~~t5`, res_cor, `se_t4~~t5`),
    cover.res_cor_56 = coverage(`est_t5~~t6`, res_cor, `se_t5~~t6`),
    cover.res_cor_67 = coverage(`est_t6~~t7`, res_cor, `se_t6~~t7`),
    cover.res_cor_78 = coverage(`est_t7~~t8`, res_cor, `se_t7~~t8`),
    .groups = "keep"
  ) %>% 
  left_join(design_matrix)


## Power ####
# Proportion of significant model comparisons (power or type-1-error rate)
results_power <- data_sim_longer %>% 
  mutate(significance = coalesce(`Pr(>Chisq)`, `P(>F)`)) %>% 
  group_by(condition, method, n_cases) %>% 
  summarize(power = mean(significance < .05)) %>% 
  left_join(design_matrix)


## Heywood you stop that ####
results_heywood <- data_sim_longer %>% 
  group_by(condition, method, model, n_cases) %>% 
  summarize(prop_trials_with_warnings = mean(warnings > 0),
            sum_warnings = sum(warnings)) %>% 
  left_join(design_matrix)

results_heywood_freqs <- data_sim_longer %>% 
  group_by(condition, method, model, n_cases) %>% 
  count(warnings) %>% 
  left_join(design_matrix)



# PLOTS ####
theme_set(theme_bw())

## Bias (plots) ####

plot_results <- function(var_name, df) {
  df %>% 
    filter(n_cases == 500) %>% 
    ggplot(aes(x = model, fill = method, color = method, group = method)) +
    scale_shape_manual(values = c(21, 22, 23)) +
    geom_linerange(aes(ymin = 0, ymax = eval(as.name(var_name))), 
                   position = position_dodge(width = .5), 
                   color = "lightgrey", linetype = 2) +
    geom_hline(yintercept = 0, alpha = .8) +
    geom_point(aes(y = eval(as.name(var_name)), shape = n_rep), 
               position = position_dodge(width = .5), 
               size = 3.5, alpha = .7) +
    facet_grid(rows = vars(mis_mech), cols = vars(res_cor)) +
    ggtitle(var_name) +
    ylab(var_name) 
}


names_bias <- names(results_bias) %>% 
  str_subset("bias.")

for (i in seq_along(names_bias)) {
  plot_results(names_bias[i], results_bias)
  ggsave(here("Plots", "Bias", 
              paste0(
                str_pad(i, 2, pad = "0"), names_bias[i], ".png")
              ), 
         height = 9, width = 16)
}


## Coverage (plots) ####

names_coverage <- names(results_coverage) %>% 
  str_subset("cover.")

for (i in seq_along(names_coverage)) {
  plot_results(names_coverage[i], results_coverage)
  ggsave(here("Plots", "Coverage", 
              paste0(
                str_pad(i, 2, pad = "0"), names_coverage[i], ".png")
              ), 
         height = 9, width = 16)
}

## Power (plots) ####
results_power %>%
  filter(n_cases == 500) %>% 
  ggplot(aes(x = as.factor(n_rep), color = method, group = method)) +
  geom_linerange(aes(ymin = 0, ymax = power), 
                 position = position_dodge(width = .5), 
                 color = "lightgrey", linetype = 2) +
  geom_point(aes(y = power), alpha = .7, position = position_dodge(width = .5)) +
  facet_grid(cols = vars(res_cor), rows = vars(mis_mech))
  

## Power (plots) ####
results_power %>%
  filter(n_cases == 500) %>% 
  ggplot(aes(x = as.factor(n_rep), fill = method, group = method)) +
  geom_col(aes(y = power), alpha = .7, position = "dodge") +
  facet_grid(cols = vars(res_cor), rows = vars(mis_mech))

ggsave(here("Plots", "Power", "power.png"), height = 9, width = 16)


## Heywood (plots) ####
results_heywood %>% 
  filter(n_cases == 500) %>% 
  ggplot(aes(x = model, fill = method, color = method, group = method)) +
  scale_shape_manual(values = c(21, 22, 23)) +
  geom_linerange(aes(ymin = 0, ymax = prop_trials_with_warnings), 
                 position = position_dodge(width = .5), 
                 color = "lightgrey", linetype = 2) +
  geom_hline(yintercept = 0, alpha = .8) +
  geom_point(aes(y = prop_trials_with_warnings, shape = n_rep), 
             position = position_dodge(width = .5), 
             size = 3.5, alpha = .7) +
  facet_grid(rows = vars(mis_mech), cols = vars(res_cor))

ggsave(here("Plots", "Heywood", "prop_trials_with_warnings.png"), height = 9, width = 16)


## Heywood (plots) ####
results_heywood %>% 
  filter(n_cases == 500) %>% 
  ggplot(aes(x = model, fill = method, color = method, group = method)) +
  scale_shape_manual(values = c(21, 22, 23)) +
  geom_linerange(aes(ymin = 0, ymax = sum_warnings), 
                 position = position_dodge(width = .5), 
                 color = "lightgrey", linetype = 2) +
  geom_hline(yintercept = 0, alpha = .8) +
  geom_point(aes(y = sum_warnings, shape = n_rep), 
             position = position_dodge(width = .5), 
             size = 3.5, alpha = .7) +
  facet_grid(rows = vars(mis_mech), cols = vars(res_cor))





