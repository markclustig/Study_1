# Study 1
# set.seed(2809)

# Measure start time
start_time <- Sys.time()

# Setup ####
library(tidyverse)
library(furrr)
library(MASS)
library(mice)
library(mitml)
library(miceadds)
library(lavaan)
library(jomo)
library(here)

# set up wd
i_am("Code/simulation_script.R")
# create output dir (if it does not exist)
if (!dir.exists(here("Data"))) {
  dir.create(here("Data"))
}

# set up parallel processing
plan(multisession)

# Define conditions ####
conditions <- list(
  n_rep = c(5, 8), # number of repeated measurements
  n_cases = c(100, 200, 500), # number of simulated participants
  res_var = 1, # residual variance (equal)
  res_cor = c(0, .1, .2), # residual covariance between rptd. measures (lag 1)
  int = 1,
  slope = .3,
  int_var = 1, # variance of latent intercept
  slope_var = .1, # variance of latent slope
  lat_cov = 0, # covariance of intercept and slope
  mis_mech = c("MAR", "MCAR") # missingness mechanism
)

design_matrix <- expand.grid(conditions) %>%
  mutate(condition = seq_len(nrow(.)))

write_csv(design_matrix, here("Data", "design_matrix.csv"))

# set number of repetitions
R <- 1000

# expand design_matrix so that it maps to every single simulated case
design_mapping <- do.call(
  "rbind",
  replicate(R, design_matrix, simplify = FALSE)
)


# Define models ####

# define model-building function
model_builder <- function(n_rep, add_res_cov = FALSE) {
  lat_factors <- data.frame(
    lhs = c("int", "slope", "int", "slope", "int"),
    op = c("~1", "~1", "~~", "~~", "~~"),
    rhs = c("", "", "int", "slope", "slope"),
    free = TRUE,
    ustart = NA
  )
  
  int_load <- data.frame(
    lhs = "int",
    op = "=~",
    rhs = paste0("t", seq(n_rep)),
    free = FALSE,
    ustart = 1
  )
  
  slope_load <- data.frame(
    lhs = "slope",
    op = "=~",
    rhs = paste0("t", seq(n_rep)),
    free = FALSE,
    ustart = seq(0, n_rep - 1)
  )
  
  rep_ints <- data.frame(
    lhs = paste0("t", seq(n_rep)),
    op = "~1",
    rhs = "",
    free = FALSE,
    ustart = 0
  )
  
  res_var <- data.frame(
    lhs = paste0("t", seq(n_rep)),
    op = "~~",
    rhs = paste0("t", seq(n_rep)),
    free = TRUE,
    ustart = NA
  )
  
  if (add_res_cov) {
    res_cov <- data.frame(
      lhs = paste0("t", seq(n_rep - 1) ),
      op = "~~",
      rhs = paste0("t", seq(n_rep - 1) + 1),
      free = TRUE,
      ustart = NA
      )
    rbind(lat_factors, int_load, slope_load, rep_ints, res_var, res_cov)
    } else rbind(lat_factors, int_load, slope_load, rep_ints, res_var)
  
}


models <- list()
models$noaut <- conditions$n_rep %>% 
  set_names() %>% 
  map(model_builder)
models$aut <- conditions$n_rep %>% 
  set_names() %>% 
  map(model_builder, add_res_cov = TRUE)

# _________________________________________________________________________ ####


# DATA SIMULATION ####

## Simulate complete data ####

# define simulation function
simulate <- function(n_rep, n_cases, res_var, res_cor, lat_cov, int, slope,
                     int_var, slope_var, mis_mech, condition) {
  # define matrix of factor loadings
  lambda <- matrix(c(rep(1, times = n_rep), seq(0, n_rep - 1)),
                   ncol = 2)

  # define mean vector of intercept and slope
  mu <- c(int, slope)

  # calculate expected values of dependent variable
  E_y <- lambda %*% mu

  # define covariance matrix of intercept and slope
  psi <- matrix(c(int_var, lat_cov, lat_cov, slope_var), ncol = 2)

  # define residual correlation matrix
  # set diagonal (residual variance) to 1
  theta_eps <- diag(x = 1, nrow = n_rep, ncol = n_rep)

  # set superdiagonal (residual correlation with next measurement)
  diag(theta_eps[-nrow(theta_eps), -1]) <- res_cor

  # set subdiagonal (residual correlation with previous measurement)
  diag(theta_eps[-1, -ncol(theta_eps)]) <- res_cor

  # (if applicable) calculate residual covariance matrix from correlation matrix
  if (res_var != 1) {
    theta_eps <- cor2cov(theta_eps, sds = rep(res_var, length(diag(theta_eps))))
  }

  # calculate covariance matrix
  sigma <- lambda %*% psi %*% t(lambda) + theta_eps

  # draw values from multivariate normal distribution
  data <- mvrnorm(n = n_cases, mu = E_y, Sigma = sigma) %>%
    as.data.frame()
  names(data) <- paste0("t", seq(n_rep))

  return(data)
}


# simulate data
data_comp <- future_pmap(design_mapping, simulate,
  .options = furrr_options(seed = TRUE)
)


## Induce missingness ####

# induce_monotone_MCAR <- function(data, p = .25) {
#   # Induces missingness such that only the last entries for each person may be
#   # missing.
#   # p: probability of missingness
#
#   # determine number of missing values per row
#   no_mis <- rbinom(n = nrow(data), size = ncol(data), prob = p)
#
#   # for each row, set the last `no_mis` values to NA
#   for (i in seq_len(nrow(data))) {
#     if (no_mis[i] > 0) {
#       data[i, (ncol(data) - no_mis[i] + 1):ncol(data)] <- NA
#     }
#   }
#   return(data)
#
#   # POTENTIAL ISSUE: THERE WILL BE (MUCH?) MORE VARIABILITY IN THE PROPORTION OF
#   # MISSING DATA BETWEEN RUNS WITH SMALL N AS COMPARED TO RUNS WITH LARGE N
# }
#


induce_binary_mis <- function(data, x = data[, 1], b = 0, p = .25,
                              frac_mis = .25, n_obs, mis_rel = TRUE) {
  # Induces missingness such that the last frac_mis * n_obs observations of some
  # subjects are missing, i.e., missingness is a "trait" of the subject.
  # x: Variable on which missingness depends. Defaults to the first col in data.
  # b: Magnitude of dependency of missingness on x. Defaults to MCAR.
  # p: Probability of missingness.
  # frac_mis: Fraction of missing observations per subject.
  # mis_rel: Whether a relative (as opposed to absolute) number of entries is
  # missing.

  # calculate missingness propensity for each participant
  z <- (x - mean(x)) / sd(x)
  e <- rnorm(n = length(z), mean = 0, sd = sqrt(1 - b^2))
  r <- qnorm(p) + z * b + e

  # calculate number of missing entries
  if (mis_rel) {
    n_mis <- round(frac_mis * ncol(data))
  } else {
    n_mis <- ncol(data) - n_obs
  }
  mis_entries <- (ncol(data) - n_mis + 1):ncol(data)

  data[which(r > 0), mis_entries] <- NA

  return(data)
}

data_mis <- data_comp %>%
  future_map_if(design_mapping$mis_mech == "MAR", induce_binary_mis,
    n_obs = 1, mis_rel = FALSE, b = .7,
    .options = furrr_options(seed = TRUE)
  ) %>%
  future_map_if(design_mapping$mis_mech == "MCAR", induce_binary_mis,
    n_obs = 1, mis_rel = FALSE, b = 0,
    .options = furrr_options(seed = TRUE)
  )


# _________________________________________________________________________ ####

# HELPER FUNCTIONS ####
## Define helper function to tidy up and extract results ####
extract_results <- function(df) {
  df %>%
    filter(str_detect(rowname, "int~|slope~|t[0-9]*~~")) %>%
    rename(any_of(
      c(
        est = "Estimate",
        se = "Std.Error"
      )
    )) %>%
    dplyr::select(rowname, est, se) %>%
    pivot_wider(names_from = rowname, values_from = c(est, se))
}

## Define replacement lavaan replacement with warnings ####

lavaan_with_warnings <- function(model, data, ...) {
  # returns list containing the lavaan result and an indicator of whether a
  # warning has occured during model fit
  
  lav_result <- lavaan(model = model, data = data, ...)
  
  lav_warning <- tryCatch(
    lavInspect(lav_result, "post.check"), 
    warning = function(w) paste0(w$message)
  )
  
  if (is.logical(lav_warning)) {
    lav_warning <- NA
    warning_happened <- 0
  } else {
    warning_happened <- 1
  }
  
  list(result = lav_result,
       warning_text = lav_warning,
       warning_happened = warning_happened)
  
}


## Define analysis function (imputed data) ####
analyze_mi <- function(mi_list, mod_autocor, mod_no_autocor) {
  
  # fit models
  mod_aut_result_with_warnings <- mi_list %>%
    lapply(., function(z) lavaan_with_warnings(data = z, model = mod_autocor))
  mod_noaut_result_with_warnings <- mi_list %>%
    lapply(., function(z) lavaan_with_warnings(data = z, model = mod_no_autocor))

  # extract model fit and lavaan warnings
  mod_aut_result <- mod_aut_result_with_warnings %>% 
    lapply(., function(z) z$result)
  mod_aut_warnings <- mod_aut_result_with_warnings %>% 
    lapply(., function(z) z$warning_happened) %>%
    do.call(sum, .)
  names(mod_aut_warnings) <- "mod-aut_warnings"
  mod_aut_warningtext <- mod_aut_result_with_warnings %>% 
    lapply(., function(z) z$warning_text) %>% 
    str_flatten_comma(na.rm = TRUE)
  
  mod_noaut_result <- mod_noaut_result_with_warnings %>% 
    lapply(., function(z) z$result)
  mod_noaut_warnings <- mod_noaut_result_with_warnings %>% 
    lapply(., function(z) z$warning_happened) %>%
    do.call(sum, .)
  names(mod_noaut_warnings) <- "mod-noaut_warnings"
  mod_aut_warningtext <- mod_noaut_result_with_warnings %>% 
    lapply(., function(z) z$warning_text) %>% 
    str_flatten_comma(na.rm = TRUE)
    
  # model comparison
  model_test <- testModels(mod_aut_result, mod_noaut_result, method = "D3") %>%
    .[["test"]] %>%
    as.data.frame()

  # model parameters
  aut_results <- mod_aut_result %>%
    testEstimates() %>%
    .$estimates %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    extract_results() %>%
    rename_with(~ paste0("mod-aut_", .x))

  noaut_results <- mod_noaut_result %>%
    testEstimates() %>%
    .$estimates %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    extract_results() %>%
    rename_with(~ paste0("mod-noaut_", .x))

  cbind(aut_results, noaut_results, model_test, mod_aut_warnings, 
        mod_noaut_warnings, mod_aut_warningtext)
}


## Define analysis function (complete data, LD, and FIML) ####
analyze_df <- function(df, ...) {
  # provide arguments missing and (if missing = "fiml") fixed.x for lavaan

  n_rep <- ncol(df)
  mod_no_autocor <- models$noaut[[paste(n_rep)]]
  mod_autocor <- models$aut[[paste(n_rep)]]
  
  mod_aut_result_with_warnings <- lavaan_with_warnings(
    data = df, 
    model = mod_autocor, 
    ...
    )
  mod_noaut_result_with_warnings <- lavaan_with_warnings(
    data = df, 
    model = mod_no_autocor, 
    ...
  )
  
  mod_aut_result <- mod_aut_result_with_warnings$result
  mod_aut_warnings <- mod_aut_result_with_warnings$warning_happened
  names(mod_aut_warnings) <- "mod-aut_warnings"
  
  mod_noaut_result <- mod_noaut_result_with_warnings$result
  mod_noaut_warnings <- mod_noaut_result_with_warnings$warning_happened
  names(mod_noaut_warnings) <- "mod-noaut_warnings"
  
  model_test <- anova(mod_aut_result, mod_noaut_result) %>%
    .[2, ] %>%
    dplyr::select(AIC, BIC, `Chisq diff`, `Df diff`, `Pr(>Chisq)`)
  rownames(model_test) <- NULL

  aut_results <- mod_aut_result %>%
    parametertable() %>%
    mutate(rowname = paste0(lhs, op, rhs)) %>%
    extract_results() %>%
    rename_with(~ paste0("mod-aut_", .x))

  noaut_results <- mod_noaut_result %>%
    parametertable() %>%
    mutate(rowname = paste0(lhs, op, rhs)) %>%
    extract_results() %>%
    rename_with(~ paste0("mod-noaut_", .x))

  cbind(aut_results, noaut_results, model_test, mod_aut_warnings, mod_noaut_warnings)
}


## Define helper functions for data transformation (long-format imputation) ####
helper_longer <- function(df) {
  df %>%
    mutate(subject = seq_len(nrow(df))) %>%
    pivot_longer(
      cols = -"subject", names_to = c(".value", "time"),
      names_pattern = "(t)(\\d+)"
    ) %>%
    mutate(time = as.numeric(time) - 1)
}


helper_wider <- function(l) {
  l %>%
    lapply(function(x) {
      x %>%
        mutate(time = time + 1) %>%
        rename(any_of(c(subject = "clus"))) %>%
        pivot_wider(
          id_cols = c("subject"), names_from = time, values_from = t,
          names_prefix = "t"
        )
    })
}


# _________________________________________________________________________ ####


# ANALYSIS ####
## Complete-data analysis ####
list_cd <- data_comp %>%
  future_map(analyze_df, .options = furrr_options(seed = TRUE))


## Listwise deletion ####
list_ld <- data_mis %>%
  future_map(analyze_df,
    missing = "listwise",
    .options = furrr_options(seed = TRUE)
  )


## FIML ####
list_fiml <- data_mis %>%
  future_map(analyze_df,
    missing = "fiml", fixed.x = FALSE,
    .options = furrr_options(seed = TRUE)
  )


## Impute in wide format ####

### FCS (wide) ####
impute_wide_fcs <- function(df) {
  
  n_rep <- ncol(df)
  mod_no_autocor <- models$noaut[[paste(n_rep)]]
  mod_autocor <- models$aut[[paste(n_rep)]]
  
  pred <- make.predictorMatrix(df)
  meth <- make.method(df)
  meth[meth != ""] <- "norm"

  mice(
    data = df, method = meth, predictorMatrix = pred,
    m = 20, maxit = 20,
    printFlag = FALSE
  ) %>%
    mids2mitml.list() %>%
    analyze_mi(mod_autocor, mod_no_autocor)
}


list_mi_wide_fcs <- data_mis %>%
  future_map(impute_wide_fcs, .options = furrr_options(seed = TRUE))


### Jomo (wide) ####
impute_wide_jomo <- function(df) {
  
  n_rep <- ncol(df)
  mod_no_autocor <- models$noaut[[paste(n_rep)]]
  mod_autocor <- models$aut[[paste(n_rep)]]
  
  Y <- df
  X <- data.frame(rep(1, nrow(df)))

  jomo(
    Y = Y, X = X,
    nburn = 500, nbetween = 50, nimp = 20,
    out.iter = NULL, output = FALSE
  ) %>%
    jomo2mitml.list() %>%
    analyze_mi(mod_autocor, mod_no_autocor)
}

list_mi_wide_jomo <- data_mis %>%
  future_map(impute_wide_jomo, .options = furrr_options(seed = TRUE))

## Impute in long format ####

### Jomo SMC ####
# mod_jomo <- t ~ time + (time | subject)
# 
# impute_long_jomo.lmer <- function(df) {
#   
#   n_rep <- ncol(df)
#   mod_no_autocor <- models$noaut[[paste(n_rep)]]
#   mod_autocor <- models$aut[[paste(n_rep)]]
#   
#   df %>%
#     helper_longer() %>%
#     jomo.lmer(
#       formula = mod_jomo,
#       nimp = 10, nburn = 1000, nbetween = 500,
#       out.iter = NULL, output = FALSE
#     ) %>%
#     jomo2mitml.list() %>%
#     helper_wider() %>%
#     analyze_mi(mod_autocor, mod_no_autocor)
# }


# list_mi_long_jomoSMC <- data_mis %>%
#   future_map(impute_long_jomo.lmer, .options = furrr_options(seed = TRUE))


### Pan mitml ####
# mod_panImpute <- t ~ time + (time | subject)
# 
# impute_long_panImpute <- function(df) {
#   
#   n_rep <- ncol(df)
#   mod_no_autocor <- models$noaut[[paste(n_rep)]]
#   mod_autocor <- models$aut[[paste(n_rep)]]
#   
#   df %>%
#     helper_longer() %>%
#     panImpute(
#       formula = mod_panImpute,
#       n.iter = 500, n.burn = 1000, m = 10,
#       silent = TRUE
#     ) %>%
#     mitmlComplete(print = "all") %>%
#     helper_wider() %>%
#     analyze_mi(mod_autocor, mod_no_autocor)
# }
# 
# list_mi_long_panImpute <- data_mis %>%
#   future_map(impute_long_panImpute, .options = furrr_options(seed = TRUE))


### Jomo (long) ####
impute_long_jomo <- function(df) {
  
  n_rep <- ncol(df)
  mod_no_autocor <- models$noaut[[paste(n_rep)]]
  mod_autocor <- models$aut[[paste(n_rep)]]
  
  df_long <- df %>%
    helper_longer() %>%
    mutate(int = 1)

  Y <- as.data.frame(df_long[, "t"])
  X <- as.data.frame(df_long[, c("int", "time")])
  Z <- as.data.frame(df_long[, c("int", "time")])
  clus <- df_long[, "subject"]

  jomo(
    Y = Y, X = X, Z = Z, clus = clus,
    nburn = 1000, nbetween = 500, nimp = 10,
    out.iter = NULL, output = FALSE
  ) %>%
    jomo2mitml.list() %>%
    helper_wider() %>%
    analyze_mi(mod_autocor, mod_no_autocor)
}


list_mi_long_jomo <- data_mis %>%
  future_map(impute_long_jomo, .options = furrr_options(seed = TRUE))


### FCS (long) ####
impute_long_fcs <- function(df) {
  
  n_rep <- ncol(df)
  mod_no_autocor <- models$noaut[[paste(n_rep)]]
  mod_autocor <- models$aut[[paste(n_rep)]]
  
  df_long <- helper_longer(df)

  pred <- make.predictorMatrix(df_long)
  pred["t", "subject"] <- -2
  pred["t", "time"] <- 2
  meth <- make.method(df_long)
  meth["t"] <- "2l.pan"

  df_long %>%
    mice(
      method = meth, predictorMatrix = pred,
      m = 10, maxit = 10, paniter = 1000,
      printFlag = FALSE
    ) %>%
    mids2mitml.list() %>%
    helper_wider() %>%
    analyze_mi(mod_autocor, mod_no_autocor)
}


list_mi_long_fcs <- data_mis %>%
  future_map(impute_long_fcs, .options = furrr_options(seed = TRUE))


# _________________________________________________________________________ ####


# OUTPUT ####
df_cd <- bind_rows(list_cd) %>%
  cbind(design_mapping, .) %>%
  mutate(method = "cd")

df_ld <- bind_rows(list_ld) %>%
  cbind(design_mapping, .) %>%
  mutate(method = "ld")

df_fiml <- bind_rows(list_fiml) %>%
  cbind(design_mapping, .) %>%
  mutate(method = "fiml")

df_mi_wide_fcs <- bind_rows(list_mi_wide_fcs) %>%
  cbind(design_mapping, .) %>%
  mutate(method = "mi_wide_fcs")

df_mi_wide_jomo <- bind_rows(list_mi_wide_jomo) %>%
  cbind(design_mapping, .) %>%
  mutate(method = "mi_wide_jomo")

# df_mi_long_jomoSMC <- do.call("rbind", list_mi_long_jomoSMC) %>%
#   cbind(design_mapping, .) %>%
#   mutate(method = "mi_long_jomoSMC")

df_mi_long_jomo <- bind_rows(list_mi_long_jomo) %>%
  cbind(design_mapping, .) %>%
  mutate(method = "mi_long_jomo")

# df_mi_long_panImpute <- do.call("rbind", list_mi_long_panImpute) %>%
#   cbind(design_mapping, .) %>%
#   mutate(method = "mi_long_panImpute")

df_mi_long_fcs <- bind_rows(list_mi_long_fcs) %>%
  cbind(design_mapping, .) %>%
  mutate(method = "mi_long_fcs")


dfs_all <- bind_rows(mget(ls(pattern = "df_")), .id = "label")

## Export CSV
write_csv(
  dfs_all,
  here(
    "Data",
    paste0(
      "simulation_results_",
      format(Sys.time(), "%y-%m-%d_%H%M"),
      ".csv"
    )
  )
)


# Measure end time
end_time <- Sys.time()

difftime(end_time, start_time, unit = "hours")