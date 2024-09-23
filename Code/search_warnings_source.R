# Search for source of warnings()

test_jomo <- data_mis[1:40] %>%
  future_map(impute_wide_jomo, .options = furrr_options(seed = TRUE))



data_mis[1] %>%
  future_map(impute_wide_jomo, .options = furrr_options(seed = TRUE))


for (i in seq(8)){
  print(paste0("Imputing dataset", i))
  impute_wide_jomo(data_mis[[i]])
  warnings()
}


test_imputations <- function(df) {
  tryCatch(expr = impute_wide_jomo(df),
           warning = function(w) print(w)
           )
}


lapply(data_mis[1:40], test_imputations)



warnings <- c()

for (i in seq(10)) {
  tryCatch(
    expr = impute_wide_jomo(data_mis[[i]]),
    warning = function(w) {
      warning(w)
      warnings <<- append(warnings, i)
    }
  )
}

warningl <- vector("logical", 160)
warningl[warnings] <- T
warning_data <- cbind(design_mapping[1:160,], warningl)
warning_data <- warning_data[c("n_rep", "res_cor", "mis_mech", "warningl")]
warning_data %>% group_by(res_cor, n_rep, mis_mech) %>% summarize(n_warn = sum(warningl))

warnings_fcs <- c()

for (i in seq(160)) {
  tryCatch(
    expr = impute_wide_fcs(data_mis[[i]]),
    warning = function(w) {
      warning(w)
      warnings_fcs <<- append(warnings_fcs, i)
    }
  )
}

warningl_fcs <- vector("logical", 160)
warningl_fcs[warnings_fcs] <- T
warning_fcs_data <- cbind(design_mapping[1:160,], warningl_fcs)
warning_fcs_data <- warning_fcs_data[c("n_rep", "res_cor", "mis_mech", "warningl_fcs")]
warning_fcs_data %>% group_by(res_cor, n_rep, mis_mech) %>% summarize(n_warn = sum(warningl_fcs))





# FIML and LD
test_fiml <- data_mis[1:40] %>% 
  future_map(analyze_df, missing = "fiml", fixed.x = F)


test_ld <- data_mis[1:40] %>% 
  future_map(analyze_df, missing = "listwise")
