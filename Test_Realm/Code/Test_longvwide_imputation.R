# Compare long and wide format imputation
library(lavaan)
library(tidyverse)
library(mice)
library(mitml)
library(miceadds)



# load data
data(Demo.growth)

# define missingness-inducing function
induce_missing <- function(vec, x, p = .25, b){
  # vec = variable in which to induce missingness
  # x = Prädiktorvariable (z.B. einzelne Variable oder Summenwert)
  # p = Wahrscheinlichkeit fehlender Werte
  # b = Effekt des Prädiktors (MAR-Stärke)
  
  z <- (x - mean(x)) / sd(x)
  e <- rnorm(n = length(z), mean = 0, sd = sqrt(1 - b^2))
  r <- qnorm(p) + z * b + e
  
  vec[which(r > 0)] <- NA
  
  return(vec)
}

# add subject variable
Demo.growth$subject <- 1:nrow(Demo.growth)

# induce missingness
Demo.growth_mis <- Demo.growth %>% 
  mutate(t1 = induce_missing(t1, x = x1, p = .25, b = 0),
         t2 = induce_missing(t1, x = x1, p = .25, b = 0),
         t3 = induce_missing(t1, x = x1, p = .25, b = 0),
         t4 = induce_missing(t1, x = x1, p = .25, b = 0))

# convert to long format
Demo.growth_mis_long <- Demo.growth_mis %>% 
  pivot_longer(cols = -c("subject", "x1", "x2"), names_to = c(".value", "time"), names_sep = 1) %>% 
  mutate(time = as.numeric(time) -1) %>% 
  rename(y = t)

# Convert complete data to long format
Demo.growth_long <- Demo.growth %>% 
  pivot_longer(cols = -c("subject", "x1", "x2"), names_to = c(".value", "time"), names_sep = 1) %>% 
  mutate(time = as.numeric(time) -1) %>% 
  rename(y = t)

# Analysis of complete data
results_complete <- lmer(y ~ time*x1 + time*x2 + (1|subject), data = Demo.growth_long)
summary(results_complete)

# Imputation ####

# Imputation in wide format

pred <- make.predictorMatrix(Demo.growth_mis)
meth <- make.method(Demo.growth_mis)
meth[c("t1", "t2", "t3", "t4")] <- "norm"


imp_wide <- mice(data = Demo.growth_mis, method = meth, predictorMatrix = pred, remove_collinear=FALSE)

imp_wide_elongated <- imp_wide %>% 
  mids2datlist() %>% 
  lapply(., function(z){
    pivot_longer(data = z, cols = -c("subject", "x1", "x2"), names_to = c(".value", "time"), names_sep = 1) %>% 
      mutate(time = as.numeric(time) -1) %>% 
      rename(y = t)
  }) 

fit_wide_elongated <- lapply(imp_wide_elongated, function(z){
    lmer(y ~ time*x1 + time*x2 + (1|subject), data = z)
    })

results_wide_elongated <- testEstimates(fit_wide_elongated)

summary(results_wide_elongated)


 
# fit_wide <- imp_wide %>% 
#   mids2datlist() %>% 
#   with.mitml.list(sem(model = ))

