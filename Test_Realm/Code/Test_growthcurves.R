# Tests for first study (MI in longitudinal data)

library(lavaan)
library(lme4)
library(tidyverse)

# load lavaan demo dataset
data(Demo.growth)
help("Demo.growth")
head(Demo.growth)

lav_model <- '
  
  # define intercept and slope as factors with fixed loadings
  int =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
  slope =~ 0*t1 + 1*t2 + 2*t3 + 3*t4
  
  # set regression equations
  int ~ x1 + x2
  slope ~ x1 + x2
  
  # constrain error variances to be equal
  t1 ~~ epsilon*t1
  t2 ~~ epsilon*t2
  t3 ~~ epsilon*t3
  t4 ~~ epsilon*t4
  
  # time-varying covariates
  # t1 ~ c1
  # t2 ~ c2
  # t3 ~ c3
  # t4 ~ c4
'

lav_fit <- growth(lav_model, data = Demo.growth)

summary(lav_fit)

# building the model using lavaan()
lav_model_test <- '
  # define intercept and slope as factors with fixed loadings
  int =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
  slope =~ 0*t1 + 1*t2 + 2*t3 + 3*t4
  
  # set regression equations
  int ~ x1 + x2
  slope ~ x1 + x2
  
  int ~ 1
  slope ~ 1
  
  # specify variances and covariances
  int ~~ int
  slope ~~ slope
  int ~~ slope
  # constrain error variances to be equal
  t1 ~~ epsilon*t1
  t2 ~~ epsilon*t2
  t3 ~~ epsilon*t3
  t4 ~~ epsilon*t4
'

# Fit the model
lav_fit_test <- lavaan(model = lav_model_test, data = Demo.growth, )

summary(lav_fit_test)


# Analyse data using a linear mixed model
# convert to long format
Demo.growth_long <- Demo.growth %>% 
  mutate(subject = 1:nrow(.)) %>% 
  pivot_longer(cols = -c("subject", "x1", "x2"), names_to = c(".value", "time"), names_sep = 1) %>% 
  mutate(time = as.numeric(time) -1)

lme_model <- lmer(t ~ time*x1 + time*x2 + (time|subject), data = Demo.growth_long, REML = F)

summary(lme_model)
