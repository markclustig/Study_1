# Autoregressive Latent Trajectory (ALT) Model (Bollen & Curran)

model <- '
  # define intercept and slope
  # model t1 as predetermined (i.e., not derived from int and slope)
  int =~ 1*t2 + 1*t3 + 1*t4 + 1*t5
  slope =~ 1*t2 + 2*t3 + 3*t4 + 5*t5
  
  int ~ 1
  slope ~ 1
  t1 ~ 1
  
  # latent variable variance
  int ~~ int
  slope ~~ slope
  t1 ~~ t1

  # latent variable covariance
  int ~~ slope
  int ~~ t1
  slope ~~ t1
  
  # autoregressive effects
  t2 ~ t1
  t3 ~ t2
  t4 ~ t3
  t5 ~ t4
  
  # error variances
  t2 ~~ t2
  t3 ~~ t3
  t4 ~~ t4
  t5 ~~ t5
'


