#############################################
# PREP FOR STUDY 1                          #
# SIMULATE GROWTH WITH CORRELATED RESIDUALS #
#############################################

library(MASS)
library(lavaan)


# Simulate from matrix expression ####

# Set some variables 
# set number of repeated measurements
n_rep <- 15
# set number of observations
n_obs <- 500
# set residual variance
res_var <- 1
# set residual covariance
res_cov <- .5
# set latent variable covariance (intercept and slope)
lat_cov <- .3



# define matrix of factor loadings
lambda <- matrix(c(seq(n_rep) - 1, rep(1, times = n_rep) ), ncol = 2)

# define mean vector of intercept and slope
mu <- c(1, .5)

# calculate expected values of dependent variable
E_y <- lambda %*% mu

# define covariance matrix of intercept and slope
psi <- matrix(c(1, lat_cov, lat_cov, 1), ncol = 2)

# define rediual covariance matrix
# set diagonal (residual variance) to 1
theta_eps <- diag(x = res_var, nrow = n_rep, ncol = n_rep)

# set superdiagonal (residual covariance with next measurement) to some value
diag(theta_eps[-nrow(theta_eps),-1]) <- res_cov

# set subdiagonals (residual covariance with previous measurement) to some value
diag(theta_eps[-1, -ncol(theta_eps)]) <- res_cov
  
# calculate covariance matrix
sigma <- lambda %*% psi %*% t(lambda) + theta_eps


# draw values from multivariate normal distribution
data <- mvrnorm(n = n_obs, mu = E_y, Sigma = sigma)
data <- as.data.frame(data)
names(data) <- paste0("t", seq(n_rep))
head(data)

# fit LCM

model <- '
  int =~ 1*t1 + 1*t2 + 1*t3 + 1*t4 +1*t5
  slope =~ 0*t1 + 1*t2 + 2*t3 + 3*t4 + 4*t5

  int ~ 1
  slope ~ 1

  # specify variances and covariances
  int ~~ int
  slope ~~ slope
  int ~~ slope

  # error variances (free to vary and covary)
  t1 ~~ t1
  t1 ~~ t2
  t2 ~~ t2
  t2 ~~ t3
  t3 ~~ t3
  t3 ~~ t4
  t4 ~~ t4
  t4 ~~ t5
  t5 ~~ t5
'
# 
# # define function to construct model syntax for different numbers of measurements
# make_model <- function(n){
#   int <- paste0("1*", paste0("t", seq(n)), collapse = " + ")
#   slope <- paste0(seq(n)-1, "*", paste0("t", seq(n)), collapse = " + ")
#   
#   
# 
#   
# }


# misspecified model (no autocorrelation)

model_mis <- '
  int =~ 1*t1 + 1*t2 + 1*t3 + 1*t4 + 1*t5
  slope =~ 0*t1 + 1*t2 + 2*t3 + 3*t4 + 1*t5

  int ~ 1
  slope ~ 1

  # specify variances and covariances
  int ~~ int
  slope ~~ slope
  int ~~ slope

  # error variances (uncorrelated)
  t1 ~~ t1
  t2 ~~ t2
  t3 ~~ t3
  t4 ~~ t4
  t5 ~~ t5
'

fit <- lavaan(model = model, data = data)
summary(fit)
