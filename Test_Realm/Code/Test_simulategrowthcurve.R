# Simulate data from Latent Curve Model
library(MASS)
library(lavaan)

# Simulate from matrix expression ####

# define matrix of factor loadings
lambda <- matrix(c(1,1,1,1,0,1,2,3), ncol = 2)

# define mean vector of intercept and slope
mu <- c(1, .5)

# calculate expected values of dependent variable
E_y <- lambda %*% mu

# define covariance matrix of intercept and slope
psi <- matrix(c(1, .3, .3, 1.5), ncol = 2)

# define (diagonal) matrix of residual variances 
theta_eps <- diag(rep(.2, 4))

# calculate covariance matrix
sigma <- lambda %*% psi %*% t(lambda) + theta_eps


# draw values from multivariate normal distribution
data <- mvrnorm(n = 10000, mu = E_y, Sigma = sigma)
data <- as.data.frame(data)
names(data) <- c("t1", "t2", "t3", "t4")
summary(data)

# fit LCM

model <- '
  int =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
  slope =~ 0*t1 + 1*t2 + 2*t3 + 3*t4
  
  int ~ 1
  slope ~ 1
  
  # specify variances and covariances
  int ~~ int
  slope ~~ slope
  int ~~ slope
  t1 ~~ t1
  t2 ~~ t2
  t3 ~~ t3
  t4 ~~ t4
  
  # constrain error variances to be equal
  # t1 ~~ epsilon*t1
  # t2 ~~ epsilon*t2
  # t3 ~~ epsilon*t3
  # t4 ~~ epsilon*t4
'

fit <- lavaan(model = model, data = data)
summary(fit)



# case-by-case approach ####

# define mean intercept (vector, entry for each participant)
n <- 10000
mu_alpha <- rep(1, n)

# define mean slope
mu_beta <- rep(.5, n)

# draw individual intercept and slope deviations from normal distribution (random intercept)
# (define covariance matrix of intercept and slope)
psi <- matrix(c(1, .3, .3, 1.5), ncol = 2)
zeta <- as.data.frame(mvrnorm(n = n, mu = c(0,0), Sigma = psi))
names(zeta) <- c("alpha", "beta")
# draw error terms 
eps_it <- matrix(rnorm(4*n, mean = 0, sqrt(.2)), ncol = 4)

# define measurement intervals
lambda <- data.frame(t1 = rep(0, n),
                     t2 = rep(1, n),
                     t3 = rep(2, n),
                     t4 = rep(3, n))

# calculate y_it
data2 <- mu_alpha + zeta$alpha + lambda * (mu_beta + zeta$beta) + eps_it

fit2 <- lavaan(model = model, data = data2)
summary(fit2)
