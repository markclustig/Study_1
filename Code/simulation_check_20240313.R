# # # # # # # # # # # # # # # # # # #
#                                   #
#         Check Simulation          #
#                                   #
# # # # # # # # # # # # # # # # # # #

# * Setup

library(MASS)
library(lavaan)
library(mice)
library(jomo)
library(mitml)

set.seed(5930)

# sample sizes
n <- 1000
n_T <- 8

# parameters
mu <- c(0, 0.3)
Psi <- diag(c(1, 0.1)) # NOTE: 0.25 still fairly large (much smaller in other sims)
theta <- 1
rho <- 0.1

# missing data
mis_beta <- 0.50
mis_perc <- 0.25

# * Simulate data

# loadings
Lambda <- cbind(rep(1, n_T), seq.int(0, n_T - 1))
#Lambda <- cbind(rep(1, n_T), seq.int(0, n_T - 1) / (n_T - 1))

# residual covariance matrix
Theta <- matrix(0, nrow = n_T, ncol = n_T)
Theta_off_diag <- subset(expand.grid(i = 1:n_T, j = 1:n_T), abs(i - j) <= 1)
Theta[as.matrix(Theta_off_diag)] <- rho
diag(Theta) <- 1
Theta <- lavaan::cor2cov(Theta, sds = rep(sqrt(theta), n_T))

# model-implied (population) means and covariances
mu_Y <- as.vector(Lambda %*% mu)
Sigma_Y <- Lambda %*% Psi %*% t(Lambda) + Theta

Sigma_Y
cov2cor(Sigma_Y)

Y <- MASS::mvrnorm(n = n, mu = mu_Y, Sigma = Sigma_Y)

dat0 <- as.data.frame(Y)
names(dat0) <- paste0("y", 1:n_T)

# colMeans(dat0)   # means: OK (mu_Y)
# cov(dat0)        # covariances: OK (Sigma_Y)
# cor(dat0)        # correlations: OK (NOTE: can get fairly large)

# * Simulate missing data

# missingness indicator
r <- qnorm(mis_perc) + mis_beta * scale(dat0$y1)[,1] + rnorm(n, mean = 0, sd = sqrt(1 - mis_beta^2))

# delete values
dat1 <- dat0
dat1[r > 0, 2:n_T] <- NA

# f_mis <- glm(is.na(y2) ~ y1, data = dat1, family = binomial(link = probit))
# var(predict(f_mis)) / (var(predict(f_mis)) + 1)   # missing-R²: OK (mis_beta^2)
# colMeans(is.na(dat1))                             # missing-%: OK (mis_perc)

# * Specify analysis model

mod <- paste0(
  "i =~ ", paste0(paste0("1", "*", names(dat1)), collapse = " + "), "\n",
  "s =~ ", paste0(paste0(seq.int(0, n_T - 1), "*", names(dat1)), collapse = " + "), "\n",
  paste0(names(dat1)[-n_T], "~~", names(dat1)[-1], collapse = "; "), "\n"
)

# cat(mod)   # OK

# * CD

fit_CD <- lavaan::growth(model = mod, data = dat0, estimator = "ML")

# summary(fit_CD, rsquare = TRUE)  # R² etc.: OK (NOTE: can get fairly large)

# * LD

fit_LD <- lavaan::growth(model = mod, data = dat1, estimator = "ML", missing = "listwise")

# * FIML

fit_FIML <- lavaan::growth(model = mod, data = dat1, estimator = "ML", missing = "FIML", fixed.x = FALSE)

# * FCS

me <- vector("character", n_T)
pm <- matrix(0, nrow = n_T, ncol = n_T)
names(me) <- colnames(pm) <- rownames(pm) <- names(dat1)

me[2:n_T] <- "norm"
pm[2:n_T, ] <- 1
diag(pm) <- 0

init <- dat1
init_ind <- sample(which(!is.na(dat1$y2)), size = sum(is.na(dat1$y2)), replace = TRUE)
init[is.na(dat1$y2), 2:n_T] <- dat1[init_ind, 2:n_T]

imp_FCS <- mice::mice(data = dat1, method = me, predictorMatrix = pm, m = 20, maxit = 20)
# plot(imp_FCS, layout = c(2, n_T-1))  # convergence: VERY POOR (even with good starting values and 50 iterations)

lst_FCS <- mids2mitml.list(imp_FCS)
fit_FCS <- with(lst_FCS, lavaan::growth(model = mod, data = dat0, estimator = "ML"), include.data = TRUE)
testEstimates(fit_FCS)

# * JM

# imp_JM <- jomo::jomo1con(Y = as.matrix(dat1), X = matrix(1, nrow = n, ncol = 1), nburn = 500, nbetween = 50, nimp = 20)
# str(imp_JM)

imp_JM <- mitml::jomoImpute(data = dat1, type = rep(1, n_T), n.burn = 500, n.iter = 50, m = 20)
# summary(imp_JM)  # convergence: OK

lst_JM <- mitml::mitmlComplete(imp_JM)
fit_JM <- with(lst_JM, lavaan::growth(model = mod, data = dat0, estimator = "ML"), include.data = TRUE)
testEstimates(fit_JM)

# * Two-Level JM

dat2 <- reshape(data = dat1, varying = names(dat1), sep = "", direction = "long")[c("id", "time", "y")]
rownames(dat2) <- NULL

# subset(dat2, id == 2)  # compare cases: OK
# dat1[2,]

# imp_JM2 <- jomo::jomo1rancon(
#   Y = as.matrix(dat2$y), X = matrix(1, nrow = nrow(dat2), ncol = 1),
#   Z = matrix(1, nrow = nrow(dat2), ncol = 1), clus = dat2$id,
#   nburn = 500, nbetween = 50, nimp = 20
# )

imp_JM2 <- mitml::jomoImpute(data = dat2, type = c(-2, 2, 1), n.burn = 500, n.iter = 50, m = 20)
# summary(imp_JM2)  # convergence: OK

lst_JM2 <- mitml::mitmlComplete(imp_JM2)
lst_wide_JM2 <- lapply(lst_JM2, function(d) {
  dw <- reshape(d, direction = "wide")[-1]
  names(dw) <- sub("\\.", "", names(dw))
  dw
})
lst_wide_JM2 <- as.mitml.list(lst_wide_JM2)
fit_JM2 <- with(lst_wide_JM2, lavaan::growth(model = mod, data = dat0, estimator = "ML"), include.data = TRUE)
testEstimates(fit_JM2)



