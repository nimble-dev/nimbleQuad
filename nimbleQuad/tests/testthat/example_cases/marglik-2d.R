library(nimbleQuad)
code <- nimbleCode({
    for(j in 1:J) {
        for(i in 1:n)
            y[i,j] ~ dnorm(mu[j], sd = sigma)
        mu[j] ~ dnorm(0, sd = tau)
    }
    sigma ~ dgamma(4, 1)
    tau ~ dgamma(1, 2)
})


set.seed(1)
n <- 10
J <- 8
mu <- rnorm(J)
sigma <- 1
mns <- rep(mu, each = n)
y <- matrix(rnorm(n*J, mns, sigma), ncol = J)
m <- nimbleModel(code, data = list(y = y), constants = list(n=n, J=J),
                 inits = list(mu = rep(0,J), tau = 1, sigma = 1), buildDerivs = TRUE)

approx <- buildNestedApprox(m, latentNodes = c('mu'), paramNodes = c('sigma','tau'))
cm <- compileNimble(m)
capprox <- compileNimble(approx, project = m)
result <- runNestedApprox(capprox) # -123.0872

latent_sample <- result$sampleLatents(100)
result  # -123.0761; could also compare to CCD-based estimate

## INLA MLL: -123.112

result$calcMarginalLogLikImproved()   # -123.0738


approx <- buildNestedApprox(m, latentNodes = c('mu'), paramNodes = c('sigma','tau'), control = list(nQuadOuter = 9))
cm <- compileNimble(m)
capprox <- compileNimble(approx, project = m)
result <- runNestedApprox(capprox) # -123.0872

latent_sample <- result$sampleLatents(100)
result  # -123.0738


M <- 5000000
py <- rep(0, M)

set.seed(1)
for(m in seq_along(py)) {
    sigma <- rgamma(1, 4,1)
    tau <- rgamma(1,1,2)
    cv <- matrix(tau^2, n, n)
    diag(cv) <- diag(cv) + sigma^2
    ch <- chol(cv)
    logpy <- sapply(1:J, function(idx)
        dmnorm_chol(y[(1+(idx-1)*n):(idx*n)], 0, ch, prec_param = FALSE, log = TRUE))
    py[m] <- exp(sum(logpy)) #  + dgamma(tau, 1,2,log=TRUE) + dgamma(sigma, 4,1,log=TRUE))
    }

log(mean(py)) # -123.0749

## reparameterize for inla comparison

library(nimbleQuad)
code <- nimbleCode({
    for(j in 1:J) {
        for(i in 1:n)
            y[i,j] ~ dnorm(mu[j], sigma)
        mu[j] ~ dnorm(0, tau)
    }
    sigma ~ dgamma(1, rate = 5e-5)
    tau ~ dgamma(1, rate = 5e-5)
})


set.seed(1)
n <- 10
J <- 8
mu <- rnorm(J)
sigma <- 1
mns <- rep(mu, each = n)
y <- matrix(rnorm(n*J, mns, sigma), ncol = J)
m <- nimbleModel(code, data = list(y = y), constants = list(n=n, J=J),
                 inits = list(mu = rep(0,J), tau = 1, sigma = 1), buildDerivs = TRUE)

approx <- buildNestedApprox(m, latentNodes = c('mu'), paramNodes = c('sigma','tau'))
cm <- compileNimble(m)
capprox <- compileNimble(approx, project = m)
result <- runNestedApprox(capprox) # -136.7465 (-136.788 w/ INLA-style calculation - see issue 84)

latent_sample <- result$sampleLatents(100)
result  # -136.7309,   -136.7267 with nQuad=11

M <- 100000
py <- rep(0, M)

set.seed(1)
for(m in seq_along(py)) {
    sigma <- rgamma(1, 1,rate=5e-5)
    tau <- rgamma(1, 1,rate=5e-5)
    cv <- matrix(1/tau, n, n)
    diag(cv) <- diag(cv) + 1/sigma
    ch <- chol(cv)
    logpy <- sapply(1:J, function(idx)
        dmnorm_chol(y[(1+(idx-1)*n):(idx*n)], 0, ch, prec_param = FALSE, log = TRUE))
    py[m] <- exp(sum(logpy)) 
}

log(mean(py)) #  -142.2013, -142.2454
## perhaps flat priors explain discrepancies?
## similarity to INLA is reassuring, plus example at top with non-flat priors indicate we get the right answer

qpts <- c(.025,.25,.5,.75,.975)
library(INLA)
group <- as.factor(rep(1:J, each = n))
yc <- c(y)
formula <- y ~ -1 + f(group, model = "iid")
fit <- inla(formula, family="normal", data=data.frame(y=yc,group=group), quantiles = qpts,
            control.compute=list(config = TRUE))
summary(fit) # -136.76

inla.priors.used(fit)

fit$summary.random

