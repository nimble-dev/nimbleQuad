qpts <- c(.025,.25,.5,.75,.975)

## d=3 non-normal case with 8 latents
## mu could also be in latents, in either centered or non-centered parameterization

## Use Gelman-style prior; this fits well both in terms of MCMC and our nested approx.
## Not sure offhand if/how to get INLA prior to match.

code <- nimbleCode({
    for(j in 1:J) {
        for(i in 1:n)
            ## This gamma likelihood is based on INLA's parameterization.
            y[i,j] ~ dgamma(mean = exp(eta[j]), sd = sqrt(exp(eta[j])^2/phi))
        eta[j] ~ dnorm(mu, sd = sigma)
    }
    mu ~ dnorm(0,.001) # INLA prior
    sigma ~ dhalfflat()  # not INLA prior
    phi ~ dgamma(1, rate = .01) # INLA prior
})

set.seed(1)
n <- 10
J <- 8
eta <- rnorm(J)
phi <- 0.5
mns <- rep(exp(eta), each = n)
sds <- rep(sqrt(exp(eta)^2/phi), each = n)
y <- matrix(rgamma(n*J, shape = mns^2/sds^2, rate = mns/sds^2), ncol = J)

m <- nimbleModel(code, data = list(y = y), constants = list(n=n, J=J),
                 inits = list(eta = rep(0,J), mu = 0, tau=1, sigma = 1), buildDerivs = TRUE)

approx <- buildNestedApprox(m, latentNodes = c('eta'), paramNodes = c('mu','sigma','phi'))
cm <- compileNimble(m)
capprox <- compileNimble(approx, project = m)
result <- runNestedApprox(capprox)
result # phi looks good but not the others

result$improveParamMarginals(c('mu','phi','sigma'), nMarginalGrid = 7)  # These look good.

latent_sample <- result$sampleLatents(100000)
apply(latent_sample, 2, quantile, qpts)  # Not terrible, but a bit off.

## MCMC

m <- nimbleModel(code, data = list(y = y), constants = list(n=n, J=J),
                 inits = list(eta = rep(0,J), mu = 0, sigma=1, phi = 1), buildDerivs = TRUE)
cm <- compileNimble(m)
conf <- configureMCMC(m, monitors = c('mu','sigma','phi','eta'), onlySlice = TRUE)
mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project=m)
out <- runMCMC(cmcmc, niter=50000, nburnin=0) ## mixing is good

apply(out[1001:50000, ], 2, quantile, qpts)

## INLA - note use of different prior for group-level variation.
## These results are generally quite far off. 
library(INLA)
group <- as.factor(rep(1:J, each = n))
yc <- c(y)
formula <- y ~ 1 + f(group, model = "iid")
fit <- inla(formula, family="gamma", data=data.frame(y=yc,group=group), quantiles = qpts,
            control.compute=list(config = TRUE))

fit$summary.random

## Now remove group-level variation

set.seed(1)
n <- 10
J <- 8
eta <- rep(0, J)
phi <- 0.5
mns <- rep(exp(eta), each = n)
sds <- rep(sqrt(exp(eta)^2/phi), each = n)
y <- matrix(rgamma(n*J, shape = mns^2/sds^2, rate = mns/sds^2), ncol = J)

m <- nimbleModel(code, data = list(y = y), constants = list(n=n, J=J),
                 inits = list(eta = rep(0,J), mu = 0, tau=1, sigma = 1), buildDerivs = TRUE)

approx <- buildNestedApprox(m, latentNodes = c('eta'), paramNodes = c('mu','sigma','phi'))
cm <- compileNimble(m)
capprox <- compileNimble(approx, project = m)
result <- runNestedApprox(capprox)
result

result$improveParamMarginals(c('mu','phi','sigma'), nMarginalGrid = 7)  # good

latent_sample <- result$sampleLatents(100000)
apply(latent_sample, 2, quantile, qpts)  # good

## MCMC

m <- nimbleModel(code, data = list(y = y), constants = list(n=n, J=J),
                 inits = list(eta = rep(0,J), mu = 0, sigma=1, phi = 1), buildDerivs = TRUE)
cm <- compileNimble(m)
conf <- configureMCMC(m, monitors = c('mu','sigma','phi','eta'), onlySlice = TRUE)
mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project=m)
out <- runMCMC(cmcmc, niter=50000, nburnin=0) # mixing is good

apply(out[1001:50000, ], 2, quantile, qpts)


## Mimic INLA's parameterization
## In general, this causes problems:
## MCMC has mixing problems and group variation concentrates near 0.
## INLA results don't seem great (and don't match MCMC).
## Nested approx gives inconsistent tau marginals as increase number of marginal points (and don't match MCMC).
## Hard to know what to conclude.

code <- nimbleCode({
    for(j in 1:J) {
        for(i in 1:n)
            y[i,j] ~ dgamma(mean = exp(eta[j]), sd = sqrt(exp(eta[j])^2/phi))
        eta[j] ~ dnorm(mu, tau)
    }
    mu ~ dnorm(0,.001)
    tau ~ dgamma(1, rate=5e-5)
    phi ~ dgamma(1, rate=.01)
})

set.seed(1)
n <- 10
J <- 8
eta <- rnorm(J)
phi <- 0.5
mns <- rep(exp(eta), each = n)
sds <- rep(sqrt(exp(eta)^2/phi), each = n)
y <- matrix(rgamma(n*J, shape = mns^2/sds^2, rate = mns/sds^2), ncol = J)

m <- nimbleModel(code, data = list(y = y), constants = list(n=n, J=J),
                 inits = list(eta = rep(0,J), mu = 0, tau=1, phi = 1), buildDerivs = TRUE)

approx <- buildNestedApprox(m, latentNodes = c('eta'), paramNodes = c('mu','tau','phi'))
cm <- compileNimble(m)
capprox <- compileNimble(approx, project = m)
result <- runNestedApprox(capprox)
result

result$improveParamMarginals(c('mu','phi','tau'), nMarginalGrid = 7)  # 'tau' is unstable and not good

latent_sample <- result$sampleLatents(100000)
apply(latent_sample, 2, quantile, qpts)

## MCMC

library(nimbleHMC)

cm <- compileNimble(m)
mcmc <- buildHMC(m, monitors = c('mu','tau','phi','eta'))
cmcmc <- compileNimble(mcmc, project=m)
out <- runMCMC(cmcmc, niter=51000, nburnin=1000) ## Gets stuck for long periods.

m <- nimbleModel(code, data = list(y = y), constants = list(n=n, J=J),
                 inits = list(eta = rep(0,J), mu = 0, tau=1, phi = 1), buildDerivs = TRUE)
cm <- compileNimble(m)
conf <- configureMCMC(m, monitors = c('mu','tau','phi','eta'), onlySlice = TRUE)
mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project=m)
out <- runMCMC(cmcmc, niter=50000, nburnin=0)  ## Mixing not good enough to use for comparison.


## INLA
library(INLA)
group <- as.factor(rep(1:J, each = n))
yc <- c(y)
formula <- y ~ 1 + f(group, model = "iid")
fit <- inla(formula, family="gamma", data=data.frame(y=yc,group=group), quantiles = qpts,
            control.compute=list(config = TRUE))

summary(fit)  # tau is quite small and not similar to MCMC

fit$summary.random


