qpts <- c(.025,.25,.5,.75,.975)

## should we explore non-centered version with a single parameter?

## d=2 case with 8 latents

code <- nimbleCode({
    for(j in 1:J) {
        for(i in 1:n)
            y[i,j] ~ dpois(exp(lambda[j]))
        lambda[j] ~ dnorm(mu, sd=sqrt(1/tau))
    }
    mu ~ dnorm(0,.001)
    tau~dgamma(1, rate=5e-5)
})


set.seed(1)
n <- 10
J <- 8
lambda <- rnorm(J)
mns <- rep(exp(lambda), each = n)
y <- matrix(rpois(n*J, mns), ncol = J)
m <- nimbleModel(code, data = list(y = y), constants = list(n=n, J=J),
                 inits = list(lambda = rep(0,J), mu = 0, tau=1), buildDerivs = TRUE)

approx <- buildNestedApprox(m, latentNodes = c('lambda'), hyperParamNodes = c('mu','tau'))
cm <- compileNimble(m)
capprox <- compileNimble(approx, project = m)
result <- runNestedApprox(capprox)

result
## MLL: -129.9792
## MLL (INLA): -130.019


# `improveMarginals` is definitely needed and now results pretty close to INLA
result$improveMarginals(c('mu','tau'), nMarginalGrid = 7) 

latent_sample <- result$sampleLatentNodes(100000)
apply(latent_sample, 2, quantile, qpts)
# lambda values not nearly as good as INLA marginals
## MLL: -129.9331

## Do better outer approximation

approx <- buildNestedApprox(m, latentNodes = c('lambda'), hyperParamNodes = c('mu','tau'),
                            control = list(nQuadOuter=5))
cm <- compileNimble(m)
capprox <- compileNimble(approx, project = m)
system.time(result <- runNestedApprox(capprox))  
latent_sample <- result$sampleLatentNodes(100000)
apply(latent_sample, 2, quantile, qpts)



## non-centered - perhaps closer to INLA's treatment of only the group precision as hyperparameter?
## This breaks conditional independence so presumably less accurate.

code <- nimbleCode({
    for(j in 1:J) {
        for(i in 1:n)
            y[i,j] ~ dpois(exp(lambda[j]))
        lambda[j] ~ dnorm(mu, sd=sqrt(1/tau))
    }
    mu ~ dnorm(0,.001)
    tau~dgamma(1, rate=5e-5)
})


set.seed(1)
n <- 10
J <- 8
lambda <- rnorm(J)
mns <- rep(exp(lambda), each = n)
y <- matrix(rpois(n*J, mns), ncol = J)
m <- nimbleModel(code, data = list(y = y), constants = list(n=n, J=J),
                 inits = list(lambda = rep(0,J), mu = 0, tau=1), buildDerivs = TRUE)

approx <- buildNestedApprox(m, latentNodes = c('mu','lambda'), hyperParamNodes = c('tau'))
cm <- compileNimble(m)
capprox <- compileNimble(approx, project = m)
result <- runNestedApprox(capprox)

result  # rather close to INLA than other parameterization/latent-hyper split
## MLL: -129.9342

result$improveMarginals('tau', nMarginalGrid = 7)  # definitely needed, now pretty close to INLA

latent_sample <- result$sampleLatentNodes(100000)
apply(latent_sample, 2, quantile, qpts)  # not all that close to INLA for `mu`, but is close to uncorrected INLA latent samples for 'mu'
## lambda values not nearly as good as INLA marginals

## -129.926

reparam <- latent_sample[,2:9] - latent_sample[,1]
apply(reparam, 2, quantile, qpts)  # These look pretty good.

## Try better approximation.
approx <- buildNestedApprox(m, latentNodes = c('mu','lambda'), hyperParamNodes = c('tau'),
                            control = list(nQuadOuter=7))
cm <- compileNimble(m)
capprox <- compileNimble(approx, project = m)
result <- runNestedApprox(capprox)

latent_sample <- result$sampleLatentNodes(100000)
apply(latent_sample, 2, quantile, qpts)  


## INLA

library(INLA)
group <- as.factor(rep(1:J, each = n))
yc <- c(y)
formula <- y ~ 1 + f(group, model = "iid")
fit <- inla(formula, family="poisson", data=data.frame(y=yc,group=group), quantiles = qpts,
            control.compute=list(config = TRUE),
            control.fixed = list(prec.intercept = .001))
summary(fit)
fit$mlik # -129.541, -129.947

fit$summary.random

sampled <- inla.posterior.sample(n = 100000, fit)
smp <- t(sapply(sampled, function(x) x$latent[81:89,1]))
apply(smp, 2, quantile, qpts) # INLA's samples rather better than ours, perhaps because of mean+skew correction

# Without correction - doesn't seem to make a big difference to do the corrections though perhaps a bit worse, but ours are still rather worse.
sampled <- inla.posterior.sample(n = 100000, fit, skew.corr = FALSE, use.improved.mean = FALSE)
smp2 <- t(sapply(sampled, function(x) x$latent[81:89,1]))
apply(smp2, 2, quantile, qpts) 


## MCMC

code <- nimbleCode({
    for(j in 1:J) {
        for(i in 1:n)
            y[i,j] ~ dpois(exp(mu + lambda[j]))
        lambda[j] ~ dnorm(0, sd=sqrt(1/tau))
    }
    mu ~ dnorm(0,.001)
    tau~dgamma(1, rate=5e-5)
})


set.seed(1)
n <- 10
J <- 8
lambda <- rnorm(J)
mns <- rep(exp(lambda), each = n)
y <- matrix(rpois(n*J, mns), ncol = J)
m <- nimbleModel(code, data = list(y = y), constants = list(n=n, J=J),
                 inits = list(lambda = rep(0,J), mu = 0, tau=1), buildDerivs = TRUE)
cm <- compileNimble(m)
mcmc <- buildHMC(m, monitors = c('mu','tau','lambda'))
cmcmc <- compileNimble(mcmc, project=m)
out <- runMCMC(cmcmc, niter=51000,nburnin=1000)

apply(out[,c('mu','tau')], 2, quantile, qpts)


tmp <- out[, grepl("lambda", colnames(out))]

apply(tmp, 2, quantile, qpts)

## centered
tmp <- tmp + out[,'mu']
apply(tmp, 2, quantile, qpts)



### Noncentered but with fixed effect in parameters.

code <- nimbleCode({
    for(j in 1:J) {
        for(i in 1:n)
            y[i,j] ~ dpois(exp(mu+ lambda[j]))
        lambda[j] ~ dnorm(0, sd=sqrt(1/tau))
    }
    mu ~ dnorm(0,.001)
    tau~dgamma(1, rate=5e-5)
})


set.seed(1)
n <- 10
J <- 8
lambda <- rnorm(J)
mns <- rep(exp(lambda), each = n)
y <- matrix(rpois(n*J, mns), ncol = J)
m <- nimbleModel(code, data = list(y = y), constants = list(n=n, J=J),
                 inits = list(lambda = rep(0,J), mu = 0, tau=1), buildDerivs = TRUE)

approx <- buildNestedApprox(m, latentNodes = c('lambda'), hyperParamNodes = c('mu','tau'),
                            control = list(nQuadOuter=7))
cm <- compileNimble(m)
capprox <- compileNimble(approx, project = m)
result <- runNestedApprox(capprox)

latent_sample <- result$sampleLatentNodes(100000)
apply(latent_sample, 2, quantile, qpts) # Not any better.
