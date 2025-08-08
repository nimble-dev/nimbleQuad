## Code for checking 1 latent, 1 parameter cases.
## These need to be transitioned into actual tests.

## d=1 case, posterior for variance, sd a bit skewed, but not near 0.
library(nimbleQuad)

dig00 <- nimbleFunction(
    run = function(x = double(0), log = logical(0, default = FALSE)) {
        returnType(double(0))
        if(log) return(-log(x)) else return(1/x)
    }, buildDerivs = TRUE)

registerDistributions(list(dig00=list(BUGSdist="dig00()",range=c(0,Inf))))

## can't use dinvgamma(0,0) as that gives density of 0

code <- nimbleCode({
    for(i in 1:n)
        y[i] ~ dnorm(mu, sd = sqrt(sigma2))
    mu ~ dflat()
    sigma2 ~ dig00() # Gelman Sec. 3.2: pi(sigma2) \propto 1/sigma2
})


set.seed(1)
n <- 30
m <- nimbleModel(code, data = list(y = rnorm(n)), constants = list(n=n),
                 inits = list(mu = 0, sigma2 = 1), buildDerivs = TRUE)

### marginal for sigma2 is IG((n-1)/2,(n-1)*s2/2)
### marginal for mu is t(ybar, s2/n)
qpts <- c(.025,.25,.5,.75,.975)

qs <- qinvgamma(qpts,(n-1)/2, scale=(n-1)*var(m$y)/2)

approx <- buildNestedApprox(m, latentNodes = 'mu', paramNodes = 'sigma2')
cm <- compileNimble(m)
capprox <- compileNimble(approx, project = m)
result <- runNestedApprox(capprox)

qs
result
result$improveParamMarginals('sigma2', nMarginalGrid = 5) # already much better
result$improveParamMarginals('sigma2', nMarginalGrid = 7)
result$improveParamMarginals('sigma2', nMarginalGrid = 9)
result$improveParamMarginals('sigma2', nMarginalGrid = 11)
result$improveParamMarginals('sigma2', nMarginalGrid = 31) # very good


sigma2_direct <- rinvgamma(10000,(n-1)/2, scale=(n-1)*var(m$y)/2)
sigma2_sample <- result$sampleParams(10000)
qqplot(sigma2_direct, sigma2_sample)


mu_sample <- result$sampleLatents(1000)
mu_direct <- rt_nonstandard(1000, n-1, mean(m$y), sqrt(var(m$y)/n))
qqplot(mu_direct, mu_sample)

mu_sample <- result$sampleLatents(1000, includeParams = TRUE)

## Check results on transformed scale.
result <- runNestedApprox(capprox, originalScale = FALSE)
log(qs)
result
result$improveParamMarginals(1, nMarginalGrid = 5)
result$improveParamMarginals(1, nMarginalGrid = 31)

logsigma2_sample <- result$sampleParams(10000)
qqplot(log(sigma2_direct), logsigma2_sample)

mu_sample <- result$sampleLatents(1000, includeParams = TRUE)


## One doesn't need inner AGHQ here given exact normality for latent, but just make sure it runs.

approx <- buildNestedApprox(m, latentNodes = 'mu', paramNodes = 'sigma2',
                            control = list(nQuadInner = 3))
cm <- compileNimble(m)
capprox <- compileNimble(approx, project = m)
result <- runNestedApprox(capprox)

mu_sample <- result$sampleLatents(100000)
quantile(mu_sample, qpts)
qt_nonstandard(qpts, n-1, mean(m$y), sqrt(var(m$y)/n))

## Use more accurate outer integration.
## This could in principle improve inference for latents,
## but we'd probably need to knock down the MC simulation
## by simulating more samples in `sampleLatents`.

approx <- buildNestedApprox(m, latentNodes = 'mu', paramNodes = 'sigma2',
                            control = list(nQuadOuter = 7))
cm <- compileNimble(m)
capprox <- compileNimble(approx, project = m)
result <- runNestedApprox(capprox)

mu_sample2 <- result$sampleLatents(100000)
quantile(mu_sample2, qpts)

## mimic INLA priors to compare to INLA

code <- nimbleCode({
    for(i in 1:n)
        y[i] ~ dnorm(mu, sd = sqrt(1/tau))
    mu ~ dnorm(0,.001)
    tau~dgamma(1, rate=5e-5)
})


set.seed(1)
n <- 30
m <- nimbleModel(code, data = list(y = rnorm(n)), constants = list(n=n),
                 inits = list(mu = 0, tau = 1), buildDerivs = TRUE)

approx <- buildNestedApprox(m, latentNodes = 'mu', paramNodes = 'tau')
cm <- compileNimble(m)
capprox <- compileNimble(approx, project = m)
result <- runNestedApprox(capprox)  # MLL=-55.09768  (-55.0941 by INLA arithmetic)

result

result$improveParamMarginals('tau', nMarginalGrid = 11)

mu_sample <- result$sampleLatents(100000)
apply(mu_sample, 2, quantile, qpts)
# MLL: -55.0941

cm <- compileNimble(m)
mcmc <- buildMCMC(m)
cmcmc <- compileNimble(mcmc, project = m)
out <- runMCMC(cmcmc,niter=1000000)

apply(out, 2, quantile, qpts)


## INLA

fit <- inla(y~1, family="gaussian", data=data.frame(y=m$y), quantiles = qpts, control.fixed = list(prec.intercept = .001))
summary(fit)

# MLL: -54.715, -55.098

## direct estimation is quite uncertain
n <- 1e6

mus <- rnorm(n,0,sd=sqrt(1/.001))
taus <- rgamma(n,1,rate=5e-5)

theta <- cbind(mus,sqrt(1/taus))
y <- m$y
logpy <- apply(theta, 1, function(x) sum(dnorm(y,x[1],x[2],log=T)))
log(mean(exp(logpy)))   # -63.8818

## Try with more constrained priors.

code <- nimbleCode({
    for(i in 1:n)
        y[i] ~ dnorm(mu, sd = sigma)
    mu ~ dnorm(0,sd=3)
    sigma~dunif(0,5)
})


set.seed(1)
n <- 30
y <- rnorm(n)
m <- nimbleModel(code, data = list(y = y), constants = list(n=n),
                 inits = list(mu = 0, sigma = 1), buildDerivs = TRUE)

approx <- buildNestedApprox(m, latentNodes = 'mu', paramNodes = 'sigma')
cm <- compileNimble(m)
capprox <- compileNimble(approx, project = m)
result <- runNestedApprox(capprox)

result
## MLL: -45.36005
## MLL INLA arithmetic: -45.3573
# MLL grid: -45.35235


## direct

n <- 1e7

mus <- rnorm(n,0,sd=3)
sigmas <- runif(n,0,5)

theta <- cbind(mus,sigmas)

set.seed(1)
logpy <- apply(theta, 1, function(x) sum(dnorm(y,x[1],x[2],log=T)))
log(mean(exp(logpy)))  # -45.34327  vs. -45.35235 for improved and -45.36005 for AG

## priors where I can use INLA too

code <- nimbleCode({
    for(i in 1:n)
        y[i] ~ dnorm(mu, tau)
    mu ~ dnorm(0, sd=3)
    tau ~ dgamma(1, 1)
})


set.seed(1)
n <- 30
m <- nimbleModel(code, data = list(y = rnorm(n)), constants = list(n=n),
                 inits = list(mu = 0, tau = 1), buildDerivs = TRUE)

approx <- buildNestedApprox(m, latentNodes = 'mu', paramNodes = 'tau')
cm <- compileNimble(m)
capprox <- compileNimble(approx, project = m)
result <- runNestedApprox(capprox)
result
## MLL: -44.04466
## INLA calc: -44.0411
## grid: -44.04109


## INLA

fit <- inla(y ~ 1,  family="gaussian", data=data.frame(y=y),quantiles = qpts,
               control.family = list(
                 hyper = list(
                   prec = list(
                       prior = "loggamma", param=c(1,1)))),
            control.fixed = list(prec.intercept = 1/9))

summary(fit)

fit$mlik # -43.66045 (integr), -44.04549 (Gaussian)

## direct

n <- 1e7

mus <- rnorm(n,0,sd=3)
sigmas <- sqrt(1/rgamma(n,1,1))

theta <- cbind(mus,sigmas)
y <- m$y
logpy <- apply(theta, 1, function(x) sum(dnorm(y,x[1],x[2],log=T)))
log(mean(exp(logpy)))  # -44.03872 (rather closer to nestedApprox than to INLA)



## Now have a case where the posterior for the parameter is not symmetric

code <- nimbleCode({
    for(i in 1:n)
        y[i] ~ dnorm(mu, sd = sqrt(sigma2))
    mu ~ dflat()
    sigma2 ~ dig00() # Gelman Sec. 3.2: pi(sigma2) \propto 1/sigma2
})


set.seed(1)
n <- 30
m <- nimbleModel(code, data = list(y = rnorm(n, 0, 0.1)), constants = list(n=n),
                 inits = list(mu = 0, sigma2 = 1), buildDerivs = TRUE)

### marginal for sigma2 is IG((n-1)/2,(n-1)*s2/2)
### marginal for mu is t(ybar, s2/n)
qpts <- c(.025,.25,.5,.75,.975)

qs <- qinvgamma(qpts,(n-1)/2, scale=(n-1)*var(m$y)/2)

approx <- buildNestedApprox(m, latentNodes = 'mu', paramNodes = 'sigma2')
cm <- compileNimble(m)
capprox <- compileNimble(approx, project = m)
result <- runNestedApprox(capprox)

qs
result
result$improveParamMarginals('sigma2', nMarginalGrid = 5) # already much better
result$improveParamMarginals('sigma2', nMarginalGrid = 7)
result$improveParamMarginals('sigma2', nMarginalGrid = 9)
result$improveParamMarginals('sigma2', nMarginalGrid = 11)
result$improveParamMarginals('sigma2', nMarginalGrid = 31) # very good


sigma2_direct <- rinvgamma(10000,(n-1)/2, scale=(n-1)*var(m$y)/2)
sigma2_sample <- result$sampleParams(10000)
qqplot(sigma2_direct, sigma2_sample)


mu_sample <- result$sampleLatents(1000)
mu_direct <- rt_nonstandard(1000, n-1, mean(m$y), sqrt(var(m$y)/n))
qqplot(mu_direct, mu_sample)
