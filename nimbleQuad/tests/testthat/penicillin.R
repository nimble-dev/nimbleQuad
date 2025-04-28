library(nimbleQuad)

## Note that INLA's default prior is somewhat informative here in terms of intercept.
## Not clear it's useful to relax that as key here is to assess accuracy
## for model as constructed.

data(penicillin, package="faraway")

qpts <- c(.025,.25,.5,.75,.975)

code <- nimbleCode({
    for(i in 1:n) {
        mu[i] <- inprod(b[1:nTreat], x[i, 1:nTreat]) + re[blend[i]]
        y[i] ~ dnorm(mu[i], tau = Tau)
    }
    # Priors to match INLA
    Tau ~ dgamma(1, 5e-05)
    Tau_re ~ dgamma(1, 5e-05)
    for( i in 1:nTreat ){ b[i] ~ dnorm(0, tau = 0.001) }
    for( i in 1:nBlend ){ re[i] ~ dnorm(0, tau = Tau_re) }
})
X <- model.matrix(~treat, data = penicillin)
data = list(y = penicillin$yield)
constants = list(nTreat = 4, nBlend = 5, n = nrow(penicillin),
                 x = X, blend = as.numeric(penicillin$blend))
inits <- list(Tau = 1, Tau_re = 1, b = c(mean(data$y), rep(0,3)), re = rep(0,5))

m <- nimbleModel(code, data = data, constants = constants,
                 inits = inits, buildDerivs = TRUE)

cm <- compileNimble(m)
approx <- buildNestedApprox(model = m, hyperParamNodes = c('Tau', 'Tau_re'), latentNodes = c('b', 're'))
capprox <- compileNimble(approx, project = m)

result <- runNestedApprox(capprox)

## Seems somewhat better than INLA.
result$improveMarginals(nodes = 'Tau_re', nMarginalGrid=7)
result$improveMarginals(nodes = 'Tau', nMarginalGrid=7)

smp <- result$sampleLatentNodes(n=100000)
apply(smp, 2, quantile, qpts)

## Use more accurate outer grid, this gets fixed effects well-aligned with MCMC.
## Random effects seem too uncertain, while INLA is too certain (but closer to MCMC).

m <- nimbleModel(code, data = data, constants = constants,
                 inits = inits, buildDerivs = TRUE)

cm <- compileNimble(m)
approx <- buildNestedApprox(model = m, hyperParamNodes = c('Tau', 'Tau_re'), latentNodes = c('b', 're'),
                            control = list(nQuadOuter=7))
capprox <- compileNimble(approx, project = m)

smp <- result$sampleLatentNodes(n=100000)
apply(smp, 2, quantile, qpts)


## MCMC
library(nimbleHMC)

cm <- compileNimble(m)
mcmc <- buildHMC(m, monitors = c('Tau_re','Tau','b','re'))
cmcmc <- compileNimble(mcmc, project=m)
out <- runMCMC(cmcmc, niter=51000,nburnin=1000)

apply(out, 2, quantile, qpts)

## INLA
library(INLA)
formula <- yield ~ treat + f(blend, model="iid")
fit <- inla(formula, family = "gaussian", data=penicillin, 
	control.compute=list(config = TRUE), quantiles = qpts)

summary(fit)

fit$summary.random
