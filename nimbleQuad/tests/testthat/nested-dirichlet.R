## INLA/BUGS-examples Dirichlet case
## https://inla.r-inla-download.org/r-inla.org/doc/likelihood/pom.pdf

library(nimbleQuad, lib.loc='/tmp/nq-fgl')

qpts <- c(.025,.25,.5,.75,.975)

library(brms)
data(inhaler)

## Mimic INLA parameterization
code <- nimbleCode({
    psi[1:K] ~ ddirch(threes[1:K])
    alpha[1] <- logit(psi[1])
    alpha[2] <- logit(psi[1]+psi[2])
    alpha[3] <- logit(psi[1]+psi[2]+psi[3])
    for (i in 1:n) {
        rating[i] ~ dcat(p[i,1:K])
        eta[i] <- beta_int + beta_treat*treat[i] + beta_period*period[i] + beta_carry*carry[i]
        for(k in 1:(K-1)) {
            gamma[i,k] <- alpha[k] - eta[i]
            F[i,k] <- expit(gamma[i,k])
        }
        p[i,1] <- F[i,1]
        p[i,2] <- F[i,2] - F[i,1]
        p[i,3] <- F[i,3] - F[i,2]
        p[i,4] <- 1-F[i,3]
    }
    beta_int ~ dflat() ## intercept has flat prior in INLA
    beta_treat ~ dnorm(0, .001)
    beta_period ~ dnorm(0, .001)
    beta_carry ~ dnorm(0, .001)
})
        

K <- 4
m <- nimbleModel(code, data = list(rating = inhaler$rating),inits = list(psi = rep(.25, 4), beta_int = 0, beta_treat = 0, beta_period = 0, beta_carry = 0),
                 constants = list(K = K, n = nrow(inhaler), period = inhaler$period, carry = inhaler$carry, treat = inhaler$treat, threes = rep(3, K)), buildDerivs = TRUE)

approx <- buildNestedApprox(m, latentNodes = c('beta_int','beta_treat','beta_period','beta_carry'), hyperParamNodes = c('psi'))
cm <- compileNimble(m)
capprox <- compileNimble(approx, project = m)
result <- runNestedApprox(capprox, originalScale = FALSE)  # Use originalScale=FALSE for comparison with INLA. Way underdispersed.
result$improveMarginals(1:3, nMarginalGrid = 7) ## Pretty good.

## AGHQ
approx2 <- buildNestedApprox(m, latentNodes = c('beta_int','beta_treat','beta_period','beta_carry'), hyperParamNodes = c('psi'), control = list(hyperGridRule='AGHQ',nQuadOuter=5))
capprox2 <- compileNimble(approx2, project = m)
result2 <- runNestedApprox(capprox2, originalScale = FALSE)  # Way underdispersed.

latent_sample <- result$sampleLatentNodes(100000)
apply(latent_sample, 2, quantile, qpts)  # beta_int is off; others seem good


result <- runNestedApprox(capprox, originalScale =TRUE)  
result$improveMarginals('psi', nMarginalGrid = 7)

## MCMC

library(nimbleHMC)
mcmc <- buildHMC(m, monitors = c('beta_int','beta_treat','beta_period','beta_carry',
                                 'psi','alpha'))
cmcmc <- compileNimble(mcmc, project = m)

out <- runMCMC(cmcmc, niter = 11000, nburnin = 1000)

## For comparison with nimble's nested approx
param <- cbind(logit(out[,'psi[1]']), logit(out[,'psi[2]']/(1-out[,'psi[1]'])),
               logit(out[,'psi[3]']/(1-out[,'psi[1]']-out[,'psi[2]'])))
apply(param, 2, quantile, qpts)
               
##  For comparison with INLA
theta <- cbind(out[,'alpha[1]'],
               log(out[,'alpha[2]']-out[,'alpha[1]']),
               log(out[,'alpha[3]']-out[,'alpha[2]']))
apply(theta, 2, quantile, qpts)

## INLA

library(INLA)
library(brms)
data(inhaler)

inla_pom <- inla(rating ~ treat + period + carry, data = inhaler, family='pom',
                 control.family=list(hyper=list(theta1=list(prior="dirichlet", param=3))))

summary(inla_pom)

inla.priors.used(inla_pom)

## inla_pom$summary.random  ## none present

