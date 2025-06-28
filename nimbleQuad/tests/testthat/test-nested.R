source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))
source(system.file(file.path('tests', 'testthat', 'AD_test_utils.R'), package = 'nimble'))
EDopt <- nimbleOptions("enableDerivs")
BMDopt <- nimbleOptions("buildModelDerivs")
nimbleOptions(enableDerivs = TRUE)
nimbleOptions(buildModelDerivs = TRUE)


test_that("Error trapping for invalid models", {

    ## Either no parameter or no latent nodes
    code <- nimbleCode({
        for(i in 1:n)
            y[i] ~ dnorm(b0+b1*x[i], 1)
        b0 ~ dnorm(0,1)
        b1 ~ dnorm(0,1)
    })
    n <- 5
    m <- nimbleModel(code, data = list(y = rnorm(n)), constants = list(n=n, x=rnorm(n)))
    
    expect_error(buildNestedApprox(m), "No latent nodes detected in model")
    expect_error(buildNestedApprox(m, paramNodes = c('b0','b1')), "No latent nodes detected in model")
    expect_error(buildNestedApprox(m, latentNodes = c('b0','b1')), "No parameter nodes detected in model")

    ## Missing priors

    code <- nimbleCode({
        for(i in 1:n)
            y[i] ~ dnorm(b0+b1*x[i], sd = sigma)
        b0 ~ dnorm(0,sd = 10)
        ## b1 ~ dnorm(0,sd = 10)  # Missing
        sigma ~ dhalfflat()
    })
    n <- 50
    m <- nimbleModel(code, data = list(y = rnorm(n)), constants = list(n=n, x=rnorm(n)))
    
    expect_error(buildNestedApprox(m, latentNodes = c('b0','b1'), paramNodes = 'sigma'),
                 "do not have prior distributions")

})


test_that("Determination of params and latents", {
    ## Fully centered
    code <- nimbleCode({
        for(i in 1:5) {
            y[i] ~ dnorm(b1 * x[i] + mu[i], tau_y)
            mu[i] ~ dnorm(mu0, tau)
        }
        mu0 ~ dflat()
        b1 ~ dflat()
        tau_y <- 1/(sigma*sigma)
        sigma ~ dhalfflat()
        tau ~ dhalfflat()
    })
    m <- nimbleModel(code, data = list(y = rnorm(5)))
    approx <- buildNestedApprox(m)
    expect_identical(c('mu0', 'sigma','tau'), approx$innerMethods$getNodeNamesVec())
    expect_identical(c('b1', paste0('mu[', 1:5, "]")), approx$innerMethods$getNodeNamesVec(returnParams = FALSE))

    ## Partially centered
    code <- nimbleCode({
        for(i in 1:5) {
            y[i] ~ dnorm(b0 + b1 * x[i] + mu[i], tau_y)
            mu[i] ~ dnorm(0, tau)
        }
        b0 ~ dflat()
        b1 ~ dflat()
        tau_y <- 1/(sigma*sigma)
        sigma ~ dhalfflat()
        tau ~ dhalfflat()
    })
    m <- nimbleModel(code, data = list(y = rnorm(5)))
    approx <- buildNestedApprox(m)
    expect_identical(c('sigma','tau'), approx$innerMethods$getNodeNamesVec())
    expect_identical(c('b0', 'b1', paste0('mu[', 1:5, "]")), approx$innerMethods$getNodeNamesVec(returnParams = FALSE))

    ## Non-normal latents
    code <- nimbleCode({
        for(i in 1:5) {
            y[i] ~ dbern(p[i])
            p[i] ~ dbeta(alpha, beta)
        }
        alpha ~ dhalfflat()
        beta ~ dhalfflat()
    })
    m <- nimbleModel(code, data = list(y = rnorm(5)))
    approx <- buildNestedApprox(m)
    expect_identical(c('alpha','beta'), approx$innerMethods$getNodeNamesVec())
    expect_identical(paste0('p[', 1:5, "]"), approx$innerMethods$getNodeNamesVec(returnParams = FALSE))
    
    ## No REs.
    code <- nimbleCode({
        for(i in 1:5) {
            y[i] ~ dnorm(b0 + b1 * x[i], tau_y)
        }
        b0 ~ dflat()
        b1 ~ dflat()
        tau_y <- 1/(sigma*sigma)
        sigma ~ dhalfflat()
    })
    m <- nimbleModel(code, data = list(y = rnorm(5)))
    ## If `latentNodes` not specified, error with no latents.
    approx <- buildNestedApprox(m, latentNodes = c('b0','b1'))
    expect_identical(c('sigma'), approx$innerMethods$getNodeNamesVec())
    expect_identical(c('b0', 'b1'), approx$innerMethods$getNodeNamesVec(returnParams = FALSE))

    ## Noncentered
    code <- nimbleCode({
        for(i in 1:5) {
            y[i] ~ dnorm(b0 + b1 * x[i] + tau*mu[i], tau_y)
            mu[i] ~ dnorm(0, 1)
        }
        b0 ~ dflat()
        b1 ~ dflat()
        tau_y <- 1/(sigma*sigma)
        sigma ~ dhalfflat()
        tau ~ dhalfflat()
    })
    m <- nimbleModel(code, data = list(y = rnorm(5)))
    ## If `latentNodes` not specified, error with no latents.
    approx <- buildNestedApprox(m, latentNodes = c('mu'))
    expect_identical(c('b0','b1','sigma','tau'), approx$innerMethods$getNodeNamesVec())
    expect_identical(paste0('mu[', 1:5, "]"), approx$innerMethods$getNodeNamesVec(returnParams = FALSE))

    ## Missing covariates
    code <- nimbleCode({
        for(i in 1:5) {
            y[i] ~ dnorm(b1 * x[i] + mu[i], tau_y)
            mu[i] ~ dnorm(mu0, tau)
            # x[i] ~ dnorm(mu_x, tau_x)
        }
        for(i in 1:nmissing)
            x[miss[i]] ~ dnorm(mu_x, tau_x)
        mu0 ~ dflat()
        b1 ~ dflat()
        tau_y <- 1/(sigma*sigma)
        sigma ~ dhalfflat()
        tau ~ dhalfflat()
        tau_x ~ dhalfflat()
        mu_x ~ dhalfflat()
    })
    x <- rnorm(5)
    miss <- c(2,4)
    x[miss] <- NA
    m <- nimbleModel(code, data = list(y = rnorm(5), x =x), constants = list(miss = miss, nmissing = 2), inits = list(x = rnorm(5)))
    approx <- buildNestedApprox(m)
    expect_identical(c('mu0', 'sigma','tau', 'tau_x', 'mu_x'), approx$innerMethods$getNodeNamesVec())
    expect_identical(c('b1', paste0('mu[', 1:5, "]"), paste0('x[', miss, ']')), approx$innerMethods$getNodeNamesVec(returnParams = FALSE))

    ## Why is mu2[1] last in the returned node names?
    code <- nimbleCode({
        for(i in 1:5) {
            y[i] ~ dnorm(b1 * x[i] + mu[i], sd = exp(b2*x[i] + mu2[i]))
            mu[i] ~ dnorm(mu0, tau)
            mu2[i] ~ dnorm(mu0, tau)
        }
        mu0 ~ dflat()
        b1 ~ dflat()
        b2 ~ dflat()
        tau_y <- 1/(sigma*sigma)
        sigma ~ dhalfflat()
        tau ~ dhalfflat()
    })
    m <- nimbleModel(code, data = list(y = rnorm(5)))
    approx <- buildNestedApprox(m)
    expect_identical(c('mu0','tau'), approx$innerMethods$getNodeNamesVec())
    expect_identical(c('b1', 'b2', paste0('mu[', 1:5, "]"), paste0('mu2[', 1:5, "]")), approx$innerMethods$getNodeNamesVec(returnParams = FALSE))
  
})





nimbleOptions(enableDerivs = EDopt)
nimbleOptions(buildModelDerivs = BMDopt)
