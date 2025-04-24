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
    expect_error(buildNestedApprox(m, hyperParamNodes = c('b0','b1')), "No latent nodes detected in model")
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
    
    expect_error(buildNestedApprox(m, latentNodes = c('b0','b1'), hyperParamNodes = 'sigma'),
                 "do not have prior distributions")

})








nimbleOptions(enableDerivs = EDopt)
nimbleOptions(buildModelDerivs = BMDopt)
