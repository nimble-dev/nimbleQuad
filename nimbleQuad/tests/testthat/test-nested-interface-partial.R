source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))
source(system.file(file.path('tests', 'testthat', 'AD_test_utils.R'), package = 'nimble'))
EDopt <- nimbleOptions("enableDerivs")
BMDopt <- nimbleOptions("buildModelDerivs")
nimbleOptions(enableDerivs = TRUE)
nimbleOptions(buildModelDerivs = TRUE)


temporarilyAssignInGlobalEnv <- function(value, replace = FALSE) {
    name <- deparse(substitute(value))
    assign(name, value, envir = .GlobalEnv)
    if(!replace) {
        rmCommand <- substitute(remove(name, envir = .GlobalEnv))
        do.call('on.exit', list(rmCommand, add = TRUE), envir = parent.frame())
    }
}

qpts <- c(.025,.25,.5,.75,.975)


test_that("Simple 1d param case - basic tests against known numerical results, including marginal distribution functions", {
    dig00 <- nimbleFunction(
        run = function(x = double(0), log = logical(0, default = FALSE)) {
            returnType(double(0))
            if(log) return(-log(x)) else return(1/x)
        }, buildDerivs = TRUE)
    
    registerDistributions(list(dig00=list(BUGSdist="dig00()",range=c(0,Inf))))
    temporarilyAssignInGlobalEnv(dig00)
    temporarilyAssignInGlobalEnv(rig00)
    
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
    
    ## marginal for sigma2 is IG((n-1)/2,(n-1)*s2/2)
    ## marginal for mu is t(ybar, s2/n)
    
    qs <- qinvgamma(qpts, (n-1)/2, scale=(n-1)*var(m$y)/2)

    approx <- buildNestedApprox(m, latentNodes = 'mu', paramNodes = 'sigma2')
    cm <- compileNimble(m)
    capprox <- compileNimble(approx, project = m)
    result <- runNestedApprox(capprox)
    exp_table <- result$expectations  # from summary table
    qs_table <- result$quantiles      # from summary table
    qs_est <- result$qmarginal('sigma2')
    expect_lt(max(abs(qs - qs_est)), 0.01)


        result$improveParamMarginals('sigma2', nMarginalGrid = 5)
        ## Expect identical results since 1d improvement done by default (and with 5 points).
        expect_identical(qs_est, result$qmarginal('sigma2')) 
        expect_identical(exp_table, result$expectations)
        expect_identical(qs_table, result$quantiles)
        
        result$improveParamMarginals('sigma2', nMarginalGrid = 11)

        qs_est_impr <- result$qmarginal('sigma2')
        expect_lt(max(abs(qs - qs_est_impr)), .0004)  
        expect_false(identical(qs_est, qs_est_impr))
        
        result$improveParamMarginals('sigma2', nMarginalGrid = 21)

        qs_est_impr2 <- result$qmarginal('sigma2')
        expect_lt(max(abs(qs - qs_est_impr2)), 1e-4) 
        expect_false(identical(qs_est_impr, qs_est_impr2))

    
    new_qpts <- c(0.3, 0.72)
    qs_est <- result$qmarginal('sigma2', new_qpts)
    qs <- qinvgamma(new_qpts, (n-1)/2, scale=(n-1)*var(m$y)/2)
    expect_lt(max(abs(qs - qs_est)), 1e-4)

    new_qpts <- 0.99
    qs_est <- result$qmarginal('sigma2', new_qpts)
    qs <- qinvgamma(new_qpts, (n-1)/2, scale=(n-1)*var(m$y)/2)
    expect_lt(abs(qs - qs_est), 3e-4)

    rtrue <- rinvgamma(1000, (n-1)/2, scale=(n-1)*var(m$y)/2)
    rapprox <- result$rmarginal('sigma2', 1000)
    expect_gt(ks.test(rtrue,rapprox)$p.value, 0.05)

    grid <- seq(.35, 2.25, len = 50)
    dtrue <- dinvgamma(grid, (n-1)/2, scale=(n-1)*var(m$y)/2, log = TRUE)
    dapprox <- result$dmarginal('sigma2', grid, log=TRUE)
    expect_lt(max(abs(dtrue - dapprox)), .0003)

    prec_approx <- result$emarginal('sigma2', function(x) 1/x)
    sd_approx <- result$emarginal('sigma2', function(x) sqrt(x))
    prob_approx <- result$emarginal('sigma2', function(x, val) x < val, 2)

    rtrue <- rinvgamma(1e5, (n-1)/2, scale=(n-1)*var(m$y)/2)
    prec_true <- mean(1/rtrue)
    sd_true <- mean(sqrt(rtrue))
    prob_true <- pinvgamma(2, (n-1)/2, scale=(n-1)*var(m$y)/2)

    expect_lt(abs(prec_true - prec_approx), 1e-3)
    expect_lt(abs(sd_true - sd_approx), 1e-3)
    expect_lt(abs(prob_true - prob_approx), 1e-4)


    ## More constrained priors for MLL calc to be valid.
    code <- nimbleCode({
        for(i in 1:n)
            y[i] ~ dnorm(mu, sd = sigma)
        mu ~ dnorm(0,sd=3)
        sigma ~ dunif(0,5)
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

    ## Basic harmonic mean estimator (high variance, so need many samples).
    if(FALSE) {
        set.seed(1)
        n <- 5e7
        mus <- rnorm(n,0,sd=3)
        sigmas <- runif(n,0,5)
        theta <- cbind(mus,sigmas)
        logpy <- apply(theta, 1, function(x) sum(dnorm(y,x[1],x[2],log=T)))
        mll <- log(mean(exp(logpy)))
    } else mll <- -45.351158

    ## WORK ON THIS ##
    
    expect_lt(abs(mll - result$marginalLogLik), .01)
    expect_lt(abs(mll - result$marginalLogLikImproved), .004)

    ## This is probably getting to the resolution of the accuracy of the harmonic mean estimator...
    result$setParamGrid(nQuad = 15)
    result$calcMarginalLogLikImproved()
    result$marginalLogLikImproved
    expect_lt(abs(mll - result$marginalLogLikImproved), .002)
})


nimbleOptions(enableDerivs = EDopt)
nimbleOptions(buildModelDerivs = BMDopt)

    
    
    

