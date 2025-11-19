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
    expect_identical(c('sigma'), approx$innerMethods$getNodeNameSingle())
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

    ## User putting RE params in params.
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
    approx <- buildNestedApprox(m, paramNodes = c('b0','b1','tau','sigma'))
    expect_identical(c('b0','b1','tau', 'sigma'), approx$innerMethods$getNodeNamesVec())
    expect_identical(paste0('mu[', 1:5, "]"), approx$innerMethods$getNodeNamesVec(returnParams = FALSE))

    ## Why is mu2[1] last in the returned node names?
    ## Bug in getConditionallyIndepSets that causes a separate Laplace for mu2[1].
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

    code <- nimbleCode({
        for(i in 1:5) {
            y[i] ~ dnorm(b1 * x[i] + mu[i], tau_y)
            mu[i] ~ dnorm(mu0, tau)
            x[i] ~ dnorm(mu_x, tau_x)
        }
        mu0 ~ dflat()
        b1 ~ dflat()
        tau_y <- 1/(sigma*sigma)
        sigma ~ dhalfflat()
        tau ~ dhalfflat()
        tau_x ~ dhalfflat()
        mu_x ~ dhalfflat()
    })
    x <- rnorm(5)
    x[2] <- NA
    m <- nimbleModel(code, data = list(y = rnorm(5), x =x), inits = list(x = rnorm(5)))
    expect_message(approx <- buildNestedApprox(m), "latent nodes: x[2], b1, mu (5 elements)", fixed = TRUE)

    x[4] <- NA
    m <- nimbleModel(code, data = list(y = rnorm(5), x =x), inits = list(x = rnorm(5)))
    expect_message(approx <- buildNestedApprox(m), "latent nodes: b1, mu (5 elements), x (2 elements)", fixed = TRUE)

})

test_that("Controlling quadrature grids", {
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
    
    expect_message(approx <- buildNestedApprox(m, latentNodes = c('eta'), paramNodes = c('mu','sigma','phi')),
                   "with CCD grid")
    expect_identical(approx$innerMethods$nQuad, 1)
    expect_identical(approx$nQuadParam, 3)
    expect_message(approx <- buildNestedApprox(m, latentNodes = c('eta'), paramNodes = c('mu','sigma','phi'),
                                               control = list(paramGridRule = "AGHQ", nQuadParam = 5)),
                   "with AGHQ grid")
    expect_identical(approx$nQuadParam, 5)
    expect_message(approx <- buildNestedApprox(m, latentNodes = c('eta'), paramNodes = c('mu','sigma','phi'),
                                               control = list(nQuadLatent = 3)),
                   "AGHQ approximation for the latent")
    expect_identical(approx$innerMethods$nQuad, 3)
})

    

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

    if(Sys.info()['sysname'] != "Windows") {  # Issue 71
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
    }
    
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

test_that("Basic interface and user input errors", {
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
    ## This takes too long to include in testing.
    # expect_output(runNestedApprox(approx), "Running an uncompiled nested approximation")

    cm <- compileNimble(m)
    capprox <- compileNimble(approx, project = m)

    ## Check use of `originalScale = FALSE`.
    result <- runNestedApprox(capprox)
    resultTrans <- runNestedApprox(capprox, originalScale = FALSE)

    expect_identical(result$expectations$mu, resultTrans$expectations$param_trans1)
    expect_identical(result$quantiles$mu, resultTrans$quantiles$param_trans1)

    expect_identical(result$quantiles$sigma, exp(resultTrans$quantiles$param_trans2))
    
    set.seed(1)
    paramSamples <- result$sampleParams(50)
    set.seed(1)
    paramSamplesTrans <- resultTrans$sampleParams(50)
    tmp <- paramSamplesTrans
    tmp[,'param_trans2'] <- exp(tmp[,'param_trans2'])
    tmp[,'param_trans3'] <- exp(tmp[,'param_trans3'])
    colnames(tmp) <- c('mu','sigma','phi')
    expect_identical(tmp, paramSamples)

    ## Check sampling options.
    paramSamples <- result$sampleParams(10000)
    paramSamplesUnmatch <- result$sampleParams(10000, matchMarginals = FALSE)
    matched_qs <- apply(paramSamples, 2, quantile, c(.025,.25,.5,.75,.975))
    unmatched_qs <- apply(paramSamplesUnmatch, 2, quantile, c(.025,.25,.5,.75,.975))
    expect_lt(max(abs(unlist(result$quantiles) - matched_qs)), .001)
    expect_lt(max(abs(unlist(result$quantiles) - unmatched_qs)), .02)

    latentSamples <- result$sampleLatents(500, includeParams = TRUE)
    expect_identical(dim(latentSamples), c(500L, 11L))
    grid <- result$approx$getParamGrid()
    expect_identical(sort(grid[,1]), sort(unique(latentSamples[,'mu'])))
    expect_identical(exp(sort(grid[,2])), sort(unique(latentSamples[,'sigma'])))

    q1 <- result$quantiles
    if(Sys.info()['sysname'] != "Windows") {  # Issue 71
        result$improveParamMarginals()
        q2 <- result$quantiles
        result$improveParamMarginals(nMarginalGrid = 11)
        q3 <- result$quantiles
        result$improveParamMarginals(nMarginalGrid = 11, nQuad = 5)
        q4 <- result$quantiles
        expect_false(identical(q1,q2))
        expect_false(identical(q2,q3))
        expect_false(identical(q3,q4))
        expect_error(result$improveParamMarginals(quadRule = 'CCD'), "Only AGHQ-based quadrature rules")
    }
    
    ## Try various interface interactions.
    expect_silent(result <- runNestedApprox(capprox))
    ll1 <- result$marginalLogLik
    
    result <- runNestedApprox(capprox, nSamplesLatents = 25, nSamplesParams = 50)
    expect_identical(dim(result$samples), c(25L,8L))
    expect_identical(dim(result$paramSamples), c(50L,3L))
    expect_identical(colnames(result$samples), m$expandNodeNames('eta'))
    expect_identical(colnames(result$paramSamples), c('mu','sigma','phi'))
    ll2 <- result$marginalLogLikImproved
    expect_false(identical(ll1, ll2))
   
    result$calcMarginalLogLikImproved()
    expect_identical(ll2, result$marginalLogLikImproved)
    set.seed(1)
    smp1 <- result$sampleLatents(20)
    set.seed(1)
    smp2 <- result$sampleLatents(20)
    expect_identical(smp1, smp2)

    result$setParamGrid(quadRule = 'AGHQ')
    result$calcMarginalLogLikImproved()
    ll3 <- result$marginalLogLikImproved
    expect_false(identical(ll2, ll3))

    set.seed(1)
    smp3 <- result$sampleLatents(20)
    expect_false(identical(smp2, smp3))

    ## This should not affect the grid as default should be 3.
    result$setParamGrid(quadRule = 'AGHQ', nQuad = 3)
    result$calcMarginalLogLikImproved()
    expect_identical(ll3, result$marginalLogLikImproved)

    result$setParamGrid(quadRule = 'AGHQ', nQuad = 5)
    result$calcMarginalLogLikImproved()
    ll4 <- result$marginalLogLikImproved
    expect_false(identical(ll2, ll3))

    set.seed(1)
    smp4 <- result$sampleLatents(20)
    expect_false(identical(smp3, smp4))

})

nimbleOptions(enableDerivs = EDopt)
nimbleOptions(buildModelDerivs = BMDopt)

    
    
    

