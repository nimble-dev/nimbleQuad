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

    expect_lt(abs(mll - result$marginalLogLik), .01)
    expect_lt(abs(mll - result$marginalLogLikImproved), .004)

    ## This is probably getting near the resolution of the accuracy of the harmonic mean estimator...
    result$setParamGrid(nQuad = 15)
    result$calcMarginalLogLikImproved()
    result$marginalLogLikImproved
    expect_lt(abs(mll - result$marginalLogLikImproved), .002)
})


test_that("Marginal log-likelihood, 2-d case", {
    code <- nimbleCode({
        for(j in 1:J) {
            for(i in 1:n)
                y[i,j] ~ dpois(mu[j])
            mu[j] ~ dgamma(a, b)
        }
        a ~ dgamma(1, 1)
        b ~ dgamma(1, 1)
    })

    set.seed(1)
    n <- 10
    J <- 8
    a <- 2
    b <- 1
    mu <- rgamma(J, a, b)
    mns <- rep(mu, each = n)
    y <- matrix(rpois(n*J, mns), ncol = J)
    m <- nimbleModel(code, data = list(y = y), constants = list(n=n, J=J),
                     inits = list(mu = rep(0,J), a=1,b=1), buildDerivs = TRUE)


    if(FALSE) {
        M <- 1e6
        py <- rep(0, M)
        
        dens <- function(idx) {
            ytmp <- y[(1+(idx-1)*n):(idx*n)]
            return(a*log(b) - lgamma(a) - sum(lgamma(ytmp+1)) + lgamma(a+sum(ytmp)) - (a+sum(ytmp)) * log(b + n))
        }
        
        set.seed(1)
        for(i in seq_along(py)) {
            a <- rgamma(1, 1,1)
            b <- rgamma(1,1,1)
            logpy <- sum(sapply(1:J, dens))
            py[i] <- exp(logpy)
        }
        
        mll <- log(mean(py))
    } else mll <- -139.80807

    approx <- buildNestedApprox(m)
    cm <- compileNimble(m)
    capprox <- compileNimble(approx, project = m)
    result <- runNestedApprox(capprox)
    
    expect_lt(abs(mll - result$marginalLogLik), 0.1)

    tmp <- result$sampleLatents(n = 10)  # Side effect of calculating logLik.
    expect_lt(abs(mll - result$marginalLogLikImproved), 0.07)

    result$setParamGrid(nQuad = 7)
    result <- runNestedApprox(capprox)
    result$calcMarginalLogLikImproved()
    expect_lt(abs(mll - result$marginalLogLikImproved), 0.05)

    approx <- buildNestedApprox(m, control = list(nQuadParam = 7, nQuadLatent = 5))
    capprox <- compileNimble(approx, project = m)
    result <- runNestedApprox(capprox) 
    result$calcMarginalLogLikImproved()
    expect_lt(abs(mll - result$marginalLogLikImproved), 0.02)

})

test_that("3-d case", {

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
                     inits = list(eta = rep(0,J), mu = 0, sigma = 1, phi = 1), buildDerivs = TRUE)

    cm <- compileNimble(m)

    if(FALSE) {
        conf <- configureMCMC(m, monitors = c('mu','sigma','phi','eta'), onlySlice = TRUE)
        mcmc <- buildMCMC(conf)
        cmcmc <- compileNimble(mcmc, project=m)
        out <- runMCMC(cmcmc, niter=51000, nburnin=1000)
        qs_mcmc <- apply(out, 2, quantile, qpts)
        save(qs_mcmc, file = 'mcmc-results1.Rda')
    } else load(system.file(file.path('tests', 'testthat', 'mcmc-results1.Rda'), package = 'nimbleQuad'))
    
    approx <- buildNestedApprox(m, latentNodes = c('eta'), paramNodes = c('mu','sigma','phi'))
    capprox <- compileNimble(approx, project = m)
    result <- runNestedApprox(capprox)  # Fairly different from MCMC for mu and sigma.

    if(Sys.info()['sysname'] != "Windows") {  # Issue 71
        result$improveParamMarginals(c('mu','sigma','phi'), nMarginalGrid = 7)
        expect_lt(max(abs(qs_mcmc[,c('mu','sigma','phi')] - unlist(result$quantiles))), .06)
    }
    
    latent_sample <- result$sampleLatents(10000)
    qs_nest <- apply(latent_sample, 2, quantile, qpts)
    ## Not as close as one might hope, but relative to true eta values, results are consistent between
    ## approximation and MCMC.
    expect_lt(max(abs(qs_mcmc[ , grep("eta", colnames(qs_mcmc))] - qs_nest)), .18) # .162
})

test_that("3-d case, no RE variation", {
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
    eta <- rep(0, J)
    phi <- 0.5
    mns <- rep(exp(eta), each = n)
    sds <- rep(sqrt(exp(eta)^2/phi), each = n)
    y <- matrix(rgamma(n*J, shape = mns^2/sds^2, rate = mns/sds^2), ncol = J)
    
    m <- nimbleModel(code, data = list(y = y), constants = list(n=n, J=J),
                     inits = list(eta = rep(0,J), mu = 0, sigma=1, phi = 1), buildDerivs = TRUE)
    cm <- compileNimble(m)

    if(FALSE) {
        conf <- configureMCMC(m, monitors = c('mu','sigma','phi','eta'), onlySlice = TRUE)
        mcmc <- buildMCMC(conf)
        cmcmc <- compileNimble(mcmc, project=m)
        system.time(out <- runMCMC(cmcmc, niter=501000, nburnin=1000, thin = 10) )
        qs_mcmc <- apply(out, 2, quantile, qpts)
        save(qs_mcmc, file = 'mcmc-results2.Rda')
    } else load(system.file(file.path('tests', 'testthat', 'mcmc-results2.Rda'), package = 'nimbleQuad'))

    approx <- buildNestedApprox(m, latentNodes = c('eta'), paramNodes = c('mu','sigma','phi'))
    capprox <- compileNimble(approx, project = m)
    result <- runNestedApprox(capprox)  # Fairly different from MCMC for mu and sigma.

    if(Sys.info()['sysname'] != "Windows") {  # Issue 71
        result$improveParamMarginals(c('mu','sigma','phi'), nMarginalGrid = 7)
        expect_lt(max(abs(qs_mcmc[,c('mu','sigma','phi')] - unlist(result$quantiles))), .01)
    }
    
    latent_sample <- result$sampleLatents(10000)
    qs_nest <- apply(latent_sample, 2, quantile, qpts)

    expect_lt(max(abs(qs_mcmc [ , grep("eta", colnames(qs_mcmc))]- qs_nest)), .4)  # Tails can be rather far off. Using 3 inner grid points has little effect. 

    ## Tried to use INLA priors, but INLA results quite far off from MCMC and from our nested approx,
    ## and our nested approx quite far off from MCMC too (ESS > 300 for params/latents) - see nested-3dgamma.R.

})

test_that("Poisson 2-d case", {
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

    ## INLA
    if(FALSE) { 
        library(INLA)
        group <- as.factor(rep(1:J, each = n))
        yc <- c(y)
        formula <- y ~ 1 + f(group, model = "iid")
        fit <- inla(formula, family="poisson", data=data.frame(y=yc,group=group), quantiles = qpts,
                    control.compute=list(config = TRUE),
                    control.fixed = list(prec.intercept = .001))
        qs_inla_param <- rbind(fit$summary.fixed[1:8], fit$summary.hyperpar)[ , c('0.025quant', '0.25quant', '0.5quant', '0.75quant', '0.975quant')]
        qs_inla <- fit$summary.random$group[ , c('0.025quant', '0.25quant', '0.5quant', '0.75quant', '0.975quant')]
        sampled <- inla.posterior.sample(n = 10000, fit)
        smp <- t(sapply(sampled, function(x) x$latent[81:89,1]))
        qs_inla_smp <- apply(smp, 2, quantile, qpts) 
    }
    
    ## MCMC
    if(FALSE) {
        m <- nimbleModel(code, data = list(y = y), constants = list(n=n, J=J),
                         inits = list(lambda = rep(0,J), mu = 0, tau=1), buildDerivs = TRUE)
        cm <- compileNimble(m)
        mcmc <- buildHMC(m, monitors = c('mu','tau','lambda'))
        cmcmc <- compileNimble(mcmc, project=m)
        out <- runMCMC(cmcmc, niter=51000,nburnin=1000)
        qs_mcmc <- apply(out, 2, quantile, qpts)
        save(qs_mcmc, file = 'mcmc-results3.Rda')
    } else load(system.file(file.path('tests', 'testthat', 'mcmc-results3.Rda'), package = 'nimbleQuad'))

    ## Treating mu as latent. Results are closer to INLA and HMC than if mu is parameter.
    m <- nimbleModel(code, data = list(y = y), constants = list(n=n, J=J),
                 inits = list(lambda = rep(0,J), mu = 0, tau=1), buildDerivs = TRUE)
    approx <- buildNestedApprox(m, latentNodes = c('lambda','mu'), paramNodes = c('tau'))
    cm <- compileNimble(m)
    capprox <- compileNimble(approx, project = m)
    result <- runNestedApprox(capprox)  
    expect_lt(max(abs(qs_mcmc[,c("tau")] - unlist(result$quantiles))), .22)  # .198, INLA is .11
    result$improveParamMarginals(c('tau'), nMarginalGrid = 7)
    expect_lt(max(abs(qs_mcmc[,c("tau")] - unlist(result$quantiles))), .13) # .11; strangely worse in right tail after 'improvement'.

    latent_sample <- result$sampleLatents(10000)
    expect_lt(max(abs(qs_mcmc[,grep("mu|lambda", colnames(qs_mcmc))] -
                      apply(latent_sample, 2, quantile, qpts)[,c(2:9,1)])), .08) # .064; .036 for INLA marginals and .028 for INLA samples; 
})

test_that("nested REs in Bernoulli GLMM", {
    set.seed(1)
    n <- 10000
    ntowns <- 262
    nstates <- 47

    gender <- sample(c("M","F"), n, replace = TRUE)
    livingarrangement <- sample(c("H","D","I"), n, prob = c(.5,.2,.3), replace = TRUE)
    race <- sample(c("AI","AS","B","H", "W"), n, prob = c(.03,.05,.1,.32,.5), replace = TRUE)

    sex <- as.numeric(gender == "F")
    live1 <- as.numeric(livingarrangement == "D")
    live2 <- as.numeric(livingarrangement == "I")

    race1 <- as.numeric(race == "AI")
    race2 <- as.numeric(race == "AS")
    race3 <- as.numeric(race == "B")
    race4 <- as.numeric(race == "H")

    town <- rcat(n, prob = rep(1/ntowns, ntowns))
    town2state <- sample(seq_len(nstates), ntowns, replace = TRUE)
    state <- town2state[town]
    sigma_state <- 0.5
    sigma_town <- 0.25

    beta0 <- 0.1
    beta_sex <- 0.2
    beta_live <- c(-.1, .05)
    beta_race <- c(-.2, .1, -.05, 0)

    u1 <- rnorm(nstates, 0, sigma_state)
    u2 <- rnorm(ntowns, 0, sigma_town)

    eta <- beta0 + beta_sex*sex + beta_live[1]*live1 + beta_live[2]*live2 +
        beta_race[1]*race1 + beta_race[2]*race2 + beta_race[3]*race3 + beta_race[4]*race4 + u1[state] + u2[town]

    p <- expit(eta)

    y <- rbinom(n, 1, p)

    if(FALSE) {
        library(INLA)
        dat <- data.frame(y = y, gender = gender, race = race, livingarrangement=livingarrangement, state = as.character(state), town = as.character(town))

        ## Need to make fixed effect priors consistent for comparison of statistical results.
        fit <- inla(
            y ~ gender + race + livingarrangement + 
                f(state,model = "iid",hyper = list(prec = list(prior = "pc.prec",param = c(.5,.5)))) +
                f(town,model = "iid",hyper = list(prec = list(prior = "pc.prec",param = c(.5,.5)))),
            data = dat, quantiles = qpts,
            family = 'binomial', control.compute=list(config=TRUE))
        ## Order of random effects has gotten shuffled...
        ord1 <- order(as.numeric(fit$summary.random$state$ID))
        ord2 <- order(as.numeric(fit$summary.random$town$ID))

        qs_inla <- t(rbind(fit$summary.random$state[ord1, paste0(as.character(qpts), "quant")],
                         fit$summary.random$town[ord2, paste0(as.character(qpts), "quant")]))
            
        sampled <- inla.posterior.sample(n = 10000, fit)
        smp <- t(sapply(sampled, function(x) x$latent[10001:10307,1]))
        smp <- cbind(smp[ , 1:(nstates-2)][ , ord1],
                     smp[ , (nstates+1-2):(nstates+ntowns-2)][ , ord2])
        qs_inla_smp <- apply(smp, 2, quantile, qpts)
    }
    
    if(FALSE) {
        ## Use noncentered for better HMC mixing.
        code <- nimbleCode({
            for(j in 1:nstates)
                u1[j] ~ dnorm(0, sd = 1)
            for(j in 1:ntowns)
                u2[j] ~ dnorm(0, sd = 1)
            sigma_state ~ dhalfflat()
            sigma_town ~ dhalfflat()

            for(i in 1:n) {
                eta[i] <- beta0 + beta_sex*sex[i] + beta_live[1]*live1[i] + beta_live[2]*live2[i] +
                    beta_race[1]*race1[i] + beta_race[2]*race2[i] + beta_race[3]*race3[i] + beta_race[4]*race4[i] +
                    sigma_state*u1[state[i]] + sigma_town*u2[town[i]]
                y[i] ~ dbern(expit(eta[i]))
            }
            beta0 ~ dflat()
            beta_sex ~ dflat()
            for(j in 1:2)
                beta_live[j] ~ dflat()
            for(j in 1:4)
                beta_race[j] ~ dflat()
        })

        set.seed(1)
        m <- nimbleModel(code, data = list(y=y), constants = list(ntowns = ntowns, nstates = nstates, n=n, state = state, town = town, sex = sex, race1=race1,race2=race2,race3=race3,race4=race4,live1=live1,live2=live2), inits = list(beta0 = 0, beta_sex = 0, beta_live = rep(0,2), beta_race = rep(0,4), sigma_state = 1, sigma_town = 1, u1 = rnorm(nstates), u2 = rnorm(ntowns)), calculate = FALSE, buildDerivs = TRUE) 

        cm <- compileNimble(m)

        library(nimbleHMC)
        mcmc <- buildHMC(m)
        cmcmc <- compileNimble(mcmc, project = m)

        system.time(out <- runMCMC(cmcmc, niter = 21000, nburnin = 1000)) #  n=1e4: 1367 sec.
        qs_mcmc <- apply(out, 2, quantile, qpts)

        u1cols <- grep("u1", colnames(out))
        u2cols <- grep("u2", colnames(out))

        ## Rescale noncentered estimates.
        out[,u1cols] <- out[,'sigma_state']*out[,u1cols]
        out[,u2cols] <- out[,'sigma_town']*out[,u2cols]
        qs_mcmc <- apply(out,2,quantile,qpts)

        save(qs_mcmc, file = 'mcmc-results4.Rda')
    } else load(system.file(file.path('tests', 'testthat', 'mcmc-results4.Rda'), package = 'nimbleQuad'))



    code <- nimbleCode({
        for(j in 1:nstates)
            u1[j] ~ dnorm(0, sd = sigma_state)
        for(j in 1:ntowns)
            u2[j] ~ dnorm(0, sd = sigma_town)
        sigma_state ~ dhalfflat()
        sigma_town ~ dhalfflat()

        for(i in 1:n) {
            eta[i] <- beta0 + beta_sex*sex[i] + beta_live[1]*live1[i] + beta_live[2]*live2[i] +
                beta_race[1]*race1[i] + beta_race[2]*race2[i] + beta_race[3]*race3[i] + beta_race[4]*race4[i] +
                u1[state[i]] + u2[town[i]]
            y[i] ~ dbern(expit(eta[i]))
        }
        beta0 ~ dflat()
        beta_sex ~ dflat()
        for(j in 1:2)
            beta_live[j] ~ dflat()
        for(j in 1:4)
            beta_race[j] ~ dflat()
    })
    set.seed(1)
    m <- nimbleModel(code, data = list(y=y), constants = list(ntowns = ntowns, nstates = nstates, n=n, state = state, town = town, sex = sex, race1=race1,race2=race2,race3=race3,race4=race4,live1=live1,live2=live2), inits = list(beta0 = 0, beta_sex = 0, beta_live = rep(0,2), beta_race = rep(0,4), sigma_state = 1, sigma_town = 1, u1 = rnorm(nstates), u2 = rnorm(ntowns)), calculate = FALSE, buildDerivs = TRUE)  
    cm <- compileNimble(m)
    
    approx <- buildNestedApprox(m, latentNodes = c('u1','u2','beta0','beta_sex','beta_race','beta_live'),
                                paramNodes = c('sigma_state','sigma_town'))
    capprox <- compileNimble(approx, project = m)
    result <- runNestedApprox(capprox)
    expect_lt(max(abs(qs_mcmc[,c('sigma_state','sigma_town')] - unlist(result$quantiles))), .015)  # 0.01; Better than INLA
    
    result$improveParamMarginals(c("sigma_state","sigma_town"), nMarginalGrid = 5)
    expect_lt(max(abs(qs_mcmc[,c('sigma_state','sigma_town')] - unlist(result$quantiles))), .009)  # 0.007; Better than INLA

    latent_sample <- result$sampleLatents(10000)
    qs_nest <- apply(latent_sample,2, quantile, qpts)
    
    expect_lt(max(abs(qs_mcmc[,1:8] - qs_nest[,c(1,3,4,5:8,2)])), .025)  # .019; Fixed effects
                                        
    ## Random effects (omit u1[c(9,28)]).
    expect_lt(max(abs(qs_mcmc[,c(11:18,20:37,39:319)] - qs_nest[,9:315])), .05) # Max is .046, while for INLA, .034 and .036 for marginal and sampled.
    ## Unlike crossed case, not much evidence of systematic error.
    
})

test_that("crossed REs in Bernoulli GLMM", {
    set.seed(1)
    n <- 10000
    ntowns <- 262
    nstates <- 47

    gender <- sample(c("M","F"), n, replace = TRUE)
    livingarrangement <- sample(c("H","D","I"), n, prob = c(.5,.2,.3), replace = TRUE)
    race <- sample(c("AI","AS","B","H", "W"), n, prob = c(.03,.05,.1,.32,.5), replace = TRUE)

    sex <- as.numeric(gender == "F")
    live1 <- as.numeric(livingarrangement == "D")
    live2 <- as.numeric(livingarrangement == "I")

    race1 <- as.numeric(race == "AI")
    race2 <- as.numeric(race == "AS")
    race3 <- as.numeric(race == "B")
    race4 <- as.numeric(race == "H")

    town <- rcat(n, prob = rep(1/ntowns, ntowns))
    if(n == 1000)
        town <- (rep(1:ntowns, 5))[1:n]
    state <- rcat(n, prob = rep(1/nstates, nstates))

    sigma_state <- 0.5
    sigma_town <- 0.25

    beta0 <- 0.1
    beta_sex <- 0.2
    beta_live <- c(-.1, .05)
    beta_race <- c(-.2, .1, -.05, 0)

    u1 <- rnorm(nstates, 0, sigma_state)
    u2 <- rnorm(ntowns, 0, sigma_town)

    eta <- beta0 + beta_sex*sex + beta_live[1]*live1 + beta_live[2]*live2 +
        beta_race[1]*race1 + beta_race[2]*race2 + beta_race[3]*race3 + beta_race[4]*race4 + u1[state] + u2[town]

    p <- expit(eta)

    y <- rbinom(n, 1, p)

    if(FALSE) {
        library(INLA)
        dat <- data.frame(y = y, gender = gender, race = race, livingarrangement=livingarrangement, state = as.character(state), town = as.character(town))

        ## Need to make fixed effect priors consistent for full comparison of statistical results,
        ## though it doesn't seem to make a difference in this case.
        fit <- inla(
            y ~ gender + race + livingarrangement + 
                f(state,model = "iid",hyper = list(prec = list(prior = "pc.prec",param = c(.5,.5)))) +
                f(town,model = "iid",hyper = list(prec = list(prior = "pc.prec",param = c(.5,.5)))),
            data = dat, quantiles = qpts,
            family = 'binomial', control.compute=list(config=TRUE))
        ## Order of random effects has gotten shuffled...
        ord1 <- order(as.numeric(fit$summary.random$state$ID))
        ord2 <- order(as.numeric(fit$summary.random$town$ID))

        qs_inla <- t(rbind(fit$summary.random$state[ord1, paste0(as.character(qpts), "quant")],
                         fit$summary.random$town[ord2, paste0(as.character(qpts), "quant")]))
            
        sampled <- inla.posterior.sample(n = 10000, fit)
        ## sampled <- inla.posterior.sample(n = 10000, fit, use.improved.mean = FALSE, skew.corr = FALSE)
        smp <- t(sapply(sampled, function(x) x$latent[10001:10309,1]))
        smp <- cbind(smp[ , 1:nstates][ , ord1],
                     smp[ , (nstates+1):(nstates+ntowns)][ , ord2])
        qs_inla_smp <- apply(smp, 2, quantile, qpts)
        
    }
    if(FALSE) {
        ## Use noncentered for better HMC mixing.
        code <- nimbleCode({
            for(j in 1:nstates)
                u1[j] ~ dnorm(0, sd = 1)
            for(j in 1:ntowns)
                u2[j] ~ dnorm(0, sd = 1)
            sigma_state ~ dhalfflat()
            sigma_town ~ dhalfflat()

            for(i in 1:n) {
                eta[i] <- beta0 + beta_sex*sex[i] + beta_live[1]*live1[i] + beta_live[2]*live2[i] +
                    beta_race[1]*race1[i] + beta_race[2]*race2[i] + beta_race[3]*race3[i] + beta_race[4]*race4[i] +
                    sigma_state*u1[state[i]] + sigma_town*u2[town[i]]
                y[i] ~ dbern(expit(eta[i]))
            }
            beta0 ~ dflat()
            beta_sex ~ dflat()
            for(j in 1:2)
                beta_live[j] ~ dflat()
            for(j in 1:4)
                beta_race[j] ~ dflat()
        })

        set.seed(1)
        m <- nimbleModel(code, data = list(y=y), constants = list(ntowns = ntowns, nstates = nstates, n=n, state = state, town = town, sex = sex, race1=race1,race2=race2,race3=race3,race4=race4,live1=live1,live2=live2), inits = list(beta0 = 0, beta_sex = 0, beta_live = rep(0,2), beta_race = rep(0,4), sigma_state = 1, sigma_town = 1, u1 = rnorm(nstates), u2 = rnorm(ntowns)), calculate = FALSE, buildDerivs = TRUE)  

        cm <- compileNimble(m)

        library(nimbleHMC)
        mcmc <- buildHMC(m)
        cmcmc <- compileNimble(mcmc, project = m)

        system.time(out <- runMCMC(cmcmc, niter = 21000, nburnin = 1000)) #  n=1e4: 1367 sec.

        u1cols <- grep("u1", colnames(out))
        u2cols <- grep("u2", colnames(out))

        ## Rescale noncentered estimates.
        out[,u1cols] <- out[,'sigma_state']*out[,u1cols]
        out[,u2cols] <- out[,'sigma_town']*out[,u2cols]
        
        qs_mcmc <- apply(out, 2, quantile, qpts)
        save(qs_mcmc, file = 'mcmc-results5.Rda')
    } else load(system.file(file.path('tests', 'testthat', 'mcmc-results5.Rda'), package = 'nimbleQuad'))
   

    code <- nimbleCode({
        for(j in 1:nstates)
            u1[j] ~ dnorm(0, sd = sigma_state)
        for(j in 1:ntowns)
            u2[j] ~ dnorm(0, sd = sigma_town)
        sigma_state ~ dhalfflat()
        sigma_town ~ dhalfflat()

        for(i in 1:n) {
            eta[i] <- beta0 + beta_sex*sex[i] + beta_live[1]*live1[i] + beta_live[2]*live2[i] +
                beta_race[1]*race1[i] + beta_race[2]*race2[i] + beta_race[3]*race3[i] + beta_race[4]*race4[i] +
                u1[state[i]] + u2[town[i]]
            y[i] ~ dbern(expit(eta[i]))
        }
        beta0 ~ dflat()
        beta_sex ~ dflat()
        for(j in 1:2)
            beta_live[j] ~ dflat()
        for(j in 1:4)
            beta_race[j] ~ dflat()
    })
    set.seed(1)
    m <- nimbleModel(code, data = list(y=y), constants = list(ntowns = ntowns, nstates = nstates, n=n, state = state, town = town, sex = sex, race1=race1,race2=race2,race3=race3,race4=race4,live1=live1,live2=live2), inits = list(beta0 = 0, beta_sex = 0, beta_live = rep(0,2), beta_race = rep(0,4), sigma_state = 1, sigma_town = 1, u1 = rnorm(nstates), u2 = rnorm(ntowns)), calculate = FALSE, buildDerivs = TRUE) 
    cm <- compileNimble(m)
    
    approx <- buildNestedApprox(m, latentNodes = c('u1','u2','beta0','beta_sex','beta_race','beta_live'),
                                paramNodes = c('sigma_state','sigma_town'))
    capprox <- compileNimble(approx, project = m)
    result <- runNestedApprox(capprox)
    expect_lt(max(abs(qs_mcmc[,c('sigma_state','sigma_town')] - unlist(result$quantiles))), .008)  # .0065
    
    result$improveParamMarginals(c("sigma_state","sigma_town"), nMarginalGrid = 5)
    expect_lt(max(abs(qs_mcmc[,c('sigma_state','sigma_town')] - unlist(result$quantiles))), .007)  # .0057

    latent_sample <- result$sampleLatents(10000)
    qs_nest <- apply(latent_sample,2, quantile, qpts)
    
    expect_lt(max(abs(qs_mcmc[,1:8] - qs_nest[,c(1,3:8,2)])), .006)  # .0047; Fixed effects
                                        
    expect_lt(max(abs(qs_mcmc[,c(11:319)] - qs_nest[,9:317])), .07)  ## Max is 0.057 for us and 0.028 for INLA samples (x w/o mean/skew correction) and 0.021 for INLA marginals.
    ## At some point I got a max of .023 here...?
})

test_that("inhaler (Dirichlet) example", {

    load(system.file(file.path('tests', 'testthat', 'inhaler.Rda'), package = 'nimbleQuad'))

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

    if(FALSE) {
        library(INLA)
        fit <- inla(rating ~ treat + period + carry, data = inhaler, family='pom', quantiles = qpts,
                    control.family=list(hyper=list(theta1=list(prior="dirichlet", param=3))))
        qs_inla_fixed <- fit$summary.fixed[ , paste0(as.character(qpts), "quant")]
        qs_inla_theta <- fit$summary.hyperpar[ , paste0(as.character(qpts), "quant")]
    }

    if(FALSE) {
        library(nimbleHMC)
        K <- 4
        m <- nimbleModel(code, data = list(rating = inhaler$rating),inits = list(psi = rep(.25, 4), beta_int = 0, beta_treat = 0, beta_period = 0, beta_carry = 0),
                         constants = list(K = K, n = nrow(inhaler), period = inhaler$period, carry = inhaler$carry, treat = inhaler$treat, threes = rep(3, K)), buildDerivs = TRUE)
        cm <- compileNimble(m)
        mcmc <- buildHMC(m, monitors = c('beta_int','beta_treat','beta_period','beta_carry',
                                         'psi','alpha'))
        cmcmc <- compileNimble(mcmc, project = m)

        set.seed(1)
        out <- runMCMC(cmcmc, niter = 26000, nburnin = 1000)

        ## For comparison with nimble's nested approx
        param <- cbind(logit(out[,'psi[1]']), logit(out[,'psi[2]']/(1-out[,'psi[1]'])),
                       logit(out[,'psi[3]']/(1-out[,'psi[1]']-out[,'psi[2]'])))
        qs_param <- apply(param, 2, quantile, qpts)

        qs_param_orig <- apply(out[ , c('psi[1]','psi[2]','psi[3]','psi[4]')], 2, quantile, qpts)
        cov_param_orig <- cov(out[ , c('psi[1]','psi[2]','psi[3]','psi[4]')])

        ##  For comparison with INLA
        theta <- cbind(out[,'alpha[1]'],
                       log(out[,'alpha[2]']-out[,'alpha[1]']),
                       log(out[,'alpha[3]']-out[,'alpha[2]']))
        qs_theta <- apply(theta, 2, quantile, qpts)
        qs_mcmc <- apply(out, 2, quantile, qpts)        
        save(qs_mcmc, qs_param, qs_theta, qs_param_orig, cov_param_orig, file = 'mcmc-results6.Rda')
    } else load(system.file(file.path('tests', 'testthat', 'mcmc-results6.Rda'), package = 'nimbleQuad'))

    
    K <- 4
    m <- nimbleModel(code, data = list(rating = inhaler$rating),inits = list(psi = rep(.25, 4), beta_int = 0, beta_treat = 0, beta_period = 0, beta_carry = 0),
                     constants = list(K = K, n = nrow(inhaler), period = inhaler$period, carry = inhaler$carry, treat = inhaler$treat, threes = rep(3, K)), buildDerivs = TRUE)
    
    approx <- buildNestedApprox(m, latentNodes = c('beta_int','beta_treat','beta_period','beta_carry'), paramNodes = c('psi'))
    cm <- compileNimble(m)
    capprox <- compileNimble(approx, project = m)
    result <- runNestedApprox(capprox)
    result_orig <- runNestedApprox(capprox, originalScale = FALSE)  # Use originalScale=FALSE for comparison with INLA. 

    expect_lt(max(abs(qs_param - unlist(result_orig$quantiles))), .05)  # .040
    result_orig$improveParamMarginals(1:3, nMarginalGrid = 7)          
    expect_lt(max(abs(qs_param - unlist(result_orig$quantiles))), .04)  # .033
    ## INLA max diff is 0.026 on theta scale.

    paramSamples <- result$sampleParams(10000)
    qs_paramSamples <- apply(paramSamples, 2, quantile, qpts)
    expect_lt(max(abs(qs_param_orig - qs_paramSamples)), .025)  # .019
    expect_lt(max(abs(cov_param_orig-cov(paramSamples))), .0008)  # .00062
    
    ## We want latents on original scale. They will be either with result or result_orig, but naming clearer with result.
    latent_sample <- result$sampleLatents(10000)

    ## `beta_int` is surprisingly far off in terms of .25,.5,.75 quantiles, but other
    ## fixed effects in INLA and NIMBLE look good.
    qs_nest <- apply(latent_sample, 2, quantile, qpts)
    fixed <- c('beta_int','beta_treat','beta_period','beta_carry')
    expect_lt(max(abs(qs_mcmc[ , fixed] - qs_nest[ , fixed])), 0.3)   # .265
    # max(abs(qs_mcmc[ , fixed] - t(qs_inla_fixed)))    # .174
    # expect_lt(max(abs(qs_nest[ , fixed] - t(qs_inla_fixed))), 0.3)    # .28
    
})

test_that("Wishart example", {

    code <- nimbleCode({
        for(j in 1:J) {
            for(i in 1:n)
                y[i,j] ~ dnorm(b[j,1] + b[j,2]*x[i,j,1] + b[j,3]*x[i,j,2], sd = sigma)
            b[j,1:3] ~ dmnorm(z[1:3], Q[1:3,1:3])
        }
        Q[1:3,1:3] ~ dwish(R=R[1:3,1:3], df = 5)
        sigma ~ dhalfflat()  
    })
    
    n <- 25
    J <- 8
    
    set.seed(1)
    x <- array(rnorm(n*J*2),c(n,J,2))
    
    Sigma <- matrix(c(1,-.2,0,-.2,1,.5,0,.5,1),3)
    L <- t(chol(Sigma))
    b <- L%*%matrix(rnorm(J*3),3,J)
    
    y <- matrix(0, n, J)
    
    for(j in 1:J) {
        y[,j] <- b[1,j] + b[2,j]*x[,j,1] + b[3,j]*x[,j,2] + rnorm(n)
    }    
        
    
    if(FALSE) {
        library(nimbleHMC)
        m <- nimbleModel(code, data=list(y=y,x=x),inits = list(z = rep(0,3), Q = diag(3), b = matrix(0,J,3),sigma=1), constants = list(R=diag(3),n=n,J=J), buildDerivs = TRUE)
        cm <- compileNimble(m)
        mcmc <- buildHMC(m, monitors = c('b','Q','sigma'))
        cmcmc <- compileNimble(mcmc, project = m)
        
        system.time(out <- runMCMC(cmcmc, niter = 11000, nburnin = 1000))
        
        trans <- parameterTransform(m, nodes = 'Q')
        trSmp <- t(apply(out[,1:9], 1, trans$transform))
        out <- cbind(log(out[,'sigma']), trSmp, out[,10:ncol(out)][ , c(1,9,17,2,10,18,3,11,19,4,12,20,5,13,21,6,14,22,7,15,23,8,16,24)])
        qs_mcmc <- apply(out, 2, quantile, qpts)
        save(qs_mcmc, file = 'mcmc-results7.Rda')
    } else load(system.file(file.path('tests', 'testthat', 'mcmc-results7.Rda'), package = 'nimbleQuad'))



    m <- nimbleModel(code, data=list(y=y,x=x),inits = list(z = rep(0,3), Q = diag(3), b = matrix(0,J,3),sigma=1), constants = list(R=diag(3),n=n,J=J), buildDerivs = TRUE)
    cm <- compileNimble(m)
    
    approx <- buildNestedApprox(m)
    capprox <- compileNimble(approx, project = m)
    result <- runNestedApprox(capprox, originalScale = FALSE)  
    expect_lt(max(abs(qs_mcmc[ , 1:7] - unlist(result$quantiles))), .1)  # .092
    result$improveParamMarginals(1:7, nMarginalGrid = 5) # Much better; fairly slow to compute with so many params and grid points.
    expect_lt(max(abs(qs_mcmc[ , 1:7] - unlist(result$quantiles))), .025) # .022

    result <- runNestedApprox(capprox)
    latent_sample <- result$sampleLatents(10000)
    qs_nest <- apply(latent_sample,2,quantile,qpts) 
    expect_lt(max(abs(qs_mcmc[ , 8:31] - qs_nest)), .08)  # .065; This was <0.04 at some point when I ran it.
})

test_that("LKJ example", {
    uppertri_mult_diag <- nimbleFunction(
        run = function(mat = double(2), vec = double(1)) {
            returnType(double(2))
            p <- length(vec)
            out <- matrix(nrow = p, ncol = p, init = FALSE)
            for(i in 1:p)
                out[ , i] <- mat[ , i] * vec[i]
            return(out)
        }, buildDerivs = list(run = list(ignore='i'))
    )
    temporarilyAssignInGlobalEnv(uppertri_mult_diag)
    
    code <- nimbleCode({
        for(j in 1:J) {
            for(i in 1:n)
                y[i,j] ~ dnorm(b[j,1] + b[j,2]*x[i,j,1] + b[j,3]*x[i,j,2], sd = sigma)
            b[j,1:3] ~ dmnorm(z[1:3], cov = C[1:3,1:3])
        }
        sigma ~ dhalfflat()
        Ustar[1:3,1:3] ~ dlkj_corr_cholesky(1.3, 3)
        U[1:3,1:3] <- uppertri_mult_diag(Ustar[1:3, 1:3], sds[1:3])
        C[1:3,1:3] <- t(U[1:3,1:3])%*%U[1:3,1:3]
        for(i in 1:3)
            sds[i] ~ dhalfflat()
    })

    n <- 25
    J <- 8

    set.seed(1)
    x <- array(rnorm(n*J*2),c(n,J,2))

    Sigma <- matrix(c(1,-.2,0,-.2,1,.5,0,.5,1),3)
    L <- t(chol(Sigma))
    b <- L%*%matrix(rnorm(J*3),3,J)

    y <- matrix(0, n, J)
    for(j in 1:J) {
        y[,j] <- b[1,j] + b[2,j]*x[,j,1] + b[3,j]*x[,j,2] + rnorm(n)
    }

    if(FALSE) {
        library(nimbleHMC)
        m <- nimbleModel(code, data=list(y=y,x=x),inits = list(z = rep(0,3), Ustar = diag(3), sds = rep(1, 3), b = matrix(0,J,3),sigma=1), constants = list(n=n,J=J), buildDerivs = TRUE)
        cm <- compileNimble(m)
        mcmc <- buildHMC(m, monitors = c('b','Ustar','sigma','sds'))
        cmcmc <- compileNimble(mcmc, project = m)
        
        system.time(out <- runMCMC(cmcmc, niter = 11000, nburnin = 1000))
        
        trans <- parameterTransform(m, nodes = 'Ustar')
        trSmp <- t(apply(out[,1:9], 1, trans$transform))

        out <- cbind(log(out[,'sigma']), trSmp, log(out[ , grep("sds", colnames(out))]), out[ , grep("b\\[", colnames(out))])
        
        ## transformer <- function(tmp) {
        ##     U <- uppertri_mult_diag(matrix(tmp[1:9],3), tmp[10:12])
        ##    return(c(t(U)%*%U))
        ## }
        ## cvHMC <- t(apply(out[ ,c(1:9, 34:36)], 1, transformer))
        ## t(apply(cvHMC, 2, quantile, qpts))
        qs_mcmc <- apply(out, 2, quantile, qpts)
        save(qs_mcmc, file = 'mcmc-results8.Rda')
        
    } else load(system.file(file.path('tests', 'testthat', 'mcmc-results8.Rda'), package = 'nimbleQuad'))

    m <- nimbleModel(code, data=list(y=y,x=x),inits = list(z = rep(0,3), Ustar = diag(3), sds = rep(1, 3), b = matrix(0,J,3),sigma=1), constants = list(n=n,J=J), buildDerivs = TRUE)
    cm <- compileNimble(m)
    
    approx <- buildNestedApprox(m)
    capprox <- compileNimble(approx, project = m)
    
    result <- runNestedApprox(capprox, originalScale = FALSE) ## some rather far off
    max(abs(qs_mcmc[,1:7] - unlist(result$quantiles)))
    
    result$improveParamMarginals(1:7, nMarginalGrid=5)  # time-consuming
    expect_lt(max(abs(qs_mcmc[,1:7] - unlist(result$quantiles))), .06)  # .056, .044 with nMarginalGrid=7

    result <- runNestedApprox(capprox)
    latent_sample <- result$sampleLatents(10000)
    
    qs_nest <- apply(latent_sample,2,quantile,qpts)
    expect_lt(max(abs(qs_mcmc[,8:31] - qs_nest[ , c(1,4,7,10,13,16,19,22,2,5,8,11,14,17,20,23,3,6,9,12,15,18,21,24)])), .04) # .035
})

## CP tried to set up a test with a spatial GLMM but was stymied by a
## combination of long run times and parameter identifiability issues
## (the latter may mostly reflect using small problem sizes).
## It was hard to get good MCMC mixing and alignment between nested approx
## and MCMC.

nimbleOptions(enableDerivs = EDopt)
nimbleOptions(buildModelDerivs = BMDopt)
