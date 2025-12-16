# Tests of Laplace approximation
EDopt <- nimbleOptions("enableDerivs")
BMDopt <- nimbleOptions("buildModelDerivs")
nimbleOptions(enableDerivs = TRUE)
nimbleOptions(buildModelDerivs = TRUE)
nimbleOptions(allowDynamicIndexing = FALSE)

test_that("simple LME case works", {
  # This test uses BFGS for inner and outer optimization method.
  # nlminb results in an outer Hessian that is not negative definite.
  set.seed(1)
  g <- rep(1:10, each = 5)
  n <- length(g)
  x <- runif(n)
  m <- nimbleModel(
    nimbleCode({
      for(i in 1:n) {
        y[i] ~ dnorm((fixed_int + random_int[g[i]]) + (fixed_slope + random_slope[g[i]])*x[i], sd = sigma_res)
      }
      for(i in 1:ng) {
        random_int[i] ~ dnorm(0, sd = sigma_int)
        random_slope[i] ~ dnorm(0, sd = sigma_slope)
      }
      sigma_int ~ dunif(0, 10)
      sigma_slope ~ dunif(0, 10)
      sigma_res ~ dunif(0, 10)
      fixed_int ~ dnorm(0, sd = 100)
      fixed_slope ~ dnorm(0, sd = 100)
    }),
    constants = list(g = g, ng = max(g), n = n, x = x),
    buildDerivs = TRUE
  )
  params <- c("fixed_int", "fixed_slope", "sigma_int", "sigma_slope", "sigma_res")
  values(m, params) <- c(10, 0.5, 3, .25, 0.2)
  m$simulate(m$getDependencies(params, self = FALSE))
  m$setData('y')
  y <- m$y
  library(lme4)
  manual_fit <- lmer(y ~ x + (1 + x || g), REML = FALSE)

  mLaplace <- buildLaplace(model = m, control=list(innerOptimMethod="BFGS",
                                                   outerOptimMethod="BFGS"))
  cm <- compileNimble(m)
  cmLaplace <- compileNimble(mLaplace, project = m)
  opt <- cmLaplace$findMLE()
  nimres <- cmLaplace$summary(opt, randomEffectsStdError = TRUE)
  lme4res <- summary(manual_fit)
  expect_equal(nimres$params$estimate[4:5], as.vector(lme4res$coefficients[,"Estimate"]), tol=1e-5)
  ## 2025-09-19: CI on MacOS is failing with tol=1e-5 
  expect_equal(nimres$params$estimate[1:3], as.data.frame(VarCorr(manual_fit))[,"sdcor"], tol = 1e-4)
  expect_equal(nimres$params$stdError[4:5], as.vector(lme4res$coefficients[,"Std. Error"]), tol=0.03)
  expect_equal(nimres$randomEffects$estimate, as.vector(t(ranef(manual_fit)$g)), tol = 1e-4)
})

test_that("simple LME with correlated intercept and slope works (and check with nQuad=3)", {
  nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions=TRUE)
  set.seed(1)
  g <- rep(1:10, each = 10)
  n <- length(g)
  x <- runif(n)
  m <- nimbleModel(
    nimbleCode({
      for(i in 1:n) {
        y[i] ~ dnorm((fixed_int + random_int_slope[g[i], 1]) + (fixed_slope + random_int_slope[g[i], 2])*x[i], sd = sigma_res)
      }
      cov[1, 1] <- sigma_int^2
      cov[2, 2] <- sigma_slope^2
      cov[1, 2] <- rho * sigma_int * sigma_slope
      cov[2, 1] <- rho * sigma_int * sigma_slope
      for(i in 1:ng) {
        random_int_slope[i, 1:2] ~ dmnorm(zeros[1:2], cov = cov[1:2, 1:2])
      }
      sigma_int ~ dunif(0, 10)
      sigma_slope ~ dunif(0, 10)
      sigma_res ~ dunif(0, 10)
      fixed_int ~ dnorm(0, sd = 100)
      fixed_slope ~ dnorm(0, sd = 100)
      rho ~ dunif(-1, 1)
    }),
    constants = list(g = g, ng = max(g), n = n, x = x, zeros = rep(0, 2)),
    buildDerivs = TRUE
  )
  params <- c("fixed_int", "fixed_slope", "sigma_int", "sigma_slope", "sigma_res", "rho")
  values(m, params) <- c(10, 0.5, 3, 0.25, 0.2, 0.45)
  m$simulate(m$getDependencies(params, self = FALSE))
  m$setData('y')
  y <- m$y
  library(lme4)
  manual_fit <- lmer(y ~ x + (1 + x | g), REML = FALSE)
  mLaplace <- buildLaplace(model = m)#, control=list(innerOptimStart="last.best"))
  cm <- compileNimble(m)
  cmLaplace <- compileNimble(mLaplace, project = m)

  params_in_order <- setupMargNodes(m)$paramNodes

  pStart <- values(m, params_in_order)

  init_llh <- cmLaplace$calcLogLik(pStart)
  init_gr_llh <- cmLaplace$gr_logLik(pStart)

  opt <- cmLaplace$findMLE()
  nimres <- cmLaplace$summary(opt, randomEffectsStdError = TRUE)
  nimsumm <- summaryLaplace(cmLaplace, opt, randomEffectsStdError = TRUE)

  lme4res <- summary(manual_fit)
  expect_equal(nimres$params$estimate[4:5], as.vector(lme4res$coefficients[,"Estimate"]), tol=1e-4)
  sdparams <- nimres$params$estimate[-c(4,5)]
  expect_equal(sdparams[c(1,2,4,3)], as.data.frame(VarCorr(manual_fit))[,"sdcor"], tol = 1e-3)
  expect_equal(nimres$params$stdError[4:5], as.vector(lme4res$coefficients[,"Std. Error"]), tol=.03)
  expect_equal(nimres$randomEffects$estimate, as.vector(t(ranef(manual_fit)$g)), tol = 5e-3)

  cmLaplace$updateSettings(nQuad = 3)
  init_llh_3 <- cmLaplace$calcLogLik(pStart)
  max_llh_3 <- cmLaplace$calcLogLik(opt$par  )
  expect_equal(init_llh, init_llh_3, tolerance = 1e-7)
  expect_equal(opt$value, max_llh_3, tolerance = 1e-4)

  for(v in m$getVarNames()) cm[[v]] <- m[[v]]
  cm$calculate()
  cmLaplace$updateSettings(useInnerCache=FALSE)
  cmLaplace$updateSettings(nQuad = 1)
  CrunLaplaceRes <- runLaplace(cmLaplace, pStart = pStart)
  expect_equal(opt$par, CrunLaplaceRes$MLE$par, tolerance = 1e-4)
  expect_equal(opt$hessian, CrunLaplaceRes$MLE$hessian, tolerance = 1e-4)
  expect_equal(nimsumm$randomEffects$estimate,
               CrunLaplaceRes$summary$randomEffects$estimate, tolerance = 1e-4)
  expect_equal(nimsumm$randomEffects$stdError,
               CrunLaplaceRes$summary$randomEffects$stdError, tolerance = 1e-4)


  for(v in m$getVarNames()) cm[[v]] <- m[[v]]
  cm$calculate()
  cmLaplace$updateSettings(useInnerCache=FALSE)
  CrunLaplaceRes <- runLaplace(cmLaplace, pStart = pStart)
  expect_equal(opt$par, CrunLaplaceRes$MLE$par, tolerance = 1e-4)
  expect_equal(opt$hessian, CrunLaplaceRes$MLE$hessian, tolerance = 1e-4)
  expect_equal(nimsumm$randomEffects$estimate,
               CrunLaplaceRes$summary$randomEffects$estimate, tolerance = 1e-4)
  expect_equal(nimsumm$randomEffects$stdError,
               CrunLaplaceRes$summary$randomEffects$stdError, tolerance = 1e-4)

})

test_that("simple LME with correlated intercept and slope works through runLaplace", {
  set.seed(1)
  g <- rep(1:10, each = 10)
  n <- length(g)
  x <- runif(n)
  m <- nimbleModel(
    nimbleCode({
      for(i in 1:n) {
        y[i] ~ dnorm((fixed_int + random_int_slope[g[i], 1]) + (fixed_slope + random_int_slope[g[i], 2])*x[i], sd = sigma_res)
      }
      cov[1, 1] <- sigma_int^2
      cov[2, 2] <- sigma_slope^2
      cov[1, 2] <- rho * sigma_int * sigma_slope
      cov[2, 1] <- rho * sigma_int * sigma_slope
      for(i in 1:ng) {
        random_int_slope[i, 1:2] ~ dmnorm(zeros[1:2], cov = cov[1:2, 1:2])
      }
      sigma_int ~ dunif(0, 10)
      sigma_slope ~ dunif(0, 10)
      sigma_res ~ dunif(0, 10)
      fixed_int ~ dnorm(0, sd = 100)
      fixed_slope ~ dnorm(0, sd = 100)
      rho ~ dunif(-1, 1)
    }),
    constants = list(g = g, ng = max(g), n = n, x = x, zeros = rep(0, 2)),
    buildDerivs = TRUE
  )
  params <- c("fixed_int", "fixed_slope", "sigma_int", "sigma_slope", "sigma_res", "rho")
  values(m, params) <- c(10, 0.5, 3, 0.25, 0.2, 0.45)
  m$simulate(m$getDependencies(params, self = FALSE))
  m$setData('y')
  y <- m$y
  library(lme4)
  manual_fit <- lmer(y ~ x + (1 + x | g), REML = FALSE)
  mLaplace <- buildLaplace(model = m)#, control=list(innerOptimStart="last.best"))
  cm <- compileNimble(m)
  cmLaplace <- compileNimble(mLaplace, project = m)

  pStart <- values(m, params)

  res <- runLaplace(cmLaplace)
  opt <- res$MLE
  nimsumm <- res$summary

  #opt <- cmLaplace$findMLE()
  #nimres <- cmLaplace$summary(opt, randomEffectsStdError = TRUE)
  #nimsumm <- summaryLaplace(cmLaplace, opt, randomEffectsStdError = TRUE)

  lme4res <- summary(manual_fit)
  expect_equal(nimsumm$params$estimate[4:5], as.vector(lme4res$coefficients[,"Estimate"]), tol=1e-4)
  sdparams <- nimsumm$params$estimate[-c(4,5)]
  expect_equal(sdparams[c(1,2,4,3)], as.data.frame(VarCorr(manual_fit))[,"sdcor"], tol = 1e-3)
  expect_equal(nimsumm$params$stdError[4:5], as.vector(lme4res$coefficients[,"Std. Error"]), tol=.03)
  expect_equal(nimsumm$randomEffects$estimate, as.vector(t(ranef(manual_fit)$g)), tol = 5e-3)
})

test_that("Laplace with non-empty calcNodesOther works", {
  m <- nimbleModel(
    nimbleCode({
      for(i in 1:3) {
        mu[i] ~ dnorm(0, sd = 10)
      }
      mu_a[1] <- mu[1] + mu[2]
      mu_a[2] <- mu[2] + mu[3]
      a[1] ~ dnorm(mu_a[1], sd = 2)
      y[1] ~ dnorm(a[1], sd = 3)
      a[2] ~ dnorm(mu_a[2], sd = 2)
      y[2] ~ dnorm(a[2], sd =3)
      y[3] ~ dnorm(mu[3], sd = 3)
    }),
    data = list(y = c(2, 3, 5)),
    inits = list(a = c(1, 2), mu = c(1, 2, 3)),
    buildDerivs = TRUE
  )

  mLaplace <- buildLaplace(model = m)
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit

  opt <- cmLaplace$findMLE()
  expect_equal(opt$par, c(4, -2, 5), tol = 1e-3)
  expect_equal(opt$value, -6.420377, tol = 1e-6)

  ## Check covariance matrix
  summ <- cmLaplace$summary(opt, jointCovariance = TRUE)
  ## TMB cpp code:
  #include <TMB.hpp>
  #template<class Type>
  #Type objective_function<Type>::operator() ()
  # {
  #   DATA_VECTOR(y);
  #   PARAMETER_VECTOR(mu);
  #   PARAMETER_VECTOR(a);
  #   int i;
  #   // Negative log-likelihood
  #   Type ans = -dnorm(a[0], mu[0]+mu[1], Type(2.0), true);
  #   ans -= dnorm(a[1], mu[1]+mu[2], Type(2.0), true);
  #   for(i = 0; i < 2; i++){
  #     ans -= dnorm(y[i], a[i], Type(3.0), true);
  #   }
  #   ans -= dnorm(y[2], mu[2], Type(3.0), true);
  #   return ans;
  # }
  ## TMB R code:
  # library(TMB)
  # compile("test.cpp")
  # dyn.load(dynlib("test"))
  # data <- list(y = m$y)
  # parameters <- list(mu = c(1, 2, 3), a = c(1, 2))
  #
  # ## Fit model
  # obj <- MakeADFun(data, parameters, random="a", DLL="test")
  # tmbres <- nlminb(obj$par, obj$fn, obj$gr)
  # tmbrep <- sdreport(obj, getJointPrecision = TRUE)
  # tmbvcov <- inverse(tmbrep$jointPrecision)

  ## Covariance matrix from TMB
  tmbvcov <- matrix(nrow = 5, ncol = 5)
  tmbvcov[1,] <- c( 35, -2.20000e+01,  9.000000e+00,  9.000000e+00, -9.000000e+00)
  tmbvcov[2,] <- c(-22,  2.20000e+01, -9.000000e+00,  8.463230e-13,  9.000000e+00)
  tmbvcov[3,] <- c( 9,  -9.00000e+00,  9.000000e+00, -3.462231e-13,  3.462231e-13)
  tmbvcov[4,] <- c( 9,   8.46323e-13, -3.462231e-13,  9.000000e+00,  3.462231e-13)
  tmbvcov[5,] <- c(-9,   9.00000e+00,  3.462231e-13,  3.462231e-13,  9.000000e+00)

  expect_equal(summ$vcov, tmbvcov, tol=1e-5)

  ## Check covariance matrix for params only
  tryResult <- try({
      summ2 <- cmLaplace$summary(opt, originalScale = TRUE, randomEffectsStdError = TRUE, jointCovariance = FALSE)
      expect_equal(summ2$vcov, tmbvcov[1:3,1:3], tol=1e-5)
  })
  if(inherits(tryResult, 'try-error')) {
      print(class(cmLaplace))
      print(cL)
  }


  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$findMLE()
  expect_equal(opt$par, optNoSplit$par, tol = 1e-2)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)
})

test_that("Laplace with 2x1D parameters (one needs transformation) and non-normal data works", {
  m <- nimbleModel(
    nimbleCode({
      mu ~ dnorm(0, sd = 10.0)
      sigma ~ dunif(0, 100)
      for (i in 1:5){
        theta[i] ~ dnorm(mu, sd = sigma)
        logit(p[i]) <- theta[i]
        y[i] ~ dbinom(10, prob = p[i])
      }
    }),
    data = list(y = c(8, 6, 5, 3, 7)),
    inits = list(mu = 1, sigma = 1, theta = rep(0, 5)),
    buildDerivs = TRUE
  )

  mLaplace <- buildLaplace(model = m)
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit

  opt <- cmLaplace$findMLE()
  ## Compare with results from TMB
  expect_equal(opt$par, c(0.330241, 0.3059177), tol = 1e-4)
  expect_equal(opt$value, -9.703857, tol = 1e-6)
  ## Check covariance matrix on the transformed scale
  summ <- cmLaplace$summary(opt, originalScale = FALSE, jointCovariance = TRUE)
  tmbvcov <- matrix(nrow = 7, ncol = 7)
  tmbvcov[1,] <- c(0.10337427,  0.04574391,  0.09719623, 0.08526807,  0.07943536,  0.06797944,  0.09118502)
  tmbvcov[2,] <- c(0.04574391,  3.21994672,  0.91522073, 0.10980129, -0.28810783, -1.07845809,  0.51064309)
  tmbvcov[3,] <- c(0.09719623,  0.91522073,  0.40584816, 0.09981763, -0.01342937, -0.23826114,  0.21393310)
  tmbvcov[4,] <- c(0.08526807,  0.10980129,  0.09981763, 0.14821768,  0.05824110,  0.03110420,  0.08580658)
  tmbvcov[5,] <- c(0.07943536, -0.28810783, -0.01342937, 0.05824110,  0.16979550,  0.16423022,  0.02255625)
  tmbvcov[6,] <- c(0.06797944, -1.07845809, -0.23826114, 0.03110420,  0.16423022,  0.50464751, -0.10296956)
  tmbvcov[7,] <- c(0.09118502,  0.51064309,  0.21393310, 0.08580658,  0.02255625, -0.10296956,  0.22602059)
  expect_equal(summ$vcov, tmbvcov, tol=1e-3)
  ## Stand error for sigma (original parameter)
  summ2 <- cmLaplace$summary(opt, originalScale = TRUE)
  expect_equal(summ2$params$stdError[2], 0.5472659, tol=1e-4)

  # Check covariance matrix for transformed params only
  summ3 <- cmLaplace$summary(opt, originalScale = FALSE, randomEffectsStdError = TRUE, jointCovariance = FALSE)
  expect_equal(summ3$vcov, tmbvcov[1:2,1:2], tol=1e-3)

  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$findMLE()
  expect_equal(opt$par, optNoSplit$par, tol = 1e-2)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)
  ## TMB cpp code:
  #include <TMB.hpp>
  #template<class Type>
  #Type objective_function<Type>::operator() ()
  # {
  #   DATA_VECTOR(y);
  #   PARAMETER(mu);
  #   PARAMETER(sigmaTrans);
  #   PARAMETER_VECTOR(theta);
  #   // Transformation for sigma
  #   Type sigma = 100 * exp(sigmaTrans) / (1 + exp(sigmaTrans));
  #   // Negative log-likelihood
  #   Type ans = 0;
  #   vector<Type> p(5);
  #   for(int i = 0; i < 5; i++){
  #     p[i] = exp(theta[i]) / (1 + exp(theta[i]));
  #     ans -= dnorm(theta[i], mu, sigma, true) + dbinom(y[i], Type(10), p[i], true);
  #   }
  #   ADREPORT(sigma);
  #   return ans;
  # }
  ## TMB R code:
  # library(TMB)
  # compile("test.cpp")
  # dyn.load(dynlib("test"))
  # data <- list(y = m$y)
  # parameters <- list(mu = m$mu, sigmaTrans = logit(m$sigma/100), theta = m$theta)
  # ## Fit model
  # obj <- MakeADFun(data, parameters, random="theta", DLL="test")
  # tmbopt <- nlminb(obj$par, obj$fn, obj$gr)
  # tmbrep <- sdreport(obj, getJointPrecision = TRUE)
  # tmbvcov <- inverse(tmbrep$jointPrecision)
})

test_that("Laplace with no random effects (simple linear regression) works", {
  set.seed(1)
  x <- rnorm(5)
  y <- sapply(-1 + x, rnorm, n = 1, sd = 1)
  m <- nimbleModel(
    nimbleCode({
      a ~ dnorm(0, sd = 10.0)
      b ~ dnorm(0, sd = 10.0)
      sigma ~ dunif(0, 100)
      for(i in 1:5){
        mu_y[i] <- a + b*x[i]
        y[i] ~ dnorm(mu_y[i], sd = sigma)
      }
    }),
    constants = list(x = x),
    data = list(y = y),
    inits = list(a = -1, b = 1, sigma = 1),
    buildDerivs = TRUE
  )

  mLaplace <- buildLaplace(model = m)
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit

  opt <- cmLaplace$findMLE()
  summ <- cmLaplace$summary(opt)
  ## Compare results with those from TMB
  expect_equal(opt$par, c(-0.8899436, 1.1940911, 0.5744841), tol = 1e-4)
  expect_equal(opt$value, -4.323288, tol = 1e-7)
  expect_equal(summ$params$stdError, c(0.2598061, 0.2988869, 0.1816661), tol = 1e-5)

  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$findMLE()
  expect_equal(opt$par, optNoSplit$par, tol = 1e-4)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)

  summL <- summaryLaplace(cmLaplace, opt, randomEffectsStdError = TRUE, jointCovariance = TRUE)
  expect_equal(nrow(summL$randomEffects), 0)
  expect_equal(nrow(summL$vcov), 3)
  ## TMB cpp code
  #include <TMB.hpp>
  #template<class Type>
  # Type objective_function<Type>::operator() ()
  # {
  #   DATA_VECTOR(y);
  #   DATA_VECTOR(x);
  #   PARAMETER(a);
  #   PARAMETER(b);
  #   PARAMETER(sigma);
  #   Type nll = -sum(dnorm(y, a+b*x, sigma, true));
  #   return nll;
  # }
  ## R code
  # compile("lm.cpp")
  # dyn.load(dynlib("lm"))
  # set.seed(1)
  # x <- rnorm(5)
  # y <- sapply(-1 + x, rnorm, n = 1, sd = 1)
  # data <- list(y=y, x=x)
  # parameters <- list(a=-1, b=1, sigma=1)
  # obj <- MakeADFun(data, parameters, DLL="lm")
  # obj$hessian <- TRUE
  # tmbres <- do.call("optim", obj)
  # tmbsumm <- summary(sdreport(obj))
})

## Possible future feature (was drafted, not completed):
##
## test_that("Laplace with no priors for unconstrained parameters works", {
##   ## Here we re-use some of tests above and remove priors for parameters
##   ## Test 1
##   m <- nimbleModel(
##     nimbleCode({
##       y ~ dnorm(a, sd = 2)
##       a ~ dnorm(mu, sd = 3)
##       # mu ~ dnorm(0, sd = 5)
##     }), data = list(y = 4), inits = list(a = -1),
##     buildDerivs = TRUE
##   )

##   mLaplace <- buildLaplace(model = m, control = list(allowNonPriors = TRUE))
##   mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE, allowNonPriors = TRUE))
##   cm <- compileNimble(m)
##   cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
##   cmLaplace <- cL$mLaplace
##   cmLaplaceNoSplit <- cL$mLaplaceNoSplit

##   opt <- cmLaplace$findMLE()
##   expect_equal(opt$par, 4, tol = 1e-4)
##   expect_equal(opt$value, dnorm(4, 4, sd = sqrt(13), log = TRUE))
##   summ <- cmLaplace$summary(opt, originalScale = TRUE, randomEffectsStdError = TRUE, jointCovariance = TRUE)
##   expect_equal(summ$randomEffects$estimate, 4, tol = 1e-5)
##   # Covariance matrix
##   vcov <- matrix(c(1/(1/4+1/9), 0, 0, 0), nrow = 2) + matrix(c(4/13, 1), ncol = 1) %*% (13) %*% t(matrix(c(4/13, 1), ncol = 1))
##   expect_equal(vcov, summ$vcov, tol = 1e-6)

##   for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
##   optNoSplit <- cmLaplaceNoSplit$findMLE()
##   expect_equal(opt$par, optNoSplit$par, tol = 1e-2)
##   expect_equal(opt$value, optNoSplit$value, tol = 1e-7)

##   ## Test 2
##   set.seed(1)
##   x <- rnorm(5)
##   y <- sapply(-1 + x, rnorm, n = 1, sd = 1)
##   m <- nimbleModel(
##     nimbleCode({
##       sigma ~ dunif(0, 100)
##       for(i in 1:5){
##         mu_y[i] <- a + b*x[i]
##         y[i] ~ dnorm(mu_y[i], sd = sigma)
##       }
##     }),
##     constants = list(x = x),
##     data = list(y = y),
##     buildDerivs = TRUE
##   )

##   mLaplace <- buildLaplace(model = m, control = list(allowNonPriors = TRUE))

##   mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE, allowNonPriors = TRUE))
##   cm <- compileNimble(m)
##   cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
##   cmLaplace <- cL$mLaplace
##   cmLaplaceNoSplit <- cL$mLaplaceNoSplit

##   opt <- cmLaplace$findMLE()
##   summ <- cmLaplace$summary(opt)
##   ## Compare results with those from TMB
##   expect_equal(opt$par, c(0.5744841, -0.8899436, 1.1940911), tol = 1e-5)
##   expect_equal(opt$value, -4.323288, tol = 1e-7)
##   expect_equal(summ$params$stdError, c(0.1816661, 0.2598061, 0.2988869), tol = 1e-5)

##   for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
##   optNoSplit <- cmLaplaceNoSplit$findMLE()
##   expect_equal(opt$par, optNoSplit$par, tol = 1e-2)
##   expect_equal(opt$value, optNoSplit$value, tol = 1e-7)

##   ## Test 3
##   set.seed(1)
##   y <- array(rnorm(8, 6, 5), dim = c(2, 2, 2))
##   cov_a <- matrix(c(2, 1.5, 1.5, 2), nrow = 2)
##   m <- nimbleModel(
##     nimbleCode({
##       # for(i in 1:2) mu[i] ~ dnorm(0, sd = 10)
##       mu_a[1] <- 0.8 * mu[1]
##       mu_a[2] <- 0.2 * mu[2]
##       for(i in 1:2) a[i, 1:2] ~ dmnorm(mu_a[1:2], cov = cov_a[1:2, 1:2])
##       for(i in 1:2) {
##         for(j in 1:2) {
##           y[1, j, i] ~ dnorm( 0.5 * a[i, 1], sd = 1.8)
##           y[2, j, i] ~ dnorm( 0.1 * a[i, 2], sd = 1.2)
##         }
##       }
##     }),
##     data = list(y = y),
##     inits = list(a = matrix(c(-2, -3, 0,  -1), nrow = 2)),
##     constants = list(cov_a = cov_a),
##     buildDerivs = TRUE
##   )

##   mLaplace <- buildLaplace(model = m, control = list(allowNonPriors = TRUE))
##   mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE, allowNonPriors = TRUE))
##   cm <- compileNimble(m)
##   cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
##   cmLaplace <- cL$mLaplace
##   cmLaplaceNoSplit <- cL$mLaplaceNoSplit

##   opt <- cmLaplace$findMLE()

##   expect_equal(opt$par, c(12.98392, 406.04878), tol = 1e-4)
##   expect_equal(opt$value, -41.86976, tol = 1e-6)
##   # Check covariance matrix
##   summ <- cmLaplace$summary(opt, jointCovariance = TRUE)
##   tmbvcov <- matrix(nrow = 6, ncol = 6)
##   tmbvcov[1,] <- c(6.625000e+00, 4.687500e+00,  4.050000e+00,  4.050000e+00, -2.693817e-11, -2.695275e-11)
##   tmbvcov[2,] <- c(4.687500e+00, 9.250000e+02,  2.965628e-11,  2.967848e-11,  1.800000e+02,  1.800000e+02)
##   tmbvcov[3,] <- c(4.050000e+00, 2.951367e-11,  3.995242e+00,  2.484758e+00,  5.596302e-01, -5.596302e-01)
##   tmbvcov[4,] <- c(4.050000e+00, 2.951367e-11,  2.484758e+00,  3.995242e+00, -5.596302e-01,  5.596302e-01)
##   tmbvcov[5,] <- c(-2.691772e-11, 1.800000e+02,  5.596302e-01, -5.596302e-01,  3.684693e+01,  3.515307e+01)
##   tmbvcov[6,] <- c(-2.691772e-11, 1.800000e+02, -5.596302e-01,  5.596302e-01,  3.515307e+01,  3.684693e+01)

##   expect_equal(summ$vcov[c(5,6,1,3,2,4), c(5,6,1,3,2,4)], tmbvcov, tol = 1e-4)

##   for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
##   optNoSplit <- cmLaplaceNoSplit$findMLE()
##   expect_equal(opt$par, optNoSplit$par, tol = 1e-4)
##   expect_equal(opt$value, optNoSplit$value, tol = 1e-7)

##   ## Test 4
##   m <- nimbleModel(
##     nimbleCode({
##       # for(i in 1:3) {
##       #   mu[i] ~ dnorm(0, sd = 10)
##       # }
##       mu_a[1] <- mu[1] + mu[2]
##       mu_a[2] <- mu[2] + mu[3]
##       a[1] ~ dnorm(mu_a[1], sd = 2)
##       y[1] ~ dnorm(a[1], sd = 3)
##       a[2] ~ dnorm(mu_a[2], sd = 2)
##       y[2] ~ dnorm(a[2], sd =3)
##       y[3] ~ dnorm(mu[3], sd = 3)
##     }),
##     data = list(y = c(2, 3, 5)),
##     inits = list(a = c(1, 2)),
##     buildDerivs = TRUE
##   )

##   mLaplace <- buildLaplace(model = m, control = list(allowNonPriors = TRUE))
##   mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE, allowNonPriors = TRUE))
##   cm <- compileNimble(m)
##   cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
##   cmLaplace <- cL$mLaplace
##   cmLaplaceNoSplit <- cL$mLaplaceNoSplit

##   opt <- cmLaplace$findMLE()
##   expect_equal(opt$par, c(4, -2, 5), tol = 1e-3)
##   expect_equal(opt$value, -6.420377, tol = 1e-6)
##   ## Check covariance matrix
##   summ <- cmLaplace$summary(opt, jointCovariance = TRUE)

##   ## Covariance matrix from TMB
##   tmbvcov <- matrix(nrow = 5, ncol = 5)
##   tmbvcov[1,] <- c( 35, -2.20000e+01,  9.000000e+00,  9.000000e+00, -9.000000e+00)
##   tmbvcov[2,] <- c(-22,  2.20000e+01, -9.000000e+00,  8.463230e-13,  9.000000e+00)
##   tmbvcov[3,] <- c( 9,  -9.00000e+00,  9.000000e+00, -3.462231e-13,  3.462231e-13)
##   tmbvcov[4,] <- c( 9,   8.46323e-13, -3.462231e-13,  9.000000e+00,  3.462231e-13)
##   tmbvcov[5,] <- c(-9,   9.00000e+00,  3.462231e-13,  3.462231e-13,  9.000000e+00)

##   expect_equal(summ$vcov[c(3:5, 1:2), c(3:5, 1:2)], tmbvcov, tol=1e-5)

##   for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
##   optNoSplit <- cmLaplaceNoSplit$findMLE()
##   expect_equal(opt$par, optNoSplit$par, tol = 1e-2)
##   expect_equal(opt$value, optNoSplit$value, tol = 1e-7)

## })

test_that("Laplace with crossed random effects works", {
  library(lme4)
  data(Penicillin)
  N <- nrow(Penicillin)
  plate <- rep(1:24, each = 6)
  np <- 24
  sample <- rep(1:6, 24)
  ns <- 6

  m <- nimbleModel(
    nimbleCode({
      ## Intercept
      beta ~ dnorm(0, sd = 100)
      ## Standard deviations
      sigma ~ dgamma(1.0, 1.0)
      sigma_p ~ dgamma(1.0, 1.0)
      sigma_s ~ dgamma(1.0, 1.0)
      ## Random effects for plate
      for(i in 1:np){
        mup[i] ~ dnorm(0, sd = sigma_p)
      }
      ## Random effects for sample
      for(i in 1:ns){
        mus[i] ~ dnorm(0, sd = sigma_s)
      }
      ## Observations
      for(i in 1:N){
        mu_y[i] <- beta + mus[sample[i]] + mup[plate[i]]
        y[i] ~ dnorm(mu_y[i], sd = sigma)
      }
    }),
    constants = list(N = N, np = np, ns = ns, plate = plate, sample = sample),
    data = list(y = Penicillin$diameter),
    inits = list(beta = 20, sigma = 1, sigma_p = 1, sigma_s = 1, mus = rep(0, ns), mup = rep(0, np)),
    buildDerivs = TRUE
  )
  mLaplace <- buildLaplace(model = m)#, control=list(innerOptimStart = "last.best"))
  cm <- compileNimble(m)
  cmLaplace <- compileNimble(mLaplace, project = m)
  ## cmLaplace$updateSettings(innerOptimMethod = "nlminb")
  opt <- cmLaplace$findMLE()
  nimres <- cmLaplace$summary(opt, randomEffectsStdError = TRUE)

  lme4_fit <- lmer(diameter ~ 1 + (1|plate) + (1|sample), data = Penicillin, REML = FALSE)
  lme4res <- summary(lme4_fit)

  expect_equal(nimres$params$estimate[1], lme4res$coefficients[,"Estimate"], tol=1e-3)
  expect_equal(nimres$params$estimate[c(3,4,2)], as.data.frame(VarCorr(lme4_fit))[,"sdcor"], tol = 5e-4)
  # Note that with innerOptimMethod "nlminb", the next check is far off, within only about 0.2
  # on Mac, and getting a NaN on ubuntu CI tests. (Also I don't know why those differ.)
  expect_equal(nimres$params$stdError[1], lme4res$coefficients[,"Std. Error"], tol=2e-3)
  expect_equal(nimres$randomEffects$estimate[25:30], as.vector(t(ranef(lme4_fit)$sample)), tol = 1e-3)
  expect_equal(nimres$randomEffects$estimate[1:24], as.vector(t(ranef(lme4_fit)$plate)), tol = 1e-4)
})

test_that("Laplace with crossed random effects works, without using normality", {
  library(lme4)
  data(Penicillin)
  N <- nrow(Penicillin)
  plate <- rep(1:24, each = 6)
  np <- 24
  sample <- rep(1:6, 24)
  ns <- 6

  m <- nimbleModel(
    nimbleCode({
      ## Intercept
      beta ~ dnorm(0, sd = 100)
      ## Standard deviations
      sigma ~ dgamma(1.0, 1.0)
      sigma_p ~ dgamma(1.0, 1.0)
      sigma_s ~ dgamma(1.0, 1.0)
      ## Random effects for plate
      for(i in 1:np){
        mup[i] ~ dnorm(0, sd = sigma_p)
      }
      ## Random effects for sample
      for(i in 1:ns){
        mus[i] ~ dnorm(0, sd = sigma_s)
      }
      ## Observations
      for(i in 1:N){
        mu_y[i] <- beta + mus[sample[i]] + mup[plate[i]]
        y[i] ~ dnorm(mu_y[i], sd = sigma)
      }
    }),
    constants = list(N = N, np = np, ns = ns, plate = plate, sample = sample),
    data = list(y = Penicillin$diameter),
    inits = list(beta = 20, sigma = 1, sigma_p = 1, sigma_s = 1, mus = rep(0, ns), mup = rep(0, np)),
    buildDerivs = TRUE
  )
  mLaplace <- buildLaplace(model = m, control = list(ADuseNormality = FALSE))#, control=list(innerOptimStart = "last.best"))
  cm <- compileNimble(m)
  cmLaplace <- compileNimble(mLaplace, project = m)
  ## cmLaplace$updateSettings(innerOptimMethod = "nlminb")
  opt <- cmLaplace$findMLE()
  nimres <- cmLaplace$summary(opt, randomEffectsStdError = TRUE)

  lme4_fit <- lmer(diameter ~ 1 + (1|plate) + (1|sample), data = Penicillin, REML = FALSE)
  lme4res <- summary(lme4_fit)

  expect_equal(nimres$params$estimate[1], lme4res$coefficients[,"Estimate"], tol=1e-3)
  expect_equal(nimres$params$estimate[c(3,4,2)], as.data.frame(VarCorr(lme4_fit))[,"sdcor"], tol = 5e-4)
  # Note that with innerOptimMethod "nlminb", the next check is far off, within only about 0.2
  # on Mac, and getting a NaN on ubuntu CI tests. (Also I don't know why those differ.)
  expect_equal(nimres$params$stdError[1], lme4res$coefficients[,"Std. Error"], tol=2e-3)
  expect_equal(nimres$randomEffects$estimate[25:30], as.vector(t(ranef(lme4_fit)$sample)), tol = 1e-3)
  expect_equal(nimres$randomEffects$estimate[1:24], as.vector(t(ranef(lme4_fit)$plate)), tol = 1e-4)
})


test_that("Laplace with nested random effects works", {
  library(lme4)
  data(Pastes)
  lme4_fit <- lmer(strength ~ 1 + (1|batch) + (1|batch:cask), data = Pastes, REML = FALSE)
  lme4res <- summary(lme4_fit)

  m <- nimbleModel(
    nimbleCode({
      ## Intercept
      beta ~ dnorm(0, sd = 100)
      ## Standard deviations
      sigma ~ dgamma(1.0, 1.0)
      sigma1 ~ dgamma(1.0, 1.0)
      sigma2 ~ dgamma(1.0, 1.0)
      ## Random effects for batch
      for(i in 1:10){
        mub[i] ~ dnorm(0, sd = sigma1)
      }
      ## Random effects for batch:cask
      for(i in 1:30){
        mubc[i] ~ dnorm(0, sd = sigma2)
      }
      ## Observations
      for(i in 1:60){
        mu_y[i] <- beta + mub[batch[i]] + mubc[cask[i]]
        y[i] ~ dnorm(mu_y[i], sd = sigma)
      }
    }),
    constants = list(batch = rep(1:10, each = 6), cask = rep(1:30, each = 2)),
    data = list(y = Pastes$strength),
    buildDerivs = TRUE
  )
  mLaplace <- buildLaplace(model = m)
  cm <- compileNimble(m)
  cmLaplace <- compileNimble(mLaplace, project = m)
  ## It seems that default start values (0, 1, 1, 1) for this example do not work well
  ## for optimisation; use c(2, 2, 2, 2) instead
  #expect_output(
    opt <- cmLaplace$findMLE(pStart = c(2,2,2,2))
  #, "optim does not converge for the inner optimization")
  nimres <- cmLaplace$summary(opt, randomEffectsStdError = TRUE)

  expect_equal(nimres$params$estimate[1], lme4res$coefficients[,"Estimate"], tol = 1e-5)
  expect_equal(nimres$params$estimate[c(4, 3, 2)], as.data.frame(VarCorr(lme4_fit))[,"sdcor"], tol = 5e-5)
  expect_equal(nimres$params$stdError[1], lme4res$coefficients[,"Std. Error"], tol = 5e-5)
  expect_equal(nimres$randomEffects$estimate[seq(1, 40, by = 4)], as.vector(t(ranef(lme4_fit)$batch)), tol = 5e-4)
  expect_equal(nimres$randomEffects$estimate[-seq(1, 40, by = 4)], as.vector(t(ranef(lme4_fit)$`batch:cask`)), tol = 5e-4)
})

test_that("Laplace error trapping of wrong-length parameters works", {
  library(nimble)
  library(testthat)

  m <- nimbleModel(
    nimbleCode({
      d[1:3] ~ ddirch(alpha[1:3]) # params
      for(i in 1:3) x[i] ~ dnorm(d[i], 1) # randomEffects
      for(i in 1:3) y[i] ~ dnorm(x[i], 1) # data
    }),
    data = list(y = rnorm(3), alpha = rep(1.1, 3)),
    inits = list(x = rnorm(3), d = c(.2, .3, .5)),
    buildDerivs = TRUE
  )
  m$calculate()
  mLaplace <- buildLaplace(model = m)
  cm <- compileNimble(m)
  cmLaplace <- compileNimble(mLaplace, project = m)

  ## cat("Eight messages beginning with [Warning] are expected:\n")

  # should work
  expect_no_error(cmLaplace$calcLogLik(c(.4, .5, .1)))
  expect_no_error(cmLaplace$calcLaplace(c(.4, .5, .1)))
  expect_no_error(cmLaplace$gr_logLik(c(.4, .5, .1)))
  expect_no_error(cmLaplace$gr_Laplace(c(.4, .5, .1)))

  # should throw errors
  expect_output(expect_error(cmLaplace$calcLogLik(c(.4, .5))), "should be length")
  expect_output(expect_error(cmLaplace$calcLaplace(c(.4, .5))), "should be length")
  expect_output(expect_error(cmLaplace$gr_logLik(c(.4, .5))), "should be length")
  expect_output(expect_error(cmLaplace$gr_Laplace(c(.4, .5))), "should be length")

  # should work
  expect_no_error(cmLaplace$calcLogLik(c(.4, .5), trans = TRUE))
  expect_no_error(cmLaplace$calcLaplace(c(.4, .5), trans = TRUE))
  expect_no_error(cmLaplace$gr_logLik(c(.4, .5), trans = TRUE))
  expect_no_error(cmLaplace$gr_Laplace(c(.4, .5), trans = TRUE))

  # should throw errors
  expect_output(expect_error(cmLaplace$calcLogLik(c(.4, .5, .1), trans = TRUE)), "should be length")
  expect_output(expect_error(cmLaplace$calcLaplace(c(.4, .5, .1), trans = TRUE)), "should be length")
  expect_output(expect_error(cmLaplace$gr_logLik(c(.4, .5, .1), trans = TRUE)), "should be length")
  expect_output(expect_error(cmLaplace$gr_Laplace(c(.4, .5, .1), trans = TRUE)), "should be length")

  ##
  output <- cmLaplace$findMLE(c(.4, .5, .1))
  expect_true(all(output$counts > 0))
  # We couldn't throw an error from a nimbleList-returning method
  # so we emit a message containing "[Warning]".
  expect_output(output <- cmLaplace$findMLE(c(.4, .5)), "should be length")
  expect_identical(output$counts, integer())
})

test_that("Laplace works with different numbers of REs in different cond. ind. sets", {
  # This checks on Issue #1312, which was really a bug with nimOptim
  # that arose from having multiple nimOptim calls share the same
  # control list.
  # This test does not check correctness of result, only that it runs.
  code <- nimbleCode({
    for(i in 1:2) {
      param[i] ~ dnorm(0, 1)
      for(j in 1:num_re[i]) {
        re[i,j] ~ dnorm(param[i], 1)
      }
      y[i] ~ dnorm(sum(re[i,1:num_re[i]]), 1)
    }
  })

  num_re <- c(3,7)   ## different numbers of REs in two conditionally independent sets
  constants <- list(num_re = num_re)
  data <- list(y = c(0,0))

  Rmodel <- nimbleModel(code, constants, data, buildDerivs = TRUE)
  Rlaplace <- buildLaplace(Rmodel, 'param', 're')

  Cmodel <- compileNimble(Rmodel)
  Claplace <- compileNimble(Rlaplace, project = Rmodel)

  expect_no_error(Claplace$findMLE(c(0,0)))
})

test_that("Laplace with N(0,1) random effects works", {
  # This test also uses dflat and dhalfflat
  set.seed(1)
  code <- nimbleCode({
    beta0 ~ dflat()
    beta1 ~ dflat()
    sigma ~ dhalfflat()
    for(i in 1:5) eps[i] ~ dnorm(0, 1)
    for(i in 1:5) sigma_eps[i] <- eps[i] * sigma
    for(i in 1:25) {
      y[i] ~ dpois(exp(beta0 + beta1*X[i] + sigma_eps[group[i]]))
    }
    for(i in 1:10) z[i] ~ dnorm(2*beta0, 1) #calcNodesOther
    foo <- step(beta0)
  })
  X <- rnorm(25)
  group <- rep(1:5, each = 5)
  eps <- rnorm(5, 0, sd = 2)
  y <- rpois(25, exp(3 + .2*X + rep(eps, each=5)))
  z <- rnorm(10, 2*3, sd = 1)
  m <- nimbleModel(code, data = list(y = y, z = z),
                   constants = list(X = X, group=group), buildDerivs=TRUE)

  # Defaults not expected to be useful
  SMN <- setupMargNodes(m)
  expect_identical(SMN$randomEffectsNodes, character())

  SMN <- setupMargNodes(m, #paramNodes = c("beta0", "beta1", "sigma"),
                        randomEffectsNodes = 'eps[1:5]')
  expect_identical(SMN$randomEffectsSets,
                   list('eps[1]','eps[2]','eps[3]','eps[4]','eps[5]'))
  expect_identical(SMN$calcNodesOther,
                   m$expandNodeNames(c('lifted_d2_times_beta0', 'z[1:10]')))
  expect_identical(SMN$paramNodes,
                   c("beta0", "beta1", "sigma"))

  mLaplace <- buildLaplace(m, SMN)
  cm <- compileNimble(m)
  cmLaplace <- compileNimble(mLaplace, project = m)
  cmLaplace$updateSettings(innerOptimMethod="nlminb") # findMLE will hang using BFGS
  res <- cmLaplace$findMLE(c(0,0,1))
  # TMB code in test_N01.cpp
##   #include <TMB.hpp>
## template<class Type>
## Type objective_function<Type>::operator() ()
## {
##   DATA_VECTOR(y);
##   DATA_VECTOR(z);
##   DATA_VECTOR(X);
##   DATA_IVECTOR(group);
##   PARAMETER_VECTOR(eps);
##   PARAMETER_VECTOR(beta);
##   PARAMETER(sigma);
##   int i;
##   // Negative log-likelihood
##   Type ans = Type(0.);
##   for(i = 0; i < 5; ++i)
##     ans -= dnorm(eps[i], Type(0.), Type(1.), true);
##   for(i = 0; i < 25; ++i)
##     ans -= dpois(y[i], exp(beta[0] + beta[1] * X[i] + sigma*eps[group[i]]), true);
##   for(i = 0; i < 10; ++i)
##     ans -= dnorm(z[i], Type(2.)*beta[0], Type(1.), true);
##   return ans;
## }
##   library(TMB)
## compile("test_N01.cpp")
## dyn.load(dynlib("test_N01"))
## data <- list(y = y, X = X, group = group-1, z = z)
## parameters <- list(beta = c(0, 0), sigma = 1, eps = rep(0, 5))
## obj <- MakeADFun(data = data, parameters = parameters, random = "eps", DLL = "test_N01")
## tmbres <- nlminb(obj$par, obj$fn, obj$gr)
## tmbrep <- sdreport(obj, getJointPrecision = TRUE)
  ## tmbvcov <- solve(tmbrep$jointPrecision)
  ##write.table(tmbvcov, file = "", sep=",",col.names = FALSE, row.names=FALSE)
  expect_equal(res$par, c(3.1276930, 0.1645356, 1.5657498), tolerance = 1e-4 )
  summ <- cmLaplace$summary(res, randomEffectsStdError=TRUE, jointCovariance=TRUE)
  ## From the write.table call just above
  ## (which is symmetric anyway, so byrow =TRUE doesn't really matter)
  TMB_vcov <- matrix(byrow = TRUE, nrow = 8, data =
    c(c(0.0153602444576517,0.0117648503870507,0.0284134827252613,0.0199805060648755,0.00318486286937286,-0.0141707177248526,-0.00040366417968837,0.018866970112233),
      c(0.0117648503870507,0.0180821401472876,0.0412180770222714,0.0268113103701797,-0.00119093828503259,-0.013523875123908,-0.000258765680997534,0.0314340905527759),
      c(0.0284134827252614,0.0412180770222714,0.36562970159108,0.167252131855179,-0.0909996094062978,0.00047823545907378,0.000125154168856971,0.28920161604805),
      c(0.0199805060648755,0.0268113103701797,0.167252131855179,0.10843325451055,-0.0439547462832083,-0.00708716995685574,-0.000521588459390236,0.154869927065379),
      c(0.00318486286937281,-0.00119093828503263,-0.0909996094062981,-0.0439547462832083,0.0453386248870613,-0.0201751621932702,-0.000656047342397189,-0.0981200179208084),
      c(-0.0141707177248526,-0.013523875123908,0.000478235459074001,-0.00708716995685565,-0.0201751621932703,0.0245443674575151,-9.27215078179135e-06,0.0144354576429851),
      c(-0.00040366417968837,-0.000258765680997534,0.00012515416885698,-0.000521588459390233,-0.000656047342397191,-9.27215078179097e-06,0.00290834901828464,0.000631332975051338),
      c(0.0188669701122331,0.031434090552776,0.28920161604805,0.154869927065379,-0.0981200179208082,0.0144354576429849,0.000631332975051331,0.283268865188007)))

  expect_equal(summ$vcov, TMB_vcov[c(6:8, 1:5), c(6:8, 1:5)], tol = 1e-4)
  # Check covariance matrix for params only
  summ2 <- cmLaplace$summary(res, originalScale = TRUE, randomEffectsStdError = TRUE, jointCovariance = FALSE)
  expect_equal(summ2$vcov, TMB_vcov[6:8,6:8], tol=1e-4)
})

## Now that innerOptim inits has controls for method and values,
## we need to check over these tests and functionality.

## test_that("Setting Different Initial Values for Inner Optim", {
##   m <- nimbleModel(
##     nimbleCode({
##       y ~ dnorm(0.2 * a, sd = 2)
##       a ~ dnorm(0.5 * mu, sd = 3)
##       mu ~ dnorm(0, sd = 5)
##     }), data = list(y = 4), inits = list(a = -1, mu = 0),
##     buildDerivs = TRUE
##   )

##   mLaplace <- buildLaplace(model = m)
##   cm <- compileNimble(m)
##   cL <- compileNimble(mLaplace, project = m)

##   cL$setInnerOptimWarning(TRUE) ## Print Errors.
##   cL$setInnerCache(FALSE) ## Recalculate inner optim to check starting values.

##   ## Test different starting values:
##   cL$setInnerOptimInits("zero")
##   expect_output(cL$calcLogLik(37), "Warning: optim did not converge for the inner optimization of AGHQuad or Laplace approximation")

##   cL$setInnerOptimInits("last.best")
##   expect_output(cL$calcLogLik(37), NA)  # Small change to actually recalculate. No warning.

##   set.seed(21)
##   cL$setInnerOptimInits("random")
##   expect_output(cL$calcLogLik(37), NA)  # Shouldn't warn.

##   values(cm, "a") <- 0  ## Bad init.
##   cL$setInnerOptimInits("model")
##   expect_output(cL$calcLogLik(37), "Warning: optim did not converge for the inner optimization of AGHQuad or Laplace approximation")

##   values(cm, "a") <- 18
##   cL$setInnerOptimInits("model")
##   expect_output(cL$calcLogLik(37), NA) ## Good init.

##   cL$setInnerOptimInits("last")  ## Last isn't great for this new value.
##   expect_output(cL$calcLogLik(15), "Warning: optim did not converge for the inner optimization of AGHQuad or Laplace approximation")

##   ## Inspect model to see if values are updated properly after running:
##   cL$setInnerCache(FALSE)
##   cL$setInnerOptimInits("random")
##   cL$calcLogLik(15)
##   old.val <- cm$a
##   cL$setModelValues(15)
##   new.val <- cm$a
##   expect_false(old.val == new.val)
## })

test_that("regression test that replacing inner/outer control list works (given NCT issue 571 bug with nimOptimDefaultControl", {
    ## Use N=4 to reduce compilation time for separate Laplaces.
    pumpCode <- nimbleCode({ 
        for (i in 1:N){
            theta[i] ~ dgamma(alpha, beta)
            lambda[i] <- theta[i] * t[i]
            x[i] ~ dpois(lambda[i])
        }
        alpha ~ dexp(1.0)
        beta ~ dgamma(0.1, 1.0)
    })
    pumpConsts <- list(N = 10, t = c(94.3, 15.7, 62.9, 126, 5.24, 31.4, 1.05, 1.05, 2.1, 10.5))
    pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22))
    pumpInits <- list(alpha = 0.1, beta = 0.1, theta = rep(0.1, pumpConsts$N))
    pump <- nimbleModel(code = pumpCode, name = "pump", constants = pumpConsts, 
                        data = pumpData, inits = pumpInits, buildDerivs = TRUE)
    
    pumpLaplace <- buildLaplace(pump)
    Cpump <- compileNimble(pump)
    CpumpLaplace <- compileNimble(pumpLaplace, project = pump)

    MLEres0 <- CpumpLaplace$findMLE()
    CpumpLaplace$updateSettings(replace_innerOptimControl=TRUE)
    Cpump$theta <- rep(0.1, 10)
    Cpump$alpha <- Cpump$beta <- 0.1
    MLEres1 <- CpumpLaplace$findMLE()
    CpumpLaplace$updateSettings(replace_outerOptimControl=TRUE)
    Cpump$theta <- rep(0.1, 10)
    Cpump$alpha <- Cpump$beta <- 0.1
    MLEres2 <- CpumpLaplace$findMLE()
    expect_equal(MLEres0$par, MLEres1$par)
    expect_equal(MLEres0$par, MLEres2$par)
})

nimbleOptions(enableDerivs = EDopt)
nimbleOptions(buildModelDerivs = BMDopt)
