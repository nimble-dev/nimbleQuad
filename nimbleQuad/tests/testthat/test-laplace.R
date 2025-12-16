# Tests of Laplace approximation
EDopt <- nimbleOptions("enableDerivs")
BMDopt <- nimbleOptions("buildModelDerivs")
nimbleOptions(enableDerivs = TRUE)
nimbleOptions(buildModelDerivs = TRUE)
nimbleOptions(allowDynamicIndexing = FALSE)


test_that("Laplace simplest 1D works", {
  m <- nimbleModel(
    nimbleCode({
      y ~ dnorm(a, sd = 2)
      a ~ dnorm(mu, sd = 3)
      mu ~ dnorm(0, sd = 5)
    }), data = list(y = 4), inits = list(a = -1, mu = 0),
    buildDerivs = TRUE
  )

  mLaplace <- buildLaplace(model = m)
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit

  opt <- cmLaplace$findMLE()
  expect_equal(opt$par, 4, tol = 1e-6) # tolerance was reduced in this test when we switched to nlminb
  # V[a] = 9
  # V[y] = 9 + 4 = 13
  # Cov[a, y] = V[a] = 9 (not needed)
  # y ~ N(mu, 13)
  expect_equal(opt$value, dnorm(4, 4, sd = sqrt(13), log = TRUE))
  # muhat = y = 4
  # ahat = (9*y+4*mu)/(9+4) = y = 4
  # Jacobian of ahat wrt mu is 4/13
  # Hessian of joint loglik wrt a: -(1/4 + 1/9)
  # Hessian of marginal loglik wrt mu is -1/13
  summ <- cmLaplace$summary(opt, originalScale = TRUE, randomEffectsStdError = TRUE, jointCovariance = TRUE)
  expect_equal(summ$randomEffects$estimate, 4, tol = 1e-5)
  # check behavior of summaryLaplace
  summ2 <- summaryLaplace(cmLaplace, opt, randomEffectsStdError = TRUE, jointCovariance = TRUE)
  expect_equal(nrow(summ2$randomEffects), 1)
  expect_equal(nrow(summ2$params), 1)
  expect_equal(row.names(summ2$randomEffects), "a")
  expect_equal(row.names(summ2$params), "mu")
  # Covariance matrix
  vcov <- matrix(c(0, 0, 0, c(1/(1/4+1/9))), nrow = 2) + matrix(c(1, 4/13), ncol = 1) %*% (13) %*% t(matrix(c(1, 4/13), ncol = 1))
  expect_equal(vcov, summ$vcov, tol = 1e-6)
  # Check covariance matrix for params only
  summ3 <- cmLaplace$summary(opt, originalScale = TRUE, randomEffectsStdError = TRUE, jointCovariance = FALSE)
  expect_equal(summ3$vcov, vcov[1,1,drop=FALSE], tol=1e-6)

  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$findMLE()
  expect_equal(opt$par, optNoSplit$par, tol = 1e-2)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)
})

test_that("Laplace simplest 1D with a constrained parameter works", {
  m <- nimbleModel(
    nimbleCode({
      y ~ dnorm(a, sd = 2)
      a ~ dnorm(mu, sd = 3)
      mu ~ dexp(1.0)
    }), data = list(y = 4), inits = list(a = -1, mu = 0),
    buildDerivs = TRUE
  )

  mLaplace <- buildLaplace(model = m)
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit

  opt <- cmLaplace$findMLE()
  # V[a] = 9
  # V[y] = 9 + 4 = 13
  # Cov[a, y] = V[a] = 9 (not needed)
  # y ~ N(mu, 13)
  # muhat = y = 4
  # ahat = (9*y+4*mu)/(9+4) = y = 4
  # Jacobian of ahat wrt transformed param log(mu) is 4/13*mu = 4*mu/13 = 16/13
  # Hessian of joint loglik wrt a: -(1/4 + 1/9)
  # Hessian of marginal loglik wrt transformed param log(mu) is (y*mu - 2*mu*mu)/13 = -4^2/13
  # Variance of transformed param is 13/16
  expect_equal(opt$par, 4, tol = 1e-4)
  expect_equal(opt$value, dnorm(4, 4, sd = sqrt(13), log = TRUE))
  expect_equal(opt$hessian[1,1], -4^2/13, tol = 1e-4)
  summ <- cmLaplace$summary(opt, originalScale = FALSE, randomEffectsStdError = TRUE, jointCovariance = TRUE)
  expect_equal(summ$randomEffects$estimate, 4, tol = 1e-4)
  expect_equal(summ$params$estimate, log(4), tol = 1e-4)
  # check summaryLaplace
  summL <- summaryLaplace(cmLaplace, opt, originalScale = FALSE, randomEffectsStdError = TRUE, jointCovariance = TRUE)
  expect_equal(summL$params['param_trans_1','estimate'], log(4), tol = 1e-4)

  # Covariance matrix on transformed scale
  vcov_transform <- matrix(c(0, 0, 0, 1/(1/4+1/9)), nrow = 2) + matrix(c(1, 16/13), ncol = 1) %*% (13/16) %*% t(matrix(c(1, 16/13), ncol = 1))
  expect_equal(vcov_transform, summ$vcov, tol = 1e-4)
  summ2 <- cmLaplace$summary(opt, originalScale = TRUE, randomEffectsStdError = TRUE, jointCovariance = TRUE)
  # Covariance matrix on original scale
  vcov <- diag(c(4, 1)) %*% vcov_transform %*% diag(c(4, 1))
  expect_equal(vcov, summ2$vcov, tol = 1e-4)
  # Check covariance matrix for params only
  tryResult <- try({
  summ3 <- cmLaplace$summary(opt, originalScale = TRUE, randomEffectsStdError = TRUE, jointCovariance = FALSE);
  expect_equal(summ3$vcov, vcov[1,1,drop=FALSE], tol=1e-5)
  summ4 <- cmLaplace$summary(opt, originalScale = FALSE, randomEffectsStdError = TRUE, jointCovariance = FALSE)
  expect_equal(summ4$vcov, vcov_transform[1,1,drop=FALSE], tol=5e-5)
  })
  if(inherits(tryResult, "try-error")) {
      print(class(cmLaplace))
      print(cmLaplace)
  }

  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$findMLE()
  expect_equal(opt$par, optNoSplit$par, tol = 1e-2)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)
})

test_that("Laplace simplest 1D (constrained) with multiple data works", {
  set.seed(1)
  m <- nimbleModel(
    nimbleCode({
      mu ~ dnorm(0, sd = 5)
      a ~ dexp(rate = exp(mu))
      for (i in 1:5){
        y[i] ~ dnorm(a, sd = 2)
      }
    }), data = list(y = rnorm(5, 1, 2)), inits = list(mu = 2, a = 1),
    buildDerivs = TRUE
  )
  mLaplace <- buildLaplace(model = m)
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit

  opt <- cmLaplace$findMLE()
  summ <- cmLaplace$summary(opt, originalScale = FALSE, jointCovariance = TRUE)
  # Results are checked using those from TMB
  # TMB cpp code:
  #include <TMB.hpp>
  #template<class Type>
  #Type objective_function<Type>::operator() ()
  # {
  #   DATA_VECTOR(y);
  #   PARAMETER(mu);
  #   PARAMETER(log_a);
  #   int n = y.size();
  #   Type a = exp(log_a); // Invserse transformation
  #   // Negative log-likelihood
  #   Type ans = -dexp(a, exp(mu), true);
  #   ans -= log_a; // logdet Jacobian of inverse transformation: exp
  #   for(int i = 0; i < n; i++){
  #     ans -= dnorm(y[i], a, Type(2), true);
  #   }
  #   return ans;
  # }
  # TMB R code:
  # library(TMB)
  # compile("test.cpp")
  # dyn.load(dynlib("test"))
  # data <- list(y = m$y)
  # parameters <- list(mu = 2, log_a = 0)
  #
  # ## Fit model
  # obj <- MakeADFun(data, parameters, random="log_a", DLL="test")
  # tmbres <- nlminb(obj$par, obj$fn, obj$gr)
  # tmbrep <- sdreport(obj, getJointPrecision = TRUE)
  # tmbvcov <- inverse(tmbrep$jointPrecision)
  expect_equal(opt$par, 0.2895238, tol = 1e-3)
  expect_equal(opt$value, -10.47905, tol = 1e-7)
  expect_equal(summ$randomEffects$estimate, -0.005608619, tol = 1e-3)
  vcov <- matrix(c(2.741033, -1.628299, -1.628299, 1.414499), nrow = 2, byrow = TRUE)
  expect_equal(summ$vcov, vcov, 2e-3)
  # Check covariance matrix for params only
  summ2 <- cmLaplace$summary(opt, originalScale = TRUE, randomEffectsStdError = TRUE, jointCovariance = FALSE)
  expect_equal(summ2$vcov, vcov[1,1,drop=FALSE], tol=1e-3)


  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$findMLE()
  expect_equal(opt$par, optNoSplit$par, tol = 1e-2)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)
})

test_that("Laplace simplest 1D (constrained) with deterministic intermediates and multiple data works", {
  set.seed(1)
  m <- nimbleModel(
    nimbleCode({
      mu ~ dnorm(0, sd = 5)
      a ~ dexp(rate = exp(0.5 * mu))
      for (i in 1:5){
        y[i] ~ dnorm(0.2 * a, sd = 2)
      }
    }), data = list(y = rnorm(5, 1, 2)), inits = list(mu = 2, a = 1),
    buildDerivs = TRUE
  )
  mLaplace <- buildLaplace(model = m)
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit

  opt <- cmLaplace$findMLE()
  summ <- cmLaplace$summary(opt, originalScale = FALSE, jointCovariance = TRUE)
  # Results are checked using those from TMB
  # TMB cpp code:
  # #include <TMB.hpp>
  # template<class Type>
  # Type objective_function<Type>::operator() ()
  # {
  #   DATA_VECTOR(y);
  #   PARAMETER(mu);
  #   PARAMETER(log_a);
  #   int n = y.size();
  #   Type a = exp(log_a); // Invserse transformation
  #   // Negative log-likelihood
  #   Type ans = -dexp(a, exp(0.5 * mu), true);
  #   ans -= log_a; // logdet Jacobian of inverse transformation: exp
  #   for(int i = 0; i < n; i++){
  #     ans -= dnorm(y[i], 0.2 * a, Type(2), true);
  #   }
  #   ADREPORT(a);
  #   return ans;
  # }
  ## R code:
  # library(TMB)
  # compile("test.cpp")
  # dyn.load(dynlib("test"))
  # data <- list(y = m$y)
  # parameters <- list(mu = 2, log_a = 0)
  #
  # ## Fit model
  # obj <- MakeADFun(data, parameters, random="log_a", DLL="test")
  # tmbres <- nlminb(obj$par, obj$fn, obj$gr)
  # tmbrep <- sdreport(obj, getJointPrecision = TRUE)
  # tmbvcov <- inverse(tmbrep$jointPrecision)
  expect_equal(opt$par, -2.639534, 2e-4)
  expect_equal(opt$value, -10.47905, tol = 1e-5)
  expect_equal(summ$randomEffects$estimate, 1.603742, tol = 1e-4)
  vcov <- matrix(c(10.967784, -3.258191, -3.258191, 1.415167), nrow = 2, byrow = TRUE)
  expect_equal(summ$vcov, vcov, 2e-3)
  # Check covariance matrix for params only
  summ2 <- cmLaplace$summary(opt, originalScale = TRUE, randomEffectsStdError = TRUE, jointCovariance = FALSE)
  expect_equal(summ2$vcov, vcov[1,1,drop=FALSE], tol=2e-3)

  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$findMLE()
  expect_equal(opt$par, optNoSplit$par, tol = 1e-2)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-4)
})

test_that("Laplace 1D with deterministic intermediates works", {
  # Note this test has some slop. In old versions there were inner optimization
  # warnings issued for this case. So we turned on warnings and checked for them.
  # Now the warnings aren't issued. As a result, we really need a new test to
  # check that warnings are correctly emitted.
  m <- nimbleModel(
    nimbleCode({
      y ~ dnorm(0.2 * a, sd = 2)
      a ~ dnorm(0.5 * mu, sd = 3)
      mu ~ dnorm(0, sd = 5)
    }), data = list(y = 4), inits = list(a = -1, mu = 0),
    buildDerivs = TRUE
  )

  mLaplace <- buildLaplace(model = m)
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit

  #expect_output(
    opt <- cmLaplace$findMLE()
  #, "Warning: inner optimzation had a non-zero convergence code\\. Use checkInnerConvergence\\(TRUE\\) to see details\\.")
  expect_equal(opt$par, 40, tol = 1e-4) # 40 = 4 * (1/.2) * (1/.5)
  # V[a] = 9
  # V[y] = 0.2^2 * 9 + 4 = 4.36
  expect_equal(opt$value, dnorm(0.1*40, 0.1*40, sd = sqrt(4.36), log = TRUE))
  # y ~ N(0.2*0.5*mu, 4.36)
  # muhat = y/(0.2*0.5) = 40
  # ahat = (9*0.2*y + 4*0.5*mu)/(4+9*0.2^2) = 20
  # Jacobian of ahat wrt mu is 4*0.5/(4+9*0.2^2) = 0.4587156
  # Hessian of joint loglik wrt a: -(0.2^2/4 + 1/9)
  # Hessian of marginal loglik wrt param mu is -(0.2*0.5)^2/4.36 = -0.002293578
  cmLaplace$updateSettings(innerOptimWarning=TRUE)
  #expect_output(
    summ <- cmLaplace$summary(opt, originalScale = TRUE, randomEffectsStdError = TRUE,
                              jointCovariance = TRUE)
  #,  "optim did not converge for the inner optimization")
  expect_equal(summ$randomEffects$estimate, 20, tol = 1e-4)
  # Covariance matrix
  vcov <- matrix(c(0, 0, 0, 1/(0.2^2/4+1/9)), nrow = 2) + matrix(c(1, 0.4587156), ncol = 1) %*% (1/0.002293578) %*% t(matrix(c(1, 0.4587156), ncol = 1))
  expect_equal(vcov, summ$vcov, tol = 1e-4)
  # Check covariance matrix for params only
  #expect_output(
    summ2 <- cmLaplace$summary(opt, originalScale = TRUE, randomEffectsStdError = TRUE, jointCovariance = FALSE)
   #,            "does not converge")
  expect_equal(summ2$vcov, vcov[1,1,drop=FALSE], tol=1e-4)

  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  #expect_output(
    optNoSplit <- cmLaplaceNoSplit$findMLE()
  #, "Warning: inner optimzation had a non-zero convergence code\\. Use checkInnerConvergence\\(TRUE\\) to see details\\.")
  expect_equal(opt$par, optNoSplit$par, tol = 1e-2)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)
})

test_that("Laplace 1D with a constrained parameter and deterministic intermediates works", {
  ## Again (see above), the innerOptimWarning and expect_output are
  ## defunct portions of this test.
  m <- nimbleModel(
    nimbleCode({
      y ~ dnorm(0.2 * a, sd = 2)
      a ~ dnorm(0.5 * mu, sd = 3)
      mu ~ dexp(1.0)
    }), data = list(y = 4), inits = list(a = -1, mu = 0),
    buildDerivs = TRUE
  )

  mLaplace <- buildLaplace(model = m)
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit

  #expect_output(
    opt <- cmLaplace$findMLE()
  #, "Warning: inner optimzation had a non-zero convergence code\\. Use checkInnerConvergence\\(TRUE\\) to see details\\.")

  # V[a] = 9
  # V[y] = 0.2^2 * 9 + 4 = 4.36
  # y ~ N(0.2*0.5*mu, 4.36)
  # muhat = y/(0.2*0.5) = 40
  # ahat = (9*0.2*y + 4*0.5*mu)/(4+9*0.2^2) = 20
  # Jacobian of ahat wrt transformed param log(mu) is 4*0.5*mu/(4+9*0.2^2) = 18.34862
  # Hessian of joint loglik wrt a: -(0.2^2/4 + 1/9)
  # Hessian of marginal loglik wrt param mu is -(0.2*0.5)^2/4.36
  # Hessian of marginal loglik wrt transformed param log(mu) is (0.2*0.5*y*mu - 2*0.1^2*mu*mu)/4.36 = -3.669725
  expect_equal(opt$par, 40, tol = 1e-4)
  expect_equal(opt$hessian[1,1], -3.669725, tol = 1e-3)
  expect_equal(opt$value, dnorm(0.1*40, 0.1*40, sd = sqrt(4.36), log = TRUE))

  cmLaplace$updateSettings(innerOptimWarning=TRUE, useInnerCache=FALSE)
  #cmLaplace$setInnerOptimWarning(TRUE)
  #cmLaplace$setInnerCache(FALSE)
  #expect_output(
    summ <- cmLaplace$summary(opt, originalScale = FALSE, randomEffectsStdError = TRUE,
                              jointCovariance = TRUE)
   #,            "optim did not converge for the inner optimization")
  expect_equal(summ$randomEffects$estimate, 20, tol = 1e-4)
  # Covariance matrix on transformed scale
  vcov_transform <- matrix(c(0, 0, 0, 1/(0.2^2/4+1/9)), nrow = 2) + matrix(c(1, 18.34862), ncol = 1) %*% (1/3.669725) %*% t(matrix(c(1, 18.34862), ncol = 1))
  expect_equal(vcov_transform, summ$vcov, tol = 1e-3)
  #expect_output(
    summ2 <- cmLaplace$summary(opt, originalScale = TRUE, randomEffectsStdError = TRUE,
                               jointCovariance = TRUE)
  #, "optim did not converge for the inner optimization")
  # Covariance matrix on original scale
  vcov <- diag(c(40,1)) %*% vcov_transform %*% diag(c(40,1))
  expect_equal(vcov, summ2$vcov, tol = 1e-3)

  # Check summary based on not recomputing random effects and hessian.
  summ.orig <- cmLaplace$summary(opt, originalScale = TRUE, randomEffectsStdError = TRUE, jointCovariance = FALSE)
  # Check covariance matrix for params only
  #expect_output(
    summ3 <- cmLaplace$summary(opt, originalScale = TRUE, randomEffectsStdError = TRUE, jointCovariance = FALSE)
  #,  "optim did not converge for the inner optimization")

  # Make sure that the recompute didn't change the values of the random effects or standard errors.
  expect_equal(summ.orig$randomEffects[["estimate"]], summ3$randomEffects[["estimate"]], tol = 1e-12)
  expect_equal(summ.orig$randomEffects[["stdError"]], summ3$randomEffects[["stdError"]], tol = 1e-12)

  expect_equal(summ3$vcov, vcov[1,1,drop=FALSE], tol=1e-3)
  #expect_output(
    summ4 <- cmLaplace$summary(opt, originalScale = FALSE, randomEffectsStdError = TRUE, jointCovariance = FALSE)
   #,             "optim did not converge for the inner optimization")
  expect_equal(summ4$vcov, vcov_transform[1,1,drop=FALSE], tol=1e-4)

  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  cmLaplace$updateSettings(innerOptimWarning=TRUE)
#  cmLaplaceNoSplit$setInnerOptimWarning(TRUE)
  #expect_output(
  optNoSplit <- cmLaplaceNoSplit$findMLE()
  #,  "optim did not converge for the inner optimization")
  expect_equal(opt$par, optNoSplit$par, tol = 1e-2)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)
})

test_that("Laplace 1D with deterministic intermediates and multiple data works", {
  m <- nimbleModel(
    nimbleCode({
      for(i in 1:n)
        y[i] ~ dnorm(mu_y, sd = 2) # larger multiplier to amplify cov terms in result below
      mu_y <- 0.8*a
      a ~ dnorm(mu_a, sd = 3)
      mu_a <- 0.5 * mu
      mu ~ dnorm(0, sd = 5)
    }),
    data = list(y = c(4, 5, 6)),
    constants = list(n = 3),
    inits = list(a = -1, mu = 0),
    buildDerivs = TRUE
  )
  mLaplace <- buildLaplace(model = m)
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit

  opt <- cmLaplace$findMLE()
  expect_equal(opt$par, 12.5, tol = 1e-4) # 12.5 = mean(y) * (1/.8) * (1/.5) where mean(y) = 5
  # V[a] = 9
  # V[y[i]] = 0.8^2 * 9 + 4 = 9.76
  # Cov[a, y[i]] = 0.8 * 9 = 7.2
  # Cov[y[i], y[j]] = 0.8^2 * 9 = 5.76
  Cov_ay1y2y3 <- matrix(nrow = 4, ncol = 4)
  Cov_ay1y2y3[1, 1:4] <- c(9, 7.2, 7.2, 7.2)
  Cov_ay1y2y3[2, 1:4] <- c(7.2, 9.76, 5.76, 5.76)
  Cov_ay1y2y3[3, 1:4] <- c(7.2, 5.76, 9.76, 5.76)
  Cov_ay1y2y3[4, 1:4] <- c(7.2, 5.76, 5.76, 9.76)
  Cov_y1y2y3 <- Cov_ay1y2y3[2:4, 2:4]
  chol_cov <- chol(Cov_y1y2y3)
  res <- dmnorm_chol(c(4, 5, 6), 0.8*0.5*12.5, cholesky = chol_cov, prec_param=FALSE, log = TRUE)
  expect_equal(opt$value, res)
  # y[i] ~ N(0.4*mu, 9.76)
  # mean(y) = 5
  # muhat = mean(y)/(0.8*0.5) = 12.5
  # ahat = (9*0.8*sum(y) + 4*0.5*mu)/(4+9*0.8^2*3) = 6.25
  # Jacobian of ahat wrt mu is 4*0.5/(4+9*0.8^2*3) = 0.09398496
  # Hessian of joint loglik wrt a: -(3*0.8^2/4 + 1/9)
  # Hessian of marginal loglik wrt mu: -0.02255639 (numerical, have not got AD work)
  summ <- cmLaplace$summary(opt, originalScale = FALSE, randomEffectsStdError = TRUE, jointCovariance = TRUE)
  expect_equal(summ$randomEffects$estimate, 6.25, tol = 1e-6)
  # Covariance matrix
  vcov <- matrix(c(0, 0, 0, 1/(0.8^2*3/4+1/9)), nrow = 2) + matrix(c(1, 0.09398496), ncol = 1) %*% (1/0.02255639) %*% t(matrix(c(1, 0.09398496), ncol = 1))
  expect_equal(vcov, summ$vcov, tol = 1e-7)
  # Check covariance matrix for params only
  summ2 <- cmLaplace$summary(opt, originalScale = TRUE, randomEffectsStdError = TRUE, jointCovariance = FALSE)
  expect_equal(summ2$vcov, vcov[1,1,drop=FALSE], tol=1e-6)

  # check that
  mLaplaceCheck <- buildLaplace(model = m, paramNodes = 'mu', randomEffectsNodes = 'a')
  nim1D <-  mLaplace$AGHQuad_nfl[[1]]
  expect_identical(nim1D$paramNodes, "mu")
  expect_identical(nim1D$paramDeps, "mu_a")
  expect_identical(nim1D$randomEffectsNodes, "a")
  expect_identical(nim1D$innerCalcNodes, c("a", "mu_y", "y[1]", "y[2]", "y[3]"))
  expect_identical(nim1D$calcNodes, c("mu_a", nim1D$innerCalcNodes))
  expect_identical(nim1D$inner_updateNodes, "mu_a")
  expect_identical(nim1D$inner_constantNodes, c("y[1]", "y[2]", "y[3]"))
  expect_identical(nim1D$joint_updateNodes, character())
  expect_identical(nim1D$joint_constantNodes, c("y[1]", "y[2]", "y[3]"))

  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$findMLE() # some warnings are ok here
  expect_equal(opt$par, optNoSplit$par, tol = 1e-2)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)
})

test_that("Laplace 1D with a constrained parameter and deterministic intermediates and multiple data works", {
  m <- nimbleModel(
    nimbleCode({
      for(i in 1:n)
        y[i] ~ dnorm(mu_y, sd = 2)
      mu_y <- 0.8*a
      a ~ dnorm(mu_a, sd = 3)
      mu_a <- 0.5 * mu
      mu ~ dexp(1.0)
    }),
    data = list(y = c(4, 5, 6)),
    constants = list(n = 3),
    inits = list(a = -1, mu = 0),
    buildDerivs = TRUE
  )
  mLaplace <- buildLaplace(model = m)
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit

  opt <- cmLaplace$findMLE()
  expect_equal(opt$par, 12.5, tol = 1e-4)
  # V[a] = 9
  # V[y[i]] = 0.8^2 * 9 + 4 = 9.76
  # Cov[a, y[i]] = 0.8 * 9 = 7.2
  # Cov[y[i], y[j]] = 0.8^2 * 9 = 5.76
  # y[i] ~ N(0.4*mu, 9.76)
  # mean(y) = 5
  # muhat = mean(y)/(0.8*0.5) = 12.5
  # ahat = (9*0.8*sum(y) + 4*0.5*mu)/(4+9*0.8^2*3) = 6.25
  # Jacobian of ahat wrt transformed param log(mu) is 4*0.5*mu/(4+9*0.8^2*3) = 1.174812
  # Hessian of joint loglik wrt a: -(3*0.8^2/4 + 1/9)
  # Hessian of marginal loglik wrt transformed param: -3.524436 (numerical, have not got AD work)
  Cov_ay1y2y3 <- matrix(nrow = 4, ncol = 4)
  Cov_ay1y2y3[1, 1:4] <- c(9, 7.2, 7.2, 7.2)
  Cov_ay1y2y3[2, 1:4] <- c(7.2, 9.76, 5.76, 5.76)
  Cov_ay1y2y3[3, 1:4] <- c(7.2, 5.76, 9.76, 5.76)
  Cov_ay1y2y3[4, 1:4] <- c(7.2, 5.76, 5.76, 9.76)
  Cov_y1y2y3 <- Cov_ay1y2y3[2:4, 2:4]
  chol_cov <- chol(Cov_y1y2y3)
  res <- dmnorm_chol(c(4, 5, 6), 0.8*0.5*12.5, cholesky = chol_cov, prec_param=FALSE, log = TRUE)
  expect_equal(opt$value, res)
  # Check covariance matrix
  summ <- cmLaplace$summary(opt, originalScale = FALSE, randomEffectsStdError = TRUE, jointCovariance = TRUE)
  expect_equal(summ$randomEffects$estimate, 6.25, tol = 1e-6)
  # Covariance matrix on transformed scale
  vcov_transform <- matrix(c(0, 0, 0, 1/(0.8^2*3/4+1/9)), nrow = 2) + matrix(c(1, 1.174812), ncol = 1) %*% (1/3.524436) %*% t(matrix(c(1, 1.174812), ncol = 1))
  expect_equal(vcov_transform, summ$vcov, tol = 1e-6)
  summ2 <- cmLaplace$summary(opt, originalScale = TRUE, randomEffectsStdError = TRUE, jointCovariance = TRUE)
  # Covariance matrix on original scale
  vcov <- diag(c(12.5, 1)) %*% vcov_transform %*% diag(c(12.5, 1))
  expect_equal(vcov, summ2$vcov, tol = 1e-5)
  # Check covariance matrix for params only
  summ3 <- cmLaplace$summary(opt, originalScale = TRUE, randomEffectsStdError = TRUE, jointCovariance = FALSE)
  expect_equal(summ3$vcov, vcov[1,1,drop=FALSE], tol=1e-5)
  summ4 <- cmLaplace$summary(opt, originalScale = FALSE, randomEffectsStdError = TRUE, jointCovariance = FALSE)
  expect_equal(summ4$vcov, vcov_transform[1,1,drop=FALSE], tol=1e-6)

  # check that
  mLaplaceCheck <- buildLaplace(model = m, paramNodes = 'mu', randomEffectsNodes = 'a')
  nim1D <-  mLaplace$AGHQuad_nfl[[1]]
  expect_identical(nim1D$paramNodes, "mu")
  expect_identical(nim1D$paramDeps, "mu_a")
  expect_identical(nim1D$randomEffectsNodes, "a")
  expect_identical(nim1D$innerCalcNodes, c("a", "mu_y", "y[1]", "y[2]", "y[3]"))
  expect_identical(nim1D$calcNodes, c("mu_a", nim1D$innerCalcNodes))
  expect_identical(nim1D$inner_updateNodes, "mu_a")
  expect_identical(nim1D$inner_constantNodes, c("y[1]", "y[2]", "y[3]"))
  expect_identical(nim1D$joint_updateNodes, character())
  expect_identical(nim1D$joint_constantNodes, c("y[1]", "y[2]", "y[3]"))

  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$findMLE() # some warnings are ok here
  expect_equal(opt$par, optNoSplit$par, tol = 1e-2)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)
})



nimbleOptions(enableDerivs = EDopt)
nimbleOptions(buildModelDerivs = BMDopt)
