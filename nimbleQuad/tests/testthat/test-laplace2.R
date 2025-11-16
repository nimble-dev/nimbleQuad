library(nimbleQuad) # to get nimbleQuad's names first
# Tests of Laplace approximation
source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))
source(system.file(file.path('tests', 'testthat', 'AD_test_utils.R'), package = 'nimble'))
EDopt <- nimbleOptions("enableDerivs")
BMDopt <- nimbleOptions("buildModelDerivs")
nimbleOptions(enableDerivs = TRUE)
nimbleOptions(buildModelDerivs = TRUE)
nimbleOptions(allowDynamicIndexing = FALSE)

test_that("Laplace simplest 2x1D works, with multiple data for each", {
  set.seed(1)
  y <- matrix(rnorm(6, 4, 5), nrow = 2)
  m <- nimbleModel(
    nimbleCode({
      for(i in 1:2) {
        mu_y[i] <- 0.8*a[i]
        for(j in 1:3)
          y[i, j] ~ dnorm(mu_y[i], sd = 2)
        a[i] ~ dnorm(mu_a, sd = 3)
      }
      mu_a <- 0.5 * mu
      mu ~ dnorm(0, sd = 5)
    }), data = list(y = y), inits = list(a = c(-2, -1), mu = 0),
    buildDerivs = TRUE
  )

  mLaplace <- buildLaplace(model = m)
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit

  opt <- cmLaplace$findMLE()
  expect_equal(opt$par, mean(y)/(0.8*0.5), tol = 1e-4) # optim's reltol is about 1e-8 but that is for the value, not param.
  # V[a] = 9
  # V[y[i]] = 0.8^2 * 9 + 4 = 9.76
  # Cov[a, y[i]] = 0.8 * 9 = 7.2
  # Cov[y[i], y[j]] = 0.8^2 * 9 = 5.76, within a group
  Cov_A_Y <- matrix(nrow = 8, ncol = 8)
  Cov_A_Y[1, 1:8] <- c(  9,    0,  7.2,  7.2,  7.2,    0,    0,    0)
  Cov_A_Y[2, 1:8] <- c(  0,    9,    0,    0,    0,  7.2,  7.2,  7.2)
  Cov_A_Y[3, 1:8] <- c(7.2,    0, 9.76, 5.76, 5.76,    0,    0,    0)
  Cov_A_Y[4, 1:8] <- c(7.2,    0, 5.76, 9.76, 5.76,    0,    0,    0)
  Cov_A_Y[5, 1:8] <- c(7.2,    0, 5.76, 5.76, 9.76,    0,    0,    0)
  Cov_A_Y[6, 1:8] <- c(  0,  7.2,    0,    0,    0, 9.76, 5.76, 5.76)
  Cov_A_Y[7, 1:8] <- c(  0,  7.2,    0,    0,    0, 5.76, 9.76, 5.76)
  Cov_A_Y[8, 1:8] <- c(  0,  7.2,    0,    0,    0, 5.76, 5.76, 9.76)
  Cov_Y <- Cov_A_Y[3:8, 3:8]
  chol_cov <- chol(Cov_Y)
  res <- dmnorm_chol(as.numeric(t(y)), mean(y), cholesky = chol_cov, prec_param=FALSE, log = TRUE)
  expect_equal(opt$value, res)
  # muhat = mean(y)/(0.8*0.5)
  # ahat[1] = (9*0.8*sum(y[1,]) + 4*0.5*mu)/(4+9*0.8^2*3)
  # ahat[2] = (9*0.8*sum(y[2,]) + 4*0.5*mu)/(4+9*0.8^2*3)
  # Jacobian of ahat[i] wrt mu is 4*0.5/(4+9*0.8^2*3) = 0.09398496
  # Hessian of joint loglik wrt a[i]a[i]: -(3*0.8^2/4 + 1/9); wrt a[i]a[j]: 0
  # Hessian of marginal loglik wrt mu: -0.04511278 (numerical, have not got AD work)
  muhat <- mean(y)/(0.8*0.5)
  ahat <- c((9*0.8*sum(y[1,]) + 4*0.5*muhat)/(4+9*0.8^2*3), (9*0.8*sum(y[2,]) + 4*0.5*muhat)/(4+9*0.8^2*3))
  summ <- cmLaplace$summary(opt, originalScale = TRUE, randomEffectsStdError = TRUE, jointCovariance = TRUE)
  expect_equal(summ$randomEffects$estimate, ahat, tol = 1e-6)
  # Covariance matrix
  vcov <- diag(c(0, rep(1/(3*0.8^2/4 + 1/9), 2))) + matrix(c(1, rep(0.09398496, 2)), ncol = 1) %*% (1/0.04511278) %*% t(matrix(c(1, rep(0.09398496, 2)), ncol = 1))
  if(Sys.info()['sysname'] != "Windows")  # Issue 66
      expect_equal(vcov, summ$vcov, tol = 1e-7)
      
  ## Check covariance matrix for params only
  tryResult <- try({
      summ2 <- cmLaplace$summary(opt, originalScale = TRUE, randomEffectsStdError = TRUE, jointCovariance = FALSE)
      expect_equal(summ2$vcov, vcov[1,1,drop=FALSE], tol=1e-6)
  })
  if(inherits(tryResult, 'try-error')) {
      print(class(cmLaplace))
      print(cL)
  }


  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$findMLE() # some warnings are ok here
  expect_equal(opt$par, optNoSplit$par, tol = 1e-2)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)
})

test_that("Laplace with 2x1D random effects needing joint integration works, without intermediate nodes", {
  set.seed(1)
  y <- matrix(rnorm(6, 4, 5), nrow = 2)
  m <- nimbleModel(
    nimbleCode({
      for(i in 1:2) {
        a[i] ~ dnorm(mu_a, sd = 3)
      }
      for(j in 1:3) # Note this is different than above.
        # These are 3 observations each of 2D
        y[1:2, j] ~ dmnorm(a[1:2], cov = cov_y[1:2, 1:2])
      mu_a <- 0.5 * mu
      mu ~ dnorm(0, sd = 5)
    }),
    data = list(y = y),
    inits = list(a = c(-2, -1), mu = 0),
    constants = list(cov_y = matrix(c(2, 1.5, 1.5, 2), nrow = 2)),
    buildDerivs = TRUE
  )

  mLaplace <- buildLaplace(model = m)
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit

  opt <- cmLaplace$findMLE()
  expect_equal(opt$par, mean(y)/(0.5), tol = 1e-4) # optim's reltol is about 1e-8 but that is for the value, not param.
  # V[a] = 9
  # V[y[1:2, i]] =  diag(2)*9 + cov_y = [ (9 + 2), 0 + 1.5; 0+1.5, (9+2)]
  # Cov[a[1], y[1,j]] = 9
  # Cov[a[1], y[2,j]] = 0
  # Cov[a[2], y[1,j]] = 0
  # Cov]a[2], y[2,j]] = 9
  # Cov[y[1,i], y[1,j]] = 9
  # Cov[y[2,i], y[2,j]] = 9
  Cov_A_Y <- matrix(nrow = 8, ncol = 8)
  Cov_A_Y[1, 1:8] <- c(  9,    0,    9,    0,    9,    0,    9,    0)
  Cov_A_Y[2, 1:8] <- c(  0,    9,    0,    9,    0,    9,    0,    9)
  Cov_A_Y[3, 1:8] <- c(  9,    0,   11,  1.5,    9,    0,    9,    0)
  Cov_A_Y[4, 1:8] <- c(  0,    9,  1.5,   11,    0,    9,    0,    9)
  Cov_A_Y[5, 1:8] <- c(  9,    0,    9,    0,   11,  1.5,    9,    0)
  Cov_A_Y[6, 1:8] <- c(  0,    9,    0,    9,  1.5,   11,    0,    9)
  Cov_A_Y[7, 1:8] <- c(  9,    0,    9,    0,    9,    0,   11,  1.5)
  Cov_A_Y[8, 1:8] <- c(  0,    9,    0,    9,    0,    9,  1.5,   11)
  Cov_Y <- Cov_A_Y[3:8, 3:8]
  chol_cov <- chol(Cov_Y)
  res <- dmnorm_chol(as.numeric(y), mean(y), cholesky = chol_cov, prec_param=FALSE, log = TRUE)
  expect_equal(opt$value, res)
  # Check covariance matrix
  summ <- cmLaplace$summary(opt, jointCovariance = TRUE)
  # Covariance matrix from TMB:
  # TMB cpp code (test.cpp) below:
  # include <TMB.hpp>
  # template<class Type>
  # Type objective_function<Type>::operator() ()
  # {
  #   DATA_MATRIX(y);
  #   DATA_MATRIX(Sigma);
  #   PARAMETER(mu);
  #   PARAMETER_VECTOR(a);
  #   int i;
  #   Type ans = 0.0;
  #   // Negative log-likelihood
  #   for(i = 0; i < 2; i++){
  #     ans -= dnorm(a[i], 0.5*mu, Type(3.0), true);
  #   }
  #   vector<Type> residual(2);
  #   using namespace density;
  #   MVNORM_t<Type> neg_log_dmvnorm(Sigma);
  #   for(i = 0; i < 3; i++)
  #   {
  #     residual = vector<Type>(y.col(i)) - a;
  #     ans += neg_log_dmvnorm(residual);
  #   }
  #   return ans;
  #   }
  # TMB R code:
  # library(TMB)
  # compile("test.cpp")
  # dyn.load(dynlib("test"))
  # data <- list(y = m$y, Sigma = m$cov_y)
  # parameters <- list(mu = 0, a = c(-2, -1))
  #
  # ## Fit model
  # obj <- MakeADFun(data, parameters, random="a", DLL="test")
  # tmbopt <- nlminb(obj$par, obj$fn, obj$gr)
  # tmbrep <- sdreport(obj, getJointPrecision = TRUE)
  # tmbvcov <- inverse(tmbrep$jointPrecision)
  tmbvcov <- matrix(nrow = 3, ncol = 3)
  tmbvcov[1,] <- c(20.333333, 1.1666667, 1.1666667)
  tmbvcov[2,] <- c(1.166667, 0.6651515, 0.5015152)
  tmbvcov[3,] <- c(1.166667, 0.5015152, 0.6651515)
  expect_equal(summ$vcov, tmbvcov, tol = 1e-4)

  # Check covariance matrix for params only
  summ2 <- cmLaplace$summary(opt, originalScale = TRUE, randomEffectsStdError = TRUE, jointCovariance = FALSE)
  expect_equal(summ2$vcov, tmbvcov[1,1,drop=FALSE], tol=1e-4)

  #cmLaplace$setInnerCache(FALSE)
  cmLaplace$updateSettings(useInnerCache=FALSE)
  summ2_recomp <- cmLaplace$summary(opt, originalScale = TRUE, randomEffectsStdError = TRUE, jointCovariance = FALSE)
  expect_equal(summ2$vcov, summ2_recomp$vcov, tol=1e-10)

  summL <- summaryLaplace(cmLaplace, opt, jointCovariance = TRUE)
  expect_equal(nrow(summL$randomEffects), 2)

  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$findMLE() # some warnings are ok here
  expect_equal(opt$par, optNoSplit$par, tol = 1e-2)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)
})

test_that("Laplace with 2x1D random effects needing joint integration works, with intermediate nodes", {
  set.seed(1)
  y <- matrix(rnorm(6, 4, 5), nrow = 2)
  m <- nimbleModel(
    nimbleCode({
      for(i in 1:2) {
        mu_y[i] <- 0.8*a[i]
        a[i] ~ dnorm(mu_a, sd = 3)
      }
      for(j in 1:3)
        y[1:2, j] ~ dmnorm(mu_y[1:2], cov = cov_y[1:2, 1:2])
      mu_a <- 0.5 * mu
      mu ~ dnorm(0, sd = 5)
    }),
    data = list(y = y),
    inits = list(a = c(-2, -1), mu = 0),
    constants = list(cov_y = matrix(c(2, 1.5, 1.5, 2), nrow = 2)),
    buildDerivs = TRUE
  )

  mLaplace <- buildLaplace(model = m)
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit

  opt <- cmLaplace$findMLE()
  expect_equal(opt$par, mean(y)/(0.8*0.5), tol = 1e-4) # optim's reltol is about 1e-8 but that is for the value, not param.
  # V[a] = 9
  # V[y[1:2, i]] =  diag(2)*0.8^2 * 9 + cov_y = [ (5.76 + 2), 0 + 1.5; 0+1.5, (5.76+2)]
  # Cov[a[1], y[1,j]] = 0.8*9 = 7.2
  # Cov[a[1], y[2,j]] = 0
  # Cov[a[2], y[1,j]] = 0
  # Cov]a[2], y[2,j]] = 0.8*9
  # Cov[y[1,i], y[1,j]] = 0.8^2*9 = 5.76
  # Cov[y[2,i], y[2,j]] = 9
  Cov_A_Y <- matrix(nrow = 8, ncol = 8)
  Cov_A_Y[1, 1:8] <- c(  9,    0,  7.2,    0,  7.2,    0,  7.2,    0)
  Cov_A_Y[2, 1:8] <- c(  0,    9,    0,  7.2,    0,  7.2,    0,  7.2)
  Cov_A_Y[3, 1:8] <- c(7.2,    0, 7.76,  1.5, 5.76,    0, 5.76,    0)
  Cov_A_Y[4, 1:8] <- c(  0,  7.2,  1.5, 7.76,    0, 5.76,    0, 5.76)
  Cov_A_Y[5, 1:8] <- c(7.2,    0, 5.76,    0, 7.76,  1.5, 5.76,    0)
  Cov_A_Y[6, 1:8] <- c(  0,  7.2,    0, 5.76,  1.5, 7.76,    0, 5.76)
  Cov_A_Y[7, 1:8] <- c(7.2,    0, 5.76,    0, 5.76,    0, 7.76,  1.5)
  Cov_A_Y[8, 1:8] <- c(  0,  7.2,    0, 5.76,    0, 5.76,  1.5, 7.76)
  Cov_Y <- Cov_A_Y[3:8, 3:8]
  chol_cov <- chol(Cov_Y)
  res <- dmnorm_chol(as.numeric(y), mean(y), cholesky = chol_cov, prec_param=FALSE, log = TRUE)
  expect_equal(opt$value, res)

  # Check covariance matrix
  summ <- cmLaplace$summary(opt, jointCovariance = TRUE)
  # Covariance matrix from TMB:
  # TMB cpp code (test.cpp) below:
  # include <TMB.hpp>
  # template<class Type>
  # Type objective_function<Type>::operator() ()
  # {
  #   DATA_MATRIX(y);
  #   DATA_MATRIX(Sigma);
  #   PARAMETER(mu);
  #   PARAMETER_VECTOR(a);
  #   int i;
  #   Type ans = 0.0;
  #   // Negative log-likelihood
  #   for(i = 0; i < 2; i++){
  #     ans -= dnorm(a[i], 0.5*mu, Type(3.0), true);
  #   }
  #   vector<Type> residual(2);
  #   using namespace density;
  #   MVNORM_t<Type> neg_log_dmvnorm(Sigma);
  #   for(i = 0; i < 3; i++)
  #   {
  #     residual = vector<Type>(y.col(i)) - 0.8 * a;
  #     ans += neg_log_dmvnorm(residual);
  #   }
  #   return ans;
  #   }
  # TMB R code:
  # library(TMB)
  # compile("test.cpp")
  # dyn.load(dynlib("test"))
  # data <- list(y = m$y, Sigma = m$cov_y)
  # parameters <- list(mu = 0, a = c(-2, -1))
  #
  # ## Fit model
  # obj <- MakeADFun(data, parameters, random="a", DLL="test")
  # tmbopt <- nlminb(obj$par, obj$fn, obj$gr)
  # tmbrep <- sdreport(obj, getJointPrecision = TRUE)
  # tmbvcov <- inverse(tmbrep$jointPrecision)
  tmbvcov <- matrix(nrow = 3, ncol = 3)
  tmbvcov[1,] <- c(21.645833, 1.8229167, 1.8229167)
  tmbvcov[2,] <- c(1.822917, 1.0380050, 0.7849117)
  tmbvcov[3,] <- c(1.822917, 0.7849117, 1.0380050)

  expect_equal(summ$vcov, tmbvcov, tol = 1e-4)
  # Check covariance matrix for params only
  summ2 <- cmLaplace$summary(opt, originalScale = TRUE, randomEffectsStdError = TRUE, jointCovariance = FALSE)
  expect_equal(summ2$vcov, tmbvcov[1,1,drop=FALSE], tol=1e-4)

  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$findMLE() # some warnings are ok here
  expect_equal(opt$par, optNoSplit$par, tol = 1e-2)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)
})

test_that("Laplace with 2x2D random effects for 1D data that are separable works, with intermediate nodes", {
  set.seed(1)
  # y[i, j] is jth datum from ith group
  y <- array(rnorm(8, 6, 5), dim = c(2, 2, 2))
  cov_a <- matrix(c(2, 1.5, 1.5, 2), nrow = 2)
  m <- nimbleModel(
    nimbleCode({
      for(i in 1:2) mu[i] ~ dnorm(0, sd = 10)
      mu_a[1] <- 0.8 * mu[1]
      mu_a[2] <- 0.2 * mu[2]
      for(i in 1:2) a[i, 1:2] ~ dmnorm(mu_a[1:2], cov = cov_a[1:2, 1:2])
      for(i in 1:2) {
        for(j in 1:2) {
          y[1, j, i] ~ dnorm( 0.5 * a[i, 1], sd = 1.8) # this ordering makes it easier below
          y[2, j, i] ~ dnorm( 0.1 * a[i, 2], sd = 1.2)
        }
      }
    }),
    data = list(y = y),
    inits = list(a = matrix(c(-2, -3, 0,  -1), nrow = 2), mu = c(0, .5)),
    constants = list(cov_a = cov_a),
    buildDerivs = TRUE
  )

  mLaplace <- buildLaplace(model = m)
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit

  opt <- cmLaplace$findMLE()

  ## Wei: I tested this using TMB instead of the code below
  # TMB cpp code:
  # #include <TMB.hpp>
  # template<class Type>
  # Type objective_function<Type>::operator() ()
  # {
  #   DATA_ARRAY(y);
  #   DATA_MATRIX(Sigma);
  #   PARAMETER_VECTOR(mu);
  #   PARAMETER_MATRIX(a);
  #   int i, j;
  #   Type ans = 0.0;
  #   vector<Type> mu_a(2);
  #   mu_a(0) = 0.8 * mu(0);
  #   mu_a(1) = 0.2 * mu(1);
  #   // Negative log-likelihood
  #   vector<Type> residual(2);
  #   using namespace density;
  #   MVNORM_t<Type> neg_log_dmvnorm(Sigma);
  #   for(i = 0; i < 2; i++)
  #   {
  #     residual = vector<Type>(a.row(i)) - mu_a;
  #     ans += neg_log_dmvnorm(residual);
  #   }
  #   for(i = 0; i < 2; i++){
  #     for(j = 0; j < 2; j++){
  #       ans -= dnorm(y(0, j, i), 0.5*a(i, 0), Type(1.8), true);
  #       ans -= dnorm(y(1, j, i), 0.1*a(i, 1), Type(1.2), true);
  #     }
  #   }
  #   return ans;
  # }
  # TMB R code:
  # library(TMB)
  # compile("test.cpp")
  # dyn.load(dynlib("test"))
  # data <- list(y = m$y, Sigma = m$cov_a)
  # parameters <- list(mu = m$mu, a = m$a)
  # ## Fit model
  # obj <- MakeADFun(data, parameters, random="a", DLL="test")
  # tmbopt <- nlminb(obj$par, obj$fn, obj$gr)
  # tmbrep <- sdreport(obj, getJointPrecision = TRUE)
  # tmbvcov <- inverse(tmbrep$jointPrecision)
  expect_equal(opt$par, c(12.98392, 406.04878), tol = 1e-4)
  expect_equal(opt$value, -41.86976, tol = 1e-6)
  # Check covariance matrix
  summ <- cmLaplace$summary(opt, jointCovariance = TRUE)
  tmbvcov <- matrix(nrow = 6, ncol = 6)
  tmbvcov[1,] <- c(6.625000e+00, 4.687500e+00,  4.050000e+00,  4.050000e+00, -2.693817e-11, -2.695275e-11)
  tmbvcov[2,] <- c(4.687500e+00, 9.250000e+02,  2.965628e-11,  2.967848e-11,  1.800000e+02,  1.800000e+02)
  tmbvcov[3,] <- c(4.050000e+00, 2.951367e-11,  3.995242e+00,  2.484758e+00,  5.596302e-01, -5.596302e-01)
  tmbvcov[4,] <- c(4.050000e+00, 2.951367e-11,  2.484758e+00,  3.995242e+00, -5.596302e-01,  5.596302e-01)
  tmbvcov[5,] <- c(-2.691772e-11, 1.800000e+02,  5.596302e-01, -5.596302e-01,  3.684693e+01,  3.515307e+01)
  tmbvcov[6,] <- c(-2.691772e-11, 1.800000e+02, -5.596302e-01,  5.596302e-01,  3.515307e+01,  3.684693e+01)

  # The ordering of a[1, 1:2] and a[2, 1:2] is flipped between nimble and TMB:
  expect_equal(summ$vcov[c(1:3, 5, 4, 6), c(1:3, 5, 4, 6)], tmbvcov, tol = 1e-4)

  # Check covariance matrix for params only
  summ2 <- cmLaplace$summary(opt, originalScale = TRUE, randomEffectsStdError = TRUE, jointCovariance = FALSE)
  expect_equal(summ2$vcov, tmbvcov[1:2,1:2], tol=1e-4)

  summL <- summaryLaplace(cmLaplace, opt, jointCovariance = TRUE)
  expect_identical(summL$randomEffects$estimate, summ$randomEffects$estimate)

  # For this case, we build up the correct answer more formulaically
  # Define A as the vector a[1, 1], a[1, 2], a[2, 1], a[2, 2]
  # cov_A <- matrix(0, nrow = 4, ncol = 4)
  # cov_A[1:2, 1:2] <- cov_a
  # cov_A[3:4, 3:4] <- cov_a
  # # Define Y as the vector y[1,1,1],y[2,1,1],y[1,2,1],y[2,2,1], then same with last index 2
  # # Define E[Y] as IA %*% A, where:
  # IA <- matrix(0, nrow = 8, ncol = 4)
  # IA[c(1, 3), 1] <- 0.5
  # IA[c(2, 4), 2] <- 0.1
  # IA[c(5, 7), 3] <- 0.5
  # IA[c(6, 8), 4] <- 0.1
  #
  # # define cov_y_given_a as the Cov[Y | A]
  # cov_y_given_a <- matrix(0, nrow = 8, ncol = 8)
  # diag(cov_y_given_a) <- rep(c(1.8^2, 1.2^2), 4)
  # # And finally get cov_Y, the marginal (over A) covariance of Y
  # cov_Y <- IA %*% cov_A %*% t(IA) + cov_y_given_a
  # chol_cov <- chol(cov_Y)
  #
  # # make a log likelihood function
  # nlogL <- function(mu) {
  #   mean_Y <- rep(c(0.8*0.5*mu[1], 0.2*0.1*mu[2]), 4)
  #   -dmnorm_chol(as.numeric(y), mean_Y, cholesky = chol_cov, prec_param=FALSE, log = TRUE)
  # }
  # # maximize it
  # opt_manual <- optim(c(20, 100), nlogL, method = "BFGS")
  # expect_equal(opt$par, opt_manual$par, tol = 1e-4)
  # expect_equal(opt$value, -opt_manual$value, tol = 1e-5)

  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$findMLE() # some warnings are ok here
  expect_equal(opt$par, optNoSplit$par, tol = 1e-4)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)
})

test_that("Laplace with 2x2D random effects for 2D data that need joint integration works, with intermediate nodes", {
  set.seed(1)
  cov_a <- matrix(c(2, 1.5, 1.5, 2), nrow = 2)
  cov_y <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  y <- rmnorm_chol(1, c(1, 1), chol(cov_y), prec_param = FALSE)
  y <- rbind(y, rmnorm_chol(1, c(1, 1), chol(cov_y), prec_param = FALSE))
  m <- nimbleModel(
    nimbleCode({
      for(i in 1:2) mu[i] ~ dnorm(0, sd = 10)
      mu_a[1] <- 0.8 * mu[1]
      mu_a[2] <- 0.2 * mu[2]
      for(i in 1:2) a[i, 1:2] ~ dmnorm(mu_a[1:2], cov = cov_a[1:2, 1:2])
      mu_y[1:2] <- 0.5*a[1, 1:2] + 0.1*a[2, 1:2]
      for(i in 1:2) {
        y[i, 1:2] ~ dmnorm(mu_y[1:2], cov = cov_y[1:2, 1:2])
      }
    }),
    data = list(y = y),
    inits = list(a = matrix(c(-2, -3, 0,  -1), nrow = 2), mu = c(0, 0.5)),
    constants = list(cov_a = cov_a, cov_y = cov_y),
    buildDerivs = TRUE
  )

  mLaplace <- buildLaplace(model = m)
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit

  opt <- cmLaplace$findMLE()
  ## Check using TMB results
  expect_equal(opt$par, c(0.5603309, 11.7064674 ), tol = 1e-4)
  expect_equal(opt$value, -4.503796, tol = 1e-7)
  # Check covariance matrix
  summ <- cmLaplace$summary(opt, jointCovariance = TRUE)
  tmbvcov <- matrix(nrow = 6, ncol = 6)
  tmbvcov[1,] <- c(4.4270833,  11.111111, 1.4583333, 3.1250000, 0.6597222,  1.9097222)
  tmbvcov[2,] <- c(11.1111111, 70.833333, 2.6388889, 7.6388889, 5.8333333, 12.5000000)
  tmbvcov[3,] <- c(1.4583333,   2.638889, 1.5000000, 0.8333333, 0.7777778,  0.2777778)
  tmbvcov[4,] <- c(3.1250000,   7.638889, 0.8333333, 4.1666667, 0.2777778,  2.7777778)
  tmbvcov[5,] <- c(0.6597222,   5.833333, 0.7777778, 0.2777778, 1.5000000,  0.8333333)
  tmbvcov[6,] <- c(1.9097222,  12.500000, 0.2777778, 2.7777778, 0.8333333,  4.1666667)
  # The ordering of a[1, 1:2] and a[2, 1:2] is flipped between nimble and TMB:
  expect_equal(summ$vcov[c(1:3, 5, 4, 6), c(1:3, 5, 4, 6)], tmbvcov, tol = 1e-4)

  # Check covariance matrix for params only
  summ2 <- cmLaplace$summary(opt, originalScale = TRUE, randomEffectsStdError = TRUE, jointCovariance = FALSE)
  expect_equal(summ2$vcov, tmbvcov[1:2,1:2], tol=1e-4)

  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$findMLE() # some warnings are ok here
  expect_equal(opt$par, optNoSplit$par, tol = 1e-4)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)

  ## TMB cpp code:
  #include <TMB.hpp>
  #template<class Type>
  #Type objective_function<Type>::operator() ()
  # {
  #   DATA_MATRIX(y);
  #   DATA_MATRIX(cov_a);
  #   DATA_MATRIX(cov_y);
  #   PARAMETER_VECTOR(mu);
  #   PARAMETER_MATRIX(a);
  #   int i;
  #   Type ans = 0.0;
  #
  #   using namespace density;
  #   // Negative log-likelihood of mv normal
  #   vector<Type> mu_a(2);
  #   mu_a(0) = 0.8 * mu(0);
  #   mu_a(1) = 0.2 * mu(1);
  #   vector<Type> residual_a(2);
  #   MVNORM_t<Type> dmvnorm_a(cov_a);
  #   for(i = 0; i < 2; i++)
  #   {
  #     residual_a = vector<Type>(a.row(i)) - mu_a;
  #     ans += dmvnorm_a(residual_a);
  #   }
  #   vector<Type> mu_y(2);
  #   mu_y(0) = 0.5*a(0, 0) + 0.1*a(1, 0);
  #   mu_y(1) = 0.5*a(0, 1) + 0.1*a(1, 1);
  #   vector<Type> residual_y(2);
  #   MVNORM_t<Type> dmvnorm_y(cov_y);
  #   for(i = 0; i < 2; i++){
  #     residual_y = vector<Type>(y.row(i)) - mu_y;
  #     ans += dmvnorm_y(residual_y);
  #   }
  #   return ans;
  # }
  # library(TMB)
  # compile("test.cpp")
  # dyn.load(dynlib("test"))
  # data <- list(y = m$y,  cov_a = m$cov_a, cov_y = m$cov_y)
  # parameters <- list(mu = m$mu, a = m$a)
  #
  # ## Fit model
  # obj <- MakeADFun(data, parameters, random="a", DLL="test")
  # tmbopt <- nlminb(obj$par, obj$fn, obj$gr)
  # tmbrep <- sdreport(obj, getJointPrecision = TRUE)
  # tmbvcov <- inverse(tmbrep$jointPrecision)
})


test_that("Laplace with 2x2D random effects for 2D data that need joint integration works, with intermediate nodes, without using normality", {
  set.seed(1)
  cov_a <- matrix(c(2, 1.5, 1.5, 2), nrow = 2)
  cov_y <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  y <- rmnorm_chol(1, c(1, 1), chol(cov_y), prec_param = FALSE)
  y <- rbind(y, rmnorm_chol(1, c(1, 1), chol(cov_y), prec_param = FALSE))
  m <- nimbleModel(
    nimbleCode({
      for(i in 1:2) mu[i] ~ dnorm(0, sd = 10)
      mu_a[1] <- 0.8 * mu[1]
      mu_a[2] <- 0.2 * mu[2]
      for(i in 1:2) a[i, 1:2] ~ dmnorm(mu_a[1:2], cov = cov_a[1:2, 1:2])
      mu_y[1:2] <- 0.5*a[1, 1:2] + 0.1*a[2, 1:2]
      for(i in 1:2) {
        y[i, 1:2] ~ dmnorm(mu_y[1:2], cov = cov_y[1:2, 1:2])
      }
    }),
    data = list(y = y),
    inits = list(a = matrix(c(-2, -3, 0,  -1), nrow = 2), mu = c(0, 0.5)),
    constants = list(cov_a = cov_a, cov_y = cov_y),
    buildDerivs = TRUE
  )

  mLaplace <- buildLaplace(model = m, control = list(ADuseNormality = FALSE))
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE,ADuseNormality = FALSE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit

  opt <- cmLaplace$findMLE()
  ## Check using TMB results
  expect_equal(opt$par, c(0.5603309, 11.7064674 ), tol = 1e-4)
  expect_equal(opt$value, -4.503796, tol = 1e-7)
  # Check covariance matrix
  summ <- cmLaplace$summary(opt, jointCovariance = TRUE)
  tmbvcov <- matrix(nrow = 6, ncol = 6)
  tmbvcov[1,] <- c(4.4270833,  11.111111, 1.4583333, 3.1250000, 0.6597222,  1.9097222)
  tmbvcov[2,] <- c(11.1111111, 70.833333, 2.6388889, 7.6388889, 5.8333333, 12.5000000)
  tmbvcov[3,] <- c(1.4583333,   2.638889, 1.5000000, 0.8333333, 0.7777778,  0.2777778)
  tmbvcov[4,] <- c(3.1250000,   7.638889, 0.8333333, 4.1666667, 0.2777778,  2.7777778)
  tmbvcov[5,] <- c(0.6597222,   5.833333, 0.7777778, 0.2777778, 1.5000000,  0.8333333)
  tmbvcov[6,] <- c(1.9097222,  12.500000, 0.2777778, 2.7777778, 0.8333333,  4.1666667)
  # The ordering of a[1, 1:2] and a[2, 1:2] is flipped between nimble and TMB:
  expect_equal(summ$vcov[c(1:3, 5, 4, 6), c(1:3, 5, 4, 6)], tmbvcov, tol = 1e-4)

  # Check covariance matrix for params only
  summ2 <- cmLaplace$summary(opt, originalScale = TRUE, randomEffectsStdError = TRUE, jointCovariance = FALSE)
  expect_equal(summ2$vcov, tmbvcov[1:2,1:2], tol=1e-4)

  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$findMLE() # some warnings are ok here
  expect_equal(opt$par, optNoSplit$par, tol = 1e-4)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)

  ## TMB cpp code:
  #include <TMB.hpp>
  #template<class Type>
  #Type objective_function<Type>::operator() ()
  # {
  #   DATA_MATRIX(y);
  #   DATA_MATRIX(cov_a);
  #   DATA_MATRIX(cov_y);
  #   PARAMETER_VECTOR(mu);
  #   PARAMETER_MATRIX(a);
  #   int i;
  #   Type ans = 0.0;
  #
  #   using namespace density;
  #   // Negative log-likelihood of mv normal
  #   vector<Type> mu_a(2);
  #   mu_a(0) = 0.8 * mu(0);
  #   mu_a(1) = 0.2 * mu(1);
  #   vector<Type> residual_a(2);
  #   MVNORM_t<Type> dmvnorm_a(cov_a);
  #   for(i = 0; i < 2; i++)
  #   {
  #     residual_a = vector<Type>(a.row(i)) - mu_a;
  #     ans += dmvnorm_a(residual_a);
  #   }
  #   vector<Type> mu_y(2);
  #   mu_y(0) = 0.5*a(0, 0) + 0.1*a(1, 0);
  #   mu_y(1) = 0.5*a(0, 1) + 0.1*a(1, 1);
  #   vector<Type> residual_y(2);
  #   MVNORM_t<Type> dmvnorm_y(cov_y);
  #   for(i = 0; i < 2; i++){
  #     residual_y = vector<Type>(y.row(i)) - mu_y;
  #     ans += dmvnorm_y(residual_y);
  #   }
  #   return ans;
  # }
  # library(TMB)
  # compile("test.cpp")
  # dyn.load(dynlib("test"))
  # data <- list(y = m$y,  cov_a = m$cov_a, cov_y = m$cov_y)
  # parameters <- list(mu = m$mu, a = m$a)
  #
  # ## Fit model
  # obj <- MakeADFun(data, parameters, random="a", DLL="test")
  # tmbopt <- nlminb(obj$par, obj$fn, obj$gr)
  # tmbrep <- sdreport(obj, getJointPrecision = TRUE)
  # tmbvcov <- inverse(tmbrep$jointPrecision)
})


nimbleOptions(enableDerivs = EDopt)
nimbleOptions(buildModelDerivs = BMDopt)
