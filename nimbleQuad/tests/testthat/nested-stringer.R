
library(aghq)
library(nimble)
library(RTMB)

qpts <- c(.025,.25,.5,.75,.975)

## Simulate Some Data:
set.seed(5212)
dat <- NULL
sigma_re <- c(0.05, 0.2)
beta <- c(1.25, 0.2, -0.1, 0.25)
sitere <- rnorm(10, 0, sigma_re[1])
nyear <- 12
nsite <- 10
nrep <- 20
coefs <- 
for( i in 1:nyear ){
  dat <- rbind(dat, data.frame(year = i, site = rep(1:nsite, each = nrep), 
                              weather = rnorm(nsite*nrep,rnorm(1),0.01), temperature = rnorm(nsite*nrep,rnorm(1),0.01),
                              year_re = rnorm(1, 0, sigma_re[2]), site_re = sitere[rep(1:nsite, each = nrep)]))
}
dat$exp_count <- exp(beta[1] + beta[2]*dat$weather + beta[3]*dat$temperature + beta[4]*dat$year + dat$year_re + dat$site_re)
dat$count <- rpois(nrow(dat), lambda = dat$exp_count)

pars <- list(log_sigma_site = 0, log_sigma_year = 0, beta = c(0,0,0,0), year_re = numeric(nyear), site_re = numeric(nsite))
negll_rtmb <- function(pars){
  getAll(pars)
  sigma_site <- exp(log_sigma_site)
  sigma_year <- exp(log_sigma_year)

  ## Add priors w/ jacobian:
  negll <- 0
  negll <- negll - (dgamma(sigma_site, shape=1, scale=1, log = TRUE) + log_sigma_site) - 
                   (dgamma(sigma_year, shape=1, scale=1, log = TRUE) + log_sigma_year) -
                    sum(dnorm(beta, 0, 1, log = TRUE))

  negll <- negll - sum(dnorm(site_re, 0, sigma_site, log = TRUE)) - sum(dnorm(year_re, 0, sigma_year, log = TRUE))
  mu <- beta[1] + beta[2]*dat$weather + beta[3]*dat$temperature + beta[4]*dat$year + site_re[dat$site] + year_re[dat$year]
  negll <- negll - sum(dpois(dat$count, exp(mu), log = TRUE))
  negll
}

## In this example going to have beta as hyperparams...
obj <- MakeADFun(negll_rtmb, pars, silent = TRUE, random = c("site_re", "year_re"))
opt <- nlminb(obj$par, obj$fn, obj$gr)


## SLOW, but I added fixed effects to marginals:
## Requires hessian from finite differencing: numhessian = TRUE or do it manually:
# obj$he <- function(x){numDeriv::jacobian(obj$gr, x)}

## aghq doesn't actually see the random effects...
# fit.aghq <- aghq(obj, 3, obj$par, control = default_control(negate = TRUE)) ## negated unnormalized log-posterior

## marginal_laplace_tmb makes use of the MakeADFun random = "" and Laplace.
fit.aghq <- marginal_laplace_tmb(obj, 3, obj$par, control = default_control(negate = TRUE, numhessian = TRUE)) ## negated unnormalized log-posterior

## Summaries:
sum.aghq <- summary(fit.aghq) ## Modes and 
get_mode(fit.aghq)

## sigma_site
pdf_cdf_1 <- compute_pdf_and_cdf(fit.aghq$marginals[[1]],interpolation = fit.aghq$control$interpolation, transformation = make_transformation('log','exp'))
plot(pdf_cdf_1$theta, pdf_cdf_1$pdf, type = 'l')
## Option as simulation:
lines(density(compute_quantiles(fit.aghq$marginals[[1]],runif(10000))), col = 'red')
plot(pdf_cdf_1$transparam, pdf_cdf_1$pdf_transparam, type = 'l')
abline(v = sigma_re[1], col = 'red')

## Sigma_year
pdf_cdf_2 <- compute_pdf_and_cdf(fit.aghq$marginals[[2]],interpolation = fit.aghq$control$interpolation, transformation = make_transformation('log','exp'))
plot(pdf_cdf_2$theta, pdf_cdf_2$pdf, type = 'l')
plot(pdf_cdf_2$transparam, pdf_cdf_2$pdf_transparam, type = 'l')
abline(v = sigma_re[2], col = 'red')

## beta[1]
pdf_cdf_3 <- compute_pdf_and_cdf(fit.aghq$marginals[[3]],interpolation = fit.aghq$control$interpolation)
plot(pdf_cdf_3$theta, pdf_cdf_3$pdf, type = 'l')
abline(v = beta[1], col = 'red')

## I couldn't find the actual random draw they must have done for the summaries...
samples <- sample_marginal(fit.aghq, 1000)
plot(density(samples$samps[1,]))
plot(density(samples$samps[2,]))

## NIMBLE VERSION:
###################################
code <- nimbleCode({
    sigma_site ~ dgamma(1,1)
    sigma_year ~ dgamma(1,1)
    for( i in 1:4) beta[i] ~ dnorm(0, 1)
    for( i in 1:nsite ) u[i] ~ dnorm(0, sd=sigma_site)
    for( i in 1:nyear ) v[i] ~ dnorm(0, sd=sigma_year)

    for(i in 1:n) {
      log(lambda[i]) <- beta[1] + beta[2]*weather[i] + beta[3]*temperature[i] + beta[4]*year[i] + u[site[i]] + v[year[i]]
      y[i] ~ dpois(lambda[i])
    }
})

m <- nimbleModel(code, data = list(y = dat$count), 
                 constants = list(n=nrow(dat), nsite=nsite, nyear=nyear, 
                                  weather = dat$weather, temperature = dat$temperature, 
                                  year=dat$year, site=dat$site),
                 inits = list(sigma_site = 0, sigma_year = 0, beta = c(0,0,0,0)), buildDerivs = TRUE)

nim_laplace <- buildLaplace(m, paramNodes = c('sigma_site','sigma_year', 'beta'), randomEffectsNodes = c('u', 'v'))
cm <- compileNimble(m)
cnim_laplace <- compileNimble(approx, project = m)

## Set up NIMBLE for AGHQ:
ff <- list(
  fn = function(x){cnim_laplace$calcLogDens(x, trans = TRUE, includeJacobian = TRUE, includePrior = TRUE),
  gr = function(x){cnim_laplace$gr_LogDens(x, trans = TRUE, includeJacobian = TRUE, includePrior = TRUE)},
  he = function(x){numDeriv::jacobian(gr_LogDens, x, trans = TRUE, includeJacobian = TRUE, includePrior = TRUE)}
)

## In theory most of the same things as with RTMB should be available... But I suspect that they cannot simulate from the Random Effects.
## Having issues with my install to run this right now... Or maybe I coded a bug.
fit.aghq <- aghq(ff, 3, values(model, c("sigma_site", "sigma_year", "beta")), control = default_control(negate = TRUE)) ## negated unnormalized log-posterior
