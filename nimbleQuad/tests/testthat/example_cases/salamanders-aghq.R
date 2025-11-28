library(aghq)
library(RTMB)
data("Salamanders",package = "glmmTMB")
qpts <- c(.025,.25,.5,.75,.975)

Salamanders$sitei <- as.numeric(factor(Salamanders$site))
X <- nimconst$X
Salamanders$z <- (Salamanders$count == 0)*1
pars <- list(logtau_re = log(1/0.2742), logitp = -0.9338, 
  beta = rnorm(ncol(X)), 
  re = rnorm(max(Salamanders$sitei), 0, 0.1))
negll <- function(pars){
  getAll(pars)

  p <- 1/(1+exp(-logitp))
  tau_re <- exp(logtau_re)
  
  negll <- 0
  negll <- negll - dnorm(logitp, -1, 1/sqrt(0.2), log = TRUE)
  negll <- negll - dgamma(tau_re, shape = 1, scale = 1/0.00005, log = TRUE) + logtau_re
  negll <- negll - sum(dnorm(beta[2:length(beta)], 0, 1/sqrt(0.001), log = TRUE))
  negll <- negll - sum(dnorm(re, 0, 1/sqrt(tau_re), log = TRUE))
  
  lam <- exp(X %*% beta + re[Salamanders$sitei])
  dp <- dpois(Salamanders$count, lam)
  
  negll <- negll - sum(log(p*Salamanders$z + (1-p)*dp))
  return(negll)
}

obj <- MakeADFun(negll, pars, random = c("re", "beta"), silent = TRUE)
fit <- nlminb(obj$par, obj$fn, obj$gr)

fit.aghq <- marginal_laplace_tmb(obj, k=11, fit$par, control = default_control(negate = TRUE, numhessian = TRUE))
aghq.logitp <- compute_pdf_and_cdf(fit.aghq$marginals[[2]],interpolation = fit.aghq$control$interpolation)
aghq.taure <- compute_pdf_and_cdf(fit.aghq$marginals[[1]],interpolation = fit.aghq$control$interpolation, transformation = make_transformation('log','exp'))

samples.aghq <- sample_marginal(fit.aghq, 10000)
aghq.b <- apply(samples.aghq$samps[rownames(samples.aghq$samps) == "beta",], 1, quantile, qpts)  
aghq.re <- apply(samples.aghq$samps[rownames(samples.aghq$samps) == "re",], 1, quantile, qpts)  
