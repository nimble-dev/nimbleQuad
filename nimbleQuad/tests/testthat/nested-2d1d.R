qpts <- c(.025,.25,.5,.75,.975)

## d=1 case with two latents

code <- nimbleCode({
    for(i in 1:n)
        y[i] ~ dnorm(b0 + b1*x[i], sd = sqrt(1/tau))
    b0 ~ dnorm(0,.001)
    b1 ~ dnorm(0,.001)
    tau~dgamma(1, rate=5e-5)
})


set.seed(1)
n <- 30
x <- rnorm(n)
y <- 0.3*x + rnorm(n)
m <- nimbleModel(code, data = list(y = y, x = x), constants = list(n=n),
                 inits = list(b0 = 0, b1 = .5, tau = 1), buildDerivs = TRUE)

approx <- buildNestedApprox(m, latentNodes = c('b0','b1'), hyperParamNodes = 'tau')
cm <- compileNimble(m)
capprox <- compileNimble(approx, project = m)
result <- runNestedApprox(capprox)

result

result$improveMarginals('tau', nMarginalGrid = 31)

mu_sample <- result$sampleLatentNodes(100000)
apply(mu_sample, 2, quantile, qpts)

## INLA

library(INLA)
fit <- inla(y~x, family="gaussian", data=data.frame(y=y,x=x), quantiles = qpts)
summary(fit)
