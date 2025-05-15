## Poi-gamma case so check effect of non-normal latents

library(nimbleQuad)
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

approx <- buildNestedApprox(m, latentNodes = c('mu'), paramNodes = c('a','b'))
cm <- compileNimble(m)
capprox <- compileNimble(approx, project = m)
result <- runNestedApprox(capprox) # -139.9048

latent_sample <- result$sampleLatents(100)
result  # -139.8688

approx <- buildNestedApprox(m, latentNodes = c('mu'), paramNodes = c('a','b'),
                            control = list(nQuadOuter = 7))
cm <- compileNimble(m)
capprox <- compileNimble(approx, project = m)
result <- runNestedApprox(capprox)  

latent_sample <- result$sampleLatents(100)
result  # -139.8579


approx <- buildNestedApprox(m, latentNodes = c('mu'), paramNodes = c('a','b'),
                            control = list(nQuadOuter = 7, nQuadInner = 3))
cm <- compileNimble(m)
capprox <- compileNimble(approx, project = m)
result <- runNestedApprox(capprox) # -139.8575 

latent_sample <- result$sampleLatents(100)
result  


M <- 500000
py <- rep(0, M)

dens <- function(idx) {
    ytmp <- y[(1+(idx-1)*n):(idx*n)]
    return(a*log(b) - lgamma(a) - sum(lgamma(ytmp+1)) + lgamma(a+sum(ytmp)) - (a+sum(ytmp)) * log(b + n))
}

set.seed(1)
for(m in seq_along(py)) {
    a <- rgamma(1, 1,1)
    b <- rgamma(1,1,1)
    logpy <- sum(sapply(1:J, dens))
    py[m] <- exp(logpy)
}

log(mean(py)) # -139.8087, -139.8075, -139.8108

## INLA can't fit this model as latents are not normal
