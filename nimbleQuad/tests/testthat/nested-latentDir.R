library(nimbleQuad)
code <- nimbleCode({
    for(j in 1:J) {
        for(i in 1:n)
            y[i,j] ~ dcat(p[j,1:3])
        p[j,1:3] ~ ddirch(alpha[1:3])
    }
    alpha[1] <- alpha0
    alpha[2] <- alpha0
    alpha[3] <- alpha0
    alpha0 ~ dunif(0, 20)
})


set.seed(1)
n <- 10
J <- 8
K <- 3
alpha0 <- 2
p <- t(sapply(seq_len(J), function(x) rdirch(1, rep(alpha0, K))))
y <- apply(p, 1, function(q) rcat(n, q))
m <- nimbleModel(code, data = list(y = y), constants = list(n=n, J=J),
                 inits = list(alpha0 = 1, p = p), buildDerivs = TRUE)

## Check that Laplace works.

cm <- compileNimble(m)
mLaplace <- buildLaplace(model = m)
cmLaplace <- compileNimble(mLaplace, project = m)
result <- runLaplace(cmLaplace, jointCovariance = TRUE)
result <- runLaplace(cmLaplace, jointCovariance = TRUE, originalScale = FALSE)
result$summary

## Does nested approx work with transformed, non-1:1 latents?

approx <- buildNestedApprox(m, latentNodes = c('p'), hyperParamNodes = c('alpha0'))
cm <- compileNimble(m)
capprox <- compileNimble(approx, project = m)
result <- runNestedApprox(capprox) 

latent_sample <- result$sampleLatentNodes(100, includeParams = TRUE)

result <- runNestedApprox(capprox, originalScale = FALSE)
latent_sample_trans <- result$sampleLatentNodes(100, includeParams = TRUE)
