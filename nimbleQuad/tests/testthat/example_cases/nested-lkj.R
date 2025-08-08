library(nimbleQuad)

qpts <- c(.025,.25,.5,.75,.975)

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

code <- nimbleCode({
    for(j in 1:J) {
        for(i in 1:n)
            y[i,j] ~ dnorm(b[j,1] + b[j,2]*x[i,j,1] + b[j,3]*x[i,j,2], sd = sigma)
        b[j,1:3] ~ dmnorm(z[1:3], cholesky = U[1:3,1:3], prec_param = 0)
    }
    sigma ~ dhalfflat()
    Ustar[1:3,1:3] ~ dlkj_corr_cholesky(1.3, 3)
    U[1:3,1:3] <- uppertri_mult_diag(Ustar[1:3, 1:3], sds[1:3])
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
for(j in 1:J)
    y[,j] <- b[1,j] + b[2,j]*x[,j,1] + b[3,j]*x[,j,2] + rnorm(n)

m <- nimbleModel(code, data=list(y=y,x=x),inits = list(z = rep(0,3), Ustar = diag(3), sds = rep(1, 3), b = matrix(0,J,3),sigma=1), constants = list(n=n,J=J), buildDerivs = TRUE)
cm <- compileNimble(m)

approx <- buildNestedApprox(m)
capprox <- compileNimble(approx, project = m)

result <- runNestedApprox(capprox) ## sigma seems good, sds a bit off
latent_sample <- result$sampleLatents(100000)
param_sample <- result$sampleParams(100000)  

result2 <- runNestedApprox(capprox, originalScale = FALSE) ## transformed Ustar seems a bit off
result2$improveParamMarginals(2:4, nMarginalGrid=7)  # compare elements 2,3,4 to trSmp - looks good!


t(apply(latent_sample,2,quantile,qpts)) # good
t(apply(param_sample,2,quantile,qpts)) # some Ustar decent, some a bit off

result$improveParamMarginals('sds', nMarginalGrid=7) # helps a lot!

transformer <- function(tmp) {
    U <- uppertri_mult_diag(matrix(tmp[1:9],3), tmp[10:12])
    return(c(t(U)%*%U))
}

tmp <- param_sample[,2:13]
cv <- t(apply(tmp, 1, transformer))

t(apply(cv, 2, quantile, qpts))  # offdiag medians look good but not diag medians or tail quantiles

library(nimbleHMC)
m <- nimbleModel(code, data=list(y=y,x=x),inits = list(z = rep(0,3), Ustar = diag(3), sds = rep(1, 3), b = matrix(0,J,3),sigma=1), constants = list(n=n,J=J), buildDerivs = TRUE)
cm <- compileNimble(m)
mcmc <- buildHMC(m, monitors = c('b','Ustar','sigma','sds'))
cmcmc <- compileNimble(mcmc, project = m)

system.time(out <- runMCMC(cmcmc, niter = 11000, nburnin = 1000))

t(apply(out, 2, quantile, qpts))

trans <- parameterTransform(m, nodes = 'Ustar')
trSmp <- t(apply(out[,1:9], 1, trans$transform))
t(apply(trSmp, 2, quantile, qpts))

cvHMC <- t(apply(out[ ,c(1:9, 34:36)], 1, transformer))
t(apply(cvHMC, 2, quantile, qpts))


cm <- compileNimble(m)
mLaplace <- buildLaplace(model = m)

cmLaplace <- compileNimble(mLaplace, project = m)
runLaplace(cmLaplace)

apply(out, 2, mean)
apply(out, 2, sd)
## latent estimates seem pretty good, as well as sigma. Ustar decent. sds seem off.
