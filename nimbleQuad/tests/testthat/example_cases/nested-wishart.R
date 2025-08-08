library(nimbleQuad)

qpts <- c(.025,.25,.5,.75,.975)

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
for(j in 1:J)
    y[,j] <- b[1,j] + b[2,j]*x[,j,1] + b[3,j]*x[,j,2] + rnorm(n)

m <- nimbleModel(code, data=list(y=y,x=x),inits = list(z = rep(0,3), Q = diag(3), b = matrix(0,J,3),sigma=1), constants = list(R=diag(3),n=n,J=J), buildDerivs = TRUE)
cm <- compileNimble(m)

approx <- buildNestedApprox(m)
capprox <- compileNimble(approx, project = m)
result <- runNestedApprox(capprox, originalScale = FALSE)  # Decent alignment with HMC.
result 
result$improveParamMarginals(1:7, nMarginalGrid = 5) # Better, still a bit off, takes a while with 3^6=792 grid points
latent_sample <- result$sampleLatents(100000)
apply(latent_sample,2,quantile,qpts)  # Pretty good.


result <- runNestedApprox(capprox, originalScale = TRUE)

library(nimbleHMC)
mcmc <- buildHMC(m, monitors = c('b','Q','sigma'))
cmcmc <- compileNimble(mcmc, project = m)

system.time(out <- runMCMC(cmcmc, niter = 11000, nburnin = 1000))

trans <- parameterTransform(m, nodes = 'Q')
trSmp <- t(apply(out[,1:9], 1, trans$transform))

apply(trSmp, 2, quantile, qpts)

apply(out[,10:ncol(out)][ , c(1,9,17,2,10,18,3,11,19,4,12,20,5,13,21,6,14,22,7,15,23,8,16,24)], 2, quantile, qpts)

quantile(log(out[,'sigma']),qpts)


cm <- compileNimble(m)

mLaplace <- buildLaplace(model = m)

cmLaplace <- compileNimble(mLaplace, project = m)

MLE <- cmLaplace$findMLE()
cmLaplace$summary(MLE) 
cmLaplace$summary(MLE, originalScale=FALSE)  # names not on transformed scale
runLaplace(cmLaplace)

runLaplace(cmLaplace,originalScale=FALSE)

apply(out, 2, quantile, qpts)  # REs and sigma seem good; Q values rather off both transformed and untransformed.

