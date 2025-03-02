## Tooling that allows us to use known derivative information for latent
## normal nodes in Laplace, rather than computing via AD.

getParam_BASE <- nimbleFunctionVirtual(
    run = function() {},
    methods = list(
        getMean = function(index = integer()) {
            returnType(double(1))
        },
        getPrecision = function(index = integer()) {
            returnType(double(2))
        },
        calcGradient = function(reTransform = double(1), index = integer(), first = integer(),
                                last = integer()) {
            returnType(double(1))
        }
    )
)

## A place holder to not take up much memory.
emptyParam <- nimbleFunction(
    contains = getParam_BASE,
    setup = function() {},
    run = function() {},
    methods = list(
        getPrecision = function(index = integer()) {
            returnType(double(2))
            return(matrix(1, nrow = 1, ncol = 1))
        },
        getMean = function(index = integer()) {
            returnType(double(1))
            return(numeric(1, length = 1))
        },
        calcGradient = function(reTransform = double(1), index = integer(), first = integer(),
                                last = integer()) {
            returnType(double(1))
            return(numeric(1, length = 1))
        }
    )
)

## Need at least one dnorm to use this.  NodeNames relate to node names in the
## model that are dmnrom distributed gNodes (length of all randomEffectsNodes)
## indicates a 1 if dmnorm, 0 o/w.  This makes it easy to get the correct
## indices when I just pass it the random-effect index in a loop.
gaussParam <- nimbleFunction(
    contains = getParam_BASE,
    setup = function(model, nodeNames, gNodes) {
    indexConvert <- cumsum(gNodes)
    },
    run = function() {},
    methods = list(
        getPrecision = function(index = integer()) {
            i <- indexConvert[index]
            Q <- matrix(model$getParam(nodeNames[i], "tau"), nrow = 1, ncol = 1)
            returnType(double(2))
            return(Q)
        },
        getMean = function(index = integer()) {
            i <- indexConvert[index]
            mu <- numeric(model$getParam(nodeNames[i], "mean"), length = 1)
            returnType(double(1))
            return(mu)
        },
        ## Avoid too much memory creation by adding this internal.
        calcGradient = function(reTransform = double(1), index = integer(), first = integer(),
                                last = integer()) {
            i <- indexConvert[index]
            ans <- -model$getParam(nodeNames[i], "tau") * (reTransform[first] -
                model$getParam(nodeNames[i], "mean"))
            returnType(double(1))
            return(numeric(value = ans, length = 1))
        }
    )
)

## Need at least one dmnorm to use this.  NodeNames relate to node names in the
## model that are dmnrom distributed gNodes (length of all randomEffectsNodes)
## indicates a 1 if dmnorm, 0 o/w.  This makes it easy to get the correct
## indices when I just pass it the random-effect index in a loop.
multiGaussParam <- nimbleFunction(
    contains = getParam_BASE,
    setup = function(model, nodeNames, gNodes) {
        indexConvert <- cumsum(gNodes)
    },
    run = function() {},
    methods = list(
        getPrecision = function(index = integer()) {
            i <- indexConvert[index]
            Q <- model$getParam(nodeNames[i], "prec")
            returnType(double(2))
            return(Q)
        },
        getMean = function(index = integer()) {
            i <- indexConvert[index]
            mu <- model$getParam(nodeNames[i], "mean")
            returnType(double(1))
            return(mu)
        },
        ## Avoid too much memory creation by adding this internal.
        calcGradient = function(reTransform = double(1), index = integer(), first = integer(),
                                last = integer()) {
            i <- indexConvert[index]
            U <- model$getParam(nodeNames[i], "cholesky")
            bstar <- (reTransform[first:last] - model$getParam(nodeNames[i], "mean"))
            
            ## If they saved it as a cholesky of precision. Let's use that for speed.
            if (model$getParam(nodeNames[i], "prec_param") == 1) {
                ans <- -(t(U) %*% (U %*% bstar))[, 1]
            } else {
                ans <- -backsolve(U, forwardsolve(t(U), bstar))
                ## ans <- - model$getParam(nodeNames[i], 'prec') %*%
                ## (reTransform[first:last] - model$getParam(nodeNames[i], 'mean'))
            }
            returnType(double(1))
            return(ans)
        }
    )
)
