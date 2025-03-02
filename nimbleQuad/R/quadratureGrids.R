QUAD_CACHE_BASE <- nimbleFunctionVirtual(
    run = function() {},
    methods = list(
        cacheQuadGrid = function(nQuad = double(), nodes = double(2), wgts = double(1),
                                 modeIndex = integer()) {
        },
        nodes = function(indx = integer(0, default = 0)) {
            returnType(double(2))
        },
        weights = function(indx = integer(0, default = 0)) {
            returnType(double(1))
        },
        modeI = function() {
            returnType(integer())
        },
        gridSize = function() {
            returnType(integer())
        },
        checkGrid = function(nQuad = double(0, default = -1), prune = double(0, default = 0)) {
            returnType(logical())
        },
        pruneGrid = function(prune = double(0, default = 0)) {
        }
    )
)

## Method for summing likelihoods on real scale with possible small values.
## Returns back on log scale.
#' @export
logSumExp = nimbleFunction(run = function(log1 = double(), log2 = double()) {
    if (log1 > log2) {
        ans <- log(1 + exp(log2 - log1)) + log1
    } else ans <- log(1 + exp(log1 - log2)) + log2
    returnType(double())
    return(ans)
}, buildDerivs = list(run = list()))

quadGridCache <- nimbleFunction(
    contains = QUAD_CACHE_BASE,
    setup = function() {
        quadGridList_internal <- quadGridListDef
        nodes_cached <- matrix(0, nrow = 1, ncol = 1)
        weights_cached <- c(0, 0)
        modeIndex_cached <- -1
        nGrid_cached <- 0L
        gridBuilt <- FALSE
        prune_ <- 0
        nQuad_ <- -1
        numError <- 1e-10  ## ***CJP More precise?
        d <- 1
    },
    run = function() {
    },
    methods = list(
        cacheQuadGrid = function(nQuad = double(), nodes = double(2), wgts = double(1),
                                 modeIndex = integer()) {
        nodes_cached <<- nodes
        weights_cached <<- wgts
        modeIndex_cached <<- modeIndex
        nGrid_cached <<- dim(nodes_cached)[1]
        gridBuilt <<- TRUE
        nQuad_ <<- nQuad
        prune_ <<- 0
        d <<- dim(nodes_cached)[2]
    },
    checkGrid = function(nQuad = double(0, default = -1), prune = double(0, default = 0)) {
        returnType(logical())
        if (!gridBuilt | ((nQuad > 0) & (nQuad != nQuad_)) | ((prune_ != prune) & (prune_ > 0)))
            return(FALSE) else return(TRUE)
    },
    ## Keep the biggest prune proportion weights.
    ## Note that our weights are for arbitrary functions.
    ## As a result, pruning will be weights adjust for a multivariate normal.
    pruneGrid = function(prune = double(0, default = 0)) {
        if (prune_ == 0) {
            if (!gridBuilt & prune > 0) {
                print("Warning: Cannot prune grid as the quadrature grid isn't built yet.")
            } else {
                ntrim <- 0
                ## Adjust weights:
                weights_adj <- numeric(value = 0, length = nGrid_cached)
                for (i in seq_along(weights_cached)) {
                    weights_adj[i] <- exp(sum(dnorm(nodes_cached[i, ], mean = 0, sd = 1, log = TRUE))) *
                        weights_cached[i]
                }
                ## Leave 3 points in total.
                while (ntrim/nGrid_cached < prune & ntrim < nGrid_cached - 3) {
                    keep <- which(weights_adj > min(weights_adj) + numError)  ## error check as weights might be equal but off by numerical.
                    if (dim(keep)[1] > 0) {
                        weights_adj <- weights_adj[keep]
                        weights_cached <<- weights_cached[keep]
                        nodes_cached <<- matrix(nodes_cached[keep, ], nrow = length(keep), ncol = d)
                        ntrim <- nGrid_cached - dim(nodes_cached)[1]
                        ## Update mode index:
                        if (modeIndex_cached > 0) {
                            modei <- which(keep == modeIndex_cached)
                            if (dim(modei)[1] > 0)
                                modeIndex_cached <<- modei[1] else modeIndex_cached <<- -1
                        }
                    } else {
                        ## Exit loop.
                        ntrim <- nGrid_cached
                    }
                }
                nGrid_cached <<- dim(nodes_cached)[1]
            }
        }
        prune_ <<- prune
    },
    nodes = function(indx = integer(0, default = 0)) {
        returnType(double(2))
        if (indx > 0)
            return(matrix(nodes_cached[indx, ], nrow = 1))
        if (indx == -1 & modeIndex_cached > 0)
            return(matrix(nodes_cached[modeIndex_cached, ], nrow = 1))
        return(nodes_cached)
    },
    weights = function(indx = integer(0, default = 0)) {
        returnType(double(1))
        if (indx > 0)
            return(numeric(value = weights_cached[indx], length = 1))
        if (indx == -1 & modeIndex_cached > 0)
            return(numeric(value = weights_cached[modeIndex_cached], length = 1))
        return(weights_cached)
    },
    modeI = function() {
        returnType(integer())
        return(modeIndex_cached)
    },
    gridSize = function() {
        returnType(integer())
        return(nGrid_cached)
    }
)
)

## Wrapper to make quadrature nodes accesible in a nimble function list. ***CJP
## check better naming convention on nQuad_.
#' @export
configureQuadGrid <- nimbleFunction(
    name = "quadGridClass",
    setup = function(d = 1, nQuad_ = 3, quadRule = "AGHQ", control = list()) {
        ## Can list all possible quad rules here and set it.
        possibleRules <- c("AGHQ", "CCD", "AGHQSPRSE", "USER")
        
        quadRules <- extractControlElement(control, "quadRules", quadRule)
        
        if (!any(quadRule == quadRules))
            quadRules <- c(quadRule, quadRules)

        if (!all(quadRules %in% possibleRules))
            stop("Error:  Only AGHQ, CCD, Sparse AGHQ, or USER suplied rules are currently implemented.")

        prune_ <- extractControlElement(control, "prune", 0)
        if (prune_ > 1 | prune_ < 0)
            stop("Can only prune a proportion of quadrature points.")

        quadGridList_internal <- quadGridListDef
        quadGridCache_nfl <- nimbleFunctionList(QUAD_CACHE_BASE)
        quadRule_nfl <- nimbleFunctionList(QUAD_RULE_BASE)

        I_AGHQ <- I_CCD <- I_USER <- I_AGHQSPRSE <- 1
        I_RULE <- which(quadRules == quadRule)[1]

        for (i in seq_along(quadRules)) {
            quadGridCache_nfl[[i]] <- quadGridCache()
            if (quadRules[i] == "AGHQ") {
                I_AGHQ <- i
                quadRule_nfl[[i]] <- quadRule_AGHQ()
            }
            if (quadRules[i] == "CCD") {
                I_CCD <- i
                quadRule_nfl[[i]] <- quadRule_CCD()
            }
            if (quadRules[i] == "AGHQSPRSE") {
                I_AGHQSPRSE <- i
                quadRule_nfl[[i]] <- quadRule_AGHQSPARSE()
            }
            if (quadRules[i] == "USER") {
                I_USER <- i
                quadRule_nfl[[i]] <- quadRule_USER()
            }
        }

        modeIndex <- -1
        nGrid <- 0
        gridBuilt <- FALSE
    },
    run = function() {
    },
    methods = list(
        ## NOCHNG means keep it as is, and nQuad = -1.
        buildGrid = function(method = character(0, default = "NOCHNG"),
                             nQuad = integer(0, default = -1)) {
            if (method != "NOCHNG") setRule(method)
            if (nQuad != -1) nQuad_ <<- nQuad
            if (!quadGridCache_nfl[[I_RULE]]$checkGrid(nQuad_, prune_) | !gridBuilt) {
                newgrid <- quadGridList_internal$new()
                newgrid <- quadRule_nfl[[I_RULE]]$buildGrid(nQuad = nQuad_, d = d)
                quadGridCache_nfl[[I_RULE]]$cacheQuadGrid(nQuad = nQuad_, nodes = newgrid$nodes,
                                                          wgts = newgrid$wgts, modeIndex = newgrid$modeIndex)
                gridBuilt <<- TRUE
            }

            modeIndex <<- quadGridCache_nfl[[I_RULE]]$modeI()
            nGrid <<- quadGridCache_nfl[[I_RULE]]$gridSize()
        },
        ## Prune grid and then cache it again.
        pruneGrid = function(prune = double(0, default = 0)) {
            if (prune > 1 | prune < 0) stop("Can only prune a proportion of quadrature points.")

            if (I_RULE == I_CCD) print("Warning:  CCD grid cannot be pruned.")

            if (I_RULE != I_CCD) {
                ## Need to rebuild the grid if pruning the grid a second time.
                if (!quadGridCache_nfl[[I_RULE]]$checkGrid(nQuad = nQuad_, prune = prune)) {
                    gridBuilt <<- FALSE
                    buildGrid()
                }
                if (prune > 0) quadGridCache_nfl[[I_RULE]]$pruneGrid(prune)
            }
            prune_ <<- prune
        },
        ## Surely there is a better way to do this...
        setRule = function(method = character(0, default = "AGHQ")) {
            if (method == "AGHQ") I_RULE <<- I_AGHQ
            if (method == "CCD") I_RULE <<- I_CCD
            if (method == "AGHQSPRSE") I_RULE <<- I_AGHQSPRSE
            if (method == "USER") I_RULE <<- I_USER
        },
        setDim = function(ndim = integer(0, default = 1)) {
            if (ndim <= 0) stop("Can't input negative dimensions") else d <<- ndim
            ## Make sure the next grid gets built.
            gridBuilt <<- FALSE
        },
        weights = function(indx = integer(0, default = 0)) {
            if (!gridBuilt) buildGrid()
            if (indx == -1 & modeIndex > 0) indx <- modeIndex
            returnType(double(1))
            return(quadGridCache_nfl[[I_RULE]]$weights(indx = indx))
        },
        nodes = function(indx = integer(0, default = 0)) {
            if (!gridBuilt) buildGrid()
            if (indx == -1 & modeIndex > 0) indx <- modeIndex
            returnType(double(2))
            return(quadGridCache_nfl[[I_RULE]]$nodes(indx = indx))
        },
        gridSize = function() {
            if (!gridBuilt) buildGrid()
            returnType(double())
            return(nGrid)
        },
        modeI = function() {
            if (!gridBuilt) buildGrid()
            returnType(double())
            return(modeIndex)
        }
    )
)  ## End of configureQuadGrid


#' Build Adaptive Gauss-Hermite Quadrature Grid
#'
#' Create quadrature grid for use in AGHQuad methods in Nimble.
#'
#' @param d Dimension of quadrature grid being requested.
#'
#' @param nQuad Number of quadrature nodes requested on build.
#'
#' @name buildAGHQGrid
#' 
#' @details
#'
#' This function is used by used by \code{buildOneAGHQuad1D}
#' and \code{buildOneAGHQuad} create the quadrature grid using
#' adaptive Gauss-Hermite quadrature. Handles single or multiple dimension 
#' grids and computes both grid locations and weights. Additionally, acts
#' as a cache system to do transformations, and return marginalized log density.
#'
#' Any of the input node vectors, when provided, will be processed using
#'   \code{nodes <- model$expandNodeNames(nodes)}, where \code{nodes} may be
#'   \code{paramNodes}, \code{randomEffectsNodes}, and so on. This step allows
#'   any of the inputs to include node-name-like syntax that might contain
#'   multiple nodes. For example, \code{paramNodes = 'beta[1:10]'} can be
#'   provided if there are actually 10 scalar parameters, 'beta[1]' through
#'   'beta[10]'. The actual node names in the model will be determined by the
#'   \code{exapndNodeNames} step.
#'
#' Available methods include
#' 
#' \itemize{
#'
#'   \item \code{buildAGHQ}. Builds a adaptive Gauss-Hermite quadrature grid in d dimensions.
#'   Calls \code{buildAGHQOne} to build the one dimensional grid and then expands in each dimension.
#'   Some numerical issues occur in Eigen decomposition making the grid weights only accurate up to 
#'   35 quadrature nodes.
#'
#'   \item Options to get internally cached values are \code{getGridSize},
#'   \code{getModeIndex} for when there are an odd number of quadrature nodes,
#'   \code{getLogDensity} for the cached values, \code{getAllNodes} for the 
#'   quadrature grids, \code{getNodes} for getting a single indexed nodes,
#'   \code{getAllNodesTransformed} for nodes transformed to the parameter scale,
#'   \code{getNodesTransformed} for a single transformed node, \code{getAllWeights} 
#'   to get all quadrature weights, \code{getWeights} single indexed weight.
#'
#'   \item \code{transformGrid(cholNegHess, inner_mode, method)} transforms 
#'   the grid using either cholesky trasnformations,
#'   as default, or spectral that makes use of the Eigen decomposition. For a single
#'   dimension \code{transformGrid1D} is used.
#'
#'   \item As the log density is evaluated externally, it is saved via \code{saveLogDens},
#'   which then is summed via \code{quadSum}.
#'
#'   \item \code{buildGrid} builds the grid the initial time and is only run once in code. After,
#'   the user must choose to \code{setGridSize} to update the grid size.
#'
#'
#'   \item \code{check}. If TRUE (default), a warning is issued if
#'         \code{paramNodes}, \code{randomEffectsNodes} and/or \code{calcNodes}
#'         are provided but seek to have missing elements or unnecessary
#'         elements based on some default inspection of the model. If
#'         unnecessary warnings are emitted, simply set \code{check=FALSE}.
#'
#'   \item \code{innerOptimControl}. A list of control parameters for the inner 
#'         optimization of Laplace approximation using \code{optim}. See 
#'         'Details' of \code{\link{optim}} for further information.
#'
#'   \item \code{innerOptimMethod}. Optimization method to be used in 
#'         \code{optim} for the inner optimization. See 'Details' of 
#'         \code{\link{optim}}. Currently \code{optim} in NIMBLE supports: 
#'         '\code{Nelder-Mead}', '\code{BFGS}', '\code{CG}', and 
#'         '\code{L-BFGS-B}'. By default, method '\code{CG}' is used when 
#'         marginalizing over a single (scalar) random effect, and '\code{BFGS}' 
#'         is used for multiple random effects being jointly marginalized over.
#'
#'   \item \code{innerOptimStart}. Choice of starting values for the inner 
#'         optimization. This could be \code{'last'}, \code{'last.best'}, or a 
#'         vector of user provided values. \code{'last'} means the most recent 
#'         random effects values left in the model will be used. When finding 
#'         the MLE, the most recent values will be the result of the most recent 
#'         inner optimization for Laplace. \code{'last.best'} means the random 
#'         effects values corresponding to the largest Laplace likelihood (from 
#'         any call to the \code{calcLaplace} or \code{calcLogLik} method, 
#'         including during an MLE search) will be used (even if it was not the 
#'         most recent Laplace likelihood). By default, the initial random 
#'         effects values will be used for inner optimization.
#'
#'   \item \code{outOptimControl}. A list of control parameters for maximizing
#'         the Laplace log-likelihood using \code{optim}. See 'Details' of
#'         \code{\link{optim}} for further information.
#' }
#'
#' @references
#'
#' Golub, G. H. and Welsch, J. H. (1969). Calculation of Gauss Quadrature Rules. 
#' Mathematics of Computation 23 (106): 221-230.
#'
#' Liu, Q. and Pierce, D. A. (1994). A Note on Gauss-Hermite Quadrature. Biometrika, 81(3) 624-629.
#'
#' Jackel, P. (2005). A note on multivariate Gauss-Hermite quadrature. London: ABN-Amro. Re.
#'
NULL


## Create a caching random effects system for simulating the posterior random
## effect distribution according to Stringer: This requires the inner mode, the
## inner cholesky, the

INNER_CACHE_BASE <- nimbleFunctionVirtual(
    run = function() {
    },
    methods = list(
        buildCache = function(nGridUpdate = integer(), nLatentNodes = integer()) {},
        cache_weights = function(weight = double(), indx = integer()) {},
        cache_inner_mode = function(mode = double(1), indx = integer()) {},
        cache_inner_negHessChol = function(negHessChol = double(2), indx = integer()) {},
        weights = function() {
            returnType(double(1))
        },
        simulate = function(n = integer()) {
            returnType(double(2))
        }
    )
)

## Things I need to add for approx posterior save inner mode, save inner
## cholesky etc. for doing simulation of random-effects save outer mode and
## negHessian.  Need wgt*density and inner mode and inner cholesky for each
## point:
inner_cache_methods = nimbleFunction(
    contains = INNER_CACHE_BASE,
    setup = function(nre = 0, nGrid = 0, condIndptSets = NULL, nCondIndptSets = 1) {
        innerMode <- matrix(0, nrow = 1, ncol = 1)
        innerNegHessChol <- array(0, c(1, 1, 1))
        wgtsDens <- c(1, -1)
        cacheBuilt <- FALSE
        if (is.null(condIndptSets)) {
            condIndptSets <- nre  ## Assuming all one set.
            nCondIndptSets <- 1  ## If NULL then this is not relevant.
        }
        if (length(condIndptSets) == 1) {
            condIndptSets <- c(condIndptSets, -1)  ##  Make sure it's a vector.
        }
    },
    run = function() {
    },
    methods = list(
        buildCache = function(nGridUpdate = integer(0, default = -1), nLatentNodes = integer()) {
            nre <<- nLatentNodes
            ## If the cond indpt sets don't match up, don't use.
            if (nre != sum(condIndptSets[1:nCondIndptSets])) {
                print("  Warning: Not able to simulate latent effects from conditionally independent sets.")
                condIndptSets <<- numeric(value = nre, length = 1)
                nCondIndptSets <<- 1
            }

            if (nGridUpdate > 0 & nGridUpdate != nGrid) {
                nGrid <<- nGridUpdate
                cacheBuilt <<- FALSE
            }
            if (!cacheBuilt) {
                nGrid <<- nGridUpdate
                wgtsDens <<- numeric(value = 0, length = nGrid)
                innerMode <<- matrix(0, nrow = nGrid, ncol = nre)
                innerNegHessChol <<- array(0, c(nGrid, nre, nre))
                cacheBuilt <<- TRUE
            }
        },
        ## Note to self, this wgt will be density*wgt, strictly for simulating.
        cache_weights = function(weight = double(), indx = integer()) {
            wgtsDens[indx] <<- weight
        },
        cache_inner_mode = function(mode = double(1), indx = integer()) {
            innerMode[indx, ] <<- mode
        },
        ## Note potentially storing a lot of zeros here. Could break it into a list of cond indpt sets.
        cache_inner_negHessChol = function(negHessChol = double(2), indx = integer()) {
            innerNegHessChol[indx, , ] <<- negHessChol
        },
        weights = function() {
            returnType(double(1))
            return(wgtsDens)
        },
        ## Adding first column to be index for theta.
        simulate = function(n = integer()) {
            val <- matrix(0, nrow = n, ncol = nre + 1)
            simwgt <- wgtsDens/sum(wgtsDens)  ## Did log sum exp when doing input.

            ## Simulate theta points first.  Seems efficient to separate to not
            ## initiate too many index vectors for cond indpt sets.
            for (i in 1:n) {
                k <- rcat(1, prob = simwgt)
                val[i, 1] <- k
                jStart <- 1
                for (j in 1:nCondIndptSets) {
                    val[i, (jStart + 1):(jStart + condIndptSets[j])] <-
                        rmnorm_chol(n = 1, mean = innerMode[k, jStart:(jStart + condIndptSets[j] - 1)],
                                    cholesky = innerNegHessChol[k, jStart:(jStart + condIndptSets[j] - 1),
                                                                jStart:(jStart + condIndptSets[j] - 1)],
                                    prec_param = TRUE)
                }
                jStart <- jStart + condIndptSets[j]
            }
            returnType(double(2))
            return(val)
        }
    )
)
