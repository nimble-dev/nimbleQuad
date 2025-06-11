## Code to approximate the posterior distribution built on top of inner
## Laplace.  Will start to build the buildApproxPosterior functionality from
## the other branch here as discussed with Chris P.  Note that we will separate
## fixed effects that are normally distributed from the hyperparameters.  *
## diff from Laplace.
#' @export
buildNestedApprox <- nimbleFunction(
    name = "NestedApprox",
    setup = function(model, hyperParamNodes, latentNodes, calcNodes, calcNodesOther, control = list()) {
        split <- extractControlElement(control, "split", TRUE)
        check <- extractControlElement(control, "check", TRUE)
        innerOptimWarning <- extractControlElement(control, "innerOptimWarning", FALSE)

        nQuadOuter <- extractControlElement(control, "nQuadOuter", 3)
        nQuadInner <- extractControlElement(control, "nQuadInner", 1)
        quadRuleMarginal <- extractControlElement(control, "marginalGridRule", "AGHQ")
        pruneMargGrid <- extractControlElement(control, "marginalGridPrune", 0)

        transformMethod <- extractControlElement(control, "quadTransform", "spectral")

        ## Default starting value for Approx Posterior set here and passed to
        ## Laplace.  zero makes sense but should test others on the AGHQ grid and
        ## see how they work.
        control$innerOptimStart <- extractControlElement(control, "innerOptimStart",
                                                         "zero")

        inferenceNodes <- model$getNodeNames(includeData = FALSE, stochOnly = TRUE)
        if(!missing(hyperParamNodes) && !all(model$expandNodeNames(hyperParamNodes) %in% inferenceNodes))
            stop("some elements of `hyperParamNodes` do not have prior distributions")
        if(!missing(latentNodes) && !all(model$expandNodeNames(latentNodes) %in% inferenceNodes))
            stop("some elements of `latentNodes` do not have prior distributions")

        
        margNodes <- splitLatents(model, hyperParamNodes, latentNodes)
        paramNodes <- margNodes$paramNodes
        latentNodes <- margNodes$randomEffectsNodes
        if(!length(paramNodes))
            stop("No parameter nodes detected in model. Please check the model structure or provide parameter nodes explicitly via `hyperParamNodes`.")
        if(!length(latentNodes))
            stop("No latent nodes detected in model. Please check the model structure or provide latent nodes explicitly via `latentNodes`.")

        ## We need `theta_length` now; this should be ok, as `paramNodes` shouldn't be changed
        ## by creation of `innerMethods`.
        paramsTransform <- parameterTransform(model, paramNodes, control = list(allowDeterm = FALSE))
        theta_length <- paramsTransform$getTransformedLength()

       
        ## Configure all grids before calling AGHQ to make sure it builds
        ## correctly.  DO NOT MOVE WHEN THIS IS CALLED
        allGridRules <- c("CCD", "AGHQ", "AGHQSPARSE")

        ## Default outer grid to CCD unless low dimensional.
        hyperGridRule <- extractControlElement(control, "hyperGridRule", "none")
        pruneHyperGrid <- extractControlElement(control, "hyperGridPrune", 0)
        if(hyperGridRule == "none")
            hyperGridRule <- ifelse(theta_length >= 3, "CCD", "AGHQ")

        messageIfVerbose("Building nested posterior approximation for the following node sets:\n",
                         "  - parameter nodes: ", makeNodeString(paramNodes, model), "\n",
                         "  - latent nodes: ", makeNodeString(latentNodes, model), "\n",
                         "  with ", hyperGridRule, " grid for the parameters and ", ifelse(nQuadInner > 1, "AGHQ", "Laplace"), " approximation for the latent nodes.")
 
        if(length(intersect(latentNodes, paramNodes)))
            stop("some nodes appear in both the parameter and latent sets")
        if (theta_length > 20)
            messageIfVerbose("  [Warning] There is a large number of parameter node elements. Computation may be slow.")

        
        if(hyperGridRule == "AGHQ" && nQuadOuter %% 2 == 0)
            messageIfVerbose("  [Note] For computational efficiency, it is recommended to use an odd number of quadrature points\n         for the parameter (outer) grid (`nQuadOuter`).")
        
        ## Default to CCD
        theta_grid <- configureQuadGrid(d = 1, levels = nQuadOuter, quadRule = hyperGridRule,
                                        control = list(quadRules = allGridRules))

        
        innerMethods <- buildAGHQ(model, nQuadInner, paramNodes, latentNodes, calcNodes,
                                  calcNodesOther, control)

        if(!identical(paramNodes, innerMethods$paramNodes))
            stop("`paramNodes` has unexpectedly changed. This should not have occurred.")

        theta_indices <- innerMethods$pTransform_indices
        
        ## Need to check this as it's is now computed in the 'buildAGHQ' function:
        nre <- innerMethods$nre

        ## Simulate from conditionally independent sets.  Do this via number of
        ## sets and the length of each.
        nInternalRESets <- length(innerMethods$AGHQuad_nfl)
        lenInternalRENodeSets <- innerMethods$lenInternalRENodeSets

        ## Outer optimization settings
        outerOptimControl_ <- nimOptimDefaultControl()
        optimControlArgNames <- c("trace", "fnscale", "parscale", "ndeps", "maxit", "abstol",
                                  "reltol", "alpha", "beta", "gamma", "REPORT", "type", "lmm", "factr", "pgtol",
                                  "temp", "tmax")
        if (!is.null(control$outerOptimControl)) {
            validNames <- intersect(names(control$outerOptimControl), optimControlArgNames)
            numValidNames <- length(validNames)
            if (numValidNames > 0) {
                for (i in 1:numValidNames) {
                    outerOptimControl_[[validNames[i]]] <- control$outerOptimControl[[validNames[i]]]
                }
            }
        }
        outerOptimControl_$fnscale <- -1

        ## Hyperparameters on real scale are named "params".
        ## On transformed scale are named "theta".
        paramNodes <- innerMethods$paramNodes
        npar <- innerMethods$npar
        paramNodesAsScalars_vec <- innerMethods$paramNodesAsScalars_vec

        ## If we use this need to add to one time fixes.
        latentNodesAsScalars_vec <- innerMethods$reNodesAsScalars_vec

        ## Set up mapping of parameter names to indices of transformed elements for
        ## 1:1 cases for determination of parameters for which approximate
        ## marginals are possible and for use when users request marginals by node
        ## name.
        paramNodesComponents <- model$expandNodeNames(paramNodes, returnScalarComponents = TRUE)
        paramNodesIndices <- seq_along(paramNodesComponents)

        if (any(paramsTransform$transformType > 9, na.rm = TRUE))
            stop("buildNestedApprox: Unknown parameter transform type: ",
                 paste0(paramsTransform$transformType[paramsTransform$transformType > 9], collapse = ", "))
        
        mapping <- paramsTransform$transformData
        for (idx in seq_len(paramsTransform$nNodes)) {
            if (paramsTransform$transformType[idx] < 7) {  # 1:1 cases; see `parameterTransform.R`.
                paramNodesIndices[mapping[idx, 1]] <- mapping[idx, 3]
            } else { # Wishart/inverse-Wishart, Dirichlet, LKJ
                paramNodesIndices[mapping[idx,1]:mapping[idx,2]] <- 0
            }
        }
        
        setupOutputs(model, paramNodesComponents, paramNodesIndices)
        

        ## Indicator for removing the redundant index -1 in theta_indices
        one_time_fixes_done <- FALSE

        ## Default calculation method for AGHQuad
        computeMethod_ <- extractControlElement(control, "computeMethod", 2)
        useInnerCache_ <- extractControlElement(control, "useInnerCache", TRUE)

        ## For compilation have to set the dimension after setup.
        theta_grid$setDim(ndim = theta_length)
        inner_grid_cache_nfl <- nimbleFunctionList(INNER_CACHE_BASE)

        ## Make sure the grids match the theta_grid numbers.  Initialize them with
        ## nre = 0 in case they aren't used to not generate too much data.
        I_GRID <- theta_grid$I_RULE
        I_CCD <- theta_grid$I_CCD
        inner_grid_cache_nfl[[I_CCD]] <- inner_cache_methods(nre = 0, nGrid = 1,
            condIndptSets = lenInternalRENodeSets, nCondIndptSets = nInternalRESets)
        I_AGHQ <- theta_grid$I_AGHQ
        inner_grid_cache_nfl[[I_AGHQ]] <- inner_cache_methods(nre = 0, nGrid = 1,
            condIndptSets = lenInternalRENodeSets, nCondIndptSets = nInternalRESets)
        ## Note that it is not valid to use a sparse grid here, but will initiate this list for ordering issues.
        I_AGHQSPARSE <- theta_grid$I_AGHQSPARSE
        inner_grid_cache_nfl[[I_AGHQSPARSE]] <- inner_cache_methods(nre = 0, nGrid = 1,
            condIndptSets = lenInternalRENodeSets, nCondIndptSets = nInternalRESets)

        I_USER <- 1
        if (any(allGridRules == "USER")) {
            I_USER <- theta_grid$I_USER
            inner_grid_cache_nfl[[I_USER]] <- inner_cache_methods(nre = 0, nGrid = 1,
                condIndptSets = lenInternalRENodeSets, nCondIndptSets = nInternalRESets)
        }

        ## Store the quadrature sums for each grid:
        marginalPostDensity <- rep(-Inf, length(allGridRules))

        theta_marg_grid <- configureQuadGrid(d = theta_length - 1, levels = 1,
                                             quadRule = quadRuleMarginal, control = list(quadRules = c("AGHQ", "AGHQSPARSE")))
        theta1_nodes <- matrix(0, nrow = 1, ncol = 2)

        ## Cached values for convenience: For marginal distributions in AGHQ over
        ## d-1.
        pTransformFix <- 0
        indexFix <- 0

        ## We will want to cache the standard deviation skew terms.  Default will
        ## be not to skew (e.g. 1)
        covTheta <- matrix(0, nrow = theta_length, ncol = theta_length)
        cholNegHess <- matrix(0, theta_length, theta_length)
        A_spectral <- matrix(0, theta_length, theta_length)
        Ainverse_spectral <- matrix(0, theta_length, theta_length)
        eigenCached <- FALSE
        cholCached <- FALSE

        Atransform <- matrix(0, theta_length, theta_length)
        AinverseTransform <- matrix(0, theta_length, theta_length)

        skewedStdDev <- matrix(1, nrow = theta_length, ncol = 2)
        logSkewedWgt <- 0

        ## Points for Asymmetric Gaussian Interpolation (integration free...?)
        ## Taken from INLA Code:
        extraPoints <- c(-15, -10, -7, 7, 10, 15, -5, -3, -2, -1, -0.5, -0.25, 0, 0.25,
                         0.5, 1, 2, 3, 5)
        aghqPoints <- pracma::gaussHermite(51)$x * sqrt(2)
        zMargGrid <- sort(unique(c(extraPoints, aghqPoints)))
        nzMargGrid <- length(zMargGrid)

        ## Some cached values for summary statistics and reporting:
        marg_P <- matrix(0, nrow = nzMargGrid, ncol = npar)
        marg_theta <- array(0, c(theta_length, nzMargGrid, 2))
        ## ***Do we want the simulations of the latent effects to be cached or just
        ## returned? @CJP?
        post_sims <- matrix(0, nrow = 3, ncol = nre)

        ## Optim info:
        modeCached <- FALSE
        thetaMode <- numeric(theta_length)
        if (theta_length == 1) {
            thetaMode <- c(thetaMode, -1)
            theta_indices <- c(theta_indices, -1)
        }
        thetaNegHess <- matrix(0, nrow = theta_length, ncol = theta_length)
        logPostProbMode <- 0
        logDetNegHessTheta <- 0

        ## Other cached values:
        skewedSDCached <- FALSE
        ## Must be cached for each grid: Up to 3 currently.
        hyperGridCached <- c(FALSE, FALSE, FALSE)

        ## Indicator for removing the redundant index -1 in theta_indices
        one_time_fixes_done <- FALSE
    },
    run = function() {},
    methods = list(
        one_time_fixes = function() {
            if (one_time_fixes_done) return()
            if (theta_length == 1) {
                theta_indices <<- numeric(length = 1, value = 1)
                thetaMode <<- numeric(length = 1, value = 0)
            }
            one_time_fixes_done <<- TRUE
        },
        ## Posterior mode for hyperparameters. findMAP
        posteriorMode = function(pStart = double(1, default = Inf),
                                 hessian = logical(0, default = TRUE),
                                 parscale = character(0, default = "transformed")) {
            optRes <- innerMethods$optimize(pStart = pStart, includePrior = TRUE,
                                            includeJacobian = TRUE,
                                            hessian = TRUE, parscale = parscale)
            dm <- dim(optRes$hessian)[1]
            if(dm != theta_length)
                stop("Posterior mode could not be found. Consider adjusting the control parameters for the optimization via the `control` argument of `buildNestedApprox`.")
            if(any_nan(c(optRes$hessian)))
                stop("While attempting to find posterior mode, invalid hessian calculated. Consider adjusting the control parameters for the optimization via the `control` argument of `buildNestedApprox`.")
            if(optRes$convergence != 0)
                print("  [Warning] In optimization over parameters to find the posterior mode as the\n",
                      "            starting point for setting up the parameter grid,\n",
                      "            `optim` has a non-zero convergence code: ", optRes$convergence, ".\n",
                      "            Approximation may not be accurate.")

            
            modeCached <<- TRUE
            thetaMode <<- optRes$par
            if(theta_length == 1) {
                thetaNegHess <<- matrix(-optRes$hessian, 1, 1)
            } else thetaNegHess <<- -optRes$hessian
            logPostProbMode <<- optRes$value
            covTheta <<- inverse(thetaNegHess)
            return(optRes)
            returnType(optimResultNimbleList())
        },
        buildHyperGrid = function(quadRule = character(0, default = "NULL"),
                                  nQuadUpdate = integer(0, default = -1),
                                  prune = double(0, default = -1)) {
            one_time_fixes()
            if(nQuadUpdate != -1)
                nQuadOuter <<- nQuadUpdate
            if(quadRule != "NULL")
                setHyperGridRule(quadRule)
            if(prune != -1) 
                pruneHyperGrid <<- prune
            theta_grid$buildGrid(method = hyperGridRule, nQuad = nQuadOuter, prune = pruneHyperGrid)
            nGrid <- theta_grid$gridSize()
            inner_grid_cache_nfl[[I_GRID]]$buildCache(nGridUpdate = nGrid, nLatentNodes = nre)
            if (!modeCached) posteriorMode(rep(Inf, npar), hessian = TRUE, parscale = "transformed")
        },
        setHyperGridRule = function(quadRule = character(0, default = "AGHQ")) {
            ## Add a rule check here to make sure it's valid.
            hyperGridRule <<- quadRule
            ## Default to AGHQ and change if requested.
            I_GRID <<- I_AGHQ
            if (quadRule == "CCD") I_GRID <<- I_CCD
            if (quadRule == "USER") I_GRID <<- I_USER
            if (quadRule == "AGHQSPARSE") I_GRID <<- I_AGHQSPARSE
        },
        calcEigen = function() {
            E <- eigen(thetaNegHess, symmetric = TRUE)  ## Should be symmetric...
            for (d in 1:theta_length) {
                A_spectral[, d] <<- E$vectors[, d]/sqrt(E$values[d])
                Ainverse_spectral[, d] <<- E$vectors[, d] * sqrt(E$values[d]) # Strictly speaking, the transpose of the inverse of A.
            }
            logDetNegHessTheta <<- sum(log(E$values))
            eigenCached <<- TRUE
        },
        calcCholesky = function() {
            cholNegHess <<- chol(thetaNegHess)
            logDetNegHessTheta <<- 2 * sum(log(diag(cholNegHess)))
            cholCached <<- TRUE
        },
        ## Need this to swap between cholesky and spectral.
        setTransformations = function(method = character(0, default = "spectral")) {
            if (method == "spectral") {
                if (!eigenCached) calcEigen()
                Atransform <<- A_spectral
                AinverseTransform <<- Ainverse_spectral
            } else {
                if (!cholCached) calcCholesky()
                Atransform <<- cholNegHess # Used with backsolve in `z_to_theta`.
                AinverseTransform <<- cholNegHess
            }
        },
        ## Transform from standard (z) to param transform (theta) scale.
        z_to_theta = function(z = double(1), postMode = double(1), A = double(2),
                              method = character(0, default = "spectral")) {
            if (method == "spectral") {
                d <- dim(z)[1]
                theta <- numeric(value = 0, length = d)
                for (i in 1:d) {  # A %*% z
                    theta[i] <- postMode[i] + sum(A[i,] * z) 
                }
            } else {
                if(method == "identity")
                  theta <- z
                else
                  theta <- postMode + backsolve(A, z)
            }
            returnType(double(1))
            return(theta)
        },
        ## Transform from param transform (theta) to standard (z) scale.
        theta_to_z = function(theta = double(1), postMode = double(1), A = double(2),
                              method = character(0, default = "spectral")) {
            if (method == "spectral") {
                ## A provided will be inverse of transpose of true A.
                d <- dim(theta)[1]
                z <- numeric(value = 0, length = d)
                theta_mean <- theta - postMode
                for (i in 1:d) {  # t(A) %*% theta_mean
                    z[i] <- sum(A[,i] * theta_mean)
                }
            } else {
                if(method == "identity")
                  z <- theta
                else
                  z <- (A %*% (theta - postMode))[, 1]
            }
            returnType(double(1))
            return(z)
        },
        calcSkewedSD = function() {
            ## Require the grid to have been built and the mode found.
            buildHyperGrid()

            setTransformations(transformMethod)
            logSkewedWgt <<- 0
            for (i in 1:theta_length) {
                z <- numeric(value = 0, length = theta_length)
                z[i] <- -sqrt(2)
                theta <- z_to_theta(z, thetaMode, Atransform, transformMethod)
                logDens2Neg <- innerMethods$calcLogDens_pTransformed(pTransform = theta)
                skewedStdDev[i, 1] <<- sqrt(2/(2 * (logPostProbMode - logDens2Neg)))  ## numerator (-sqrt(2)) ^2
                z[i] <- sqrt(2)
                theta <- z_to_theta(z, thetaMode, Atransform, transformMethod)
                logDens2Pos <- innerMethods$calcLogDens_pTransformed(pTransform = theta)
                skewedStdDev[i, 2] <<- sqrt(2/(2 * (logPostProbMode - logDens2Pos)))  ## numerator (-sqrt(2)) ^2
                logSkewedWgt <<- logSkewedWgt + log(sum(skewedStdDev[i, ]/2))
                if(any(skewedStdDev[i,] < 0.3) | any(skewedStdDev[i,] > 3.333))
                    nimCat("  [Warning] Skewness in posterior of the hyperparameters in dimension ", i, " is large and a potential sign of an issue for these approximations.\n")
            }
            skewedSDCached <<- TRUE
        },
        getSkewedStdDev = function() {
            returnType(double(2))
            return(skewedStdDev)
        },
        ## INLA like function for approx marginal likelihood (based on skewed normal).
        calcMarginalLogLikApprox = function() {
            if (!skewedSDCached) calcSkewedSD()
            ## Line 2748 in r-inla/blob/devel/gmrflib/approx-inference.c Commit #
            ## ef4eb20 marg <- logPostProbMode + 0.5*theta_length*log(2*pi) -
            ## 0.5*(logDetNegHessTheta) - sum(log(skewedStdDev[,1] * skewedStdDev[,2]))
            ## *** What Paul thinks it should be. ***
            marg <- logPostProbMode + 0.5 * theta_length * log(2 * pi) - 0.5 * (logDetNegHessTheta) +
                logSkewedWgt  # sum(log((skewedStdDev[,1] + skewedStdDev[,2])/2))
            returnType(double())
            return(marg)
        },
        ## This is the meat and potatoes for being able to make inference on the latent nodes.
        ## Calculate theta on the quadrature grid points. AGHQ or CCD.
        ## Stores all values we need for simulation inference on the latent nodes.
        calcHyperGrid = function(skew = logical(0, default = TRUE)) {
            if(I_GRID == I_AGHQSPARSE)
                print("  [Note] Sparse grids cannot be used to simulate latent effects which is the main reason to compute posterior on the hyper grid.")

            buildHyperGrid()
            setTransformations(transformMethod)
            nGrid <- theta_grid$gridSize()
            
            if (!skewedSDCached & skew) calcSkewedSD()
            ans <- 0
            ## Now fill in the grid values.
            nimCat("Calculating inner AGHQ/Laplace approximation at ", nGrid, " outer (parameter) grid points (one dot per point): ")
            for (i in 1:nGrid) {
                nimCat(".")
                ## Operations at the mode:
                if (i == theta_grid$modeIndex()) {
                    wgt <- theta_grid$weights(indx = i)[1]
                    inner_grid_cache_nfl[[I_GRID]]$cache_inner_mode(
                        mode = innerMethods$get_inner_mode(atOuterMode = 1), indx = i)
                    inner_grid_cache_nfl[[I_GRID]]$cache_inner_negHessChol(
                        negHessChol = innerMethods$get_inner_cholesky(atOuterMode = 1), indx = i)
                    inner_grid_cache_nfl[[I_GRID]]$cache_weights(weight = wgt, indx = i)
                    ans <- ans + wgt
                } else {
                    wgt <- theta_grid$weights(indx = i)[1]
                    node <- theta_grid$nodes(indx = i)[1, ]

                    ## Skew the CCD values:
                    if (skew) {
                        for (d in 1:theta_length) {
                            node[d] <- node[d] * skewedStdDev[d, step(node[d]) + 1]  ## negative skew column 1, positive skew column 2
                        }
                    }
                    ## Transform to theta scale:
                    node <- z_to_theta(node, thetaMode, Atransform, transformMethod)
                    thetaLogPostDens <- innerMethods$calcLogDens_pTransformed(node)
                    wgt_dens <- wgt * exp(thetaLogPostDens - logPostProbMode)
                    ## Marginal sum:
                    ans <- ans + wgt_dens

                    ## Cache everything for simulation:
                    inner_grid_cache_nfl[[I_GRID]]$cache_inner_mode(
                        mode = innerMethods$get_inner_mode(atOuterMode = 0), indx = i)
                    inner_grid_cache_nfl[[I_GRID]]$cache_inner_negHessChol(
                        negHessChol = innerMethods$get_inner_cholesky(atOuterMode = 0), indx = i)
                    inner_grid_cache_nfl[[I_GRID]]$cache_weights(weight = wgt_dens, indx = i)
                }
                ## *** Add a convergence check?
            }
            nimCat("\n")
            if (skew) adjLogWgt <- logSkewedWgt else adjLogWgt <- 0

            ## Marginal log posterior density, a normalizing constant for other
            ## methods.
            marginalPostDensity[I_GRID] <<- log(ans) + logPostProbMode - 0.5 * logDetNegHessTheta +
                adjLogWgt

            hyperGridCached[I_GRID] <<- TRUE
        },
        ## Quadrature based marginal log-likelihood
        ## Probably not particularly accurate for CCD.
        calcMarginalLogLikQuad = function() {
            if (I_GRID == I_CCD)
                print("  [Note]: Estimating marginal log-likelihood based on CCD grid.\n           Estimation based on an AGHQ grid may be more accurate (but more computationally expensive).")
            if(!hyperGridCached[I_GRID])
                calcHyperGrid()
            returnType(double())
            return(marginalPostDensity[I_GRID])
        },
        ## Marginals AGHQ from Stringer et al.
        ## *** Investigate pruning for AGHQ.
        ## *** Is this a good name? Tooooo long.
        ## *** Need to make this for theta 1D as well. No AGHQ needed in that case.
        findMarginalPosteriorDensity = function(pIndex = integer(),
                                                nPts = integer(0, default = 5),
                                                nQuad = integer(0, default = 3),
                                                gridTransformMethod = character(0, default = "spectral"),
                                                quadRule = character(0, default = "NULL"),
                                                prune = double(0, default = -1)) {
            one_time_fixes()
                                  
            ## Build the quadrature grid points:
            if (dim(theta1_nodes)[1] != nPts) theta1_nodes <<- quadGH(levels = nPts, type = "GHe")

            if(theta_length == 1) {
                ## d-1 dimensional AGHQ not needed; just evaluate inner approximation.
                res <- matrix(0, nrow = nPts, ncol = 2)
                stdDev <- sqrt(covTheta[1, 1])
                for(i in 1:nPts) {
                    res[i, 1] <- theta1_nodes[i, 2] * stdDev + thetaMode[pIndex]
                    res[i, 2] <- innerMethods$calcLogDens_pTransformed(c(res[i, 1]))
                }
                return(res)
            }

            if(nQuad %% 2 == 0)
                cat("  [Note] For computational efficiency, it is recommended to use an odd number of quadrature points\n         (via argument `nQuad`) for marginalizing over the parameter (outer) grid.\n")

            ## Grid for additional theta.
            ## Build marginal AGHQ grid to compute the hyperparameter marginals
            ## (integrate over pT-1 theta values).      
            if( pruneMargGrid != -1)
              pruneMargGrid <<- prune
            theta_marg_grid$buildGrid(method = quadRule, nQuad = nQuad, prune = pruneMargGrid)  

            nQuadGrid <- theta_marg_grid$gridSize()

            if (!modeCached) posteriorMode(rep(Inf, npar), hessian = TRUE, parscale = "transformed")  ## *** default is now nlminb
            
            ## 1D quadrature to evaluate the theta on.
            stdDev <- sqrt(covTheta[pIndex, pIndex])

            ## Initialize optimization at theta mode.
            Atransform_i <- matrix(0, nrow = theta_length - 1, ncol = theta_length - 1)

            ## Column 1 is chosen theta values, Column 2 is marginalized values, normalized based on AGHQ.
            res <- matrix(0, nrow = nPts, ncol = 2)
            thetaj <- thetaMode
            other_theta_indices <- theta_indices[theta_indices != pIndex]

            ## For each value of thetai, we need to do AGHQ which means finding the
            ## mode of the other parameters, transforming and computing.
            nimCat("Calculating inner AGHQ/Laplace approximation at (", nPts, ") marginal points with ", nQuadGrid, " quadrature grid points (one dot per grid point): ")
            for (i in 1:nPts) {
                res[i, 1] <- theta1_nodes[i, 2] * stdDev + thetaMode[pIndex]
                thetaj[pIndex] <- res[i, 1]

                ## If this is the mode then we know optim already:
                if (theta1_nodes[i, 2] == 0) {
                    theta_iMode <- thetaMode[other_theta_indices]
                    subsetNegHess <- thetaNegHess[other_theta_indices, other_theta_indices]
                    maxPostDensi <- logPostProbMode
                } else {
                    optRes <- innerMethods$findMax_fixedp(pStartTransform = thetaMode, pTransformIndex = pIndex,
                        pTransformValue = res[i, 1], includePrior = TRUE, includeJacobian = TRUE,
                        hessian = TRUE)
                    subsetNegHess <- -optRes$hessian
                    theta_iMode <- optRes$par
                    maxPostDensi <- optRes$value
                }
                
                if (gridTransformMethod == "spectral") {
                    E <- eigen(subsetNegHess, symmetric = TRUE)
                    for (d in 1:(theta_length-1)) {
                        Atransform_i[, d] <- E$vectors[, d]/sqrt(E$values[d])
                    }
                    logDetNegHessThetai <- sum(log(E$values))
                } else {
                    Atransform_i <- chol(subsetNegHess)
                    logDetNegHessThetai <- 2 * sum(log(diag(Atransform_i)))
                }

                logDensi <- 0
                nimCat("(", i, ")")
                for (j in 1:nQuadGrid) {
                    nimCat(".")
                    if (j != theta_marg_grid$modeIndex()) {
                        nodej <- theta_marg_grid$nodes(indx = j)[1, ]
                        theta_tmp <- z_to_theta(z = nodej, postMode = theta_iMode, A = Atransform_i,
                                                method = gridTransformMethod)
                        thetaj[other_theta_indices] <- theta_tmp
                        postLogDensij <- innerMethods$calcLogDens_pTransformed(pTransform = thetaj)
                        logDensi <- logDensi + exp(postLogDensij - maxPostDensi) * theta_marg_grid$weights(indx = j)[1]
                    } else {
                        logDensi <- logDensi + theta_marg_grid$weights(indx = j)[1]
                    }
                }
                res[i, 2] <- log(logDensi) + maxPostDensi - 0.5 * logDetNegHessThetai
            }
            nimCat("\n")
            ## Because thetai values are AGHQ, we can normalize to get the proper
            ## posterior density.  This lets us get the marginal posterior via spline
            ## without any more normalizing (but note that in `fitMarginalSpline`
            ## we do also normalize.
            ## Note that this is a 1-d quadrature,
            ## normalizing P(thetai,Y) to get P(Y) rather than the expensive
            ## calculation of denominator in (8) in Bilodeau et al.
            margi <- sum(exp(res[, 2] - logPostProbMode) * theta1_nodes[, 1])
            lognormconst <- log(margi) + logPostProbMode + log(stdDev)
            res[, 2] <- res[, 2] - lognormconst
            ## *** Should I cache this?
            returnType(double(2))
            return(res)
        },
        ## This can't be until I've built the CCD grid.
        ## so that we have covTheta.
        ## Should also ensure that if they plan to skew the grid that is also done.
        findMarginalHyperIntFree = function(pIndex = integer()) {
            ## Error Trapping:
            if(pIndex <= 0 | pIndex > theta_length)
                stop("Transformed parameter index requested is larger than available.")
                
            ## Requires running `calcSkewedSD()` first.
            if (!skewedSDCached) calcSkewedSD()
            
            stdDev <- sqrt(covTheta[pIndex, pIndex])
            thetai <- numeric(value = 0, length = theta_length)
            setTransformations(transformMethod)
            for (i in 1:nzMargGrid) {
                                        # Known fixed # of points
                thetai[pIndex] <- thetaMode[pIndex] + zMargGrid[i] * stdDev
                marg_theta[pIndex, i, 1] <<- thetai[pIndex]
                ## Find the conditional mean:
                for (j in 1:theta_length) {
                    if (j != pIndex) {
                        thetai[j] <- thetaMode[j] + covTheta[pIndex, j] / covTheta[pIndex, pIndex] *
                            (thetai[pIndex] - thetaMode[pIndex])
                    }
                }
                ## Calculate asymmetric Gaussian:
                zi <- theta_to_z(thetai, thetaMode, AinverseTransform, transformMethod)

                ## logDens = sum log(exp(-z^2/sigma_(+/-))) *Not normalized.
                ## Can we normalize analytically? ***CJP?
                logDens <- 0
                for (j in 1:theta_length) {
                    side <- 2
                    if (zi[j] <= 0) side <- 1
                    logDens <- logDens - 0.5 * (zi[j]/skewedStdDev[j, side])^2
                }
                marg_theta[pIndex, i, 2] <<- logDens
            }
            returnType(double(2))
            return(marg_theta[pIndex, , ])
        },
        ## marginalTransformedSplineDensity = function(pIndex = integer()) {
        ## returnType(double(2))
        ## return(marginalSplineR(marg_theta[pIndex, , 1], marg_theta[pIndex, , 2]))
        ## },
        simulateLatentEffects = function(n = integer()) {
            if(n < 0)
              stop("Cannot simulate less than 1 values.")
            if(I_GRID == I_AGHQSPARSE)
              stop("Sparse grids can have negative weights and are not valid for simulating the latent effects.")
            if (!hyperGridCached[I_GRID]) calcHyperGrid()

            sims <- inner_grid_cache_nfl[[I_GRID]]$simulate(n)
            returnType(double(2))
            return(sims)
        },
        ## Simulation method for theta marginal on the skewed multivariate normal.
         simulateHyperParams = function(n = integer()) {
            if(n < 0)
              stop("Cannot simulate less than 1 values.")

            sims <- matrix(0, nrow = n, ncol = theta_length)
            if (!skewedSDCached) calcSkewedSD()

            setTransformations(transformMethod)

            prob <- skewedStdDev[, 2]/(skewedStdDev[, 1] + skewedStdDev[, 2])

            for (i in 1:n) {
                ## simulate z on the base scale
                z <- abs(rnorm(theta_length, 0, 1))
                for (j in 1:theta_length) {
                    dir <- rbinom(1, 1, prob[j])  ## Skew z pos if 1, neg if 0.
                    if (dir == 1) { z[j] <- skewedStdDev[j, 2] * z[j]
                    } else z[j] <- -skewedStdDev[j,1] * z[j]
                }
                ## Scale it based on method
                sims[i, ] <- z_to_theta(z, thetaMode, Atransform, transformMethod)
            }
            returnType(double(2))
            return(sims)
        },
        getParamGrid = function() {
            return(theta_grid$nodes())
            returnType(double(2))
        }
    )
)
