## Code to approximate the posterior distribution built on top of inner Laplace.
#' @export
buildNestedApprox <- nimbleFunction(
    name = "nestedApprox",
    setup = function(model, paramNodes, latentNodes, calcNodes, calcNodesOther, control = list()) {
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
        if(!missing(paramNodes) && !all(model$expandNodeNames(paramNodes) %in% inferenceNodes))
            stop("some elements of `paramNodes` do not have prior distributions")
        if(!missing(latentNodes) && !all(model$expandNodeNames(latentNodes) %in% inferenceNodes))
            stop("some elements of `latentNodes` do not have prior distributions")

        
        margNodes <- splitLatents(model, paramNodes, latentNodes)
        paramNodes <- margNodes$paramNodes
        latentNodes <- margNodes$randomEffectsNodes
        if(!length(paramNodes))
            stop("No parameter nodes detected in model. Please check the model structure or provide parameter nodes explicitly via `paramNodes`.")
        if(!length(latentNodes))
            stop("No latent nodes detected in model. Please check the model structure or provide latent nodes explicitly via `latentNodes`.")

        ## We need `nParamTrans` now; this should be ok, as `paramNodes` shouldn't be changed
        ## by creation of `innerMethods`.
        paramsTransform <- parameterTransform(model, paramNodes, control = list(allowDeterm = FALSE))
        nParamTrans <- paramsTransform$getTransformedLength()

       
        ## Configure all grids before calling AGHQ to make sure it builds
        ## correctly.  DO NOT MOVE WHEN THIS IS CALLED
        allGridRules <- c("CCD", "AGHQ", "AGHQSPARSE")

        ## Default outer grid to CCD unless low dimensional.
        paramGridRule <- extractControlElement(control, "paramGridRule", "none")
        pruneParamGrid <- extractControlElement(control, "paramGridPrune", 0)
        if(paramGridRule == "none")
            paramGridRule <- ifelse(nParamTrans >= 3, "CCD", "AGHQ")

        messageIfVerbose("Building nested posterior approximation for the following node sets:\n",
                         "  - parameter nodes: ", makeNodeString(paramNodes, model), "\n",
                         "  - latent nodes: ", makeNodeString(latentNodes, model), "\n",
                         "  with ", paramGridRule, " grid for the parameters and ", ifelse(nQuadInner > 1, "AGHQ", "Laplace"), " approximation for the latent nodes.")
 
        if(length(intersect(latentNodes, paramNodes)))
            stop("some nodes appear in both the parameter and latent sets")
        if (nParamTrans > 20)
            messageIfVerbose("  [Warning] There is a large number of parameter node elements. Computation may be slow.")

        
        if(paramGridRule == "AGHQ" && nQuadOuter %% 2 == 0)
            messageIfVerbose("  [Note] For computational efficiency, it is recommended to use an odd number of quadrature points\n         for the parameter (outer) grid (`nQuadOuter`).")
        
        ## Default to CCD (in which case `nQuadOuter` is ignored).
        paramGrid <- configureQuadGrid(d = 1, levels = nQuadOuter, quadRule = paramGridRule, control = list(quadRules = allGridRules))

        
        innerMethods <- buildAGHQ(model, nQuadInner, paramNodes, latentNodes, calcNodes,
                                  calcNodesOther, control)

        if(!identical(paramNodes, innerMethods$paramNodes))
            stop("`paramNodes` has unexpectedly changed. This should not have occurred.")

        paramTrans_indices <- innerMethods$pTransform_indices
        
        nreTrans <- innerMethods$nreTrans

        ## Simulate from conditionally independent sets.  Do this via number of
        ## sets and the length of each.
        lenInternalRENodeSets <- innerMethods$getREtransLength()
        nInternalRESets <- length(lenInternalRENodeSets)
            
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
        ## On transformed scale are named "paramTrans".
        paramNodes <- innerMethods$paramNodes
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
        

        ## Indicator for removing the redundant index -1 in paramTrans_indices
        one_time_fixes_done <- FALSE

        ## Default calculation method for AGHQuad
        computeMethod_ <- extractControlElement(control, "computeMethod", 2)
        useInnerCache_ <- extractControlElement(control, "useInnerCache", TRUE)

        ## For compilation have to set the dimension after setup.
        paramGrid$setDim(ndim = nParamTrans)
        inner_grid_cache_nfl <- nimbleFunctionList(INNER_CACHE_BASE)

        ## Make sure the grids match the paramGrid numbers.  Initialize them with
        ## nre = 0 in case they aren't used to not generate too much data.
        I_GRID <- paramGrid$I_RULE
        I_CCD <- paramGrid$I_CCD
        inner_grid_cache_nfl[[I_CCD]] <- inner_cache_methods(nre = 0, nGrid = 1,
            condIndepSets = lenInternalRENodeSets, nCondIndepSets = nInternalRESets)
        I_AGHQ <- paramGrid$I_AGHQ
        inner_grid_cache_nfl[[I_AGHQ]] <- inner_cache_methods(nre = 0, nGrid = 1,
            condIndepSets = lenInternalRENodeSets, nCondIndepSets = nInternalRESets)
        ## Note that it is not valid to use a sparse grid here, but will initiate this list for ordering issues.
        I_AGHQSPARSE <- paramGrid$I_AGHQSPARSE
        inner_grid_cache_nfl[[I_AGHQSPARSE]] <- inner_cache_methods(nre = 0, nGrid = 1,
            condIndepSets = lenInternalRENodeSets, nCondIndepSets = nInternalRESets)

        I_USER <- 1
        if (any(allGridRules == "USER")) {
            I_USER <- paramGrid$I_USER
            inner_grid_cache_nfl[[I_USER]] <- inner_cache_methods(nre = 0, nGrid = 1,
                condIndepSets = lenInternalRENodeSets, nCondIndepSets = nInternalRESets)
        }

        ## Store the quadrature sums for each grid:
        marginalPostDensity <- rep(-Inf, length(allGridRules))

        paramMargGrid <- configureQuadGrid(d = nParamTrans - 1, levels = 1,
                                             quadRule = quadRuleMarginal, control = list(quadRules = c("AGHQ", "AGHQSPARSE")))
        paramTrans1_nodes <- matrix(0, nrow = 1, ncol = 2)

        ## Cached values for convenience: For marginal distributions in AGHQ over
        ## d-1.
        pTransformFix <- 0
        indexFix <- 0

        ## We will want to cache the standard deviation skew terms.  Default will
        ## be not to skew (e.g. 1)
        covParamTrans <- matrix(0, nrow = nParamTrans, ncol = nParamTrans)
        cholNegHess <- matrix(0, nParamTrans, nParamTrans)
        A_spectral <- matrix(0, nParamTrans, nParamTrans)
        AinverseTrans_spectral <- matrix(0, nParamTrans, nParamTrans)
        eigenCached <- FALSE
        cholCached <- FALSE

        Atransform <- matrix(0, nParamTrans, nParamTrans)
        AinverseTransform <- matrix(0, nParamTrans, nParamTrans)

        skewedStdDev <- matrix(1, nrow = nParamTrans, ncol = 2)
        logSkewedWgt <- 0

        ## Points for Asymmetric Gaussian Interpolation (integration free...?)
        ## Taken from INLA Code:
        extraPoints <- c(-15, -10, -7, 7, 10, 15, -5, -3, -2, -1, -0.5, -0.25, 0, 0.25,
                         0.5, 1, 2, 3, 5)
        aghqPoints <- pracma::gaussHermite(51)$x * sqrt(2)
        zMargGrid <- sort(unique(c(extraPoints, aghqPoints)))
        nzMargGrid <- length(zMargGrid)

        ## Some cached values for summary statistics and reporting:
        marg_paramTrans <- array(0, c(nParamTrans, nzMargGrid, 2))

        ## Optim info:
        modeCached <- FALSE
        paramTransMode <- numeric(nParamTrans)
        if (nParamTrans == 1) {
            paramTransMode <- c(paramTransMode, -1)
            paramTrans_indices <- c(paramTrans_indices, -1)
        }
        paramTransNegHess <- matrix(0, nrow = nParamTrans, ncol = nParamTrans)
        logPostProbMode <- 0
        logDetNegHessParamTrans <- 0

        ## Other cached values:
        skewedSDCached <- FALSE
        ## Must be cached for each grid: Up to 4 currently.
        paramGridCached <- c(FALSE, FALSE, FALSE, FALSE)

        ## Indicator for removing the redundant index -1 in paramTrans_indices
        one_time_fixes_done <- FALSE
    },
    run = function() {},
    methods = list(
        one_time_fixes = function() {
            if (one_time_fixes_done) return()
            if (nParamTrans == 1) {
                paramTrans_indices <<- numeric(length = 1, value = 1)
                paramTransMode <<- numeric(length = 1, value = 0)
            }
            one_time_fixes_done <<- TRUE
        },
        ## Posterior mode for hyperparameters. findMAP
        findMode = function(pStart = double(1, default = Inf),
                                 hessian = logical(0, default = TRUE),
                                 parscale = character(0, default = "transformed")) {
            optRes <- innerMethods$optimize(pStart = pStart, includePrior = TRUE,
                                            includeJacobian = TRUE,
                                            hessian = TRUE, parscale = parscale)
            dm <- dim(optRes$hessian)[1]
            if(dm != nParamTrans)
                stop("Posterior mode could not be found. Consider adjusting the control parameters for the optimization via the `control` argument of `buildNestedApprox`.")
            if(any_nan(c(optRes$hessian)))
                stop("While attempting to find posterior mode, invalid hessian calculated. Consider adjusting the control parameters for the optimization via the `control` argument of `buildNestedApprox`.")
            if(optRes$convergence != 0)
                print("  [Warning] In optimization over parameters to find the posterior mode as the\n",
                      "            starting point for setting up the parameter grid,\n",
                      "            `optim` has a non-zero convergence code: ", optRes$convergence, ".\n",
                      "            Approximation may not be accurate.")

            
            modeCached <<- TRUE
            paramTransMode <<- optRes$par
            if(nParamTrans == 1) {
                paramTransNegHess <<- matrix(-optRes$hessian, 1, 1)
            } else paramTransNegHess <<- -optRes$hessian
            logPostProbMode <<- optRes$value
            covParamTrans <<- inverse(paramTransNegHess)
            return(optRes)
            returnType(optimResultNimbleList())
        },
        buildParamGrid = function(quadRule = character(0, default = "NULL"),
                                  nQuadUpdate = integer(0, default = -1),
                                  prune = double(0, default = -1)) {
            one_time_fixes()
            if(nQuadUpdate != -1)
                nQuadOuter <<- nQuadUpdate
            if(quadRule != "NULL" )
              setParamGridRule(quadRule)
            if(prune != -1)
                pruneParamGrid <<- prune
            paramGrid$buildGrid(method = paramGridRule, nQuad = nQuadOuter, prune = pruneParamGrid)
            nGrid <- paramGrid$gridSize()
            nCache <- inner_grid_cache_nfl[[I_GRID]]$gridSize()
            inner_grid_cache_nfl[[I_GRID]]$buildCache(nGridUpdate = nGrid, nLatents = nreTrans)
            ## If grid changed, then need to update here that we have to calcParamGrid again too.
            if(nGrid != nCache) paramGridCached[I_GRID] <<- FALSE
            if (!modeCached) findMode(rep(Inf, nParamTrans), hessian = TRUE, parscale = "transformed")
        },
        setParamGridRule = function(quadRule = character(0, default = "NULL")) {
            valid_rule <- TRUE
            if(quadRule == "AGHQ"){
              I_GRID <<- I_AGHQ
            }else if(quadRule == "CCD"){
              I_GRID <<- I_CCD
            }else if(quadRule == "USER"){
              I_GRID <<- I_USER
            }else if(quadRule == "AGHQSPARSE"){
              I_GRID <<- I_AGHQSPARSE
            }else{
              valid_rule <- FALSE
            }
            if(valid_rule)
              paramGridRule <<- quadRule
        },
        calcEigen = function() {
            E <- eigen(paramTransNegHess, symmetric = TRUE)  ## Should be symmetric...
            for (d in 1:nParamTrans) {
                A_spectral[, d] <<- E$vectors[, d]/sqrt(E$values[d])
                AinverseTrans_spectral[, d] <<- E$vectors[, d] * sqrt(E$values[d]) 
            }
            logDetNegHessParamTrans <<- sum(log(E$values))
            eigenCached <<- TRUE
        },
        calcCholesky = function() {
            cholNegHess <<- chol(paramTransNegHess)
            logDetNegHessParamTrans <<- 2 * sum(log(diag(cholNegHess)))
            cholCached <<- TRUE
        },
        ## Need this to swap between cholesky and spectral.
        setTransformations = function(method = character(0, default = "spectral")) {
            if (method == "spectral") {
                if (!eigenCached) calcEigen()
                Atransform <<- A_spectral
                AinverseTransform <<- AinverseTrans_spectral
            } else {
                if (!cholCached) calcCholesky()
                Atransform <<- cholNegHess # Used with backsolve in `z_to_paramTrans`.
                AinverseTransform <<- cholNegHess
            }
        },
        ## Transform from standard (z) to param transform (paramTrans) scale.
        z_to_paramTrans = function(z = double(1), postMode = double(1), A = double(2),
                              method = character(0, default = "spectral")) {
            if (method == "spectral") {
                d <- dim(z)[1]
                paramTrans <- numeric(value = 0, length = d)
                for (i in 1:d) {  # A %*% z
                    paramTrans[i] <- postMode[i] + sum(A[i,] * z) 
                }
            } else {
                if(method == "identity")
                  paramTrans <- z
                else
                  paramTrans <- postMode + backsolve(A, z)
            }
            returnType(double(1))
            return(paramTrans)
        },
        ## Transform from param transform (paramTrans) to standard (z) scale.
        paramTrans_to_z = function(paramTrans = double(1), postMode = double(1), A = double(2),
                              method = character(0, default = "spectral")) {
            if (method == "spectral") {
                ## 'A' provided needs to be transpose of inverse of true A.
                d <- dim(paramTrans)[1]
                z <- numeric(value = 0, length = d)
                paramTrans_mean <- paramTrans - postMode
                for (i in 1:d) {  # t(A) %*% paramTrans_mean
                    z[i] <- sum(A[,i] * paramTrans_mean)
                }
            } else {
                if(method == "identity")
                  z <- paramTrans
                else
                  z <- (A %*% (paramTrans - postMode))[, 1]
            }
            returnType(double(1))
            return(z)
        },
        calcSkewedSD = function() {
            ## Require the grid to have been built and the mode found.
            buildParamGrid()

            setTransformations(transformMethod)
            logSkewedWgt <<- 0
            for (i in 1:nParamTrans) {
                z <- numeric(value = 0, length = nParamTrans)
                z[i] <- -sqrt(2)
                paramTrans <- z_to_paramTrans(z, paramTransMode, Atransform, transformMethod)
                logDens2Neg <- innerMethods$calcLogDens_pTransformed(pTransform = paramTrans)
                skewedStdDev[i, 1] <<- sqrt(2/(2 * (logPostProbMode - logDens2Neg)))  ## numerator (-sqrt(2)) ^2
                z[i] <- sqrt(2)
                paramTrans <- z_to_paramTrans(z, paramTransMode, Atransform, transformMethod)
                logDens2Pos <- innerMethods$calcLogDens_pTransformed(pTransform = paramTrans)
                skewedStdDev[i, 2] <<- sqrt(2/(2 * (logPostProbMode - logDens2Pos)))  ## numerator (-sqrt(2)) ^2
                logSkewedWgt <<- logSkewedWgt + log(sum(skewedStdDev[i, ]/2))
                if(any(skewedStdDev[i,] < 0.3) | any(skewedStdDev[i,] > 3.333))
                    nimCat("  [Warning] Skewness in posterior of the (hyper)parameters in dimension ", i, " is large and a potential sign of an issue for these approximations.\n")
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
            ## ef4eb20 marg <- logPostProbMode + 0.5*nParamTrans*log(2*pi) -
            ## 0.5*(logDetNegHessParamTrans) - sum(log(skewedStdDev[,1] * skewedStdDev[,2]))
            ## *** What Paul thinks it should be. ***
            marg <- logPostProbMode + 0.5 * nParamTrans * log(2 * pi) - 0.5 * (logDetNegHessParamTrans) +
                logSkewedWgt  # sum(log((skewedStdDev[,1] + skewedStdDev[,2])/2))
            returnType(double())
            return(marg)
        },
        ## This is the meat and potatoes for being able to make inference on the latent nodes.
        ## Calculate paramTrans on the quadrature grid points. AGHQ or CCD.
        ## Stores all values we need for simulation inference on the latent nodes.
        calcParamGrid = function(skew = logical(0, default = TRUE)) {
            if(I_GRID == I_AGHQSPARSE)
                print("  [Note] Sparse grids cannot be used to simulate latent effects which is the main reason to compute posterior on the hyper grid.")

            buildParamGrid()
            setTransformations(transformMethod)
            nGrid <- paramGrid$gridSize()
            
            ## Check to make sure cache is working correctly:
            if(length(inner_grid_cache_nfl[[I_GRID]]$weights()) != nGrid)
                stop("The grid cache system does not match the quadrature grid being calculated.")
            
            if (!skewedSDCached & skew) calcSkewedSD()
            ans <- 0
            ## Now fill in the grid values.
            nimCat("Calculating inner AGHQ/Laplace approximation at ", nGrid, " outer (parameter) grid points (one dot per point): ")
            for (i in 1:nGrid) {
                nimCat(".")
                ## Operations at the mode:
                if (i == paramGrid$modeIndex()) {
                    wgt <- paramGrid$weights(idx = i)[1]
                    inner_grid_cache_nfl[[I_GRID]]$cache_inner_mode(
                        mode = innerMethods$get_inner_mode(atOuterMode = 1), idx = i)
                    inner_grid_cache_nfl[[I_GRID]]$cache_inner_negHessChol(
                        negHessChol = innerMethods$get_inner_cholesky(atOuterMode = 1), idx = i)
                    inner_grid_cache_nfl[[I_GRID]]$cache_weights(weight = wgt, idx = i)
                    ans <- ans + wgt
                } else {
                    wgt <- paramGrid$weights(idx = i)[1]
                    node <- paramGrid$nodes(idx = i)[1, ]

                    ## Skew the CCD values:
                    if (skew) {
                        for (d in 1:nParamTrans) {
                            node[d] <- node[d] * skewedStdDev[d, step(node[d]) + 1]  ## negative skew column 1, positive skew column 2
                        }
                    }
                    ## Transform to paramTrans scale:
                    node <- z_to_paramTrans(node, paramTransMode, Atransform, transformMethod)
                    paramTransLogPostDens <- innerMethods$calcLogDens_pTransformed(node)
                    wgt_dens <- wgt * exp(paramTransLogPostDens - logPostProbMode)
                    ## Marginal sum:
                    ans <- ans + wgt_dens

                    ## Cache everything for simulation:
                    inner_grid_cache_nfl[[I_GRID]]$cache_inner_mode(
                        mode = innerMethods$get_inner_mode(atOuterMode = 0), idx = i)
                    inner_grid_cache_nfl[[I_GRID]]$cache_inner_negHessChol(
                        negHessChol = innerMethods$get_inner_cholesky(atOuterMode = 0), idx = i)
                    inner_grid_cache_nfl[[I_GRID]]$cache_weights(weight = wgt_dens, idx = i)
                }
                ## *** Add a convergence check?
            }
            nimCat("\n")
            if (skew) adjLogWgt <- logSkewedWgt else adjLogWgt <- 0

            ## Marginal log posterior density, a normalizing constant for other
            ## methods.
            marginalPostDensity[I_GRID] <<- log(ans) + logPostProbMode - 0.5 * logDetNegHessParamTrans +
                adjLogWgt

            paramGridCached[I_GRID] <<- TRUE
        },
        ## Quadrature based marginal log-likelihood
        ## Probably not particularly accurate for CCD.
        calcMarginalLogLikQuad = function() {
            if (I_GRID == I_CCD)
                print("  [Note]: Estimating marginal log-likelihood based on CCD grid.\n           Estimation based on an AGHQ grid may be more accurate (but more computationally expensive).")
            if(!paramGridCached[I_GRID])
                calcParamGrid()
            returnType(double())
            return(marginalPostDensity[I_GRID])
        },
        ## Marginals AGHQ from Stringer et al.
        calcMarginalParamQuad = function(pIndex = integer(),
                                                nPts = integer(0, default = 5),
                                                nQuad = integer(0, default = 3),
                                                gridTransformMethod = character(0, default = "spectral"),
                                                quadRule = character(0, default = "NULL"),
                                                prune = double(0, default = -1)) {
            one_time_fixes()
                                  
            ## Build the quadrature grid points:
            if (dim(paramTrans1_nodes)[1] != nPts) paramTrans1_nodes <<- quadGH(levels = nPts, type = "GHe")

            if(nParamTrans == 1) {
                ## d-1 dimensional AGHQ not needed; just evaluate inner approximation.
                res <- matrix(0, nrow = nPts, ncol = 2)
                stdDev <- sqrt(covParamTrans[1, 1])
                for(i in 1:nPts) {
                    res[i, 1] <- paramTrans1_nodes[i, 2] * stdDev + paramTransMode[pIndex]
                    res[i, 2] <- innerMethods$calcLogDens_pTransformed(c(res[i, 1]))
                }
                return(res)
            }

            if(nQuad %% 2 == 0)
                cat("  [Note] For computational efficiency, it is recommended to use an odd number of quadrature points\n         (via argument `nQuad`) for marginalizing over the parameter (outer) grid.\n")

            ## Grid for additional paramTrans.
            ## Build marginal AGHQ grid to compute the hyperparameter marginals
            ## (integrate over pT-1 paramTrans values).      
            if( pruneMargGrid != -1)
              pruneMargGrid <<- prune
            paramMargGrid$buildGrid(method = quadRule, nQuad = nQuad, prune = pruneMargGrid)  

            nQuadGrid <- paramMargGrid$gridSize()

            if (!modeCached) findMode(rep(Inf, npar), hessian = TRUE, parscale = "transformed")  ## *** default is now nlminb
            
            ## 1D quadrature to evaluate the paramTrans on.
            stdDev <- sqrt(covParamTrans[pIndex, pIndex])

            ## Initialize optimization at paramTrans mode.
            Atransform_i <- matrix(0, nrow = nParamTrans - 1, ncol = nParamTrans - 1)

            ## Column 1 is chosen paramTrans values, Column 2 is marginalized values, normalized based on AGHQ.
            res <- matrix(0, nrow = nPts, ncol = 2)
            paramTrans_j <- paramTransMode
            other_paramTrans_indices <- paramTrans_indices[paramTrans_indices != pIndex]

            ## For each value of paramTrans_i, we need to do AGHQ which means finding the
            ## mode of the other parameters, transforming and computing.
            nimCat("Calculating inner AGHQ/Laplace approximation at (", nPts, ") marginal points with ", nQuadGrid, " quadrature grid points (one dot per grid point): ")
            for (i in 1:nPts) {
                res[i, 1] <- paramTrans1_nodes[i, 2] * stdDev + paramTransMode[pIndex]
                paramTrans_j[pIndex] <- res[i, 1]

                ## If this is the mode then we know optim already:
                if (paramTrans1_nodes[i, 2] == 0) {
                    paramTrans_iMode <- paramTransMode[other_paramTrans_indices]
                    subsetNegHess <- paramTransNegHess[other_paramTrans_indices, other_paramTrans_indices]
                    maxLogDensity_i <- logPostProbMode
                } else {
                    optRes <- innerMethods$findMax_fixedp(pStartTransform = paramTransMode, pTransformIndex = pIndex,
                        pTransformValue = res[i, 1], includePrior = TRUE, includeJacobian = TRUE,
                        hessian = TRUE)
                    subsetNegHess <- -optRes$hessian
                    paramTrans_iMode <- optRes$par
                    maxLogDensity_i <- optRes$value
                }
                
                if (gridTransformMethod == "spectral") {
                    E <- eigen(subsetNegHess, symmetric = TRUE)
                    for (d in 1:(nParamTrans-1)) {
                        Atransform_i[, d] <- E$vectors[, d]/sqrt(E$values[d])
                    }
                    logDetNegHessParamTrans_i <- sum(log(E$values))
                } else {
                    Atransform_i <- chol(subsetNegHess)
                    logDetNegHessParamTrans_i <- 2 * sum(log(diag(Atransform_i)))
                }

                density_i <- 0
                nimCat("(", i, ")")
                for (j in 1:nQuadGrid) {
                    nimCat(".")
                    if (j != paramMargGrid$modeIndex()) {
                        nodej <- paramMargGrid$nodes(idx = j)[1, ]
                        paramTrans_tmp <- z_to_paramTrans(z = nodej, postMode = paramTrans_iMode, A = Atransform_i,
                                                method = gridTransformMethod)
                        paramTrans_j[other_paramTrans_indices] <- paramTrans_tmp
                        postLogDensij <- innerMethods$calcLogDens_pTransformed(pTransform = paramTrans_j)
                        density_i <- density_i + exp(postLogDensij - maxLogDensity_i) * paramMargGrid$weights(idx = j)[1]
                    } else {
                        density_i <- density_i + paramMargGrid$weights(idx = j)[1]
                    }
                }
                res[i, 2] <- log(density_i) + maxLogDensity_i - 0.5 * logDetNegHessParamTrans_i
            }
            nimCat("\n")
            ## Because paramTrans_i values are AGHQ, we can normalize to get the proper
            ## posterior density.  This lets us get the marginal posterior via spline
            ## without any more normalizing (but note that in `fitMarginalSpline`
            ## we do also normalize.
            ## Note that this is a 1-d quadrature,
            ## normalizing P(paramTrans_i,Y) to get P(Y) rather than the expensive
            ## calculation of denominator in (8) in Bilodeau et al.
            margi <- sum(exp(res[, 2] - logPostProbMode) * paramTrans1_nodes[, 1])
            lognormconst <- log(margi) + logPostProbMode + log(stdDev)
            res[, 2] <- res[, 2] - lognormconst
            ## *** Should I cache this?
            returnType(double(2))
            return(res)
        },
        ## This can't be until the CCD grid is build.
        ## so that we have covParamTrans.
        ## Should also ensure that if they plan to skew the grid that is also done.
        calcMarginalParamIntegFree = function(pIndex = integer()) {
            ## Error Trapping:
            if(pIndex <= 0 | pIndex > nParamTrans)
                stop("Transformed parameter index requested is larger than available.")
                
            ## Requires running `calcSkewedSD()` first.
            if (!skewedSDCached) calcSkewedSD()
            
            stdDev <- sqrt(covParamTrans[pIndex, pIndex])
            paramTrans_i <- numeric(value = 0, length = nParamTrans)
            setTransformations(transformMethod)
            for (i in 1:nzMargGrid) {
                                        # Known fixed # of points
                paramTrans_i[pIndex] <- paramTransMode[pIndex] + zMargGrid[i] * stdDev
                marg_paramTrans[pIndex, i, 1] <<- paramTrans_i[pIndex]
                ## Find the conditional mean:
                for (j in 1:nParamTrans) {
                    if (j != pIndex) {
                        paramTrans_i[j] <- paramTransMode[j] + covParamTrans[pIndex, j] / covParamTrans[pIndex, pIndex] *
                            (paramTrans_i[pIndex] - paramTransMode[pIndex])
                    }
                }
                ## Calculate asymmetric Gaussian:
                zi <- paramTrans_to_z(paramTrans_i, paramTransMode, AinverseTransform, transformMethod)

                ## logDens = sum log(exp(-z^2/sigma_(+/-))) *Not normalized.
                ## Can we normalize analytically? ***CJP?
                logDens <- 0
                for (j in 1:nParamTrans) {
                    side <- 2
                    if (zi[j] <= 0) side <- 1
                    logDens <- logDens - 0.5 * (zi[j]/skewedStdDev[j, side])^2
                }
                marg_paramTrans[pIndex, i, 2] <<- logDens
            }
            returnType(double(2))
            return(marg_paramTrans[pIndex, , ])
        },
        ## marginalTransformedSplineDensity = function(pIndex = integer()) {
        ## returnType(double(2))
        ## return(marginalSplineR(marg_paramTrans[pIndex, , 1], marg_paramTrans[pIndex, , 2]))
        ## },
        simulateLatents = function(n = integer()) {
            if(n < 0)
              stop("number `n` of simulated values must be at least one")
            if(I_GRID == I_AGHQSPARSE)
              stop("Sparse grids can have negative weights and are not valid for simulating the latent effects.")
            if (!paramGridCached[I_GRID]) calcParamGrid()

            sims <- inner_grid_cache_nfl[[I_GRID]]$simulate(n)
            returnType(double(2))
            return(sims)
        },
        ## Simulation method for paramTrans marginal on the skewed multivariate normal.
        simulateParams = function(n = integer()) {
            if(n < 0)
              stop("number `n` of simulated values must be at least one")

            sims <- matrix(0, nrow = n, ncol = nParamTrans)
            if (!skewedSDCached) calcSkewedSD()

            setTransformations(transformMethod)

            prob <- skewedStdDev[, 2]/(skewedStdDev[, 1] + skewedStdDev[, 2])

            for (i in 1:n) {
                ## simulate z on the base scale
                z <- abs(rnorm(nParamTrans, 0, 1))
                for (j in 1:nParamTrans) {
                    dir <- rbinom(1, 1, prob[j])  ## Skew z pos if 1, neg if 0.
                    if (dir == 1) { z[j] <- skewedStdDev[j, 2] * z[j]
                    } else z[j] <- -skewedStdDev[j,1] * z[j]
                }
                ## Scale it based on method
                sims[i, ] <- z_to_paramTrans(z, paramTransMode, Atransform, transformMethod)
            }
            returnType(double(2))
            return(sims)
        },
        getParamGrid = function() {
            return(paramGrid$nodes())
            returnType(double(2))
        }
    )
)
