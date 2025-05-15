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
        allGridRules <- c("CCD", "AGHQ", "AGHQSPARSE", "USER")

        ## Default outer grid to CCD unless low dimensional.
        paramGridRule <- extractControlElement(control, "paramGridRule", "none")
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
        
        ## Default to CCD
        trans_grid <- configureQuadGrid(d = 1, nQuad_ = nQuadOuter, quadRule = paramGridRule,
                                        control = list(quadRules = allGridRules))

        
        innerMethods <- buildAGHQ(model, nQuadInner, paramNodes, latentNodes, calcNodes,
                                  calcNodesOther, control)

        if(!identical(paramNodes, innerMethods$paramNodes))
            stop("`paramNodes` has unexpectedly changed. This should not have occurred.")

        trans_indices <- innerMethods$pTransform_indices
        
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
        ## On transformed scale are named "trans".
        paramNodes <- innerMethods$paramNodes
        nParam <- innerMethods$npar
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
        

        ## Indicator for removing the redundant index -1 in trans_indices
        one_time_fixes_done <- FALSE

        ## Default calculation method for AGHQuad
        computeMethod_ <- extractControlElement(control, "computeMethod", 2)
        useInnerCache_ <- extractControlElement(control, "useInnerCache", TRUE)

        ## For compilation have to set the dimension after setup.
        trans_grid$setDim(ndim = nParamTrans)
        inner_grid_cache_nfl <- nimbleFunctionList(INNER_CACHE_BASE)

        ## Make sure the grids match the trans_grid numbers.  Initialize them with
        ## nre = 0 in case they aren't used to not generate too much data.
        I_GRID <- trans_grid$I_RULE
        I_CCD <- trans_grid$I_CCD
        inner_grid_cache_nfl[[I_CCD]] <- inner_cache_methods(nre = 0, nGrid = 1,
            condIndepSets = lenInternalRENodeSets, nCondIndepSets = nInternalRESets)
        I_AGHQ <- trans_grid$I_AGHQ
        inner_grid_cache_nfl[[I_AGHQ]] <- inner_cache_methods(nre = 0, nGrid = 1,
            condIndepSets = lenInternalRENodeSets, nCondIndepSets = nInternalRESets)
        I_AGHQSPARSE <- trans_grid$I_AGHQSPARSE
        inner_grid_cache_nfl[[I_AGHQSPARSE]] <- inner_cache_methods(nre = 0, nGrid = 1,
            condIndepSets = lenInternalRENodeSets, nCondIndepSets = nInternalRESets)

        I_USER <- 1
        if (any(allGridRules == "USER")) {
            I_USER <- trans_grid$I_USER
            inner_grid_cache_nfl[[I_USER]] <- inner_cache_methods(nre = 0, nGrid = 1,
                condIndepSets = lenInternalRENodeSets, nCondIndepSets = nInternalRESets)
        }

        ## Store the quadrature sums for each grid:
        marginalPostDensity <- rep(-Inf, length(allGridRules))

        trans_marg_grid <- configureQuadGrid(d = nParamTrans - 1, nQuad_ = 1,
                                             quadRule = quadRuleMarginal)
        trans1_nodes <- matrix(0, nrow = 1, ncol = 2)

        ## Cached values for convenience: For marginal distributions in AGHQ over
        ## d-1.
        pTransformFix <- 0
        indexFix <- 0

        ## We will want to cache the standard deviation skew terms.  Default will
        ## be not to skew (e.g. 1)
        covTrans <- matrix(0, nrow = nParamTrans, ncol = nParamTrans)
        cholNegHess <- matrix(0, nParamTrans, nParamTrans)
        A_spectral <- matrix(0, nParamTrans, nParamTrans)
        Ainverse_spectral <- matrix(0, nParamTrans, nParamTrans)
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
        marg_P <- matrix(0, nrow = nzMargGrid, ncol = nParam)
        marg_trans <- array(0, c(nParamTrans, nzMargGrid, 2))
        ## ***Do we want the simulations of the latent effects to be cached or just
        ## returned? @CJP?
        post_sims <- matrix(0, nrow = 3, ncol = nre)

        ## Optim info:
        modeCached <- FALSE
        transMode <- numeric(nParamTrans)
        if (nParamTrans == 1) {
            transMode <- c(transMode, -1)
            trans_indices <- c(trans_indices, -1)
        }
        transNegHess <- matrix(0, nrow = nParamTrans, ncol = nParamTrans)
        logPostProbMode <- 0
        logDetNegHessTrans <- 0

        ## Other cached values:
        skewedSDCached <- FALSE
        ## Must be cached for each grid: Up to 3 currently.
        paramGridCached <- c(FALSE, FALSE, FALSE)

        ## Indicator for removing the redundant index -1 in trans_indices
        one_time_fixes_done <- FALSE
    },
    run = function() {},
    methods = list(
        one_time_fixes = function() {
            if (one_time_fixes_done) return()
            if (nParamTrans == 1) {
                trans_indices <<- numeric(length = 1, value = 1)
                transMode <<- numeric(length = 1, value = 0)
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
            transMode <<- optRes$par
            if(nParamTrans == 1) {
                transNegHess <<- matrix(-optRes$hessian, 1, 1)
            } else transNegHess <<- -optRes$hessian
            logPostProbMode <<- optRes$value
            covTrans <<- inverse(transNegHess)
            return(optRes)
            returnType(optimResultNimbleList())
        },
        buildParamGrid = function(quadRule = character(0, default = "NULL"),
                                  nQuadUpdate = integer(0, default = -1)) {
            one_time_fixes()
            if(nQuadUpdate != -1)
                nQuadOuter <<- nQuadUpdate
            if(quadRule != "NULL")
                setParamGridRule(quadRule)
            trans_grid$buildGrid(method = paramGridRule, nQuad = nQuadOuter)
            nGrid <- trans_grid$gridSize()
            inner_grid_cache_nfl[[I_GRID]]$buildCache(nGridUpdate = nGrid, nLatentNodes = nre)
            if (!modeCached) findMode(rep(Inf, nParam), hessian = TRUE, parscale = "transformed")
        },
        setParamGridRule = function(quadRule = character(0, default = "AGHQ")) {
            ## Add a rule check here to make sure it's valid.
            paramGridRule <<- quadRule
            ## Default to AGHQ and change if requested.
            I_GRID <<- I_AGHQ
            if (quadRule == "CCD") I_GRID <<- I_CCD
            if (quadRule == "AGHQSPARSE") I_GRID <<- I_AGHQSPARSE
        },
        calcEigen = function() {
            E <- eigen(transNegHess, symmetric = TRUE)  ## Should be symmetric...
            for (d in 1:nParamTrans) {
                A_spectral[, d] <<- E$vectors[, d]/sqrt(E$values[d])
                AinverseTrans_spectral[, d] <<- E$vectors[, d] * sqrt(E$values[d]) 
            }
            logDetNegHessTrans <<- sum(log(E$values))
            eigenCached <<- TRUE
        },
        calcCholesky = function() {
            cholNegHess <<- chol(transNegHess)
            logDetNegHessTrans <<- 2 * sum(log(diag(cholNegHess)))
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
                Atransform <<- cholNegHess # Used with backsolve in `z_to_trans`.
                AinverseTransform <<- cholNegHess
            }
        },
        ## Transform from standard (z) to param transform (trans) scale.
        z_to_trans = function(z = double(1), postMode = double(1), A = double(2),
                              method = character(0, default = "spectral")) {
            if (method == "spectral") {
                d <- dim(z)[1]
                trans <- numeric(value = 0, length = d)
                for (i in 1:d) {  # A %*% z
                    trans[i] <- postMode[i] + sum(A[i,] * z) 
                }
            } else {
                trans <- postMode + backsolve(A, z)
            }
            returnType(double(1))
            return(trans)
        },
        ## Transform from param transform (trans) to standard (z) scale.
        trans_to_z = function(trans = double(1), postMode = double(1), A = double(2),
                              method = character(0, default = "spectral")) {
            if (method == "spectral") {
                ## 'A' provided needs to be transpose of inverse of true A.
                d <- dim(trans)[1]
                z <- numeric(value = 0, length = d)
                trans_mean <- trans - postMode
                for (i in 1:d) {  # t(A) %*% trans_mean
                    z[i] <- sum(A[,i] * trans_mean)
                }
            } else {
                z <- (A %*% (trans - postMode))[, 1]
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
                trans <- z_to_trans(z, transMode, Atransform, transformMethod)
                logDens2Neg <- innerMethods$calcLogDens_pTransformed(pTransform = trans)
                skewedStdDev[i, 1] <<- sqrt(2/(2 * (logPostProbMode - logDens2Neg)))  ## numerator (-sqrt(2)) ^2
                z[i] <- sqrt(2)
                trans <- z_to_trans(z, transMode, Atransform, transformMethod)
                logDens2Pos <- innerMethods$calcLogDens_pTransformed(pTransform = trans)
                skewedStdDev[i, 2] <<- sqrt(2/(2 * (logPostProbMode - logDens2Pos)))  ## numerator (-sqrt(2)) ^2
                logSkewedWgt <<- logSkewedWgt + log(sum(skewedStdDev[i, ]/2))
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
            ## 0.5*(logDetNegHessTrans) - sum(log(skewedStdDev[,1] * skewedStdDev[,2]))
            ## *** What Paul thinks it should be. ***
            marg <- logPostProbMode + 0.5 * nParamTrans * log(2 * pi) - 0.5 * (logDetNegHessTrans) +
                logSkewedWgt  # sum(log((skewedStdDev[,1] + skewedStdDev[,2])/2))
            returnType(double())
            return(marg)
        },
        ## This is the meat and potatoes for being able to make inference on the latent nodes.
        ## Calculate trans on the quadrature grid points. AGHQ or CCD.
        ## Stores all values we need for simulation inference on the latent nodes.
        calcParamGrid = function(skew = logical(0, default = TRUE)) {
            buildParamGrid()
            setTransformations(transformMethod)
            nGrid <- trans_grid$gridSize()

            if (!skewedSDCached & skew) calcSkewedSD()
            ans <- 0
            ## Now fill in the grid values.
            nimCat("Calculating inner AGHQ/Laplace approximation at ", nGrid, " outer (parameter) grid points (one dot per point): ")
            for (i in 1:nGrid) {
                nimCat(".")
                ## Operations at the mode:
                if (i == trans_grid$modeI()) {
                    wgt <- trans_grid$weights(indx = i)[1]
                    inner_grid_cache_nfl[[I_GRID]]$cache_inner_mode(
                        mode = innerMethods$get_inner_mode(atOuterMode = 1), indx = i)
                    inner_grid_cache_nfl[[I_GRID]]$cache_inner_negHessChol(
                        negHessChol = innerMethods$get_inner_cholesky(atOuterMode = 1), indx = i)
                    inner_grid_cache_nfl[[I_GRID]]$cache_weights(weight = wgt, indx = i)
                    ans <- ans + wgt
                } else {
                    wgt <- trans_grid$weights(indx = i)[1]
                    node <- trans_grid$nodes(indx = i)[1, ]

                    ## Skew the CCD values:
                    if (skew) {
                        for (d in 1:nParamTrans) {
                            node[d] <- node[d] * skewedStdDev[d, step(node[d]) + 1]  ## negative skew column 1, positive skew column 2
                        }
                    }
                    ## Transform to trans scale:
                    node <- z_to_trans(node, transMode, Atransform, transformMethod)
                    transLogPostDens <- innerMethods$calcLogDens_pTransformed(node)
                    wgt_dens <- wgt * exp(transLogPostDens - logPostProbMode)
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
            marginalPostDensity[I_GRID] <<- log(ans) + logPostProbMode - 0.5 * logDetNegHessTrans +
                adjLogWgt

            paramGridCached[I_GRID] <<- TRUE
        },
        ## Quadrature based marginal log-likelihood
        ## Probably not particularly accurate for CCD.
        calcMarginalLogLikQuad = function() {
            if (I_GRID == I_CCD)
                print("  [Note]: Estimating marginal log-likelihood based on CCD grid. Estimation based on an AGHQ grid may be more accurate (but more computationally expensive).")
            if(!paramGridCached[I_GRID])
                calcParamGrid()
            returnType(double())
            return(marginalPostDensity[I_GRID])
        },
        ## Marginals AGHQ from Stringer et al.
        calcMarginalParamQuad = function(pIndex = integer(),
                                                nPts = integer(0, default = 5),
                                                nQuad = integer(0, default = 3),
                                                gridTransformMethod = character(0, default = "spectral")) {
            one_time_fixes()
            ## Build the quadrature grid points:
            if (dim(trans1_nodes)[1] != nPts) trans1_nodes <<- AGHQ1D(nQuad = nPts)

            if(nParamTrans == 1) {
                ## d-1 dimensional AGHQ not needed; just evaluate inner approximation.
                res <- matrix(0, nrow = nPts, ncol = 2)
                stdDev <- sqrt(covTrans[1, 1])
                for(i in 1:nPts) {
                    res[i, 1] <- trans1_nodes[i, 2] * stdDev + transMode[pIndex]
                    res[i, 2] <- innerMethods$calcLogDens_pTransformed(c(res[i, 1]))
                }
                return(res)
            }

            if(nQuad %% 2 == 0)
                cat("  [Note] For computational efficiency, it is recommended to use an odd number of quadrature points\n         (via argument `nQuad`) for marginalizing over the parameter (outer) grid.\n")

            ## Grid for additional trans.
            ## Build marginal AGHQ grid to compute the hyperparameter marginals
            ## (integrate over pT-1 trans values).            
            trans_marg_grid$buildGrid(nQuad = nQuad)  

            nQuadGrid <- trans_marg_grid$gridSize()

            if (!modeCached) findMode(rep(Inf, nParam), hessian = TRUE, parscale = "transformed")  ## *** default is now nlminb
            
            ## 1D quadrature to evaluate the trans on.
            stdDev <- sqrt(covTrans[pIndex, pIndex])

            ## Initialize optimization at trans mode.
            Atransform_i <- matrix(0, nrow = nParamTrans - 1, ncol = nParamTrans - 1)

            ## Column 1 is chosen trans values, Column 2 is marginalized values, normalized based on AGHQ.
            res <- matrix(0, nrow = nPts, ncol = 2)
            transj <- transMode
            other_trans_indices <- trans_indices[trans_indices != pIndex]

            ## For each value of transi, we need to do AGHQ which means finding the
            ## mode of the other parameters, transforming and computing.
            nimCat("Calculating inner AGHQ/Laplace approximation at (", nPts, ") marginal points with ", nQuadGrid, " quadrature grid points (one dot per grid point): ")
            for (i in 1:nPts) {
                res[i, 1] <- trans1_nodes[i, 2] * stdDev + transMode[pIndex]
                transj[pIndex] <- res[i, 1]

                ## If this is the mode then we know optim already:
                if (trans1_nodes[i, 2] == 0) {
                    trans_iMode <- transMode[other_trans_indices]
                    subsetNegHess <- transNegHess[other_trans_indices, other_trans_indices]
                    maxLogDensity_i <- logPostProbMode
                } else {
                    optRes <- innerMethods$findMax_fixedp(pStartTransform = transMode, pTransformIndex = pIndex,
                        pTransformValue = res[i, 1], includePrior = TRUE, includeJacobian = TRUE,
                        hessian = TRUE)
                    subsetNegHess <- -optRes$hessian
                    trans_iMode <- optRes$par
                    maxLogDensity_i <- optRes$value
                }
                
                if (gridTransformMethod == "spectral") {
                    E <- eigen(subsetNegHess, symmetric = TRUE)
                    for (d in 1:(nParamTrans-1)) {
                        Atransform_i[, d] <- E$vectors[, d]/sqrt(E$values[d])
                    }
                    logDetNegHessTrans_i <- sum(log(E$values))
                } else {
                    Atransform_i <- chol(subsetNegHess)
                    logDetNegHessTrans_i <- 2 * sum(log(diag(Atransform_i)))
                }

                density_i <- 0
                nimCat("(", i, ")")
                for (j in 1:nQuadGrid) {
                    nimCat(".")
                    if (j != trans_marg_grid$modeI()) {
                        nodej <- trans_marg_grid$nodes(indx = j)[1, ]
                        trans_tmp <- z_to_trans(z = nodej, postMode = trans_iMode, A = Atransform_i,
                                                method = gridTransformMethod)
                        transj[other_trans_indices] <- trans_tmp
                        postLogDensij <- innerMethods$calcLogDens_pTransformed(pTransform = transj)
                        density_i <- density_i + exp(postLogDensij - maxLogDensity_i) * trans_marg_grid$weights(indx = j)[1]
                    } else {
                        density_i <- density_i + trans_marg_grid$weights(indx = j)[1]
                    }
                }
                res[i, 2] <- log(density_i) + maxLogDensity_i - 0.5 * logDetNegHessTrans_i
            }
            nimCat("\n")
            ## Because transi values are AGHQ, we can normalize to get the proper
            ## posterior density.  This lets us get the marginal posterior via spline
            ## without any more normalizing (but note that in `fitMarginalSpline`
            ## we do also normalize.
            ## Note that this is a 1-d quadrature,
            ## normalizing P(transi,Y) to get P(Y) rather than the expensive
            ## calculation of denominator in (8) in Bilodeau et al.
            margi <- sum(exp(res[, 2] - logPostProbMode) * trans1_nodes[, 1])
            lognormconst <- log(margi) + logPostProbMode + log(stdDev)
            res[, 2] <- res[, 2] - lognormconst
            ## *** Should I cache this?
            returnType(double(2))
            return(res)
        },
        ## This can't be until I've built the CCD grid.
        ## so that we have covTrans.
        ## Should also ensure that if they plan to skew the grid that is also done.
        calcMarginalParamIntegFree = function(pIndex = integer()) {
            ## Requires running `calcSkewedSD()` first.
            if (!skewedSDCached) calcSkewedSD()

            stdDev <- sqrt(covTrans[pIndex, pIndex])
            transi <- numeric(value = 0, length = nParamTrans)
            setTransformations(transformMethod)
            for (i in 1:nzMargGrid) {
                                        # Known fixed # of points
                transi[pIndex] <- transMode[pIndex] + zMargGrid[i] * stdDev
                marg_trans[pIndex, i, 1] <<- transi[pIndex]
                ## Find the conditional mean:
                for (j in 1:nParamTrans) {
                    if (j != pIndex) {
                        transi[j] <- transMode[j] + covTrans[pIndex, j] / covTrans[pIndex, pIndex] *
                            (transi[pIndex] - transMode[pIndex])
                    }
                }
                ## Calculate asymmetric Gaussian:
                zi <- trans_to_z(transi, transMode, AinverseTransform, transformMethod)

                ## logDens = sum log(exp(-z^2/sigma_(+/-))) *Not normalized.
                ## Can we normalize analytically? ***CJP?
                logDens <- 0
                for (j in 1:nParamTrans) {
                    side <- 2
                    if (zi[j] <= 0) side <- 1
                    logDens <- logDens - 0.5 * (zi[j]/skewedStdDev[j, side])^2
                }
                marg_trans[pIndex, i, 2] <<- logDens
            }
            returnType(double(2))
            return(marg_trans[pIndex, , ])
        },
        ## marginalTransformedSplineDensity = function(pIndex = integer()) {
        ## returnType(double(2))
        ## return(marginalSplineR(marg_trans[pIndex, , 1], marg_trans[pIndex, , 2]))
        ## },
        simulateLatents = function(n = integer()) {
            if (!paramGridCached[I_GRID]) calcParamGrid()

            sims <- inner_grid_cache_nfl[[I_GRID]]$simulate(n)
            returnType(double(2))
            return(sims)
        },
        ## Simulation method for trans marginal on the skewed multivariate normal.
         simulateParams = function(n = integer()) {
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
                sims[i, ] <- z_to_trans(z, transMode, Atransform, transformMethod)
            }
            returnType(double(2))
            return(sims)
        },
        getParamGrid = function() {
            return(trans_grid$nodes())
            returnType(double(2))
        }
    )
)
