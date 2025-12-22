#' Build Nested Bayesian Approximation Using Quadrature-based Methods
#'
#' Build a nested approximation for a given NIMBLE model, providing parameter estimation and 
#' sampling for latent nodes. The approximation uses an inner Laplace (or AGHQ) approximation 
#' to approximately marginalize over latent nodes and an outer quadrature grid on the
#' parameters.
#'
#' @param model a NIMBLE model object created by calling \code{nimbleModel}. 
#' The model must have automatic derivatives (AD) turned on, e.g., by using
#'   \code{buildDerivs=TRUE} in \code{nimbleModel}.
#' @param paramNodes optional character vector of (hyper)parameter nodes in the model. If missing, this
#' will be the stochastic non-data nodes with no parent stochastic nodes.
#' @param latentNodes optional character vector of latent nodes (e.g., random and fixed
#' effects) in the model. If missing this will be the stochastic non-data nodes 
#' that are not determined to be parameter nodes.
#' @param calcNodes optional character vector of names of nodes for calculating the
#'   integrand for Laplace/AGHQ approximation over the latent nodes; defaults are provided by
#'   \code{\link[nimble]{setupMargNodes}}. Note that users will generally not need to provide this. 
#'   There may be deterministic nodes between
#'   \code{paramNodes} and \code{calcNodes}. These will be included in
#'   calculations automatically and thus do not need to be included in
#'   \code{calcNodes} (but there is no problem if they are).
#' @param calcNodesOther optional character vector of names of nodes for calculating
#'   terms in the log-likelihood that do not depend on any
#'   \code{randomEffectsNodes}, and thus are not part of the marginalization,
#'   but should be included for purposes of finding the approximation. 
#'   Note that users will generally not need to provide this. This defaults to
#'   stochastic nodes that depend on \code{paramNodes} but are not part of and
#'   do not depend on \code{latentNodes}. There may be deterministic
#'   nodes between \code{paramNodes} and \code{calcNodesOther}. These will be
#'   included in calculations automatically and thus do not need to be included
#'   in \code{calcNodesOther} (but there is no problem if they are).
#' @param control a named list for providing additional settings for the approximation,
#' including settings for the inner Laplace/AGHQ approximation. See \code{control} section below. 
#'
#' @details
#' 
#' This function builds a nested Bayesian approximation for the provided model. NIMBLE's nested
#' approximation provides approximate posterior inference using methodology similar to the 
#' well-known INLA approach (Rue et al. 2009), implemented in the 
#' R-INLA package and to the related methods for extended Gaussian latent models (EGLMs)
#' of Stringer et al. (2023), implemented in the `aghq` R package. For more details on the 
#' nested approximation algorithms, see the NIMBLE User Manual.
#' 
#' Unlike Laplace/AGHQ approximation, the nested approximation is Bayesian, requiring
#' prior distributions for all parameters and providing functionality to estimate 
#' the marginal posterior distributions of the individual parameters and to sample from the 
#' marginal joint posterior distribution of the latent nodes. Similarly to Laplace,
#' the nested approximation uses Laplace/AGHQ to approximately marginalize over the latent nodes.
#' However, instead of then maximizing the approximate marginalized posterior density,
#' the nested approximation uses a quadrature grid on the parameters to perform approximate Bayesian
#' inference on the parameters and latent nodes.
#' 
#' The recommended way to use the nested approximation once it is built is to call
#' \code{runNestedApprox} on the returned object, and then to call additional approximation
#' functions on the output of \code{runNestedApprox} as needed. For details see \code{runNestedApprox}.
#' However, for more granular control, one can also call the internal methods of the nested approximation,
#' discussed briefly below.
#' 
#' In general, the theory that underpins these approximations assumes that the latent nodes (fixed and random effects) 
#' are Gaussian. \code{buildNestedApprox} makes no such assumptions, allowing the user to extend these approximations to any imaginable
#' set of models in NIMBLE. However, the accuracy of the approximation is then not supported theoretically, and it is up to the user to 
#' determine whether or not the posterior approximation is valid.
#'
#' @section Computational considerations:
#' 
#' The computational cost of the nested approximation can vary depending on what 
#' nodes are considered as parameter nodes and what nodes as latent nodes, as well as
#' by the number of quadrature points (for both the latent and parameter nodes) and 
#' type of grid used for the parameter nodes. Some details are provided in the User Manual.
#' 
#' \code{buildNestedApprox} will by default (unless changed
#' manually by specifying sets of nodes) determine from the model which latent nodes
#' can be integrated over (marginalized) independently. For example, in a GLMM
#' with a grouping factor and an independent random effect intercept for each
#' group (and no fixed effects), the random effects can be marginalized as a set of univariate
#' approximations rather than one multivariate approximation. On the other hand,
#' correlated or nested random effects would require multivariate marginalization,
#' as would the presence of fixed effects (since they affect all the observations).
#' Independent marginalizations result in lower-dimensional calculations (essentially 
#' exploiting sparsity in the covariance structure of the latent nodes) and therefore
#' improve computational efficiency. Note that at this time, the nested approximation
#' cannot otherwise take advantage of sparsity in the covariance structure of the latent nodes.
#' 
#' @section How input nodes are processed:
#' 
#' In many cases, the selection of parameter and latent nodes will be handled automatically in a 
#' reasonable fashion. However, random effects can be
#' written in models in multiple equivalent ways, and customized use cases may
#' call for integrating over chosen parts of a model. Hence, one can take full
#' charge of how different parts of the model will be used, specifying explicitly the
#' \code{paramNodes} and \code{latentNodes}. The User Manual provides more details on situations
#' in which one may want to specify these nodes explicitly.
#'
#' Any of the input node vectors, when provided, will be processed using
#'   \code{nodes <- model$expandNodeNames(nodes)}, where \code{nodes} may be
#'   \code{paramNodes}, \code{latentNodes}, and so on. This step allows
#'   any of the inputs to include node-name-like syntax that might contain
#'   multiple nodes. For example, \code{paramNodes = 'beta[1:10]'} can be
#'   provided if there are actually 10 scalar parameters, 'beta[1]' through
#'   'beta[10]'. The actual node names in the model will be determined by the
#'   \code{expandNodeNames} step.
#'
#' In many (but not all) cases, one only needs to provide a NIMBLE model object
#'   and then the function will construct reasonable defaults necessary for
#'   Laplace approximation to marginalize over all continuous latent nodes
#'   (both random and fixed effects) in a model. 
#'
#' \code{buildNestedApprox} uses \code{\link[nimble]{setupMargNodes}} (in a multi-step process)
#'   to try to give sensible defaults from
#'   any combination of \code{paramNodes}, \code{latentNodes},
#'   \code{calcNodes}, and \code{calcNodesOther} that are provided. 
#'
#' \code{\link[nimble]{setupMargNodes}} also determines which integration dimensions are
#' conditionally independent, i.e., which can be done separately from each
#' other. For example, when possible, 10 univariate random effects will be split
#' into 10 univariate integration problems rather than one 10-dimensional
#' integration problem. Note that models that include fixed effects as latent
#' nodes often prevent this splitting into conditionally independent sets.
#'
#' The defaults make general assumptions such as that
#'   \code{latentNodes} have \code{paramNodes} as parents or (for fixed effects)
#' are also components of a linear predictor expression. However, the
#'   steps for determining defaults are not simple, and it is possible that they
#'   will be refined in the future. It is also possible that they simply don't
#'   give what you want for a particular model. One example where they will not
#'   give desired results can occur when random effects have no prior
#'   parameters, such as `N(0,1)` nodes that will be multiplied by a scale
#'   factor (e.g., `sigma``) and added to other explanatory terms in a model. Such
#'   nodes look like top-level parameters in terms of model structure, so
#'   you must provide a \code{latentNodes} argument to indicate which
#'   they are.
#'
#' 
#' @section \code{control} list arguments:
#'
#' The \code{control} list allows additional settings to be made using named
#' elements of the list. Any elements in addition to those below are passed along
#' as the \code{control} list for the inner Laplace/AGHQ approximation (see \code{buildAGHQ}).
#' Below, `d` refers to the number of parameter nodes.
#' Supported elements include:
#'
#' \itemize{
#'     \item \code{nQuadOuter}. Number of outer quadrature points in each dimension (for parameter nodes). 
#'           Default is 3 for d > 1, 5 for d = 1. Not used with CCD grid.
#'     \item \code{nQuadInner}. Number of inner quadrature points in each dimension (for latent nodes). 
#'           Default is 1, corresponding to Laplace approximation. 
#'     \item \code{paramGridRule}. Quadrature rule for the parameter grid. Defaults to \code{"CCD"} for
#'          d > 2 and to \code{"AGHQ"} otherwise. Can also be \code{"AGHQSPARSE"} or (for user-defined grids)
#'          a user-defined nimbleFunction generator (created by calling `nimbleFunction`) with an appropriate
#'          `buildGrid` method that has arguments \code{levels} and \code{d} and that returns a matrix.
#'     \item \code{paramGridRule_userType}. If \code{paramGridRule} is a user-defined rule, this optional
#'          element can be used to indicate that the provided rule constructs a univariate rule rather
#'          than directly constructing a multivariate rule and that a multivariate rule should be constructed
#'          from the univariate rule as either a product rule (by specifying "PRODUCT") or a sparse rule
#'          (by specifying "SPARSE").
#'     \item \code{innerOptimWarning}. Whether to show inner optimization warnings. Default is \code{FALSE}.
#'     \item \code{marginalGridRule}. Rule for the grid for parameter marginalization. Default is \code{"AGHQ"}.
#'          Can also be \code{"AGHQSPARSE"}. At present, user-defined grids are not allowed.
#'     \item \code{marginalGridPrune}. Pruning parameter for marginal grid. Default is 0, corresponding to no pruning.
#'     \item \code{quadTransform}. Quadrature transformation method. Default is \code{"spectral"}, with
#'           \code{"cholesky"} as the other option.
#'   }
#' 
#' @section Parameter transformations used internally:
#'
#' If any \code{paramNodes} (parameters) or \code{latentNodes} have constraints on the range of valid values
#'   (because of the distribution they follow), they will be used on a
#'   transformed scale determined by \code{parameterTransform}. This means that internally the
#'   Laplace/AGHQ approximation itself will be done on the transformed scale for
#'   latent nodes and that the grid-based computation on the parameters will be done on the transformed scale.  
#' 
#' @section Available methods for advanced/development use:
#' 
#' Additional methods to access or control the Laplace/AGHQ
#' approximation directly (as an alternative to the recommended use of \code{runNestedApprox}) 
#' include the following, described only briefly:
#'
#' \itemize{
#' \item \code{findMode}: Find the posterior mode for hyperparameters.
#' \item \code{buildParamGrid}: Build the parameter grid using specified quadrature rule and settings.
#' \item \code{setParamGridRule}: Set the quadrature rule for the parameter grid (AGHQ, CCD, USER, AGHQSPARSE).
#' \item \code{calcEigen}: Calculate eigendecomposition of the negative Hessian for spectral transformations. 
#' \item \code{calcCholesky}: Calculate Cholesky decomposition of the negative Hessian for Cholesky transformations.
#' \item \code{setTransformations}: Set transformation method between spectral and Cholesky approaches.
#' \item \code{z_to_paramTrans}: Transform from standard (z) scale to parameter transform scale.
#' \item \code{paramTrans_to_z}: Transform from parameter transform scale to standard (z) scale.
#' \item \code{calcSkewedSD}: Calculate skewed standard deviations for asymmetric Gaussian approximations.
#' \item \code{getSkewedStdDev}: Retrieve the calculated skewed standard deviations.
#' \item \code{calcMarginalLogLikApprox}: Calculate INLA-like approximation of marginal log-likelihood
#'             using the approximate Gaussian.
#' \item \code{calcParamGrid}: Calculate inner approximation at parameter grid values and cache results.
#'        Required only for latent node simulation and quadrature-based marginal log-likelihood.
#' \item \code{calcMarginalLogLikQuad}: Calculate quadrature-based marginal log-likelihood.
#' \item \code{calcMarginalParamQuad}: Calculate univariate marginal parameter distributions using selected quadrature rule.
#' \item \code{calcMarginalParamIntegFree}: Calculate univariate marginal parameter distributions using INLA-like
#'         integration-free method based on approximate Gaussian approximations.
#' \item \code{simulateLatents}: Simulate from the posterior distribution of (transformed) latent nodes.
#' \item \code{simulateParams}: Simulate from the marginal posterior of (transformed) parameters using skewed normal.
#' \item \code{getParamGrid}: Retrieve the parameter grid points on the transformed scale.
#' }
#' 
#' 
#' @aliases nestedApprox INLA nested
#' 
#' @export 
#' 
#' @author Paul van Dam-Bates, Christopher Paciorek, Perry de Valpine
#' 
#' @references 
#' 
#' Rue, H., Martino, S., and Chopin, N. (2009). Approximate Bayesian inference for 
#' latent Gaussian models by using integrated nested Laplace approximations. 
#' Journal of the Royal Statistical Society: Series B (Statistical Methodology), 
#' 71(2), 319-392.
#' 
#' Stringer, A., Brown, P., and Stafford, J. (2023). Fast, scalable approximations to 
#' posterior distributions in extended latent Gaussian models. 
#' Journal of Computational and Graphical Statistics, 32(1), 84-98.
#' 
#' @examples 
#' data(penicillin, package="faraway")
#' code <- nimbleCode({
#'     for(i in 1:n) {
#'         mu[i] <- inprod(b[1:nTreat], x[i, 1:nTreat]) + re[blend[i]]
#'         y[i] ~ dnorm(mu[i], tau = Tau)
#'     }
#'     # Priors corresponding simply to INLA defaults and not being recommended.
#'     # Instead consider uniform or half-t distributions on the standard deviation scale
#'     # or penalized complexity priors.
#'     Tau ~ dgamma(1, 5e-05)
#'     Tau_re ~ dgamma(1, 5e-05)
#'     for( i in 1:nTreat ){ b[i] ~ dnorm(0, tau = 0.001) }
#'     for( i in 1:nBlend ){ re[i] ~ dnorm(0, tau = Tau_re) }
#' })
#' X <- model.matrix(~treat, data = penicillin)
#' data = list(y = penicillin$yield)
#' constants = list(nTreat = 4, nBlend = 5, n = nrow(penicillin),
#'                  x = X, blend = as.numeric(penicillin$blend))
#' inits <- list(Tau = 1, Tau_re = 1, b = c(mean(data$y), rep(0,3)), re = rep(0,5))
#' 
#' model <- nimbleModel(code, data = data, constants = constants,
#'                  inits = inits, buildDerivs = TRUE)
#' approx <- buildNestedApprox(model = model)
#' 
buildNestedApprox <- nimbleFunction(
    name = "nestedApprox",
    setup = function(model, paramNodes, latentNodes, calcNodes, calcNodesOther, control = list()) {
        innerOptimWarning <- extractControlElement(control, "innerOptimWarning", FALSE)
        verbose <- ifTRUE(nimble::getNimbleOption('verbose'))

        nQuadLatent <- extractControlElement(control, "nQuadLatent", 1)
        quadRuleMarginal <- extractControlElement(control, "marginalGridRule", "AGHQ")
        pruneMargGrid <- extractControlElement(control, "marginalGridPrune", 0)

        transformMethod <- extractControlElement(control, "quadTransform", "spectral")

        ## Default starting value for nested approximation set here and passed to
        ## Laplace.  Zero makes sense but should test others on the AGHQ grid and
        ## see how they work.
        control$innerOptimStart <- extractControlElement(control, "innerOptimStart",
                                                         "zero")
        if(!is.Rmodel(model))
            stop("`model` must be an R model, created by calling `nimbleModel`")
        
        inferenceNodes <- model$getNodeNames(includeData = FALSE, stochOnly = TRUE)
        if(!missing(paramNodes) && !all(model$expandNodeNames(paramNodes) %in% inferenceNodes))
            stop("some elements of `paramNodes` do not have prior distributions")
        if(!missing(latentNodes) && !all(model$expandNodeNames(latentNodes) %in% inferenceNodes))
            stop("some elements of `latentNodes` do not have prior distributions")

        ## Allow users to provide only one of the two sets and assume the rest are in the other.
        ## If one really wants to exclude a node, one would have to explicitly set
        ## both `paramNodes` and `latentNodes`.
        if(missing(paramNodes) && !missing(latentNodes))
            paramNodes <- setdiff(inferenceNodes, model$expandNodeNames(latentNodes))
        if(missing(latentNodes) && !missing(paramNodes))
            latentNodes <- setdiff(inferenceNodes, model$expandNodeNames(paramNodes))

        margNodes <- splitLatents(model, paramNodes, latentNodes)
        paramNodes <- margNodes$paramNodes
        latentNodes <- margNodes$randomEffectsNodes
        if(!length(paramNodes))
            stop("No parameter nodes detected in model. Check the model structure or provide parameter nodes explicitly via `paramNodes`.")
        if(!length(latentNodes))
            stop("No latent nodes detected in model. Check the model structure or provide latent nodes explicitly via `latentNodes`. Note that this can occur in a model with random effects that do not depend on any (hyper)parameters.")

        ## We need `nParamTrans` now; this should be ok, as `paramNodes` shouldn't be changed
        ## by creation of `innerMethods`.
        paramsTransform <- parameterTransform(model, paramNodes, control = list(allowDeterm = FALSE))
        nParamTrans <- paramsTransform$getTransformedLength()

        nQuadParam <- extractControlElement(control, "nQuadParam", ifelse(nParamTrans == 1, 5, 3))
        
        ## Configure all grids before calling AGHQ to make sure it builds
        ## correctly.  DO NOT MOVE WHEN THIS IS CALLED
        allGridRules <- c("CCD", "AGHQ", "AGHQSPARSE")

        ## Default param (outer) grid to CCD unless low dimensional.
        paramGridRule <- extractControlElement(control, "paramGridRule", "none")
        pruneParamGrid <- extractControlElement(control, "paramGridPrune", 0)
        paramGridRule_userType <- extractControlElement(control, "paramGridRule_userType", "MULTI")

        if(is.character(paramGridRule) && paramGridRule == "none")
            paramGridRule <- ifelse(nParamTrans >= 3, "CCD", "AGHQ")
    

        paramGridRuleName <- paramGridRule
        if(is.function(paramGridRule)) {
            paramGridRuleName <- environment(paramGridRule)$name
        } 

        messageIfVerbose("Building nested posterior approximation for the following node sets:\n",
                         "  - parameter nodes: ", makeNodeString(paramNodes, model), "\n",
                         "  - latent nodes: ", makeNodeString(latentNodes, model), "\n",
                         "  with ", paramGridRuleName, " grid for the parameters and ", ifelse(nQuadLatent > 1, "AGHQ", "Laplace"), " approximation for the latent nodes.")
 
        if(length(intersect(latentNodes, paramNodes)))
            stop("some nodes appear in both the parameter and latent sets")
        if (nParamTrans > 20)
            messageIfVerbose("  [Warning] There is a large number of parameter node elements. Computation may be slow.")

        
        if(is.character(paramGridRule) && paramGridRule == "AGHQ" && nQuadParam %% 2 == 0)
            messageIfVerbose("  [Note] For computational efficiency, it is recommended to use an odd number of quadrature points\n         for the parameter (outer) grid (`nQuadParam`).")
        
        ## Default to CCD (in which case `nQuadParam` is ignored).
        paramGrid <- configureQuadGrid(d = 1, levels = nQuadParam, quadRule = paramGridRule,
                                       control = list(quadRules = allGridRules, userConstruction = paramGridRule_userType))

        if(is.function(paramGridRule)) {
            paramGridRule <- "USER"
            allGridRules <- c(allGridRules, paramGridRule)
        }
        
        innerMethods <- buildAGHQ(model, nQuadLatent, paramNodes, latentNodes, calcNodes,
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
                                           quadRule = quadRuleMarginal)
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

        paramGridSkewed <- TRUE

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
            if(verbose) nimCat("Finding posterior mode for parameter(s).\n")
            optRes <- innerMethods$optimize(pStart = pStart, includePrior = TRUE,
                                            includeJacobian = TRUE,
                                            hessian = hessian, parscale = parscale)
            dm <- dim(optRes$hessian)[1]
            if(dm != nParamTrans)
                stop("Posterior mode could not be found. Consider adjusting the control parameters for the optimization via the `control` argument of `buildNestedApprox`.")
            if(any_nan(c(optRes$hessian)))
                stop("While attempting to find posterior mode, invalid hessian calculated. Consider adjusting the control parameters for the optimization via the `control` argument of `buildNestedApprox`.")
            if(optRes$convergence != 0 & verbose)
                print("  [Warning] In optimization over parameter(s) to find the posterior mode as the\n",
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
                nQuadParam <<- nQuadUpdate
            if(quadRule != "NULL" )
              setParamGridRule(quadRule)
            if(prune != -1)
                pruneParamGrid <<- prune
            paramGrid$buildGrid(method = paramGridRule, nQuad = nQuadParam, prune = pruneParamGrid)
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
                    if(verbose) nimCat("  [Warning] Skewness in posterior of the (hyper)parameters in dimension ", i, " is large and\na potential sign of an issue for these approximations.\n")
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
            marg <- logPostProbMode + 0.5 * nParamTrans * log(2 * pi) - 0.5 * (logDetNegHessParamTrans) +
                logSkewedWgt 
            returnType(double())
            return(marg)
        },
        ## This is the meat and potatoes for being able to make inference on the latent nodes.
        ## Calculate paramTrans on the quadrature grid points. AGHQ or CCD.
        ## Stores all values we need for simulation inference on the latent nodes.
        calcParamGrid = function(skew = logical(0, default = TRUE)) {
            if(I_GRID == I_AGHQSPARSE & verbose)
                print("  [Note] Sparse grids cannot be used to simulate latent effects (the main reason to compute the posterior on the parameter grid).")

            buildParamGrid()
            setTransformations(transformMethod)
            nGrid <- paramGrid$gridSize()
            
            ## Check to make sure cache is working correctly:
            if(length(inner_grid_cache_nfl[[I_GRID]]$weights()) != nGrid)
                stop("The grid cache system does not match the quadrature grid being calculated.")
            
            if (!skewedSDCached & skew) calcSkewedSD()
            ans <- 0
            ## Now fill in the grid values.
            if(verbose) nimCat("Calculating inner AGHQ/Laplace approximation at ", nGrid, " parameter (outer)\n  grid points (one dot per point): ")
            for (i in 1:nGrid) {
                if(verbose) nimCat(".")
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
            if(verbose) nimCat("\n")
            if (skew) adjLogWgt <- logSkewedWgt else adjLogWgt <- 0

            ## Marginal log posterior density, a normalizing constant for other
            ## methods.
            marginalPostDensity[I_GRID] <<- log(ans) + logPostProbMode - 0.5 * logDetNegHessParamTrans +
                adjLogWgt

            paramGridCached[I_GRID] <<- TRUE
            if(skew)
                paramGridSkewed <<- TRUE else paramGridSkewed <<- FALSE
        },
        ## Quadrature based marginal log-likelihood
        ## Probably not particularly accurate for CCD.
        calcMarginalLogLikQuad = function() {
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

            if(pIndex <= 0 | pIndex > nParamTrans)
                stop("calcMarginalParamQuad: Transformed parameter index, `pIndex`, requested is invalid.")
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

            if(nQuad %% 2 == 0 & verbose)
                cat("  [Note] For computational efficiency, it is recommended to use an odd number of quadrature points\n         (via argument `nQuad`) for marginalizing over the parameter (outer) grid.\n")

            ## Grid for additional paramTrans.
            ## Build marginal AGHQ grid to compute the hyperparameter marginals
            ## (integrate over pT-1 paramTrans values).      
            if( pruneMargGrid != -1)
              pruneMargGrid <<- prune
            paramMargGrid$buildGrid(method = quadRule, nQuad = nQuad, prune = pruneMargGrid)  

            nQuadGrid <- paramMargGrid$gridSize()

            if (!modeCached) findMode(rep(Inf, nParamTrans), hessian = TRUE, parscale = "transformed")  ## *** default is now nlminb
            
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
            if(verbose) nimCat("  - calculating inner AGHQ/Laplace approximation at (", nPts, ") marginal points\n    with ", nQuadGrid, " quadrature grid points (one dot per grid point): ")
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
                if(verbose) nimCat("(", i, ")")
                for (j in 1:nQuadGrid) {
                    if(verbose) nimCat(".")
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
            if(verbose) nimCat("\n")
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
                stop("calcMarginalParamQuad: Transformed parameter index, `pIndex`, requested is invalid.")
                
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
            vals <- paramGrid$nodes()  # z
            for(i in 1:dim(vals)[1]) {
                if (paramGridSkewed) {
                    for (d in 1:nParamTrans) {
                        vals[i,d] <- vals[i,d] * skewedStdDev[d, step(vals[i,d]) + 1]  
                    }
                }
                vals[i,] <- z_to_paramTrans(vals[i,], paramTransMode, Atransform, transformMethod)
            }
            returnType(double(2))
            return(vals)
        }
    )
)


makeNodeString <- function(nodes, model) {
    if (!length(nodes))
        return("")
    elements <- model$expandNodeNames(nodes, returnScalarComponents = TRUE)
    vars <- sapply(strsplit(elements, "[", fixed = TRUE), `[[`, 1)    
    nodesCount <- table(vars)
    items <- elements[nodesCount[vars] == 1]
    multiples <- names(nodesCount[nodesCount > 1])
    if(length(multiples)) 
        items <- c(items, paste0(multiples, " (", nodesCount[nodesCount > 1], " elements)"))
    return(paste0(items, collapse = ", "))
}
