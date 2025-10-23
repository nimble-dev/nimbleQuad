#' Main user interface for NIMBLE's nested approximation
#'
#' This file provides the main user-facing functions and helpers for running nested approximations
#' and summarizing results in the NIMBLE framework. It includes the main summary class, the main
#' wrapper for running approximations, and utilities for improving marginals, sampling, and more.
#'

### Example workflow
## Rapprox <- buildNestedApprox(model)
## capprox <- compileNimble(Rapprox, project = model)
## result <- runNestedApprox(capprox)
## improveParamMarginals(result, nodes = 'sigma')
## paramSamples <- sampleParams(result, n=1000)
## samples <- sampleLatents(result, n=1000)

## alternatively, with functions converted to class methods:
## result$improveParamMarginals(nodes = 'sigma')
## result$sampleParams(n=1000)

## Class for holding nestedApprox object and various outputs/summaries computed from it
#' @importFrom R6 R6Class
approxSummary <- R6Class("approxSummary",
    public = list(
        initialize = function(approx, quantiles, expectations, marginalsApprox,
                              marginalsRaw, indivParamTransforms,
                              originalScale, marginalLogLik,
                              marginalLogLikImproved, samples, paramSamples) {
            self$approx <- approx
            self$quantiles <- quantiles
            self$expectations <- expectations
            self$marginalsApprox <- marginalsApprox
            self$marginalsRaw <- marginalsRaw
            self$indivParamTransforms <- indivParamTransforms
            self$originalScale <- originalScale
            self$marginalLogLik <- marginalLogLik
            self$marginalLogLikImproved <- marginalLogLikImproved
            self$samples <- samples
            self$paramSamples <- paramSamples
        },
        generateParamsMatrix = function() {
            if(is(self$approx, "nestedApprox")) 
                Rapprox <- self$approx else Rapprox <- self$approx$Robject

            first <- which(!sapply(self$quantiles, is.null))[1]
            qs <- self$quantiles[[first]]
            first <- which(!sapply(self$expectations, is.null))[1]
            exps <- self$expectations[[first]]
            
            ## Form tabular info on expectations and quantiles as a dataframe (or could
            ## be a matrix as with INLA output), which should print nicely.
            params <- list()
            for (i in seq_along(exps)) {
                tmp <- sapply(self$expectations, `[`, i)
                names(tmp) <- NULL
                params[[names(exps)[i]]] <- tmp
            }
            for (i in seq_along(qs)) {
                tmp <- sapply(self$quantiles, `[`, i)
                names(tmp) <- NULL
                params[[names(qs)[i]]] <- tmp
            }
            nms <- names(params)
            params <- as.data.frame(params)
            names(params) <- nms  # Deals with `25%` formatting.
            row.names(params) <- names(self$quantiles)
            self$params <- params
        },
        print = function() {
            cat("Model (hyper)parameters: \n")
            if (is.null(self$params)) self$params <- self$generateParamsMatrix()
            if(length(self$params)) {
                print(self$params)
            } else cat("  No analytic marginals available for non-1:1 transformations; use `sampleParams`.\n")
            cat("\nMarginal log-likelihood (asymmetric Gaussian approximation): ",
                self$marginalLogLik, "(*)\n")
            if (!is.na(self$marginalLogLikImproved))
                cat("Marginal log-likelihood (grid-based): ", self$marginalLogLikImproved, "(*)\n")
            cat("(*) Marginal log-likelihood is invalid for improper priors and may not be useful\nfor non-informative priors.\n")
            
            invisible(self)
        },
        setParamGrid = function(summary, quadRule = "NULL", nQuad = -1, prune = -1){
            setParamGrid(self, quadRule, nQuad, prune)
        },
        improveParamMarginals = function(nodes, nMarginalGrid = 5, nQuad = 3, quadRule = "NULL", prune = -1, transform = "spectral") {
            improveParamMarginals(self, nodes, nMarginalGrid, nQuad, quadRule, prune, transform)
        },
        calcMarginalLogLikImproved = function() {
            calcMarginalLogLikImproved(self)
        },
        sampleParams = function(n = 1000, matchMarginals = TRUE) {
            sampleParams(self, n, matchMarginals)
        },
        sampleLatents = function(n = 1000, includeParams = FALSE) {
            sampleLatents(self, n, includeParams)
        },
        qmarginal = function(node, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975)) {
            qmarginal(self, node, quantiles)
        },
        rmarginal = function(node, n = 1000) {
            rmarginal(self, node, n)
        },
        dmarginal = function(node, x, log = FALSE) {
            dmarginal(self, node, x, log)
        },
        emarginal = function(node, functional, ...) {
            emarginal(self, node, functional, ...)
        },
        plotMarginal = function(node, log = FALSE, add = FALSE, ...){
            plotMarginal(self, node, log, add = add, ...)
        },
        approx = NULL,
        quantiles = NULL,
        expectations = NULL,
        marginalsApprox = NULL,
        marginalsRaw = NULL,
        indivParamTransforms = NULL,
        originalScale = NULL,
        marginalLogLik = NULL,
        marginalLogLikImproved = NULL,
        samples = NULL,
        paramSamples = NULL,
        params = NULL
    )
)


#' Run a nested approximation, returning a summary object with default inference
#'
#' Uses a nested approximation (compiled or uncompiled) returned from \code{buildNestedApprox}) 
#' to do default inference and return a summary object that can be used for additional tailored inference.
#' It estimates marginal distributions for parameters (quantiles and expectations), and can 
#' optionally return posterior samples for the latent nodes and parameters
#'
#' @param approx a compiled or uncompiled nestedApprox object created by \code{buildNestedApprox}.
#' @param quantiles numeric vector of quantiles to estimate for each parameter. Default is \code{c(0.025, 0.25, 0.5, 0.75, 0.975)}.
#' @param originalScale logical; if \code{TRUE}, report results on the original scales of the parameters and latent nodes. 
#'   Default is \code{TRUE}.
#' @param improve1d logical; if \code{TRUE} and there is a single parameter, improve the estimate of the estimate of marginal
#' distribution for the marginal by directly using the Laplace/AGHQ approximate marginal distribution rather than the asymmetric
#' Gaussian approximation. Default is \code{TRUE}.
#' @param nSamplesLatents number of samples of the latent nodes to draw. Default is 0.
#' @param nSamplesParams number of samples of the parameter nodes to draw. Default is 0.
#'
#' @details 
#' 
#' This is the main user interface for running a nested approximation. It carries out
#' default inference and then returns a summary object that can be used for further inference
#' by calling methods on the summary object, as seen in the examples (or running the equivalent 
#' function calls with the first argument being the summary object).
#' 
#' @section Methods available for object of class \code{approxSummary}
#'
#' Once the default inference has been run, inference can then be improved by calling different available methods within the returned object.
#' Each method is explained in detail in their documentation, but the user may choose
#' 
#' \itemize{
#'     \item \code{setParamGrid}. Allows the user to change the parameter grid used in the nested approximation.         
#'     \item \code{improveParamMarginals}. Improve univariate parameter marginals using grid-based quadrature.
#'     \item \code{calcMarginalLogLikImproved}. Calculate improved marginal log-likelihood using grid-based quadrature.
#'     \item \code{sampleParams}. Sample from the parameter posterior distribution.
#'     \item \code{sampleLatents}. Sample from the posterior distribution of the latent nodes.
#'     \item \code{qmarginal}. Compute quantiles for a parameter.
#'     \item \code{rmarginal}. Draw random samples from the marginal posterior of a parameter.
#'     \item \code{emarginal}. Compute the expectation of a function of a parameter under the marginal posterior distribution.
#'     \item \code{plotMarginal}. Plot the marginal posterior for a parameter.
#'   }
#'
#' @return An object of class \code{approxSummary} containing initial results that can be used to carry out further inference.
#' 
#' @author Christopher Paciorek
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
#' \dontrun{
#' comp_model <- compileNimble(model)
#' comp_approx <- compileNimble(approx, project = model)
#' result <- runNestedApprox(comp_approx)
#' # Improve marginals for a parameter node using AGHQ.
#' result$improveParamMarginals(nodes = 'Tau_re', nMarginalGrid = 9)
#' # Specify other quantiles of interest.
#' result$qmarginal('Tau_re', quantiles = c(.05, .95))
#' # Compute other expectations of interest, here the mean on the standard deviation scale
#' result$emarginal('Tau_re', function(x) 1/sqrt(x))
#' 
#' # Sample from the approximate posterior for the latent nodes.
#' latent_sample <- result$sampleLatents(n = 1000)
#' # For joint inference on parameters, sample from the approximate posterior for the parameters.
#' param_sample <- result$sampleParams(n = 1000)
#'
#' @export
#' 
runNestedApprox <- function(approx, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
                            originalScale = TRUE, improve1d = TRUE,
                            nSamplesLatents = 0, nSamplesParams = 0) {
    if(is(approx, "nestedApprox")) 
        Rapprox <- approx else Rapprox <- approx$Robject

    nParamTrans <- Rapprox$nParamTrans

    marginalsRaw <- list()
    length(marginalsRaw) <- nParamTrans
    marginalsApprox <- list()
    length(marginalsApprox) <- nParamTrans
    indivParamTransforms <- list()

    quantileEsts <- expectations <- list()

    ## Estimate marginals on transformed scale first.
    for (i in seq_len(nParamTrans)) {
        marginalsRaw[[i]] <- approx$calcMarginalParamIntegFree(i)
        marginalsApprox[[i]] <- fitMarginalSpline(marginalsRaw[[i]])
        if(!originalScale) {
            length(quantileEsts) <- length(expectations) <- length(indivParamTransforms) <- nParamTrans
            quantileEsts[[i]] <- estimateQuantiles(marginalsApprox[[i]], NULL, quantiles)
            expectations[[i]] <- estimateExpectations(marginalsApprox[[i]], NULL)
            names(indivParamTransforms) <- names(quantileEsts) <- names(expectations) <-
                paste0("param_trans", seq_len(nParamTrans))
        }
    }

    if(originalScale) {
        length(quantileEsts) <- length(expectations) <- length(indivParamTransforms) <- sum(Rapprox$paramNodesIndices > 0)
        cnt <- 0
        for (i in seq_along(Rapprox$paramNodesComponents)) {
            idx <- Rapprox$paramNodesIndices[i]
            if(idx > 0) {  # 1:1 case
                cnt <- cnt + 1
                indivParamTransforms[[cnt]] <- parameterTransform(Rapprox$model, Rapprox$paramNodesComponents[i])
                quantileEsts[[cnt]] <- estimateQuantiles(marginalsApprox[[idx]], indivParamTransforms[[cnt]],
                                                         quantiles)
                expectations[[cnt]] <- estimateExpectations(marginalsApprox[[idx]], indivParamTransforms[[cnt]])
            }
        }
        names(quantileEsts) <- names(expectations) <- names(indivParamTransforms) <-
            Rapprox$paramNodesComponents[Rapprox$paramNodesIndices > 0]
    }

    marginalLogLik <- approx$calcMarginalLogLikApprox()

    summary <- approxSummary$new(approx, quantileEsts, expectations, marginalsApprox,
        marginalsRaw, indivParamTransforms, originalScale, marginalLogLik, NA, NULL,
        NULL)

    ## With only one parameter, computations will not be slow unless number of latents is large,
    ## so go ahead and use better marginal and logLik estimates. 
    if (nParamTrans == 1 && improve1d) {
        improveParamMarginals(summary, ifelse(originalScale, Rapprox$paramNodesComponents[1], 1))
        summary$marginalLogLikImproved <- approx$calcMarginalLogLikQuad()
    }


    ## This is expensive. It is provided to enable an all-in-one experience.
    ## Avoid if user only needs parameter inference.
    ## Also, how do we have user tell us whether to `includeParams`?
    ## For now they must use more manual workflow if they need that.
    if (nSamplesLatents) 
        sampleLatents(summary, n = nSamplesLatents, includeParams = FALSE)

    if (nSamplesParams) 
        sampleParams(summary, n = nSamplesParams)

    return(summary)
}


#' Set the parameter grid for the nested approximation
#'
#' Allows the user to change the parameter grid used in the nested approximation.
#'
#' @param summary an approxSummary object, returned by \code{runNestedApprox}.
#' @param quadRule quadrature rule to use for the parameter grid. Can be any of
#'        \code{"CCD"}, \code{"AGHQ"}, \code{"AGHQSPARSE"} or \code{"USER"},
#'          the latter for user-defined grids.
#' @param nQuad number of quadrature points (not used for \code{"CCD"}.
#' @param prune pruning parameter for removing AGHQ points at low-density points.
#' 
#' @details If the chosen quadrature rule is \code{"AGHQSPARSE"}, then \code{sampleLatents} will no longer work
#' as it requires that all the quadrature weights are non-negative, which is no longer true for sparse AGHQ. To use a 
#' \code{"USER"} supplied quadrature rule, then this will need to be defined prior to calling \code{buildNestedApprox}.
#'
#' @return None. Modifies the summary object in place.
#' @export
setParamGrid <- function(summary, quadRule = "NULL", nQuad = -1, prune = -1){
  summary$approx$buildParamGrid(quadRule, nQuad, prune)
}


## Get the index of a parameter node.
getNodeIndex <- function(node, Rapprox) {
    if (is.character(node)) {
        mtch <- which(node == Rapprox$paramNodesComponents)
        if(length(mtch) != 1)
            stop("node `", node, "` is not a parameter element")
        idx <- Rapprox$paramNodesIndices[mtch]
        if(idx == 0)
            stop("node `", node, "`is not involved in a 1:1 parameter transformation, so marginals cannot be estimated by analytic approximation. Use `sampleParams` for inference.")
                 
    } else {
        if(node > Rapprox$innerMethods$nparTrans)  # AGHQ should have a method providing this.
            stop("Numeric index value ", node, " exceeds number of transformed parameters")
        idx <- node
    }
    return(idx)
}


##
#' Improve univariate parameter marginals using grid-based quadrature
#'
#' Uses d-1 dimensional quadrature (by default AGHQ) to get improved univariate marginal estimates for parameters.
#' Users can select to apply to specific nodes of interest to limit computation. 
#' 
#' @param summary an approxSummary object, returned by \code{runNestedApprox}.
#' @param nodes parameter nodes to improve inference for. Specified as character (when using original scale) 
#' or integer (when using transformed scale), where the scale was specified in \code{runNestedApprox}.
#' @param nMarginalGrid number of grid points for marginal calculations. Default is 5.
#' @param nQuad number of AGHQ quadrature points. Default is 5 if d=2 and 3 otherwise.
#' @param quadRule quadrature rule to use for the parameter grid. Can be any of
#'        \code{"AGHQ"}, \code{"CCD"}, \code{"AGHQSPARSE"} or \code{"USER"},
#'          the latter for user-defined grids, but standard use will be of \code{"AGHQ"}.
#' @param prune pruning parameter for removing AGHQ points at low-density points.
#' @param transform grid transformation method for internal AGHQ. Default is \code{"spectral"}.
#'
#' @return The modified \code{approxSummary} object with improved marginals.
#'
#' @details
#'
#' See \code{runNestedApprox} for example usage.
#' 
#' @export 
#' 
improveParamMarginals <- function(summary, nodes, nMarginalGrid = 5, nQuad, quadRule = "NULL", prune = -1, transform = "spectral") {
    Rapprox <- summary$approx$Robject

    originalScale <- summary$originalScale

    if(missing(nodes))
        if(originalScale) {
            nodes <- Rapprox$innerMethods$paramNodes
        } else nodes <- seq_len(Rapprox$nParamTrans)
    
    if(originalScale) {
        if(!is.character(nodes))
            stop("Results are being reported on the original scale as specified in the model. `nodes` must contain model node(s) or variable(s).")
    }
    if(!originalScale && is.character(nodes))
        stop("Results are being reported on the transformed (unconstrained) scale. `nodes` must contain one or more integer values indicating the transformed parameters.")
    
    if(is.character(nodes)) 
        nodes <- Rapprox$model$expandNodeNames(nodes, returnScalarComponents = TRUE)

    if(missing(nQuad))
        nQuad <- ifelse(Rapprox$innerMethods$nparTrans == 2, 5, 3)

    for (i in seq_along(nodes)) {
        ## Improve marginal and insert into raw and summary objects.
        idx <- getNodeIndex(nodes[i], Rapprox)
        if(!originalScale || idx > 0) {
            if(is.character(nodes[i])) paramName <- nodes[i] else paramName <- paste0("param_trans", nodes[i])
        
            summary$marginalsRaw[[idx]] <- summary$approx$calcMarginalParamQuad(idx,
                                                      nPts = nMarginalGrid, nQuad = nQuad, gridTransformMethod = transform, 
                                                      quadRule = quadRule, prune = prune)
            summary$marginalsApprox[[idx]] <- fitMarginalSpline(summary$marginalsRaw[[idx]])
            
            summary$quantiles[[paramName]] <- estimateQuantiles(summary$marginalsApprox[[idx]],
                                                                summary$indivParamTransforms[[paramName]])
            summary$expectations[[paramName]] <- estimateExpectations(summary$marginalsApprox[[idx]],
                                                                      summary$indivParamTransforms[[paramName]])
        }
    }
    summary$generateParamsMatrix()
    return(summary)
}

##
#' Calculate improved marginal log-likelihood using grid-based quadrature
#' 
#' Uses quadrature (by default AGHQ) to get an improved estimate of the marginal log-likelihood.
#'
#' @param summary an approxSummary object, returned by \code{runNestedApprox}.
#'
#' @return The improved marginal log-likelihood.
#' 
#' @details 
#' 
#' Users will not generally need to call this function directly, as it is called
#' automatically when sampling from the posterior of the latent nodes, since
#' its computation comes for free in that case.
#' 
#' Warning: the marginal log-likelihood is invalid for improper priors and may not be useful
#' for non-informative priors, because it averages the log-likelihood (approximately marginalized 
#' with respect to the latent nodes) over the prior distribution, thereby
#' including log-likelihood values corresponding to parameter values that are inconsistent with the data.
#' 
#' @export 
#' 
calcMarginalLogLikImproved <- function(summary) {
    summary$marginalLogLikImproved <- summary$approx$calcMarginalLogLikQuad()
    invisible(summary$marginalLogLikImproved)
}


##
#' Sample from the parameter posterior distribution
#'
#' Draws samples from the parameter posterior using the asymmetric Gaussian approximation. 
#' Optionally uses a copula approach to match the univariate marginals to the 
#' currently-available marginal distributions (based on either the initial asymmetric 
#' Gaussian approximation or improved marginals from 
#' calling \code{improveParamMarginals}).
#'
#' @param summary an approxSummary object, returned by \code{runNestedApprox}.
#' @param n number of samples to draw. Default is 1000.
#' @param matchMarginals logical; if \code{TRUE} (the default), match marginals using copula approach.
#'
#' @return Matrix of parameter samples.
#' 
#' @details 
#' 
#' Draws samples from the joint parameter posterior distribution (marginalized with respect to the latent nodes)
#' using the asymmetric Gaussian approximation.
#'
#' This is useful for joint inference on the parameters, including inference on functions of more than one parameter.
#'
#' See \code{runNestedApprox} for example usage.
#' 
#' @export
#'
sampleParams <- function(summary, n = 1000, matchMarginals = TRUE) {
    Rapprox <- summary$approx$Robject
    originalScale <- summary$originalScale

    samplesTrans <- summary$approx$simulateParams(n)

    if (matchMarginals) {
        ## Copula approach based on `inla.hperpar.sample` via
        ## `improve.marginals`: F^{-1}(F(x)) where F is ecdf of samples and
        ## F^{-1} is quantile fxn from approx's marginal.
        for (i in seq_len(ncol(samplesTrans))) {
                empirQuantiles <- ecdf(samplesTrans[, i])(samplesTrans[, i])
                quantiles <- estimateQuantiles(summary$marginalsApprox[[i]],
                                               NULL,
                                               empirQuantiles)
                samplesTrans[, i] <- quantiles
        }
    }

    ## Should we compile such that we can used compiled $paramsTransform?
    if (originalScale) {
        samples <- t(apply(samplesTrans, 1, Rapprox$innerMethods$paramsTransform$inverseTransform))
        if(Rapprox$nParamTrans == 1)
            samples <- matrix(samples, ncol = 1)
        colnames(samples) <- Rapprox$model$expandNodeNames(Rapprox$innerMethods$paramNodes,
                                                           returnScalarComponents = TRUE)
    } else {
        samples <- samplesTrans
        colnames(samples) <- paste0("param_trans", seq_len(ncol(samples)))
    }

    ## TODO: check that INLA's improve.marginals does nothing for non 1:1 cases.
    summary$paramSamples <- samples
    invisible(samples)
}

##
#' Sample from the posterior distribution of the latent nodes
#'
#' Draws samples from the posterior distribution of the latent nodes. 
#' Optionally includes parameter values corresponding to each sample.
#'
#' @param summary an approxSummary object, returned by \code{runNestedApprox}.
#' @param n Number of samples to draw (default: 1000).
#' @param includeParams logical; if \code{TRUE}, include parameter values corresponding to each sample. 
#' Default is \code{FALSE}.
#'
#' @return Matrix of latent samples.
#' 
#' @details 
#' 
#' The sampling approach uses stratified sampling from a weighted mixture of multivariate normals, 
#' where the weights are based on the
#' approximate marginal density at each grid point in the parameter grid. For each point, the multivariate normal
#' is based on Laplace approximation, using the maximum for the mean and the inverse Hessian for the covariance
#' matrix. This approach is not valid for sparse AGHQ due to negative quadrature weights. This can be updated by
#' \code{setParamGrid} and any quadrature rule other than \code{AGHQSPARSE}.
#' 
#' The parameter values corresponding to the samples can be requested via \code{includeParams}.
#' 
#' Note that NIMBLE's nested approximation framework does not provide marginals for the latent nodes
#' based on analytic approximation, so both joint and univariate inference on the latent nodes
#' is from sampling.
#'
#' See \code{runNestedApprox} for example usage.
#' 
#' @export
#' 
sampleLatents <- function(summary, n = 1000, includeParams = FALSE) {
    Rapprox <- summary$approx$Robject
    originalScale <- summary$originalScale

    samples <- summary$approx$simulateLatents(n)

    ## Grid-based marginal log-likelihood comes "for free" if simulate parameters.
    summary$marginalLogLikImproved <- summary$approx$calcMarginalLogLikQuad()

    if(originalScale) {
        nms <- Rapprox$innerMethods$reNodesAsScalars_vec
    } else nms <- paste0("latent_trans", seq_len(Rapprox$nreTrans))
    if(dim(samples)[2] == 2) 
        nms <- nms[1]

    if(originalScale && !all(Rapprox$innerMethods$reTransform$transformType == 1, na.rm = TRUE)) {
        samplesTrans <- t(apply(samples[ , -1], 1, Rapprox$innerMethods$reTransform$inverseTransform))
        samples <- cbind(samples[ , 1], samplesTrans)
    }
    colnames(samples) <- c("index", nms)
    
    if (includeParams) {
        paramValues <- summary$approx$getParamGrid()
        if(summary$originalScale) {
            paramValues <- apply(paramValues, 1, Rapprox$innerMethods$paramsTransform$inverseTransform)
            if(is.null(dim(paramValues)))
                paramValues <- matrix(paramValues, ncol = 1) else paramValues <- t(paramValues)
        }
        paramSamples <- paramValues[samples[, "index"], , drop = FALSE]
        if(summary$originalScale) {
            colnames(paramSamples) <- Rapprox$paramNodesComponents
        } else colnames(paramSamples) <- paste0('param_trans', seq_len(Rapprox$nParamTrans))
        samples <- cbind(samples, paramSamples)
    }
    summary$samples <- samples[, -1, drop = FALSE]
    invisible(summary$samples)
}

##
#' Compute quantiles for a parameter
#' 
#' Quantile estimation for univariate parameter marginals.
#'
#' @param summary an approxSummary object, returned by \code{runNestedApprox}.
#' @param node parameter node of interest. Specified as character (when using original scale) 
#' or integer (when using transformed scale), where the scale was specified in \code{runNestedApprox}.
#' @param quantiles numeric vector of quantiles to compute. Default is \code{c(0.025, 0.25, 0.5, 0.75, 0.975)}.
#'
#' @return Named vector of quantile estimates.
#' 
#' @details
#' 
#' Uses a spline approximation to the quantile function of the marginal posterior distribution,
#' based on a cached approximation of the marginal density on a fine grid.
#'
#' See \code{runNestedApprox} for example usage.
#' 
#' @export 
#' 
qmarginal <- function(summary, node, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975)) {
    Rapprox <- summary$approx$Robject
    idx <- getNodeIndex(node, Rapprox)
    if(is.character(node)) {
        paramTransform  <- parameterTransform(Rapprox$model, node) 
    } else paramTransform <- NULL
    quantileEsts <- estimateQuantiles(summary$marginalsApprox[[idx]], paramTransform, quantiles)
    names(quantileEsts) <- quantiles
    return(quantileEsts)
}

#' Draw random samples from the marginal posterior of a parameter
#' 
#' Random sampling for univariate parameter marginals.
#'
#' @param summary an approxSummary object, returned by \code{runNestedApprox}.
#' @param node parameter node of interest. Specified as character (when using original scale) 
#' or integer (when using transformed scale), where the scale was specified in \code{runNestedApprox}.
#' @param n number of samples to draw. Default is 1000.
#'
#' @return Numeric vector of samples.
#' 
#' @details Uses the inverse CDF method applied to the quantile function of the marginal posterior distribution.
#' 
#' @export
rmarginal <- function(summary, node, n = 1000) {
    samples <- qmarginal(summary, node, runif(n))
    names(samples) <- NULL
    return(samples)
}

#' Evaluate the marginal posterior density for a parameter.
#' 
#' Density evaluation for univariate parameter marginals.
#'
#' @param summary an approxSummary object, returned by \code{runNestedApprox}.
#' @param node parameter node of interest. Specified as character (when using original scale) 
#' or integer (when using transformed scale), where the scale was specified in \code{runNestedApprox}.
#' @param x numeric vector of values at which to evaluate the density.
#' @param log logical; if \code{TRUE}, return log-density. Default is \code{FALSE}.
#'
#' @return Numeric vector of (log-)density values.
#' 
#' @details 
#' 
#' Uses a spline approximation to the log-density of the marginal posterior distribution,
#' based on a cached approximation of the marginal density on a fine grid.
#' 
#' @export 
#' 
dmarginal <- function(summary, node, x, log = FALSE) {
    Rapprox <- summary$approx$Robject
    logDetJac <- 0
    idx <- getNodeIndex(node, Rapprox)
    if(is.character(node)) {
        paramTransform  <- parameterTransform(Rapprox$model, node)
        x <- sapply(x, paramTransform$transform)
        logDetJac <- sapply(x, paramTransform$logDetJacobian)
    }
    logPDF <- fitMarginalSpline(summary$marginalsRaw[[idx]], xnew = x, refine = FALSE) - logDetJac

    if(log) return(logPDF) else return(exp(logPDF))
}

##
#' Compute the expectation of a function of a parameter under the marginal posterior distribution
#' 
#' Posterior expectations for univariate parameter marginals.
#'
#' @param summary an approxSummary object, returned by \code{runNestedApprox}.
#' @param node parameter node of interest. Specified as character (when using original scale) 
#' or integer (when using transformed scale), where the scale was specified in \code{runNestedApprox}.
#' @param functional function to compute the expectation of.
#' @param ... Additional arguments passed to the function.
#'
#' @return Numeric value of the expectation.
#' 
#' @details 
#' 
#' Estimate the expectation of a function of a parameter using univariate numerical
#' integration based on a cached approximation of the marginal density on a fine grid.
#'
#' See \code{runNestedApprox} for example usage.
#' 
#' @export 
#' 
emarginal <- function(summary, node, functional, ...) {
    Rapprox <- summary$approx$Robject
    if(is.character(node))
        paramTransform  <- parameterTransform(Rapprox$model, node) else paramTransform <- NULL
    idx <- getNodeIndex(node, Rapprox)
    expectation <- estimateExpectations(summary$marginalsApprox[[idx]], paramTransform, functional, ...)
    return(expectation)
}

##
#' Plot the marginal posterior for a parameter
#' 
#' Univariate marginal posterior plotting for parameters.
#'
#' @param summary an approxSummary object, returned by \code{runNestedApprox}.
#' @param node parameter node of interest. Specified as character (when using original scale) 
#' or integer (when using transformed scale), where the scale was specified in \code{runNestedApprox}.
#' @param log logical; if \code{TRUE}, plot log-density. Default is \code{FALSE}.
#' @param add logical; if \code{TRUE}, add to existing plot. Default is \code{FALSE}.
#' @param ... Additional arguments passed to plotting function.
#'
#' @return None. Produces a plot.
plotMarginal <- function(summary, node, log = FALSE, add = FALSE, ...){
    minmax <- summary$qmarginal(node, c(.001, 0.999))
    x <- seq(minmax[1], minmax[2], length = 200)
    y <- summary$dmarginal(node, x, log)
    if(log) 
      ylab <- "Log Posterior Density"
    else
      ylab <- "Posterior Density"
    if(!add)
      plot(x, y, type = 'l', xlab = node, ylab = ylab,...)
    else
      lines(x, y, xlab = xlab, ylab = ylab,...)
}
