## Code for main user interface for NIMBLE's nested approximation.

## NOTE: none of this code has yet been run on any examples and there will be
## various bugs.

## Uses core algorithm code in `approxPosterior.R` and code for marginal
## summaries in `approxSummaries.R`.

### Example workflow
## Rapprox <- buildNestedApprox(model)
## capprox <- compileNimble(Rapprox, project = model)
## result <- runNestedApprox(capprox)
## improveMarginals(result, nodes = 'sigma')
## paramSamples <- sampleParamNodes(result, n=1000)
## samples <- sampleLatentNodes(result, n=1000)

## alternatively, with functions converted to class methods:
## result$improveMarginals(nodes = 'sigma')
## result$sampleParamNodes(n=1000)

## Class for holding nestedApprox object and various outputs/summaries computed
## from it when running `runNestedApprox` or individual functions that
## manipulate the approximation.
#' @importFrom R6 R6Class
approxSummary <- R6Class("approxSummary",
    public = list(
        initialize = function(approx, quantiles, expectations, marginalsApprox,
                              marginalsRaw, indivParamTransforms,
                              originalScale, marginalLogLik,
                              marginalLogLik_improved, samples, paramSamples) {
            self$approx <- approx
            self$quantiles <- quantiles
            self$expectations <- expectations
            self$marginalsApprox <- marginalsApprox
            self$marginalsRaw <- marginalsRaw
            self$indivParamTransforms <- indivParamTransforms
            self$originalScale <- originalScale
            self$marginalLogLik <- marginalLogLik
            self$marginalLogLik_improved <- marginalLogLik_improved
            self$samples <- samples
            self$paramSamples <- paramSamples
        },
        generateParamsMatrix = function() {
            if(is(approx, "NestedApprox")) 
                Rapprox <- approx else Rapprox <- approx$Robject

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
            ## params <- data.frame(mean = sapply(self$marginalsSummary, `[[`,
            ## 'mean'), sd = sapply(self$marginalsSummary, `[[`, 'sd'), row.names =
            ## paramNames)
            for (i in seq_along(qs)) {
                tmp <- sapply(self$quantiles, `[`, i)
                names(tmp) <- NULL
                params[[names(qs)[i]]] <- tmp
            }
            params <- as.data.frame(params)
            row.names(params) <- names(self$quantiles)
            self$params <- params
        },
        print = function() {
            cat("Model (hyper)parameters: \n")
            if (is.null(self$params)) self$params <- self$generateParamsMatrix()
            print(self$params)
            
            cat("\nMarginal log-likelihood (asymmetric Gaussian approximation): ",
                self$marginalLogLik, "\n")
            ## Careful with ref to AGHQ here as we do at the moment allow 'improved'
            ## calc under CCD.
            if (!is.na(self$marginalLogLik_improved))
                cat("Marginal log-likelihood (grid-based): ", self$marginalLogLik_improved, "\n")
            
            invisible(self)
        },
        improveMarginals = function(nodes, nMarginalGrid = 3, nQuad = 3) {
            improveMarginals(self, nodes, nMarginalGrid, nQuad)
        },
        calcMarginalLogLikImproved = function() {
            calcMarginalLogLikImproved(self)
        },
        sampleParamNodes = function(n = 1000, matchMarginals = TRUE) {
            sampleParamNodes(self, n, matchMarginals)
        },
        sampleLatentNodes = function(n = 1000, includeParams = FALSE) {
            sampleLatentNodes(self, n, includeParams)
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
        approx = NULL,
        quantiles = NULL,
        expectations = NULL,
        marginalsApprox = NULL,
        marginalsRaw = NULL,
        indivParamTransforms = NULL,
        originalScale = NULL,
        marginalLogLik = NULL,
        marginalLogLik_improved = NULL,
        samples = NULL,
        paramSamples = NULL,
        params = NULL
    )
)


## Main user-facing function for running a nested approximation and getting a
## results summary.  
runNestedApprox <- function(approx, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
                            originalScale = TRUE, nSamplesLatents = 0, nSamplesParams = 0) {
    if(is(approx, "NestedApprox")) 
        Rapprox <- approx else Rapprox <- approx$Robject

    nParamTrans <- Rapprox$npar

    marginalsRaw <- list()
    length(marginalsRaw) <- nParamTrans
    marginalsApprox <- list()
    length(marginalsApprox) <- nParamTrans
    indivParamTransforms <- list()

    quantileEsts <- expectations <- list()

    ## Estimate marginals on transformed scale first.
    for (i in seq_len(nParamTrans)) {
        marginalsRaw[[i]] <- approx$findMarginalHyperIntFree(i)
        marginalsApprox[[i]] <- fitMarginalSpline(marginalsRaw[[i]])
        if(!originalScale) {
            length(quantileEsts) <- length(expectations) <- length(indivParamTransforms) <- nParamTrans
            quantileEsts[[i]] <- estimateQuantiles(marginalsApprox[[i]], NULL, quantiles)
            expectations[[i]] <- estimateExpectations(marginalsApprox[[i]], NULL)
            names(indivParamTransforms) <- names(quantileEsts) <- names(expectations) <-
                paste0("param", seq_len(nParamTrans))
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
                quantileEsts[[cnt]] <- estimateQuantiles(marginalsApprox[[idx]], indivParamTransforms[[i]],
                                                         quantiles)
                expectations[[cnt]] <- estimateExpectations(marginalsApprox[[idx]], indivParamTransforms[[i]])
            }
        }
        names(quantileEsts) <- names(expectations) <- names(indivParamTransforms) <-
            Rapprox$paramNodesComponents[Rapprox$paramNodesIndices > 0]
    }

    marginalLogLik <- approx$calcMarginalLogLikApprox()

    summary <- approxSummary$new(approx, quantileEsts, expectations, marginalsApprox,
        marginalsRaw, indivParamTransforms, originalScale, marginalLogLik, NA, NULL,
        NULL)

    ## With only one parameter, computations should not be so slow, so go ahead
    ## and use better marginal and logLik estimates.
    if (length(Rapprox$paramNodesComponents) == 1) {
        improveMarginals(summary, ifelse(originalScale, Rapprox$paramNodesComponents[1], 1))
        ## TODO: make sure that `calcMarginalLogLikQuad()` works if
        ## `calcHyperGrid` has not yet been called.
        summary$marginalLogLik_improved <- approx$calcMarginalLogLikQuad()
    }

    ## Should we embed sampling in `runNestedApprox`?
    ## OTOH it would provide an all-in-one experience.
    ## OTOH, it complicates things, including the args, and would be cleaner to
    ## have users requests samples afterwards.

    ## This is expensive. Avoid if user only needs parameter inference.  How do
    ## we have user tell us whether to `includeParams`?  Perhaps tell them to
    ## use more manual workflow if they need that.
    if (nSamplesLatents) 
        sampleLatentNodes(summary, n = nSamplesLatents, includeParams = FALSE)

    if (nSamplesParams) 
        sampleParamNodes(summary, n = nSamplesParams)

    return(summary)
}

## Helper function that takes either a character string for an original node element
## and returns transformed parameter index (if in a 1:1 transformation)
## or simply checks that the index of the transformed parameter is valid.
getNodeIndex <- function(node, Rapprox) {
    if (is.character(node)) {
        mtch <- which(node == Rapprox$paramNodesComponents)
        if(length(mtch) != 1)
            stop("node `", nodes[i], "` is not a parameter element or is not involved in a 1:1 parameter transformation, so marginals cannot be estimated by analytic approximation. In the latter case, use `sampleParamNodes` for inference.")
        idx <- Rapprox$paramNodesIndices[mtch]
    } else {
        if(node > Rapprox$innerMethods$pTransform_length)
            stop("Numeric index value ", node, " exceeds number of transformed parameters")
        idx <- node
    }
    return(idx)
}


## This uses d-1 dimensional AGHQ to get improved univariate marginal estimates
## for parameters.
## Should it be called `improveParamMarginals`?
improveMarginals <- function(summary, nodes, nMarginalGrid = 3, nQuad = 3) {
    Rapprox <- summary$approx$Robject

    originalScale <- summary$originalScale
    if(originalScale && !is.character(nodes))
        stop("Results are being reported on the original scale as specified in the model. `nodes` must contain model node(s) or variable(s).")
    if(!originalScale && is.character(nodes))
        stop("Results are being reported on the transformed (unconstrained) scale. `nodes` must contain one or more integer values indicating the transformed parameters.")
    
    if(is.character(nodes)) 
        nodes <- Rapprox$model$expandNodeNames(nodes, returnScalarComponents = TRUE)

    for (i in seq_along(nodes)) {
        ## Improve marginal and insert into raw and summary objects.
        idx <- getNodeIndex(nodes[i], Rapprox)
        if(is.character(nodes[i])) paramName <- nodes[i] else paramName <- paste0("param", nodes[i])
        
        summary$marginalsRaw[[idx]] <- summary$approx$findMarginalPosteriorDensity(idx,
                                            nPts = nMarginalGrid, nQuad = nQuad)
        summary$marginalsApprox[[idx]] <- fitMarginalSpline(summary$marginalsRaw[[idx]])

        summary$quantiles[[paramName]] <- estimateQuantiles(summary$marginalsApprox[[idx]],
            summary$indivParamTransforms[[paramName]])
        summary$expectations[[paramName]] <- estimateExpectations(summary$marginalsApprox[[idx]],
            summary$indivParamTransforms[[paramName]])
    }
    summary$generateParamsMatrix()
    return(summary)
}

calcMarginalLogLikImproved <- function(summary) {
    summary$marginalLogLik_improved <- summary$approx$calcMarginalLogLikQuad()
    invisible(summary$marginalLogLik_improved)
}


## This uses whatever marginals (asymm Gaussian approx or improved) are in
## `summary`.  Potentially called from runNestedApprox or independently.
sampleParamNodes <- function(summary, n = 1000, matchMarginals = TRUE) {
    Rapprox <- summary$approx$Robject
    originalScale <- summary$originalScale

    samplesTrans <- summary$approx$simulateHyperParams(n)

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
        if(Rapprox$npar == 1)
            samples <- matrix(samples, ncol = 1)
        colnames(samples) <- Rapprox$model$expandNodeNames(Rapprox$innerMethods$paramNodes,
                                                           returnScalarComponents = TRUE)
    } else {
        samples <- samplesTrans
        colnames(samples) <- paste0("param", seq_len(ncol(samples)))
    }

    ## TODO: check that INLA's improve.marginals does nothing for non 1:1
    ## cases.
    summary$paramSamples <- samples
    invisible(samples)
}

## Potentially called from runNestedApprox or independently.
sampleLatentNodes <- function(summary, n = 1000, includeParams = FALSE) {
    Rapprox <- summary$approx$Robject
    samples <- summary$approx$simulateLatentEffects(n)

    ## Grid-based marginal log-likelihood comes "for free" if simulate parameters.
    summary$marginalLogLik_improved <- summary$approx$calcMarginalLogLikQuad()

    colnames(samples) <- c("index", Rapprox$innerMethods$reNodesAsScalars_vec)
    if (includeParams) {
        paramValues <- getParamGrid()  ## `getParamGrid` needs to be written as a nestedApprox nf method, returning `theta_grid$nodes`.
        paramValues <- t(apply(paramValues, 1, Rapprox$innerMethods$paramsTransform))
        paramSamples <- paramValues[samples[, "index"], ]
        colnames(paramSamples) <- Rapprox$paramNodesComponents
        samples <- cbind(samples, paramSamples)
    }
    summary$samples <- samples[, -1]
    invisible(summary$samples)
}

qmarginal <- function(summary, node, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975)) {
    Rapprox <- summary$approx$Robject
    if(is.character(node))
        paramTransform  <- parameterTransform(Rapprox$model, node) else paramTransform <- NULL
    idx <- getNodeIndex(node, Rapprox)
    quantileEsts <- estimateQuantiles(summary$marginalsApprox[[idx]], paramTransform, quantiles)
    names(quantileEsts) <- quantiles
    return(quantileEsts)
}

rmarginal <- function(summary, node, n = 1000) {
    samples <- qmarginal(summary, node, runif(n))
    names(samples) <- NULL
    return(samples)
}

dmarginal <- function(summary, node, x, log = FALSE) {
    Rapprox <- summary$approx$Robject
    if(is.character(node)) {
        paramTransform  <- parameterTransform(Rapprox$model, node)
        x <- sapply(x, paramTransform$transform)
    }
    idx <- getNodeIndex(node, Rapprox)

    logPDF <- fitMarginalSpline(summary$marginalsRaw[[idx]], xnew = x)
    if(log) return(logPDF) else return(exp(logPDF))
}

emarginal <- function(summary, node, functional, ...) {
    Rapprox <- summary$approx$Robject
    if(is.character(node))
        paramTransform  <- parameterTransform(Rapprox$model, node) else paramTransform <- NULL
    idx <- getNodeIndex(node, Rapprox)
    expectation <- estimateExpectations(summary$marginalsApprox[[idx]], paramTransform, functional, ...)
    return(expectation)
}
