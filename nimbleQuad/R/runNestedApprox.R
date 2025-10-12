## Code for main user interface for NIMBLE's nested approximation.

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

## Class for holding nestedApprox object and various outputs/summaries computed
## from it when running `runNestedApprox` or individual functions that
## manipulate the approximation.
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
                self$marginalLogLik, "(*)\n", sep = "")
            if(self$approx$paramGridRule == "CCD")
                extra <- "(**)" else extra <- ""
            if (!is.na(self$marginalLogLikImproved))
                cat("Marginal log-likelihood (grid-based", extra, "): ", self$marginalLogLikImproved, "(*)\n", sep = "")
            cat("(*) Marginal log-likelihood is invalid for improper priors and may not be useful\nfor non-informative priors.\n")
            if(!is.na(self$marginalLogLikImproved) && self$approx$paramGridRule == "CCD")
                cat("(**) Estimated using CCD grid. Estimation based on an AGHQ grid may be more\naccurate (but more computationally expensive).\n")
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
        marginalsApprox = NULL,  # perhaps make private
        marginalsRaw = NULL,     # perhaps make private
        indivParamTransforms = NULL,
        originalScale = NULL,
        marginalLogLik = NULL,
        marginalLogLikImproved = NULL,
        samples = NULL,
        paramSamples = NULL,
        params = NULL
    )
)


## Main user-facing function for running a nested approximation and getting a
## results summary.  
runNestedApprox <- function(approx, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
                            originalScale = TRUE, improve1d = TRUE,
                            nSamplesLatents = 0, nSamplesParams = 0) {
    if(is(approx, "nestedApprox")) {
        Rapprox <- approx
        messageIfVerbose('  [Warning] Running an uncompiled nested approximation.  Use compileNimble() for faster execution.')
        tmp <- Rapprox$innerMethods$gr_logDens_pTransformed
        tmp <- Rapprox$innerMethods$calcLogDens_pTransformed
        for(i in seq_along(Rapprox$innerMethods$AGHQuad_nfl)) {
            tmp <- Rapprox$innerMethods$AGHQuad_nfl[[i]]$gr_inner_logLik
            tmp <- Rapprox$innerMethods$AGHQuad_nfl[[i]]$he_inner_logLik
        }
    } else Rapprox <- approx$Robject

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

    ## Should we embed sampling in `runNestedApprox`?
    ## OTOH it would provide an all-in-one experience.
    ## OTOH, it complicates things, including the args, and would be cleaner to
    ## have users requests samples afterwards.

    ## This is expensive. Avoid if user only needs parameter inference.
    ## Also, how do we have user tell us whether to `includeParams`?
    ## For now they must use more manual workflow if they need that.
    if (nSamplesLatents) 
        sampleLatents(summary, n = nSamplesLatents, includeParams = FALSE)

    if (nSamplesParams) 
        sampleParams(summary, n = nSamplesParams)

    return(summary)
}

## Add option for the user to change the parameter grid in the wrapper.
setParamGrid <- function(summary, quadRule = "NULL", nQuad = -1, prune = -1){
  summary$approx$buildParamGrid(quadRule, nQuad, prune)
}

## Helper function that takes either a character string for an original node element
## and returns transformed parameter index (if in a 1:1 transformation)
## or simply checks that the index of the transformed parameter is valid.
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


## This uses d-1 dimensional AGHQ to get improved univariate marginal estimates
## for parameters.
## Note that quadRule = "NULL" makes sure the default is orginal user choice.
improveParamMarginals <- function(summary, nodes, nMarginalGrid = 5, nQuad, quadRule = "NULL", prune = -1, transform = "spectral") {
    Rapprox <- summary$approx$Robject

    originalScale <- summary$originalScale

    if(missing(nodes))
        if(originalScale) {
            nodes <- Rapprox$innerMethods$paramNodes
        } else nodes <- seq_len(Rapprox$nParamTrans)

    if(!quadRule %in% c("AGHQ", "AGHQSPARSE"))
        stop("Only AGHQ-based quadrature rules are available for integration-based estimation of marginals.")
    
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

calcMarginalLogLikImproved <- function(summary) {
    summary$marginalLogLikImproved <- summary$approx$calcMarginalLogLikQuad()
    invisible(summary$marginalLogLikImproved)
}


## This uses whatever marginals (asymm Gaussian approx or improved) are in
## `summary`.  Potentially called from runNestedApprox or independently.
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

## Potentially called from runNestedApprox or independently.
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

rmarginal <- function(summary, node, n = 1000) {
    samples <- qmarginal(summary, node, runif(n))
    names(samples) <- NULL
    return(samples)
}

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

emarginal <- function(summary, node, functional, ...) {
    Rapprox <- summary$approx$Robject
    if(is.character(node))
        paramTransform  <- parameterTransform(Rapprox$model, node) else paramTransform <- NULL
    idx <- getNodeIndex(node, Rapprox)
    expectation <- estimateExpectations(summary$marginalsApprox[[idx]], paramTransform, functional, ...)
    return(expectation)
}

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
