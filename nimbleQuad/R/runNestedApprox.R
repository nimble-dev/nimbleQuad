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

            if (self$originalScale) {
                paramNames <- paste0("param", Rapprox$npar)
            } else {
                paramNames <- names(self$quantiles)
            }
            first <- which(!sapply(self$quantiles, is.null))[1]
            qs <- self$quantiles[[first]]$quantiles
            first <- which(!sapply(self$expectations, is.null))[1]
            exps <- self$expectations[[first]]
            
            ## Form tabular info on expectations and quantiles as a dataframe (or could
            ## be a matrix as with INLA output), which should print nicely.
            params <- data.frame()
            for (i in seq_along(exps))
                params[[names(exps)[i]]] <- sapply(self$expectations, `[`, i)
            ## params <- data.frame(mean = sapply(self$marginalsSummary, `[[`,
            ## 'mean'), sd = sapply(self$marginalsSummary, `[[`, 'sd'), row.names =
            ## paramNames)
            for (i in seq_along(qs))
                params[[names(qs)[i]]] <- sapply(self$quantiles, `[`, i)
            row.names(params) <- paramNames
            self$params <- params
        },
        print = function() {
            cat("Model (hyper)parameters: \n")
            if (is.null(self$params)) self$params <- self$generateParamsMatrix()
            print(self$params)
            
            cat("Marginal log-likelihood (asymmetric Gaussian approximation): ",
                self$marginalLogLik, "\n")
            ## Careful with ref to AGHQ here as we do at the moment allow 'improved'
            ## calc under CCD.
            if (!is.na(self$marginalLogLik_improved))
                cat("Marginal log-likelihood (AGHQ): ", self$marginalLogLik_improved, "\n")
            
            invisible(self)
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
## results summary.  Note in roxygen that `functionalsScale` will have no
## effect if `originalScale=FALSE`, in which case functional will compute on
## the transformed parameter values.
runNestedApprox <- function(approx, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
                            originalScale = TRUE, nSamplesLatents = 0, nSamplesParams = 0,
                            functionals = NULL, functionalsArgs = NULL, functionalsScale = "original") {
    if(is(approx, "NestedApprox")) 
        Rapprox <- approx else Rapprox <- approx$Robject

    nParamTrans <- Rapprox$npar

    marginalsRaw <- list()
    length(marginalsRaw) <- nParamTrans
    marginalsApprox <- list()
    length(marginalsApprox) <- nParamTrans
    indivParamTransforms <- list()
    length(indivParamTransforms) <- nParamTrans

    nMarginalsReport <- 0

    for (i in seq_len(nParamTrans)) {
        ## Save computation of marginals, except non-1:1 transformed parameters.
        if (!originalScale || i %in% Rapprox$paramNodesIndices) {
            marginalsRaw[[i]] <- approx$findMarginalHyperIntFree(Rapprox$paramNodesIndices[i])
            nMarginalsReport <- nMarginalsReport + 1
        }
    }

    cnt <- 1
    quantileEsts <- expectations <- list()
    length(quantileEsts) <- length(expectations) <- nMarginalsReport

    cnt <- 0
    for (i in seq_len(nParamTrans)) {
        if (!is.null(marginalsRaw[[i]])) {
            cnt <- cnt + 1
            if (originalScale)
                indivParamTransforms[[i]] <- parameterTransform(Rapprox$model, Rapprox$paramNodesComponents[i])
            marginalsApprox[[i]] <- fitMarginalSpline(marginalsRaw[[i]])
            quantileEsts[[cnt]] <- estimateQuantiles(marginalsApprox[[i]], indivParamTransforms[[i]],
                quantiles)
            expectations[[cnt]] <- estimateExpectations(marginalsApprox[[i]], indivParamTransforms[[i]],
                functionals = functionals, functionalsArgs = functionalsArgs, scale = functionalsScale)
        }
    }

    if (originalScale) {
        names(quantileEsts) <- Rapprox$paramNodesComponents[Rapprox$paramNodesIndices > 0]
        names(expectations) <- Rapprox$paramNodesComponents[Rapprox$paramNodesIndices > 0]
    }

    marginalLogLik <- approx$calcMarginalLogLikApprox()

    summary <- approxSummary$new(approx, quantileEsts, expectations, marginalsApprox,
        marginalsRaw, indivParamTransforms, originalScale, marginalLogLik, NA, NULL,
        NULL)

    ## With only one parameter, computations should not be so slow, so go ahead
    ## and use better marginal and logLik estimates.
    if (length(Rapprox$paramNodesComponents) == 1) {
        improveMarginals(summary, Rapprox$paramNodesComponents[1], nMarginalGrid = 3,
            functionals = functionals, functionalsArgs = functionalsArgs, functionalsScale = functionalsScale)
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
    if (nSamplesLatents) {
        summary$samples <- sampleLatentNodes(summary, n = nSamplesLatents, includeParams = FALSE)
        if (Rapprox$hyperGridRule == "AGHQ")
            summary$marginalLogLik_improved <- approx$calcMarginalLogLikQuad()  ## TODO: should be named 'get'?
    }

    if (nSamplesParams) {
        summary$paramSamples <- sampleParamNodes(summary, n = nSamplesParams)
    }

    return(summary)
}



## This uses d-1 dimensional AGHQ to get improved univariate marginal estimates
## for parameters Should it be called `improveParamMarginals`?
improveMarginals <- function(summary, nodes, nMarginalGrid = 3, functionals = NULL,
                             functionalsArgs = NULL, functionalsScale = "original") {
    Rapprox <- summary$approx$Robject
    originalScale <- summary$originalScale
    nodes <- Rapprox$model$expandNodeNames(nodes, returnScalarComponents = TRUE)
    if (!all(nodes %in% Rapprox$paramNodesComponents)) {
        badNodes <- nodes[!nodes %in% Rapprox$paramNodesComponents]
        stop("improveMarginals: nodes `", paste0(badNodes, sep = "`, `"), "` are not parameter elements or are not involved in 1:1 parameter transformation, so marginals cannot be estimated by analytic approximation. In the latter case, use `sampleParamNodes` for inference.")
    }
    marginalsRaw <- list()
    length(marginalsRaw) <- Rapprox$innerMethods$npar
    for (i in seq_len(nodes)) {
        ## Improve marginal and insert into raw and summary objects.
        if (originalScale) {
            idx = Rapprox$paramNodesIndices[which(nodes[i] == paramNodesComponents)]
        } else idx <- i
        ## TODO: what should `nQuad` and should it be what is set in
        ## `buildNestedApprox` setup. And make sure it is available via
        ## `approx$nQuadMarginal`.
        summary$marginalsRaw[[idx]] <- findMarginalPosteriorDensity(idx, nPts = nMarginalGrid,
            nQuad = summary$approx$nQuadMarginal)
        if (originalScale)
            paramTransform <- parameterTransform(Rapprox$model, nodes[i]) else paramTransform <- NULL
        summary$marginalsApprox[[idx]] <- fitMarginalSpline(summary$marginalsRaw[[idx]])
        ## TODO: check that nodes[i] name will correspond to element in
        ## `marginalsSummary`
        if (originalScale)
            idx2 <- nodes[i] else idx2 <- i
        summary$quantileEsts[[idx2]] <- estimateQuantiles(summary$marginalsApprox[[idx]],
            summary$indivParamTransforms[[idx]])
        summary$expectations[[idx2]] <- estimateExpectations(summary$marginalsApprox[[idx]],
            summary$indivParamTransforms[[idx]], functionals = functionals, functionalsArgs = functionalsArgs,
            scale = functionalsScale)
    }
    summary$generateParamsMatrix()
    return(summary)
}

calcMarginalLogLikImproved <- function(summary) {
    summary$marginalLogLik_improved <- approx$calcMarginalLogLikQuad()
    invisible(summary$marginalLogLik_improved)
}


## This uses whatever marginals (asymm Gaussian approx or improved) are in
## `summary`.  Potentially called from runNestedApprox or independently.
sampleParamNodes <- function(summary, n = 1000, matchMarginals = TRUE) {
    Rapprox <- summary$approx$Robject
    originalScale <- summary$originalScale

    samplesTrans <- summary$approx$simulateHyperParams(n)
    if (originalScale) {
        samples <- t(apply(samplesTrans, 1, Rapprox$innerMethods$paramsTransform$inverseTransform))
    } else samples <- samplesTrans

    if (matchMarginals) {
        ## Copula approac based on `inla.hperpar.sample` via
        ## `improve.marginals`: F^{-1}(F(x)) where F is ecdf of samples and
        ## F^{-1} is quantile fxn from approx's marginal.
        for (i in seq_len(ncol(samples))) {
            if (!originalScale || Rapprox$paramNodesIndices[i] > 0) {
                empirQuantiles <- ecdf(samples[, i])(samples[, i])
                idx <- ifelse(originalScale, Rapprox$paramNodesIndices[i], i)

                quantiles <- estimateQuantiles(summary$marginalsApprox[[idx]],
                                               summary$paramTransforms[[idx]],
                                               empirQuantiles)
                samples[, i] <- quantiles
            }
        }
    }

    if (originalScale)
        colnames(samples) <- Rapprox$model$expandNodeNames(Rapprox$innerMethods$paramNodes,
                                                           returnScalarComponents = TRUE)
    ## TODO: check that INLA's improve.marginals does nothing for non 1:1
    ## cases.
    summary$paramSamples <- samples
    return(samples)
}

## Potentially called from runNestedApprox or independently.
sampleLatentNodes <- function(summary, n = 1000, includeParams = FALSE) {
    Rapprox <- summary$approx$Robject
    samples <- summary$approx$simulateLatentEffects(n)
    colnames(samples) <- c("index", Rapprox$innerMethods$reNodesAsScalars_vec)
    if (includeParams) {
        paramValues <- getParamGrid()  ## `getParamGrid` needs to be written as a nestedApprox nf method, returning `theta_grid$nodes`.
        paramValues <- t(apply(paramValues, 1, Rapprox$innerMethods$paramsTransform))
        paramSamples <- paramValues[samples[, "index"], ]
        colnames(paramSamples) <- Rapprox$paramNodesComponents
        samples <- cbind(samples, paramSamples)
    }
    summary$samples <- samples[, -1]
    return(summary$samples)
}



## This allows user to provide different quantiles or expectations of interest
## after running `runNestedApprox`.
updateMarginalSummaries <- function(summary, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
                                    functionals = NULL, functionalsArgs = NULL, functionalsScale = "original") {
    cnt <- 0
    for (i in seq_along(summary$marginalsRaw)) {
        if (!is.null(summary$marginalsRaw[[i]])) {
            cnt <- cnt + 1
            summary$quantileEsts[[cnt]] <- estimateQuantiles(summary$marginalsApprox[[i]],
                summary$indivParamTransforms[[i]], quantiles)
            expectations[[cnt]] <- estimateExpectations(marginalsApprox[[i]], paramTransforms[[i]],
                functionals = functionals, functionalsArgs = functionalsArgs, scale = functionalsScale)
        }
    }
    invisible(summary)
}
