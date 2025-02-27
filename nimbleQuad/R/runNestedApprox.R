## Code for main user interface for NIMBLE's nested approximation.

## NOTE: none of this code has yet been run on any examples and
## there will be various bugs.

## Uses core algorithm code in `approxPosterior.R` and code
## for marginal summaries in `approxSummaries.R`.

### Example workflow
## Rapprox <- buildNestedApprox(model)
## capprox <- compileNimble(Rapprox, project = model)
## result <- runNestedApprox(capprox)
## improveMarginals(result, nodes = 'sigma')
## paramSamples = sampleParamNodes(result, n=1000)
## samples = sampleLatentNodes(result, n=1000)
## 

## alternatively, with functions converted to class methods:
## result$improveMarginals(nodes = 'sigma')
## result$sampleParamNodes(n=1000)

## Class for holding nestedApprox object and various outputs/summaries computed from it
## when running `runNestedApprox` or individual functions that manipulate the approximation.
#' @importFrom R6 R6Class
approxSummary <- R6Class("approxSummary",
   public = list(
       initialize = function(approx, quantiles, expectations,
                             marginalsApprox, marginalsRaw,
                             indivParamTransforms, originalScale,
                             marginalLogLik, marginalLogLik_improved,
                             samples, paramSamples) {
           private$approx <- approx
           private$quantiles <- quantiles
           private$expectations <- expectations
           private$marginalsApprox <- marginalsApprox
           private$marginalsRaw <- marginalsRaw
           private$indivParamTransforms <- indivParamTransforms
           private$originalScale <- originalScale
           private$marginalLogLik <- marginalLogLik
           private$marginalLogLik_improved <- marginalLogLik_improved
           private$samples <- samples
           private$paramSamples <- paramSamples
       },
       generateParamsMatrix = function() {
           ## TODO: do I need `private$`?
           Rapprox <- ifelse(is(private$approx, "NestedApprox"), private$approx, private$approx$Robject)
           if(private$originalScale) {
               paramNames <- paste0("param", Rapprox$innerMethods$nPar)
           } else {
               paramNames <- names(private$quantiles)
           }
           first <- which(!sapply(private$quantiles, is.null))[1]
           qs <- quantiles[[first]]$quantiles
           first <- which(!sapply(private$expectations, is.null))[1]
           exps <- expectations[[first]]
           
           ## Form tabular info on expectations and quantiles as a dataframe
           ## (or could be a matrix as with INLA output),
           ## which should print nicely.
           params <- data.frame()
           for(i in seq_along(exps)) 
               params[[names(exps)[i]]] <- sapply(private$expectations, `[`, i)
                                        # params <- data.frame(mean = sapply(private$marginalsSummary, `[[`, 'mean'),
                                        #                     sd = sapply(private$marginalsSummary, `[[`, 'sd'),
                                        #                     row.names = paramNames)
            for(i in seq_along(qs)) 
                params[[names(qs)[i]]] <- sapply(private$quantiles, `[`, i)
            row.names(params) <- paramNames
            private$params <- params
        },
        print = function() {
            cat("Model (hyper)parameters: \n")
            if(is.null(private$params))
                private$params <- generateParamsMatrix()
            print(private$params)

            cat("Marginal log-likelihood (asymmetric Gaussian approximation): ", private$marginalLogLik, "\n")
            if(!is.na(private$marginalLogLik_improved))
                cat("Marginal log-likelihood (AGHQ): ", private$marginalLogLik_improved, "\n")
                
            invisible(self)
        }
    ),

    ## private methods and functions not accessible externally
    private = list(
        approx = NULL,
        marginalsSummary = NULL,
        marginalsApprox = NULL,
        marginalsRaw = NULL,
        indivParamTransforms = NULL,
        marginalLogLik = NULL,
        samples = NULL,
        params = NULL
    )
)   


## This contains initial steps not in Paul's version to handle node determination
## and mapping of node element names to indices for 1:1 marginalization work.
## Need to integrate this into his full function.
buildNestedApprox_v2 <- function(model, paramNodes, latentNodes, control = list()) {

    hyperGridRule <- extractControlElement(control, 'hyperGridRule', 'CCD')     ## Default rule for outer grid.

    ## TODO: check if handling of no `paramNodes` or `latentNodes` is correct.
    margNodes <- splitLatents(model, paramNodes, latentNodes)
    paramNodes <- margNodes$paramNodes
    latentNodes <- margNodes$randomEffectsNodes

    ## This is done later in Paul's code
    paramsTransform <- parameterTransform(model, paramNodes, control = list(allowDeterm = FALSE))
    
    ## Set up mapping of parameter names to indices of transformed elements for 1:1 cases for
    ## determination of parameters for which approximate marginals are possible and for
    ## use when users request marginals by node name.
    paramNodesComponents <- model$expandNodeNames(paramNodes, returnScalarComponents = TRUE)
    paramNodesIndices <- seq_along(paramNodesComponents)

    if(any(paramsTransform$transformType > 9, na.rm = TRUE)) 
        stop("buildNestedApprox: Unknown parameter transform type: ",
             paste0(paramsTransform$transformType[paramsTransform$transformType > 9], collapse = ', '))

    mapping <- paramsTransform$transformData
    for(idx in seq_len(paramsTransform$nNodes)) {
        if(paramsTransform$transformType < 7) {  
            paramNodeIndices[mapping[idx,1]] <- mapping[idx,3]
        } else {
            paramNodeIndices[mapping[idx,1:2]] <- 0
        }
    }

    setupOutputs(paramNodesComponents, paramNodesIndices)

    
    ## Default outer grid to CCD unless low dimensional.
    if(!'hyperGridRule' %in% names(control)) 
        hyperGridRule <- ifelse(length(paramNodes) >= 3, 'CCD', 'AGHQ')

    ## Presumably we allow CCD if user requests for nParam = 2.
    ## Do we allow CCD if user requests for nParam = 1?

    ## skewOuterGrid=TRUE
    theta_grid <- configureQuadGrid(d = 1, nQuad_ = nQuadOuter, quadRule = hyperGridRule, 
                                    control = list(quadRules = allGridRules))

    ## call `buildAGHQ`.
    innerMethods <- buildAGHQ(model, nQuadInner, paramNodes, latentNodes, 
                          margNodes$calcNodes, margNodes$calcNodesOther, control)
  
    ## TODO: Merge in with Paul's existing code.
}

## Main user-facing function for running a nested approximation and getting a results summary.
## Note in roxygen that `functionalsScale` will have no effect if `originalScale=FALSE`, in
## which case functional will compute on the transformed parameter values.
runNestedApprox <- function(approx, quantiles = c(.025,.25,.5,.75,.975), originalScale = TRUE,
                            nSamplesLatents = 0, nSamplesParams = 0, functionals = NULL, functionalsArgs = NULL, functionalsScale = "original") { 
    Rapprox <- ifelse(is(approx, "NestedApprox"), approx, approx$Robject)

    ## CHECK: do we already have this info in `Rapprox`?
    nParamTrans <- Rapprox$innerMethods$npar

    marginalsRaw <- list(); length(marginalsRaw) <- nParamTrans
    marginalsApprox <- list(); length(marginalsApprox) <- nParamTrans
    indivParamsTransform <- list(); length(indivParamsTransform) <- nParamTrans
    
    nMarginalsReport <- 1

    for(i in seq_along(nParamTrans)) {
        ## Save computation of marginals, except non-1:1 transformed parameters.
        if(!originalScale || i %in% Rapprox$paramNodesIndices) {
            marginalsRaw[[i]] <- approx$findMarginalHyperIntFree(Rapprox$paramNodesIndices[i])
            nMarginalsReport <- nMarginalsReport + 1
        }
    }

    cnt <- 1
    marginalsApprox <- quantiles <- indivParamTransforms <- list()
    length(marginalsSummary) <- nMarginalsReport

    cnt <- 0
    for(i in seq_along(nParamTrans)) {
        if(!is.null(marginalsRaw[[i]])) {
            cnt <- cnt+1
            if(originalScale)
                indivParamTransforms[[i]] <- parameterTransform(Rapprox$model, Rapprox$paramNodesComponents[i])
            marginalsApprox[[i]] <- fitMarginalSpline(marginalsRaw[[i]])
            quantiles[[cnt]] <- estimateQuantiles(marginalsApprox[[i]], paramTransforms[[i]], quantiles)
            expectations[[cnt]] <- estimateExpectations(marginalsApprox[[i]], paramTransforms[[i]], 
                                                        functionals = functionals, functionalsArgs = functionalsArgs, scale = functionalsScale)
        }
    }

    if(originalScale) {
        names(quantiles) <- Rapprox$paramNodesComponents[Rapprox$paramNodesIndices > 0]
        names(expectations) <- Rapprox$paramNodesComponents[Rapprox$paramNodesIndices > 0]
    }

    marginalLogLik <- approx$calcMarginalLogLikApprox()

    summary <- approxSummary$new(approx, quantiles, expectations, marginalsApprox, marginalsRaw,
                indivParamTransforms, originalScale,
                marginalLogLik, NA, NULL, NULL)

    ## With only one parameter, computations should not be so slow, so
    ## go ahead and use better marginal and logLik estimates.
    if(length(Rapprox$paramNodesComponents) == 1) {
        improveMarginals(summary, Rapprox$paramNodesComponents[1], nMarginalGrid = 3, 
                         functionals = functionals, functionalsArgs = functionalsArgs,
                         functionalsScale = functionalsScale)
        ## TODO: make sure that `calcMarginalLogLikQuad()` works if
        ## `calcHyperGrid` has not yet been called.
        summary$marginalLogLik_improved <-calcMarginalLogLikQuad()   
    }

    ## This is expensive. Avoid if user only needs parameter inference.
    ## How do we have user tell us whether to `includeParams`?
    ## Perhaps tell them to use more manual workflow if they need that.
    if(nSamplesLatents) {
        summary$samples <- sampleLatentNodes(approx, n = nSamplesLatents, includeParams = FALSE)
        if(Rapprox$hyperGridRule == 'AGHQ')
            marginalLogLik_improved <- approx$calcMarginalLogLikQuad()  ## TODO: should be named "get"?
    } 

    if(nSamplesParams) {
        summary$paramSamples <- sampleParamNodes(approx, n = nSamplesParams)
    } 
        
    return(summary)
}



## This uses d-1 dimensional AGHQ to get improved univariate marginal estimates for
## parameters
## Should it be called `improveParamMarginals`?
improveMarginals  <- function(summary, nodes, nMarginalGrid = 3, 
                              functionals = NULL, functionalsArgs = NULL, functionalsScale = "original") {
    Rapprox <- ifelse(is(summary$approx, "NestedApprox"), summary$approx, summary$approx$Robject)
    originalScale <- summary$originalScale
    nodes <- Rapprox$model$expandNodeNames(nodes, returnScalarComponents = TRUE)
    if(!all(nodes %in% Rapprox$paramNodesComponents)) {
        badNodes <- nodes[!nodes %in% Rapprox$paramNodesComponents]
        stop("improveMarginals: nodes `", paste0(badNodes, sep = "`, `"), "` are not parameter elements or are not involved in 1:1 parameter transformation, so marginals cannot be estimated by analytic approximation. In the latter case, use `sampleParamNodes` for inference.")
    }
    marginalsRaw <- list()
    length(marginalsRaw) <- Rapprox$innerMethods$npar
    for(i in seq_len(nodes)) {
        ## Improve marginal and insert into raw and summary objects.
        if(originalScale) {
            idx = Rapprox$paramNodesIndices[which(nodes[i] == paramNodesComponents)]
        } else idx <- i
        ## TODO: what should `nQuad` and should it be what is set in `buildNestedApprox` setup. And make sure it is available via `approx$nQuadMarginal`.
        summary$marginalsRaw[[idx]] <- findMarginalPosteriorDensity(idx, nPts = nMarginalGrid, nQuad = summary$approx$nQuadMarginal)
        if(originalScale)
            paramTransform <- parameterTransform(Rapprox$model, nodes[i]) else paramTransform <- NULL
        summary$marginalsApprox[[idx]] <- fitMarginalSpline(summary$marginalsRaw[[idx]])
        ## TODO: check that nodes[i] name will correspond to element in `marginalsSummary`
        if(originalScale)
            idx2 <- nodes[i] else idx2 <- i
        summary$quantiles[[idx2]] <- estimateQuantiles(summary$marginalsApprox[[idx]], summary$indivParamTransforms[[idx]])
        summary$expectations[[idx2]] <- estimateExpectations(summary$marginalsApprox[[idx]], summary$indivParamTransforms[[idx]],
                                                             functionals = functionals, functionalsArgs = functionalsArgs, scale = functionalsScale)
    }
    summary$generateParamsMatrix()
    return(summary)  
}


## This uses whatever marginals (asymm Gaussian approx or improved) are in `summary`.
## Potentially called from runNestedApprox or independently.
sampleParamNodes <- function(summary, n = 1000, matchMarginals = TRUE) {
    Rapprox <- ifelse(is(summary$approx, "NestedApprox"), summary$approx, summary$approx$Robject)
    originalScale <- summary$originalScale
    
    samplesTrans <- summary$approx$simulateHyperParams(n)
    if(originalScale) {
        samples <- t(apply(samplesTrans, 1, Rapprox$parameterTransform$inverseTransform))
    } else samples <- samplesTrans

    if(matchMarginals) {
        ## Copula approac based on `inla.hperpar.sample` via `improve.marginals`:
        ## F^{-1}(F(x)) where F is ecdf of samples and F^{-1} is quantile fxn from approx's marginal.
        for(i in seq_len(ncol(samples))) {
            if(!originalScale || Rapprox$paramNodesIndices[i] > 0) {
                empirQuantiles <- ecdf(samples[,i])(samples[,i])
                idx <- ifelse(originalScale, Rapprox$paramNodesIndices[i], i)
                
                quantiles <- estimateQuantiles(summary$marginalsApprox[[idx]], summary$paramTransforms[[idx]], empirQuantiles)
                samples[ , i] <- quantiles
            }
        }
    }

    if(originalScale)
        colnames(samples) <- Rapprox$model$expandNodeNames(Rapprox$innerMethods$paramNodes, returnScalarComponents = TRUE)
    ## TODO: check that INLA's improve.marginals does nothing for non 1:1 cases.
    summary$paramSamples <- samples
    return(samples)
}

## Potentially called from runNestedApprox or independently.
sampleLatentNodes <- function(summary, n = 1000, includeParams = FALSE) {
    Rapprox <- ifelse(is(summary$approx, "NestedApprox"), summary$approx, summary$approx$Robject)
    samples <- summary$approx$sampleLatentNodes(n)
    colnames(samples) <- c('index', Rapprox$innerMethods$reNodesAsScalars_vec)
    if(includeParams) {
        paramsTransform <- parameterTransform(model, Rapprox$paramNodes, control = list(allowDeterm = FALSE))
        paramValues <- getParamGrid()  ## `getParamGrid` needs to be written as a nestedApprox nf method, returning `theta_grid$nodes`.
        paramValues <- t(apply(paramValues, 1, paramsTransform$inverseTransform))
        paramSamples <- paramValues[samples[ , 'index'], ]
        colnames(paramSamples) <- Rapprox$paramNodesComponents
        samples <- cbind(samples, paramSamples)
    }
    summary$samples <- samples[ , -1]
    return(summary$samples)
}



## This allows user to provide different quantiles or expectations of interest after running `runNestedApprox`.
updateMarginalSummaries <- function(summary, quantiles = c(.025,.25,.5,.75,.975),
                                    functionals = NULL, functionalsArgs = NULL, functionalsScale = "original") {
    cnt <- 0
    for(i in seq_along(summary$marginalsRaw)) {
        if(!is.null(summary$marginalsRaw[[i]])) {
            cnt <- cnt+1
            summary$quantiles[[cnt]] <- estimateQuantiles(summary$marginalsApprox[[i]], summary$indivParamTransforms[[i]], quantiles)
            expectations[[cnt]] <- estimateExpectations(marginalsApprox[[i]], paramTransforms[[i]], 
                                                        functionals = functionals, functionalsAargs = functionalsArgs, scale = functionalsScale)
        }
    }
    return(summary)    
}
