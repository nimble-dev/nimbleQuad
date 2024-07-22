

#' Combine steps of running Laplace or adaptive Gauss-Hermite quadrature approximation
#'
#' Use an approximation (compiled or uncompiled) returned from
#' `buildLaplace` or `buildAGHQ` to find the maximum likelihood estimate and return it
#' with random effects estimates and/or standard errors.
#'
#' @aliases runAGHQ runLaplace
#'
#' @param laplace A (compiled or uncompiled) nimble laplace approximation object
#'   returned from `buildLaplace` or `buildAGHQ`. These return the same type of
#'   approximation algorithm object. `buildLaplace` is simply `buildAGHQ`
#'   with `nQuad=1`.
#'
#' @param AGHQ Same as \code{laplace}.
#'
#' @param method Optimization method for outer optimization. See \code{method}
#'   argument to \code{findMLE} method in \code{\link{buildLaplace}}.
#'
#' @param pStart Initial values for parameters to begin optimization search for
#'   the maximum likelihood estimates. If omitted, the values currently in the
#'   (compiled or uncompiled) model object will be used.
#'
#'
#' @param originalScale If \code{TRUE}, return all results on the original scale
#'   of the parameters and/or random effects as written in the model. Otherwise,
#'   return all results on potentially unconstrained transformed scales that are
#'   used in the actual computations. Transformed scales (parameterizations) are
#'   used if any parameter or random effect has contraint(s) on its support
#'   (range of allowed values). Default = \code{TRUE}.
#'
#' @param randomEffectsStdError If \code{TRUE}, include standard errors for the
#'   random effects estimates. Default = \code{TRUE}.
#'
#' @param jointCovariance If \code{TRUE}, return the full joint covariance
#'   matrix (inverse of the Hessian) of parameters and random effects. Default =
#'   \code{FALSE}.
#'
#' @details
#'
#' Adaptive Gauss-Hermite quadrature is a generalization of Laplace
#' approximation. \code{runLaplace} simply calles \code{runAGHQ} and provides a
#' convenient name.
#'
#' These functions manage the steps of calling the `findMLE` method to obtain
#' the maximum likelihood estimate of the parameters and then the
#' `summaryLaplace` function to obtain standard errors, (optionally) random
#' effects estimates (conditional modes), their standard errors, and the full
#' parameter-random effects covariance matrix.
#'
#' Note that for `nQuad > 1` (see \code{\link{buildAGHQ}}), i.e., AGHQ with
#' higher order than Laplace approximation, maximum likelihood estimation is
#' available only if all random effects integrations are univariate. With
#' multivariate random effects integrations, one can use `nQuad > 1` only to
#' calculate marginal log likelihoods at given parameter values. This is useful
#' for checking the accuracy of the log likelihood at the MLE obtained for
#' Laplace approximation (`nQuad == 1`). `nQuad` can be changed using the
#' `updateSettings` method of the approximation object.
#'
#' See \code{\link{summaryLaplace}}, which is called for the summary components.
#'
#' @return
#'
#' A list with elements \code{MLE} and \code{summary}.
#'
#' \code{MLE} is the result of the \code{findMLE} method, which contains the
#' parameter estimates and Hessian matrix. This is considered raw output, and
#' one should normally use instead the contents of \code{summary}. (For example
#' not that the Hessian matrix in \code{MLE} may not correspond to the same
#' scale as the parameter estimates if a transformation was used to operate in
#' an unconstrained parameter space.)
#'
#' \code{summary} is the result of \code{summaryLaplace} (or equivalently
#' \code{summaryAGHQ}), which contains parameter estimates and standard errors,
#' and optionally other requested components. All results in this object will be
#' on the same scale (parameterization), either original or transformed, as
#' requested.
#'
#' @export
runLaplace <- function(laplace, pStart, method = "BFGS",
                       originalScale = TRUE,
                       randomEffectsStdError = TRUE,
                       jointCovariance = FALSE) {
    if(missing(pStart)) pStart <- Inf # code to use values in model
    runAGHQ(AGHQ = laplace, pStart, method, originalScale, randomEffectsStdError,
            jointCovariance)
}

#' @rdname runLaplace
#' @export
runAGHQ <- function(AGHQ, pStart, method = "BFGS",
                    originalScale = TRUE,
                    randomEffectsStdError = TRUE,
                    jointCovariance = FALSE) {
    if(missing(AGHQ)) stop('must provide a NIMBLE Laplace or AGHQ algorithm')
    if(!identical(nfGetDefVar(AGHQ, 'name'), 'AGHQ'))
        stop('AGHQ or laplace argument must be a NIMBLE Laplace or AGHQ algorithm (compiled or uncompiled) from buildLaplace or buildAGHQ.')
    if(!is.Cnf(AGHQ))
        messageIfVerbose('  [Warning] Running an uncompiled Laplace or AGHQ algorithm.',
                         ' Use compileNimble() for faster execution.')

    if(missing(pStart)) pStart <- Inf # code to use values in the model

    opt <- try(AGHQ$findMLE(pStart = pStart, method = method, hessian = TRUE))
    if(inherits(opt, "try-error"))
        stop("method findMLE had an error.")

    summary  <- try(summaryLaplace(laplace=AGHQ, MLEoutput=opt,
                                   originalScale=originalScale,
                                   randomEffectsStdError=randomEffectsStdError,
                                   jointCovariance=jointCovariance))
    if(inherits(summary, "try-error")) {
        warning("summaryLaplace had an error. Only the MLE result will be returned.")
        summary <- NULL
    }
    list(MLE = opt, summary=summary)
}

