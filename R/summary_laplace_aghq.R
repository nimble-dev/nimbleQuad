

#' Summarize results from Laplace or adaptive Gauss-Hermite quadrature approximation
#'
#' Process the results of the `findMLE` method of a nimble Laplace or AGHQ approximation
#' into a more useful format.
#'
#' @param laplace The Laplace approximation object, typically the compiled one.
#'   This would be the result of compiling an object returned from
#'   `buildLaplace`.
#'
#' @param AGHQ Same as \code{laplace}. Note that `buildLaplace` and
#'   `buildAGHQ` create the same kind of algorithm object that can be used
#'   interchangeably. `buildLaplace` simply sets the number of quadrature points
#'   (`nQuad`) to 1 to achieve Laplace approximation as a special case of AGHQ.
#'
#' @param MLEoutput The maximum likelihood estimate using Laplace or AGHQ,
#'   returned from e.g. `approx$findMLE(...)`, where \code{approx} is the
#'   algorithm object returned by `buildLaplace` or `buildAGHQ`, or (more
#'   typically) the result of compiling that object with `compileNimble`. See
#'   `help(buildLaplace)` for more information.
#'
#' @param originalScale Should results be returned using the original
#'   parameterization in the model code (TRUE) or the potentially transformed
#'   parameterization used internally by the Laplace approximation (FALSE).
#'   Transformations are used for any parameters and/or random effects that have
#'   constrained ranges of valid values, so that in the transformed parameter
#'   space there are no constraints. 
#'
#' @param randomEffectsStdError If TRUE, calculate the standard error of the
#'   estimates of random effects values.
#'
#' @param jointCovariance If TRUE, calculate the joint covariance matrix of
#'   the parameters and random effects together. If FALSE, calculate the 
#'   covariance matrix of the parameters.
#'
#' @details
#'
#' The numbers obtained by this function can be obtained more directly by
#' `approx$summary(...)`. The added benefit of `summaryLaplace` is to arrange
#' the results into data frames (for parameters and random effects), with row
#' names for the model nodes, and also adding row and column names to the
#' covariance matrix.
#'
#' @return
#'
#' A list with data frames `params` and `randomEffects`, each with columns for
#' `estimate` and (possibly) `se` (standard error) and row names for model
#' nodes, a matrix `vcov` with the covariance matrix with row and column names,
#' and `originalScale` with the input value of `originalScale` so it is recorded
#' for later use if wanted.
#'
#' @aliases summaryAGHQ
#'
#' @name summaryLaplace
#'
#' @export
summaryLaplace <- function(laplace, MLEoutput,
                           originalScale =TRUE,
                           randomEffectsStdError = FALSE,
                           jointCovariance = FALSE) {
    summary <- laplace$summary(MLEoutput, originalScale = originalScale,
                               randomEffectsStdError = randomEffectsStdError,
                               jointCovariance = jointCovariance)
    paramNames <- summary$params$names
    paramEsts <- summary$params$estimates
    if(length(paramEsts) < length(paramNames)) paramNames <- paramNames[1:(length(paramNames)-1)]
    names(paramEsts) <- paramNames
    stdErrParams <- summary$params$stdErrors
    paramsDF <- data.frame(estimate = paramEsts, se = stdErrParams, row.names = paramNames)

    REnames <- summary$randomEffects$names
    REests <- summary$randomEffects$estimates
    if(length(REests) < length(REnames)) REnames <- REnames[1:(length(REnames)-1)]
    REstdErrs <- summary$randomEffects$stdErrors
    if(length(REstdErrs))
        REDF <- data.frame(estimate = REests, se = REstdErrs, row.names = REnames)
    else
        REDF <- data.frame(estimate = REests, row.names = REnames)

    vcov <- summary$vcov
    if (dim(vcov)[1] == length(paramNames)) {
        colnames(vcov) <- rownames(vcov) <- c(paramNames)
    } else {
        colnames(vcov) <- rownames(vcov) <- c(paramNames, REnames)
    }
    list(params = paramsDF,
         randomEffects = REDF,
         vcov = vcov,
         originalScale = originalScale)
}

#' @rdname summaryLaplace
#' @export
summaryAGHQ <- function(AGHQ, MLEoutput,
                        originalScale =TRUE,
                        randomEffectsStdError = FALSE,
                        jointCovariance = FALSE) {
    summaryLaplace(AGHQ, MLEoutput, originalScale, randomEffectsStdError, jointCovariance)
}
