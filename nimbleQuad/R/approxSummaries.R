marginalSplineR <- nimbleRcall(function(theta = double(1), logdens = double(1)) {},
                               Rfun = "marginalSpline", returnType = double(2))

fitMarginalSpline <- function(gridded, normalize = TRUE) {
    theta <- gridded[, 1]
    logdens <- gridded[, 2]
    n <- length(theta)
    rn <- range(theta)
    rnl <- diff(rn)
    thetarange <- c(min(rn) - rnl/2, max(rn) + rnl/2)
    finegrid <- seq(thetarange[1], thetarange[2], length.out = 1000)  ## This is based off of Stringer for a fine grid.
    if (n <= 3) {
        log_pdf <- as.function(polynom::poly.calc(x = theta, y = logdens))
        logPDF <- log_pdf(finegrid)
    } else {
        ss <- splines::interpSpline(theta, logdens, bSpline = TRUE, sparse = FALSE)
        if (isS4(co <- ss[["coefficients"]]))
            ss[["coefficients"]] <- as.vector(co)  ## Not sure why this might be necessary.
        logPDF <- as.numeric(stats::predict(ss, finegrid)$y)
    }
    ## Normalize the PDF:
    pdf <- exp(logPDF)
    trapezoids <- diff(finegrid) * (pdf[-n] + pdf[-1])/2  # trapezoidal rule (could use Simpson as (2M+T)/3
    if (normalize)
        norm <- sum(trapezoids) else norm <- 1
    pdf <- pdf/norm
    cdf <- c(0, cumsum(trapezoids)/norm)
    ## Return the gridded distribution information.
    return(cbind(finegrid, pdf, cdf))
}
## TODO: should we refine the grid based on excluding portions with negligible
## mass?



estimateQuantiles <- function(marginalApprox, transform = NULL,
                              quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975)) {
    finegridTrans <- marginalApprox[, "finegrid"]
    cdf <- marginalApprox[, "cdf"]

    ## For now use smoothing spline on quantile function on transformed (theta)
    ## scale.
    used <- cdf > 0.001 & cdf < 0.999  # Need to avoid cdf values too close together or spline fit will fail.
    ss <- splines::interpSpline(cdf[used], finegridTrans[used], bSpline = TRUE, sparse = FALSE)
    quantsTrans <- stats::predict(ss, quantiles)$y
    if (!is.null(transform)) {
        quants <- sapply(quantsTrans, transform$inverseTransform)
    } else quants <- quantsTrans
    names(quants) <- paste0("q", quantiles)

    return(quants)
}

## Trapezoidal rule for estimating expectation.
runTrapezRule <- function(grid, pdf, functional) {
    n <- length(grid)
    return(sum(diff(grid) * (pdf[-n] * functional[-n] + pdf[-1] * functional[-1])/2))
}


## functionals should be a named list of functions.  if user wants any of their
## functionals to take additional args, all of them must take ...
estimateExpectations <- function(marginalApprox, transform = NULL, functionals = NULL,
    functionalsArgs = NULL, scale = "original") {

    finegridTrans <- marginalApprox[, "finegrid"]
    pdfTrans <- marginalApprox[, "pdf"]
    if (!is.null(transform)) {
        finegrid <- sapply(finegridTrans, transform$inverseTransform)
    } else {
        finegrid <- finegridTrans
    }

    ## Posterior expectations: defaults
    ## Use pdf on transformed (theta) scale.
    ## If user set `originalScale = FALSE` then there will be no transform, 
    ## and this will give mean and sd on transformed scale.
    functionalVals <- finegrid
    postMean <- runTrapezRule(finegridTrans, pdfTrans, functionalVals)
    functionalVals <- (finegrid - postMean)^2
    postVar <- runTrapezRule(finegridTrans, pdfTrans, functionalVals)
    postSD <- sqrt(postVar)

    expectations <- c(mean = postMean, sd = postSD)

    ## Expectations calculated directly on original scale.
    ## postMean <-   ## sum(diff(finegrid)*(pdf[-n]*finegrid[-n]+pdf[-1]*finegrid[-1])/2)
    ## functional <- (finegrid - postMean)^2
    ## postVar <- sum(diff(finegrid)*(pdf[-n]*functional[-n]+pdf[-1]*functional[-1])/2)

    ## User-defined expectations.  if `scale = 'original'` then
    ## `functionals[[i]]` is assumed to be a function of parameters on original
    ## scale.
    ## If `originalScale=FALSE` then `transform` will be NULL, so this
    ## will take functionals to operate on transformed values, regardless of
    ## `functionalsScale`.
    if (length(functionals)) {
        userExps <- rep(0, length(functionals))
        names(userExps) <- names(functionals)
        for (i in seq_along(functionals)) {
            if (scale == "original")
                input <- finegrid else input <- finegridTrans
            if (is.null(functionalsArgs)) {
                functionalVals <- functionals[[i]](input)
            } else {
                argList <- list()
                length(argList) <- length(functionalArgs[[i]]) + 1
                argList[[1]] <- input
                argList[2:length(argList)] <- functionalArgs[[i]]
                functionalVals <- do.call(functionals[[i]], argList)
            }
            userExps[i] <- runTrapezRule(finegridTrans, pdfTrans, functionalVals)
        }
        expectations <- c(expectations, userExps)
    }

    return(expectations)
}



## Original version. May not be used.
summarizeMarginalOriginal <- function(marginalApprox, transform = inverseTransform,
    logDetJacobian = logDetJacobian, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
    functionals = NULL, ...) {
    n <- nrow(marginalApprox)

    finegridTrans <- marginalApprox[, "finegrid"]
    pdfTrans <- marginalApprox[, "pdf"]
    cdf <- marginalApprox[, "cdf"]  # same for both scales

    finegrid <- transform(finegridTrans)
    pdf <- pdfTrans[, "pdf"] * logDetJac(finegrid)
    # sum(diff(finegrid)*(pdf[-n]+pdf[-1])/2) # trapezoidal rule showing
    # normalization is preserved

    ## For now use smoothing spline on quantile function on transformed (theta)
    ## scale.  How do INLA and Stringer/aghq get their quantiles? I think Paul
    ## and Chris discussed this.
    used <- cdf > 0.001 & cdf < 0.999
    ss <- splines::interpSpline(cdf[used], finegridTrans[used], bSpline = TRUE, sparse = FALSE)
    quantsTrans <- stats::predict(ss, quantiles)$y
    quants <- transform(quantsTrans)

    ## Posterior expectations: defaults Use pdf on transformed (theta) scale.
    postMean <- estimateExpectation(finegridTrans, pdfTrans, finegrid, n)
    functional <- (finegrid - postMean)^2
    postVar <- estimateExpectation(finegridTrans, pdfTrans, functional, n)
    postSD <- sqrt(postVar)

    ## transformed scale postMean <-
    ## sum(diff(finegrid)*(pdf[-n]*finegrid[-n]+pdf[-1]*finegrid[-1])/2)
    ## functional <- (finegrid - postMean)^2 postVar <-
    ## sum(diff(finegrid)*(pdf[-n]*functional[-n]+pdf[-1]*functional[-1])/2)

    ## Posterior expectations: user-defined Not sure about use of `...`.
    if (!is.null(functionals)) {
        userExpectations <- sapply(functionals, function(fun) {
            functional <- fun(finegrid, ...)
            return(estimateExpectation(finegridTrans, pdfTrans, functional, n))
        })
    } else userExpectations <- NULL
    ## Also decide what to return in terms of pdf information. (See
    ## INLA/Stringer.)
    return(list(quantiles = quants, postMean = postMean, postVar = postVar, postSD = postSD,
        userExpectations = userExpectations))
}
