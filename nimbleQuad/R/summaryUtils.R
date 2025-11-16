marginalSplineR <- nimbleRcall(function(theta = double(1), logdens = double(1)) {},
                               Rfun = "marginalSpline", returnType = double(2))

fitMarginalSpline <- function(gridded, normalize = TRUE, xnew = NULL, refine = TRUE, extend = TRUE) {
    theta <- gridded[, 1]
    logdens <- gridded[, 2]
    n <- length(theta)
    rn <- range(theta)
    rnl <- diff(rn)
    if(!extend)
        rnl <- 0   # Don't be cautious if refining the interval.
    thetarange <- c(min(rn) - rnl/2, max(rn) + rnl/2)
    finegrid <- seq(thetarange[1], thetarange[2], length.out = 1000)  ## This is based off of Stringer for a fine grid.
    if (n <= 3) {
        log_pdf <- as.function(polynom::poly.calc(x = theta, y = logdens))
        logPDF <- log_pdf(finegrid)
        if(!is.null(xnew))
            logPDFnew <- log_pdf(xnew)
    } else {
        ss <- splines::interpSpline(theta, logdens, bSpline = TRUE, sparse = FALSE)
        if (isS4(co <- ss[["coefficients"]]))
            ss[["coefficients"]] <- as.vector(co)  ## Not sure why this might be necessary.
        logPDF <- as.numeric(stats::predict(ss, finegrid)$y)
        if(!is.null(xnew))
            logPDFnew <- as.numeric(stats::predict(ss, xnew)$y)
    }
    ## Normalize the PDF.
    pdf <- exp(logPDF)
    # Trapezoidal rule: could use Simpson as (2M+T)/3, but it would require additional spline evaluation.
    trapezoids <- diff(finegrid) * (pdf[-length(pdf)] + pdf[-1])/2  
    if (normalize) {
        norm <- sum(trapezoids) 
        pdf <- pdf/norm
        logPDF <- logPDF - log(norm)
    } else norm <- 1
    cdf <- c(0, cumsum(trapezoids)/norm)

    if(refine) {
        core <- which(cdf > 0.00001 & cdf < 0.99999)
        return(fitMarginalSpline(cbind(finegrid[core], logPDF[core]), xnew = xnew, refine = FALSE, extend = FALSE))
    }
    
    ## Return the gridded distribution information or evaluation at provided points.
    if(is.null(xnew))
        return(cbind(finegrid, pdf, cdf)) else return(logPDFnew - log(norm))
}


estimateQuantiles <- function(marginalApprox, transform = NULL,
                              quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975)) {
    if(any(quantiles < 0 | quantiles > 1))
        stop("quantile probabilities are outside [0,1]")
    finegridTrans <- marginalApprox[, "finegrid"]
    cdf <- marginalApprox[, "cdf"]

    ## For now use smoothing spline on quantile function on transformed (theta) scale.
    ## Need to avoid cdf values too close together or spline fit will fail.
    used <- cdf > 0.001 & cdf < 0.999  
    ss <- splines::interpSpline(cdf[used], finegridTrans[used], bSpline = TRUE, sparse = FALSE)
    quantsTrans <- stats::predict(ss, quantiles)$y
    if (!is.null(transform)) {
        quants <- sapply(quantsTrans, transform$inverseTransform)
    } else quants <- quantsTrans
    names(quants) <- paste0(100*quantiles, "%")

    return(quants)
}

## Trapezoidal rule for estimating expectation.
runTrapezRule <- function(grid, pdf, functional) {
    n <- length(grid)
    return(sum(diff(grid) * (pdf[-n] * functional[-n] + pdf[-1] * functional[-1])/2))
}


## `functionals` should be a named list of functions.  If user wants any of their
## functionals to take additional args, all of them must take ...
estimateExpectations <- function(marginalApprox, transform = NULL, functional = NULL, ...) {

    finegridTrans <- marginalApprox[, "finegrid"]
    pdfTrans <- marginalApprox[, "pdf"]
    if (!is.null(transform)) {
        finegrid <- sapply(finegridTrans, transform$inverseTransform)
    } else {
        finegrid <- finegridTrans
    }

    ## Posterior expectations: defaults
    ## Use pdf on transformed (theta) scale.
    ## If user set `originalScale = FALSE` then there will be no transform
    ## provided from `runNestedApprox`, 
    ## and this will give mean and sd on transformed scale.
    if(is.null(functional)) {
        functionalVals <- finegrid
        postMean <- runTrapezRule(finegridTrans, pdfTrans, functionalVals)
        functionalVals <- (finegrid - postMean)^2
        postVar <- runTrapezRule(finegridTrans, pdfTrans, functionalVals)
        postSD <- sqrt(postVar)

        expectations <- c(mean = postMean, sd = postSD)
    } else { ## User-defined expectation.
        functionalVals <- functional(finegrid, ...)
        expectations <- runTrapezRule(finegridTrans, pdfTrans, functionalVals)
    }
    
    ## Expectations calculated directly on original scale.
    ## postMean <-   ## sum(diff(finegrid)*(pdf[-n]*finegrid[-n]+pdf[-1]*finegrid[-1])/2)
    ## functional <- (finegrid - postMean)^2
    ## postVar <- sum(diff(finegrid)*(pdf[-n]*functional[-n]+pdf[-1]*functional[-1])/2)

    return(expectations)
}


