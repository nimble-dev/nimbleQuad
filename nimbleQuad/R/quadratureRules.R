## Quadrature Rules + Grids - Rules are choosing between AGHQ and CCD etc.
## Grids are implementations of the rule to generate the actual nodes and
## weights.

## Base class for nimble function list quadrature rules.
QUAD_RULE_BASE <- nimbleFunctionVirtual(
    name = "QUAD_RULE_BASE",
    run = function() {
    },
    methods = list(
        buildGrid = function(levels = integer(0, default = 0), d = integer(0, default = 1)) {
            returnType(double(2))
        }
    )
)


#' Nimble List Quadrature Data type
#'
#' Creates a quadrature nimble list type to be used internally and for making new custom
#' quadrature rules to marginalize random effects and for posterior approximations.  
#'
#' @details
#' 
#'
#' List is generated with three data types. An integer that is the mode index `modeIndex` that indicates
#' which quadrature node is the mode, values are all zero. A numeric vector, `wgts`, that is a weight for each
#' quadrature node. A matrix, `nodes`, that are the quadrature nodes made by the rule that are of dimension `nQ` rows
#' and `d` dimension columns.
#'
#' @author Paul van Dam-Bates
#' @export
quadGridListDef <- nimbleList(
    modeIndex = integer(0),
    wgts = double(1),
    nodes = double(2),
    name = "quadGridList")

#' Gauss-Hermite Quadrature Points in one dimension
#'
#' Generates GH quadrature weights and nodes for integrating a general univariate function from -Inf to Inf.
#' 
#' @param levels How many quadrature points to generate.
#'
#' @details
#' This function generates Gauss-Hermite points and returns a matrix with the first column as weights and
#' the second nodes. Some numerical issues occur in Eigen decomposition making the grid weights only accurate up to 35 quadrature nodes.
#'
#' @author Paul van Dam-Bates
#'
#' @references
#'
#' Golub, G. H. and Welsch, J. H. (1969). Calculation of Gauss Quadrature Rules. 
#' Mathematics of Computation 23 (106): 221-230.
#'
#' Liu, Q. and Pierce, D. A. (1994). A Note on Gauss-Hermite Quadrature. Biometrika, 81(3) 624-629.
#'
#' Jackel, P. (2005). A note on multivariate Gauss-Hermite quadrature. London: ABN-Amro. Re.
#'
#' @export
quadGH <- nimbleFunction(run = function(levels = integer(0, default = 1), type = character(0,default = "GHe")) {
    odd <- TRUE
    if (levels %% 2 == 0)
        odd <- FALSE

    res <- matrix(0, nrow = levels, ncol = 2)
    if (levels == 1) {
        ## Laplace Approximation:
        res[1, 2] <- 0
        res[1, 1] <- 1
    } else {
        i <- 1:(levels - 1)
        dv <- sqrt(i/2)
        ## Recreate pracma::Diag for this problem.
        if (levels == 2)
            fill_diag <- matrix(dv, 1, 1) else fill_diag <- diag(dv)

        y <- matrix(0, nrow = levels, ncol = levels)
        y[1:(levels - 1), 1:(levels - 1) + 1] <- fill_diag
        y[1:(levels - 1) + 1, 1:(levels - 1)] <- fill_diag
        E <- eigen(y, symmetric = TRUE)
        L <- E$values  # Always biggest to smallest.
        V <- E$vectors
        inds <- numeric(value = 0, length = levels)
        for (j in seq_along(L)) inds[j] <- levels - j + 1  ## Is this an efficient way to do it?
        x <- L[inds]
        ## Make mode hard zero. We know nQ is odd and > 1.
        if (odd)
            x[ceiling(levels/2)] <- 0
        V <- t(V[, inds])

        w <- V[, 1]^2
        
        res[, 1] <- w
        res[, 2] <- x
    }
    ## For GHN
    ## Update nodes and weights in terms of z = x/sqrt(2) and include
    ## Gaussian kernel in weight to integrate an arbitrary function. (i.e. excludes normal distr)
    if(type == "GHe"){
      res[,1] <- res[,1] * sqrt(2 * pi) * exp(res[,2]^2)
      res[,2] <- res[,2] * sqrt(2)
    }

    returnType(double(2))
    return(res)
})

#' Gauss-Hermite Quadrature Rule for Laplace and Approx Posteriors
#'
#' Generate a 1 dimension GHQ grid via a nimble function list.
#' 
#' @param levels Length of GHQ
#'
#' @details
#' This function a 1D Gauss-Hermite Quadrature Grid (nodes and weights). When choosing `type = "GHe"`, 
#' the nodes are adjusted to integrate a general function, adjusting the points by
#' the sqrt(2) and the weights by sqrt(2*pi) * exp(x^2). It cannot be compiled without being
#' included within a virtual nimble list "QUAD_RULE_BASE".
#'
#' @author Paul van Dam-Bates
#'
#' @references
#'
#' Jackel, P. (2005). A note on multivariate Gauss-Hermite quadrature. London: ABN-Amro. Re.
#'
#' @export
quadRule_GH = nimbleFunction(
    contains = QUAD_RULE_BASE,
    name = "quadRule_GH",
    setup = function(type = "GHe") {},
    run = function() {},
    methods = list(
        buildGrid = function(levels = integer(0, default = 0), d = integer(0, default = 1)) {
            returnType(double(2))

            if (levels > 35) {
                print("Warning:  More than 35 quadrature nodes per dimension is not supported. Setting levels to 35.")
                levels <- 35
            }
            if (levels == 0) {
                print("Warning:  No default number of quadrature points given. Assuming levels = 3 per dimension.")
                levels <- 3
            }
            nodes <- quadGH(levels, type)
                
            return(nodes)
        }
      )
    )


#' Drop Algorithm to generate permutations of dimension d with a fixed sum.
#'
#' Generates a matrix of all permutations of 'd' cols that sum to 'order' with no zeros.
#' 
#' @param d Number of columns (dimensions)
#' @param order Row sum to permute over.
#'
#' @details
#' This function generates permutation matrix in order to be used for sparse grid quadrature building.
#' It is adapted from the library `mvQuad`.
#'
#' @author Paul van Dam-Bates
#'
#' @export
drop_algorithm <- nimbleFunction(run = function(d = double(), order = double()) {
    if (d > order)
        stop("Drop algorithm requires order > dim.")
    k <- numeric(d)
    a <- order - d
    k[1] <- a
    nc <- factorial(order - 1)/(factorial(d - 1) * factorial(order - d))
    fs <- matrix(0, nrow = nc, ncol = d)
    fs[1, ] <- k
    q <- 1
    q.seq <- 1

    while (k[d] < a) {
        if (q == d) {
            i <- q
            while (i > 0) {
                i <- i - 1
                q <- i
                if (k[i] != 0)
                  i <- 0
            }
        }
        k[q] <- k[q] - 1
        q <- q + 1
        k[q] <- a - sum(k[1:(q - 1)])
        if (q < d) {
            k[(q + 1):d] <- rep(0, d - q)
        }
        q.seq <- q.seq + 1
        fs[q.seq, ] <- k
    }

    fs <- fs + 1

    returnType(double(2))
    return(fs)
})

#' Central Composite Design (CCD) used for approximate posterior distributions.
#'
#' Generate a d dimension CCD grid via a nimble function list.
#' 
#' @param d Number of dimensions.
#' @param nQuad Ignored.
#'
#' @details
#' This function generates a CCD grid to be used in approximate posteriors. It cannot be compiled without being
#' included within a virtual nimble list "QUAD_RULE_BASE".
#'
#' @author Paul van Dam-Bates
#'
#' @references
#'
#' Rue, H., Martino, S., and Chopin, N. (2009). Approximate Bayesian Inference for Latent Gaussian Models by Using 
#' Integrated Nested Laplace Approximations. Journal of the Royal Statistical Society, Series B 71 (2): 319–92.
#'
#'
#' @export
quadRule_CCD <- nimbleFunction(
    contains = QUAD_RULE_BASE,
    name = "quadRule_CCD",
    setup = function(f0 = 1.1) {
        ## Walsh Index Assignments for Resolution V Fractional Factorials
        index <- c(1, 2, 4, 8, 15, 16, 32, 51, 64, 85, 106, 128, 150, 171, 219, 237,
            247, 256, 279, 297, 455, 512, 537, 557, 594, 643, 803, 863, 998, 1024,
            1051, 1070, 1112, 1169, 1333, 1345, 1620, 1866, 2048, 2076, 2085, 2185,
            2372, 2456, 2618, 2800, 2873, 3127, 3284, 3483, 3557, 3763, 4096, 4125,
            4135, 4174, 4435, 4459, 4469, 4497, 4752, 5255, 5732, 5804, 5915, 6100,
            6369, 6907, 7069, 8192, 8263, 8351, 8422, 8458, 8571, 8750, 8858, 9124,
            9314, 9500, 10026, 10455, 10556, 11778, 11885, 11984, 13548, 14007, 14514,
            14965, 15125, 15554, 16384, 16457, 16517, 16609, 16771, 16853, 17022,
            17453, 17891, 18073, 18562, 18980, 19030, 19932, 20075, 20745, 21544,
            22633, 23200, 24167, 25700, 26360, 26591, 26776, 28443, 28905, 29577,
            32705)

    },
    run = function() {
    },
    methods = list(
        ## Taken from Simon Wood's mgcv package.
        ## https://github.com/cran/mgcv/blob/master/R/inla.r
        ## However, we do scaled design following INLA such that z*zT = 1
        ## from https://github.com/hrue/r-inla/blob/devel/gmrflib/design.c
        ## Can't update nQuad here but makes it general.
        buildGrid = function(levels = integer(0, default = 0), d = integer(0, default = 1)) {
            if ((d > 120 | d < 1)) stop("Dimension of Theta must be in [1,120]")
            
            ## Number of grid points for different dimensions of theta.
            nCCD <- index
            p <- 1
            for (i in seq_along(index)) {
                if (index[i] >= p) p <- p * 2
                nCCD[i] <- p
            }
            nC <- nCCD[d]  ## minimum 2. If 1, choose points c(0,-1,1) but they don't make sense.
            nQ <- nC + 2 * d + 1

            ## First point is mode.,
            design <- matrix(0, nQ, d)

            if (d > 1) {
                for (i in 1:d) {
                    design[index[i] + 2, i] <- 1
                    design[2:(nC + 1), i] <- fwt(x = design[2:(nC + 1), i], n = nC)
                }
                design <- design/sqrt(d)
                ## Next are the star points on the axes. (scaled)
                design[(nC + 2):(nC + d + 1), 1:d] <- diag(d) * 1
                design[(nC + d + 2):(nC + 2 * d + 1), 1:d] <- diag(d) * -1
            } else {
                design <- matrix(c(0, -1, 1), nrow = 3, ncol = 1)
                nQ <- 3
            }

            ## Weights as defined by Rue 2009.  Note that the paper weights are
            ## incorrect:
            ## https://groups.google.com/g/r-inla-discussion-group/c/sy2xYin7YJA
            ## See
            ## https://github.com/hrue/r-inla/blob/devel/gmrflib/approx-inference.c#L1894
            ## w = 1.0 / ((design->nexperiments - 1.0) * (1.0 + exp(-0.5 * SQR(f)) * (SQR(f) / nhyper - 1.0)));
            # f0 <- 1.1
            ## From INLA: z_local[i] = f * design->experiment[k][i] where f = f0*sqrt(d)
            design <- design * sqrt(d) * f0

            ## Weights that actually make sense: Including making the points at
            ## distance f0*sqrt(m) on the sphere: ***This part does not match INLA
            ## but the theory***
            wgts <- 1/((nQ - 1) * f0^2 * (2 * pi)^(-d/2) * exp(-d * f0^2/2))
            wgt0 <- (2 * pi)^(d/2) * (1 - f0^-2)
            ## INLA Weights wgts <- 1 / ((nQ - 1 ) * ( 1 + exp(- (d * f0^2)/2) *
            ## (f0^2 - 1 )) ) wgt0 <- 1 - (nQ-1)*wgts

            ## One time fixes for scalar / vector changes.
            wgt <- numeric(value = 0, length = nQ)
            wgt[1] <- wgt0
            wgt[2:nQ] <- rep(wgts, nQ - 1)

            returnType(double(2))
            output <- matrix(0, nrow = nQ, ncol = d+1)
            output[,1] <- wgt
            output[,2:(d+1)] <- design
            return(output)
        },
        ## fast Walsh transform taken from Wood MGCV inla.
        fwt = function(x = double(1), n = integer()) {
            lag <- 1
            while (lag < n) {
                offset <- lag * 2
                ngroups <- length(x)/offset
                for (group in 0:(ngroups - 1)) {
                    ## vectorized
                    j <- 1:lag + group * offset
                    k <- j + lag
                    xj <- x[j]
                    xk <- x[k]
                    x[j] <- xj + xk
                    x[k] <- xj - xk
                }
                lag <- offset
            }  ## while lag
            returnType(double(1))
            return(x)
        }
    )
)

#' User supplied quadrature grid
#'
#' Generate a d dimension custom grid via a nimble function list.
#' 
#' @param d Number of dimensions.
#' @param nQuad Ignored.
#'
#' @details
#' This function is a placeholder for a user supplied grid. 
#'
#' @author Paul van Dam-Bates
#'
#' @export
# quadRule_USER <- nimbleFunction(
    # contains = QUAD_RULE_BASE,
    # name = "quadRule_USER",
    # setup = function() {
    # },
    # run = function() {},
    # methods = list(
        # buildGrid = function(levels = integer(0, default = 0), d = integer(0, default = 1)) {
            # This will be a place holder for something others may choose to add.
            # Can look for quadRule_Custom and check if it's implemented. If it is
            # will try and use it...
            # output <- matrix(0, nrow = levels,  ncol = d+1)
            # output[,1] <- numeric(value = 0, length = levels)
            # output[,2:(d+1)] <- matrix(0, nrow = levels, d)
            # returnType(double(2))
            # return(output)
        # }
    # )
# )

