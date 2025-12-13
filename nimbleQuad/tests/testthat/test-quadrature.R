# Tests of Quadrature Rules and Grids for numerical integration:

if(!exists(runFailingWindowsTests)) unFailingWindowsTests <- FALSE

test_that("Check Basic Quad Rules work.", {
  ## Basic Virtual List Wrapper:
  rule_wrapper <- nimbleFunction(
    setup = function(){
      quadRule_nfl <- nimbleFunctionList(QUAD_RULE_BASE)
      quadRule_nfl[[1]] <- quadRule_GH(type="GHe")
      quadRule_nfl[[2]] <- quadRule_CCD(f0 = 1.1)
      quadRule_nfl[[3]] <- quadRule_GH(type="GHN")
    },
    run = function(){},
    methods = list(
      ccd = function(d = double()){
         out <- quadRule_nfl[[2]]$buildGrid(levels = 1, d = d)
         returnType(double(2))
         return(out)
      },
      ghe = function(levels = double()){
         out <- quadRule_nfl[[1]]$buildGrid(levels = levels, d = 1)
         returnType(double(2))
         return(out)
      },
      ghn = function(levels = double()){
         out <- quadRule_nfl[[3]]$buildGrid(levels = levels, d = 1)
         returnType(double(2))
         return(out)
      }      
    )
  )
  rules <- rule_wrapper()
  rulesc <- compileNimble(rules)
  ghe <- rulesc$ghe(3)
  ccd <- rulesc$ccd(2)
  ghn <- rulesc$ghn(3)

  ## Should integrate standard normal to 1.
  F1 <- sum(ghe[,1]*dnorm(ghe[,2]))
  expect_equal(F1, 1, 1e-15)  # MacOS gives F1 as 1+mach.epsilon
  ## CCD2 integrates bivariate normal
  F2 <- sum(ccd[,1]*dnorm(ccd[,2])*dnorm(ccd[,3]))
  expect_equal(F2, 1, 1e-16)
  
  ccd21 <- rulesc$ccd(21)
  F21 <- sum(ccd21[,1]*exp(colSums(apply(ccd21[,-1], 1, dnorm, log = TRUE))))
  expect_equal(F21, 1, 1e-14)

  ## pracma check on AGHQ
  gh21 <- rulesc$ghe(21)
  pgh <- pracma::gaussHermite(21)
  expect_equal(pgh$x*sqrt(2), gh21[,2], 1e-14)
  expect_equal(pgh$w*sqrt(2)*exp(pgh$x^2), gh21[,1], 1e-8)
  
  ## Basic Laplace Test:
  lap <- rulesc$ghe(1)
  expect_equal(lap[1,], c(sqrt(2*pi),0), 1e-16)
  lap2 <- rulesc$ghn(1)
  expect_equal(lap2[1,], c(1,0), 1e-16)
  
  ## Test GHN vs mvQuad:
  nw <- mvQuad::createNIGrid(dim=1, type="GHN", level=5)
  ghn5 <- rulesc$ghn(5)
  expect_equal(ghn5[,1], nw$weights[,1], 1e-14)
  expect_equal(ghn5[,2], nw$nodes[,1], 1e-14)

})

test_that("Returns Laplace when nQuad = 1, and one rule is passed.", {
  ## 1D Case
  quadGrid <- configureQuadGrid(d=1, levels=1, quadRule = "AGHQ")
  cquadGrid <- compileNimble(quadGrid)
  cquadGrid$buildGrid()
  expect_equal(cquadGrid$modeIndex(), 1)
  nodes <- cquadGrid$nodes()
  expect_equal(nodes[1,1], 0)
  wgts <- cquadGrid$weights()
  expect_equal(wgts[1], sqrt(2*pi))
  ans <- sum(dnorm(nodes)*wgts)
  expect_equal(ans, 1, tol = 1e-14)
  expect_equal(cquadGrid$gridSize(), length(wgts))

  ## Increase number of nodes:
  cquadGrid$buildGrid(nQuad = 11)
  nodes <- cquadGrid$nodes()
  wgts <- cquadGrid$weights()
  nw <- mvQuad::createNIGrid(dim=1, type="GHe", level=11)
  expect_equal(wgts, nw$weights[,1], tol = 1e-12)
  expect_equal(nodes, nw$nodes, tol = 1e-12)
  expect_equal(cquadGrid$gridSize(), length(wgts))
  expect_equal(length(wgts), nrow(nodes))
})

test_that("Quadrature Grid Configures Correctly", {
  
  ## 1D Case
  quadGrid1 <- configureQuadGrid(d=1, levels=3, quadRule = "AGHQ", control = list(quadRules = c("AGHQ", "CCD", "AGHQSPARSE")))
  cquadGrid1 <- compileNimble(quadGrid1)
  cquadGrid1$buildGrid()
  nodes <- cquadGrid1$nodes()
  wgts <- cquadGrid1$weights()
  ans <- sum(dnorm(nodes)*wgts)
  expect_equal(ans, 1, tol = 1e-14)
  expect_equal(cquadGrid1$gridSize(), length(wgts))
  expect_equal(length(wgts), nrow(nodes))

  ## Check against mvQuad
  nw <- mvQuad::createNIGrid(dim=1, type="GHe", level=3)
  expect_equal(wgts, nw$weights[,1], tol = 1e-14)
  expect_equal(nodes, nw$nodes, tol = 1e-14)

  ## Now change the dimension:
  cquadGrid1$setDim(2)
  cquadGrid1$buildGrid()
  nodes <- cquadGrid1$nodes()
  wgts <- cquadGrid1$weights()
  ans <- sum(dnorm(nodes[,1])*dnorm(nodes[,2])*wgts) ## Weights and nodes should sum to 1 on mvnorm.
  expect_equal(ans, 1, tol = 1e-14)
  expect_equal(cquadGrid1$gridSize(), length(wgts))
  expect_equal(length(wgts), nrow(nodes))

  ## Check against mvQuad
  nw <- mvQuad::createNIGrid(dim=2, type="GHe", level=3)
  expect_equal(wgts, nw$weights[,1], tol = 1e-14)
  expect_equal(nodes, nw$nodes, tol = 1e-14)

  ## Increase number of nodes:
  cquadGrid1$buildGrid(nQuad = 11)
  nodes <- cquadGrid1$nodes()
  wgts <- cquadGrid1$weights()
  nw <- mvQuad::createNIGrid(dim=2, type="GHe", level=11)
  expect_equal(wgts, nw$weights[,1], tol = 1e-12)
  expect_equal(nodes, nw$nodes, tol = 1e-12)
  expect_equal(cquadGrid1$gridSize(), length(wgts))
  expect_equal(length(wgts), nrow(nodes))
  modei <- cquadGrid1$modeIndex()
  expect_equal(nodes[modei,], c(0,0))

  ## Check even number:
  cquadGrid1$buildGrid(nQuad = 6)
  nodes <- cquadGrid1$nodes()
  wgts <- cquadGrid1$weights()
  nw <- mvQuad::createNIGrid(dim=2, type="GHe", level=6)
  expect_equal(wgts, nw$weights[,1], tol = 1e-12)
  expect_equal(nodes, nw$nodes, tol = 1e-12)
  expect_equal(cquadGrid1$gridSize(), length(wgts))  ## Make sure cache is correct:
  expect_equal(length(wgts), nrow(nodes))
  modei <- cquadGrid1$modeIndex()
  expect_equal(modei, -1)

  ## Check against a sparse grid:
  cquadGrid1$buildGrid(method = "AGHQSPARSE", nQuad = 3)
  nodes <- cquadGrid1$nodes()
  wgts <- cquadGrid1$weights()
  nw <- mvQuad::createNIGrid(dim=2, type="GHe", level=3, ndConstruction = "sparse")
  zeros <- which(rowSums(abs(nw$nodes)) < 1e-15)
  wgts_s <- nw$weights[,1]
  nodes_s <- nw$nodes
  wgts_s[zeros[1]] <- sum(nw$weights[zeros])
  wgts_s <- wgts_s[-zeros[-1]]
  nodes_s <- nodes_s[-zeros[-1],]
  nodes_s[zeros[1],] <- numeric(ncol(nodes_s))
  
  ord1 <- do.call(order, data.frame(nodes))
  ord2 <- do.call(order, data.frame(nodes_s))
  expect_equal(wgts[ord1], wgts_s[ord2], tol = 1e-12)
  expect_equal(matrix(nodes[ord1,]), matrix(nodes_s[ord2,]), tol = 1e-12)

  ## Check CCD:
  ccddesign1 <- t(matrix(c(0, 0, 1.414210000000, 0, -1.414210000000, 0,
                          0, 1.414210000000, 0, -1.414210000000, -1, 1, -1, -1, 1, 1, 1, -1), nrow = 2))
  f0 <- 1.1
  d <- 2
  ccddesign1 <- ccddesign1*f0
  cquadGrid1$buildGrid(method = "CCD")
  nodes <- cquadGrid1$nodes()
  wgts <- cquadGrid1$weights()
  ## Manually caclulate weights:
  nQ <- nrow(nodes)
  wgts_ <- c(0, rep(1/((nQ - 1) * f0^2 * (2 * pi)^(-d/2) * exp(-d * f0^2/2)), nQ-1))
  wgts_[1] <- (2 * pi)^(d/2) * (1 - f0^-2)
  expect_equal(wgts, wgts_, tol = 1e-14)
  ord1 <- do.call(order, data.frame(ccddesign1))
  ord2 <- do.call(order, data.frame(nodes))
  expect_equal(ccddesign1[ord1,], nodes[ord2,], tol = 1e-5) ## Not totally accurate as INLA just saves in text.

  expect_equal(cquadGrid1$gridSize(), length(wgts))  ## Make sure cache is correct:
  expect_equal(length(wgts), nrow(nodes))
  modei <- cquadGrid1$modeIndex()
  expect_equal(modei, 1)  ## CCD mode should be first value I think.
})

## Should separate CCD test and add a higher dimension example.

if(Sys.info()['sysname'] != "Windows"|| runFailingWindowsTests) {  # Issue 65
test_that("AGHQ Pruning works.", {

  quadGrid2 <- configureQuadGrid(d=3, levels=11, quadRule = "AGHQ", control = list(quadRules = c("AGHQ", "CCD")))
  cquadGrid2 <- compileNimble(quadGrid2)
  cquadGrid2$buildGrid(prune=0)
  nodes <- cquadGrid2$nodes()
  wgts <- cquadGrid2$weights()
  nQ <- cquadGrid2$gridSize()  
  expect_equal(nQ, nrow(nodes))
  expect_equal(nodes[cquadGrid2$modeIndex(),], numeric(3))

  expect_equal(cquadGrid2$gridSize(), length(wgts))  ## Make sure cache is correct:
  expect_equal(length(wgts), nrow(nodes))
  modei <- cquadGrid2$modeIndex()
  expect_equal(nodes[modei,], c(0,0,0))

  cquadGrid2$buildGrid(prune = 0.2)
  nodes.p <- cquadGrid2$nodes()
  wgts.p <- cquadGrid2$weights()
  expect_equal(nodes.p[cquadGrid2$modeIndex(),], numeric(3))
  nQp <- cquadGrid2$gridSize()
  expect_equal(nQp, nrow(nodes.p))
  ## Should be thinned based on a quantile function / but using bubble sort.
  wgts.adj <- apply(nodes, 1, FUN = function(x)exp(sum(dnorm(x, log = TRUE)))) * wgts
  q <- quantile(wgts.adj, 0.2) + 1e-16
  expect_equal(sum(wgts.adj > q)/nQ, nQp/nQ, tol = 1e-14)
  expect_equal(wgts[wgts.adj > q], wgts.p)

  expect_equal(cquadGrid2$gridSize(), length(wgts.p))  ## Make sure cache is correct:
  expect_equal(length(wgts.p), nrow(nodes.p))
  modei <- cquadGrid2$modeIndex()
  expect_equal(nodes.p[modei,], c(0,0,0))
 
  cquadGrid2$buildGrid(prune = 0)
  nodes.up <- cquadGrid2$nodes()
  wgts.up <- cquadGrid2$weights()
  ## Ensures that the pruning is removed.
  expect_equal(nodes, nodes.up, tol = 1e-16)
  expect_equal(wgts, wgts.up, tol = 1e-16)

  expect_equal(cquadGrid2$gridSize(), length(wgts))  ## Make sure cache is correct:
  expect_equal(length(wgts), nrow(nodes.up))
  modei <- cquadGrid2$modeIndex()
  expect_equal(nodes.up[modei,], c(0,0,0))

  ## Error checks:
  expect_error(cquadGrid2$buildGrid(prune = 0.999), "Will not prune to less than 3 quadrature points. Choose another pruning proportion or switch to Laplace, one quadrature node.")
  expect_error(cquadGrid2$setRule(method = "AGHQSPARSE"), "Quadrature Rule being requested was either not created or is invalid. Choose a valid quadrature rule.")
  expect_error(cquadGrid2$setRule(method = "WHATEVER"), "Quadrature Rule being requested was either not created or is invalid. Choose a valid quadrature rule.")
})
}
## Test a user provided quadrature rule. ***Note that `QUAD_RULE_BASE` needs to be exported which requires a new install of nimbleQuad.
test_that("User provided quadrature rule.", {
  # Try to include a user defined quadrature rule: In this case GLe.
  .GlobalEnv$RmvQuad <- function(levels, d) {
    out <- mvQuad::createNIGrid(dim=d, type = "GLe", level=levels)
    cbind(out$weights, out$nodes)
  }
  .GlobalEnv$nimMVQuad <- nimbleRcall(function(levels = double(), d = double()){}, Rfun = "RmvQuad", returnType = double(2))
  .GlobalEnv$myQuadRule <- nimbleFunction(
      contains = QUAD_RULE_BASE,
      name = "quadRule_USER",
      setup = function() {},
      run = function() {},
      methods = list(
          buildGrid = function(levels = integer(0, default = 0), d = integer(0, default = 1)) {
              output <- nimMVQuad(levels, d)
              returnType(double(2))
              return(output)
          }
      )
  )
  
  ## Case 1: User passes MULTI type grid.
  quadGrid_user <- configureQuadGrid(d=2, levels=3, quadRule = myQuadRule, control = list(quadRules = c("AGHQ", "CCD", "AGHQSPARSE"), userConstruction = "MULTI"))
  cquadGrid_user <- compileNimble(quadGrid_user)
  cquadGrid_user$buildGrid(method = "USER")
  nodes <- cquadGrid_user$nodes()
  wgts <- cquadGrid_user$weights()
  nw <- mvQuad::createNIGrid(dim=2, type="GLe", level=3, ndConstruction = "product")
  ord1 <- do.call(order, data.frame(nodes))
  ord2 <- do.call(order, data.frame(nw$nodes))
  expect_equal(wgts[ord1], nw$weights[ord2,1], tol = 1e-12)
  expect_equal(matrix(nodes[ord1,]), matrix(nw$nodes[ord2,]), tol = 1e-12)

  ## Case 2: User passes a univariate rule.
  quadGrid_user <- configureQuadGrid(d=2, levels=3, quadRule = myQuadRule, control = list(quadRules = c("AGHQ", "CCD", "AGHQSPARSE"), userConstruction = "PRODUCT"))
  cquadGrid_user <- compileNimble(quadGrid_user)
  cquadGrid_user$buildGrid(method = "USER")
  nodes <- cquadGrid_user$nodes()
  wgts <- cquadGrid_user$weights()
  nw <- mvQuad::createNIGrid(dim=2, type="GLe", level=3, ndConstruction = "product")
  ord1 <- do.call(order, data.frame(nodes))
  ord2 <- do.call(order, data.frame(nw$nodes))
  expect_equal(wgts[ord1], nw$weights[ord2,1], tol = 1e-12)
  expect_equal(matrix(nodes[ord1,]), matrix(nw$nodes[ord2,]), tol = 1e-12)

  ## Case 3: User passes a univariate rule for SPARSE construction.
  quadGrid_user <- configureQuadGrid(d=2, levels=3, quadRule = myQuadRule, control = list(quadRules = c("AGHQ", "CCD", "AGHQSPARSE"), userConstruction = "SPARSE"))
  cquadGrid_user <- compileNimble(quadGrid_user)
  cquadGrid_user$buildGrid(method = "USER")
  nodes <- cquadGrid_user$nodes()
  wgts <- cquadGrid_user$weights()
  nw <- mvQuad::createNIGrid(dim=2, type="GLe", level=3, ndConstruction = "sparse")
  ord1 <- do.call(order, data.frame(nodes))
  ord2 <- do.call(order, data.frame(nw$nodes))
  expect_equal(wgts[ord1], nw$weights[ord2,1], tol = 1e-12)
  expect_equal(matrix(nodes[ord1,]), matrix(nw$nodes[ord2,]), tol = 1e-12)

  ## Make sure we can still swap back accurately.
  cquadGrid_user$buildGrid(method = "AGHQ")
  nodes2 <- cquadGrid_user$nodes()
  wgts2 <- cquadGrid_user$weights()
  nw <- mvQuad::createNIGrid(dim=2, type="GHe", level=3, ndConstruction = "product")  
  expect_equal(nw$weights[,1], wgts2, tol = 1e-14)
  expect_equal(nw$nodes, nodes2, tol = 1e-15)
  
  ## Check in on sparse.
  cquadGrid_user$buildGrid(method = "AGHQSPARSE") 
  nodes3 <- cquadGrid_user$nodes()
  wgts3 <- cquadGrid_user$weights()
  ## Ensure this still integrates gaussian perfectly:
  expect_equal( sum(wgts3*dnorm(nodes3[,1])*dnorm(nodes3[,2])), 1, 1e-16)

  ## Compare with mvQuad: Need to shift their modes manually.
  nw <- mvQuad::createNIGrid(dim=2, type="GHe", level=3, ndConstruction = "sparse")
  zeros <- which(rowSums(abs(nw$nodes)) < 1e-15)
  wgts4 <- nw$weights[,1]
  nodes4 <- nw$nodes
  wgts4[zeros[1]] <- sum(nw$weights[zeros])
  wgts4 <- wgts4[-zeros[-1]]
  nodes4 <- nodes4[-zeros[-1],]
  nodes4[zeros[1],] <- numeric(ncol(nodes4))
  ord <- do.call(order, data.frame(nodes3))
  ord2 <- do.call(order, data.frame(nodes4))
  expect_equal(wgts3[ord], wgts4[ord2], tol = 1e-14)
  expect_equal(as.numeric(nodes3[ord,]), as.numeric(nodes4[ord2,]), tol = 1e-14)
})


## Write a buildGridTest

## Need to write a test for the "inner_cache_methods"

## Make some plots:
# quadGrid <- configureQuadGrid(d=2, levels=3, quadRule = "AGHQ", control = list(quadRules = c("AGHQ", "AGHQSPARSE", "CCD")))
# quadGrid$buildGrid()
# nodes_aghq <- quadGrid$nodes()
# quadGrid$buildGrid(nQuad = 5)
# nodes_aghq2 <- quadGrid$nodes()

# quadGrid$buildGrid(method = "CCD")
# nodes_ccd <- quadGrid$nodes()
# quadGrid$buildGrid(method = "AGHQSPARSE", nQuad = 2)
# nodes_sparse <- quadGrid$nodes()
# quadGrid$buildGrid(method = "AGHQSPARSE", nQuad = 3)
# nodes_sparse2 <- quadGrid$nodes()

# par(mfrow = c(1,3))
# plot(nodes_aghq, xlab = expression(z[1]), ylab = expression(z[2]), pch = 16, col = 'black', main = "AGHQ", xlim = c(-3,3), ylim = c(-3,3))
# points(nodes_aghq2, pch = 3, col = 'red')

# plot(nodes_sparse, xlab = expression(z[1]), ylab = expression(z[2]), pch = 16, col = 'black', main = "AGHQ Sparse", xlim = c(-3,3), ylim = c(-3,3))
# points(nodes_sparse2, pch = 3, col = 'red')

# plot(nodes_ccd, xlab = expression(z[1]), ylab = expression(z[2]), pch = 16, col = 'black', main = "CCD", xlim = c(-3,3), ylim = c(-3,3))

# C <- rbind(c(1.5, 0.75), c(0.75, 1.75))
# eC <- eigen(C)
# x     <- seq(-5, 5, 0.25)
# y     <- seq(-5, 5, 0.25)
# cholC <- chol(C)
# f     <- function(x, y) apply(cbind(x,y), 1, FUN = function(z) dmnorm_chol(z, c(0,0), cholesky = cholC, prec_param = FALSE) )
# fz     <- outer(x, y, f)

# z <- nodes
# zChol <- t(cholC %*% t(z))
# zEigen <- t(eC$vectors %*% diag(sqrt(eC$values)) %*% t(z))

# ycond <- C[1,2]/C[1,1]*x
# xcond <- C[2,1]/C[2,2]*y

# contour(x, y, fz, main = "", xlab = expression(theta[1]), ylab = expression(theta[2]))
# points(zChol, pch = 16, col = 'blue')
# points(zEigen, pch = 16, col = 'red')
# legend('bottomright', legend = c('Spectral Transformation', 'Cholesky Transformation'), 
	# col = c('red', 'blue'), pch = 16) 

# par(mfrow=c(1,2))
# contour(x, y, fz, main = "Cholesky Transformation", xlab = expression(theta[1]), ylab = expression(theta[2]))
# points(zChol, pch = 16, col = 'red')
# contour(x, y, fz, main = "Spectral Transformation", xlab = expression(theta[1]), ylab = expression(theta[2]))
# points(zEigen, pch = 16, col = 'red')


# z <- nodes_ccd
# zChol <- t(cholC %*% t(z))
# contour(x, y, fz, xlab = expression(theta[1]), ylab = expression(theta[2]))
# points(zChol, pch = 16, col = 'red')
