library(nimble)
library(testthat)
# Tests of Quadrature Rules and Grids for numerical integration:
source("../../R/quadratureRules.R")
source("../../R/quadratureGrids.R")

test_that("Quadrature Grid Configures Correctly", {
  
  ## 1D Case
  quadGrid <- configureQuadGrid(d=1, levels=3, quadRule = "AGHQ", control = list(quadRules = c("AGHQ", "CCD", "AGHQSPARSE")))
  cquadGrid <- compileNimble(quadGrid)
  cquadGrid$buildGrid()
  nodes <- cquadGrid$nodes()
  wgts <- cquadGrid$weights()
  ans <- sum(dnorm(nodes)*wgts)
  expect_equal(ans, 1, tol = 1e-14)

  ## Check against mvQuad
  nw <- mvQuad::createNIGrid(dim=1, type="GHe", level=3)
  expect_equal(wgts, nw$weights[,1], tol = 1e-14)
  expect_equal(nodes, nw$nodes, tol = 1e-14)

  ## Now change the dimension:
  cquadGrid$setDim(2)
  cquadGrid$buildGrid()
  nodes <- cquadGrid$nodes()
  wgts <- cquadGrid$weights()
  ans <- sum(dnorm(nodes[,1])*dnorm(nodes[,2])*wgts) ## Weights and nodes should sum to 1 on mvnorm.
  expect_equal(ans, 1, tol = 1e-14)

  ## Check against mvQuad
  nw <- mvQuad::createNIGrid(dim=2, type="GHe", level=3)
  expect_equal(wgts, nw$weights[,1], tol = 1e-14)
  expect_equal(nodes, nw$nodes, tol = 1e-14)

  ## Increase number of nodes:
  cquadGrid$buildGrid(nQuad = 11)
  nodes <- cquadGrid$nodes()
  wgts <- cquadGrid$weights()
  nw <- mvQuad::createNIGrid(dim=2, type="GHe", level=11)
  expect_equal(wgts, nw$weights[,1], tol = 1e-12)
  expect_equal(nodes, nw$nodes, tol = 1e-12)

  ## Check against a sparse grid:
  cquadGrid$buildGrid(method = "AGHQSPARSE", nQuad = 3)
  nodes <- cquadGrid$nodes()
  wgts <- cquadGrid$weights()
  nw <- mvQuad::createNIGrid(dim=2, type="GHe", level=3, ndConstruction = "sparse")
  nw$nodes[abs(nw$nodes) < 1e-15] <- 0  ## Make some hard zeros to align, otherwise not in same order.
  ord1 <- do.call(order, data.frame(nodes))
  ord2 <- do.call(order, data.frame(nw$nodes))
  expect_equal(wgts[ord1], nw$weights[ord2,1], tol = 1e-12)
  expect_equal(matrix(nodes[ord1,]), matrix(nw$nodes[ord2,]), tol = 1e-12)

  ## Check CCD:
  ccddesign1 <- t(matrix(c(0, 0, 1.414210000000, 0, -1.414210000000, 0,
                          0, 1.414210000000, 0, -1.414210000000, -1, 1, -1, -1, 1, 1, 1, -1), nrow = 2))
  f0 <- 1.1
  d <- 2
  ccddesign1 <- ccddesign1*f0
  cquadGrid$buildGrid(method = "CCD")
  nodes <- cquadGrid$nodes()
  wgts <- cquadGrid$weights()
  ## Manually cacl weights:
  nQ <- nrow(nodes)
  wgts_ <- c(0, rep(1/((nQ - 1) * f0^2 * (2 * pi)^(-d/2) * exp(-d * f0^2/2)), nQ-1))
  wgts_[1] <- (2 * pi)^(d/2) * (1 - f0^-2)
  expect_equal(wgts, wgts_, tol = 1e-14)
  ord1 <- do.call(order, data.frame(ccddesign1))
  ord2 <- do.call(order, data.frame(nodes))
  expect_equal(ccddesign1[ord1,], nodes[ord2,], tol = 1e-5) ## Not totally accurate as INLA just saves in text.
})

## Should separate CCD test and add a higher dimension example.

test_that("AGHQ Pruning works.", {

  quadGrid2 <- configureQuadGrid(d=3, levels=11, quadRule = "AGHQ", control = list(quadRules = c("AGHQ", "CCD")))
  cquadGrid2 <- compileNimble(quadGrid2)
  cquadGrid2$buildGrid(prune=0)
  nodes <- cquadGrid2$nodes()
  wgts <- cquadGrid2$weights()
  nQ <- cquadGrid2$gridSize()  
  expect_equal(nQ, nrow(nodes))
  expect_equal(nodes[cquadGrid2$modeIndex(),], numeric(3))

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
 
  cquadGrid2$buildGrid(prune = 0)
  nodes.up <- cquadGrid2$nodes()
  wgts.up <- cquadGrid2$weights()
  ## Ensures that the pruning is removed.
  expect_equal(nodes, nodes.up, tol = 1e-16)
  expect_equal(wgts, wgts.up, tol = 1e-16)

  ## Error checks:
  expect_error(cquadGrid2$buildGrid(prune = 0.999), "Will not prune to less than 3 quadrature points. Choose another pruning proportion or switch to Laplace, one quadrature node.")
  expect_error(cquadGrid2$setRule(method = "AGHQSPARSE"), "Quadrature Rule being requested was either not created or is invalid. Choose a valid quadrature rule.")
  expect_error(cquadGrid2$setRule(method = "WHATEVER"), "Quadrature Rule being requested was either not created or is invalid. Choose a valid quadrature rule.")
})

## *** NEED TO ASK CHRIS about how to get this function recognized by the configureQuadGrid function in testthat environment.
# test_that("User provided quadratuture rule.", {
  ## Try to include a user defined quadrature rule:
  RmvQuad <- function(levels, d) {
    out <- mvQuad::createNIGrid(dim=d, type = "GLe", level=levels)
    out <- cbind(out$weights, out$nodes)
  }
  nimMVQuad <- nimbleRcall(function(levels = double(), d = double()){}, Rfun = "RmvQuad", returnType = double(2))
  quadRule_USER <- nimbleFunction(
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

  quadGrid_user <- configureQuadGrid(d=2, levels=3, quadRule = "USER", control = list(quadRules = c("USER", "USERMULTI", "USERSPARSE")))
  cquadGrid_user <- compileNimble(quadGrid_user)
  cquadGrid_user$buildGrid(method = "USER")
  nodes <- cquadGrid_user$nodes()
  wgts <- cquadGrid_user$weights()
  nw <- mvQuad::createNIGrid(dim=2, type="GLe", level=3, ndConstruction = "product")
  ord1 <- do.call(order, data.frame(nodes))
  ord2 <- do.call(order, data.frame(nw$nodes))
  expect_equal(wgts[ord1], nw$weights[ord2,1], tol = 1e-12)
  expect_equal(matrix(nodes[ord1,]), matrix(nw$nodes[ord2,]), tol = 1e-12)

  cquadGrid_user$buildGrid(method = "USERMULTI")
  nodes2 <- cquadGrid_user$nodes()
  wgts2 <- cquadGrid_user$weights()
  expect_equal(wgts, wgts2, tol = 1e-15)
  expect_equal(nodes, nodes2, tol = 1e-15)
  
  cquadGrid_user$buildGrid(method = "USERSPARSE") 
  nodes <- cquadGrid_user$nodes()
  wgts <- cquadGrid_user$weights()
  dup <- duplicated(nodes)
  nodes <- nodes[!dup,]
  wgts <- wgts[!dup]
  nw <- mvQuad::createNIGrid(dim=2, type="GLe", level=3, ndConstruction = "sparse")
  ord1 <- do.call(order, data.frame(nodes))
  ord2 <- do.call(order, data.frame(nw$nodes))
  # expect_equal(wgts[ord1], nw$weights[ord2,1], tol = 1e-12) ## Inefficient combination of repeated values possible. May need to check duplicates... Is that faster? I don't know.
  expect_equal(matrix(nodes[ord1,]), matrix(nw$nodes[ord2,]), tol = 1e-12)
# })

## Need to write a test for the "inner_cache_methods"
