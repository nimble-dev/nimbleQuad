# nimbleQuad
[![Build Status](https://github.com/nimble-dev/nimbleQuad/actions/workflows/tests.yaml/badge.svg?branch=devel)](https://github.com/nimble-dev/nimbleQuad/actions/workflows/tests.yaml)
[![CRAN](http://www.r-pkg.org/badges/version/nimbleQuad)](https://CRAN.R-project.org/package=nimbleQuad)

Laplace approximation, quadrature and nested approximation methods for NIMBLE.

Please see the [NIMBLE](https://r-nimble.org/) website and [main NIMBLE repository](https://github.com/nimble-dev/nimble) for more information.

### Package Requirements

`nimbleQuad` must be used with version `1.4.0` or higher of the `nimble` package, because it makes use of the automatic differentiation (AD) feature of `nimble`released in `nimble` version 1.0.0 and extended in version 1.4.0.

For using the methods of `nimbleQuad` on a model, derivative calculations need to be built into for the model object.  This is accomplished using the `buildDerivs = TRUE` argument in the call to `nimbleModel` as:
```
nimbleModel(code, constants, data, inits, buildDerivs = TRUE)
```






