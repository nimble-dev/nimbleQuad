#              CHANGES IN VERSION 1.4.0 (December 2025)

These change notes reflect changes in functionality relative to the 
implementation of Laplace/AGHQ approximation in the `nimble` package, 
version 1.3.0.

## USER LEVEL CHANGES

- Move Laplace/AGHQ approximation from `nimble` to `nimbleQuad` package.

- Add a nested approximation method that provides approximate posterior
  inference using methodology similar to the well-known INLA approach 
  (implemented in the R-INLA package) and to the related methods for extended
  Gaussian latent models (EGLMs) of Stringer et al. (2023), implemented in 
  the `aghq` R package.
  
- Greatly improve efficiency of Laplace/AGHQ approximation through
  major revision of derivative calculations, including removal of "computeMethod"
  options (PR 52).
  
- Greatly improve efficiency of Laplace/AGHQ approximation when using the
  multivariate normal distribution, via a new `dmnormAD` version of the
  distribution customized for efficient computation with NIMBLE's AD 
  system.
  
- Greatly improve efficiency of building Laplace/AGHQ when determining parameter 
  dependencies in `setup_OneAGHQuad` (PR 48).
    
- Greatly improve efficiency of Laplace/AGHQ summary computation of Hessian.
  
- Automatically make use of analytic gradient and Hessian with respect
  to the random variables for univariate and multivariate normal distributions
  rather than AD in Laplace/AGHQ approximation. This also causes the
  "outer" Laplace/AGHQ optimization to not use the AD-based gradient.
  
- Greatly improve the efficiency of building overall Laplace approximation when
  using many individual Laplace approximations (PR #59).
  
- Inform user when overall Laplace approximation is decomposed into multiple
  Laplace approximations.
  
- Improve messaging when Laplace/AGHQ approximation does not converge.

- Provide a variety of quadrature grids, extending from standard AGHQ, 
  including grid pruning, sparse grids, and allowing user-provided grids.
  
- Automatically remove random effect nodes with no dependent data from Laplace
  (Issue #56).
  
- Add ability to find posterior mode (as well as prior-penalized MLE, as used in
  fisheries) via Laplace approximation. 
  
## BUG FIXES

- Fix handling of `fnscale` when `updateSettings` is used.

- Fix `buildAGHQ` for case with no random effects.

- Fix length of latent nodes in Laplace/AGHQ under non-1:1 transformation 
  (PR #41).
  
- Fix seg fault when "model" option passed as `innerOptimStart` (Issue #62).


## DEVELOPER LEVEL CHANGES
 
- Completely refactor and generalize the quadrature grids implementation,
  including making consistent with `mvQuad`.


