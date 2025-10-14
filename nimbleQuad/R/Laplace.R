## NIMBLE Laplace approximation
## AGHQuad/Laplace base class
AGHQuad_BASE <- nimbleFunctionVirtual(
  run = function() {},
  methods = list(
    reset = function(gr_RE = logical(0, default = TRUE),
                     he_RE = logical(0, default = TRUE),
                     gr_P_RE = logical(0, default = TRUE),
                     gr_P_RE_wrt_RE = logical(0, default = TRUE),
                     he_P_RE_wrt_RE2_uptri = logical(0, default = TRUE)){
    },
    calcLogLik2 = function(p = double(1)){
      returnType(double())
    },
    gr_logLik2 = function(p = double(1)){
      returnType(double(1))
    },
    negHess = function(p = double(1), reTransform = double(1), forceReset = logical(0, default = FALSE)){
      returnType(double(2))
    },
    update_max_logLik_RE = function(p = double(1)){
      returnType(double(1))
    },
    he_P_RE_wrt_RE_wrt_P_b = function(p = double(1), reTransform = double(1),
                                 forceUpdate = logical(0, default = FALSE), forceReset = logical(0, default = FALSE)){
      returnType(double(2))
    },
    jac_gr_P_RE_wrt_RE_inDir= function(p = double(1), reTransform = double(1),
                                         inDir = double(1),
                                         forceUpdate = logical(0, default = FALSE),
                                         forceReset = logical(0, default = FALSE)) {
      returnType(double(2))
    },
    reset_outer_logLik = function(){},
    get_reTransLength = function(){returnType(double(0))},
    save_outer_logLik = function(logLikVal = double()){},
    get_param_value = function(atOuterMode = integer(0, default = 0)){
      returnType(double(1))
    },
    get_inner_mode = function(atOuterMode = integer(0, default = 0)){
      returnType(double(1))
    },
    get_inner_negHessian = function(atOuterMode = integer(0, default = 0)){
      returnType(double(2))
    },
    get_inner_negHessian_chol = function(atOuterMode = integer(0, default = 0)){
      returnType(double(2))
    },
    check_convergence = function(){
      returnType(double())
    },
    updateSettings = function(optimMethod = character(0, default="NULL"),
                              optimStart = character(0, default="NULL"),
                              optimStartValues = double(1, default=Inf),
                              optimWarning = integer(0, default = -1),
                              useInnerCache = integer(0, default=-1),
                              nQuad = integer(0, default=-1),
                              quadTransform = character(0, default="NULL"),
                              optimControl = optimControlNimbleList(default=nimOptimDefaultControl()),
                              replace_optimControl = logical(0, default=FALSE)) {
    },
    ## set_nQuad = function(nQUpdate = integer()){},
    ## set_transformation = function(transformation = character()){},
    ## set_warning = function(warn = logical()){},
    ## set_reInitMethod = function(method = character(), value=double(1)){},
    set_randomeffect_values = function(p = double(1)){}
    ## set_inner_cache = function(cache = logical(0, default = TRUE)){}
  )
)

setup_OneAGHQuad <- function(model, paramNodes, randomEffectsNodes, calcNodes, paramDeps, control) {
  # common setup steps for 1D and >1D cases
  optimControl_ <- extractControlElement(control, 'optimControl', nimOptimDefaultControl())
  optimMethod_ <- extractControlElement(control, 'optimMethod', 'nlminb')
  optimStart_ <- extractControlElement(control, 'optimStart', 'last.best')
  optimStartValues_ <- extractControlElement(control, 'optimStartValues', 0)
  nre  <- length(model$expandNodeNames(randomEffectsNodes, returnScalarComponents = TRUE))

  if(missing(paramDeps))
    paramDeps <- model$getDependencies(paramNodes, determOnly = TRUE, self=FALSE)
  if(length(paramDeps)) {
    calcNodesParents <- model$getParents(calcNodes, determOnly = TRUE)
    paramDeps <- paramDeps[!paramDeps %in% calcNodes & paramDeps %in% calcNodesParents]
  }
  innerCalcNodes <- calcNodes
  if(length(paramDeps))
    calcNodes <- model$expandNodeNames(c(paramDeps, calcNodes), sort = TRUE)
  wrtNodes <- c(paramNodes, randomEffectsNodes)
  reTrans <- parameterTransform(model, randomEffectsNodes)
  npar <- length(model$expandNodeNames(paramNodes, returnScalarComponents = TRUE))

  if(npar > 1) p_indices <- as.numeric(1:npar)
  else p_indices <- as.numeric(c(1, -1))

  quadRule <- extractControlElement(control, 'innerQuadRule', 'AGHQ') ##***CJP latentQuadRule?

  list(optimControl_=optimControl_,
       optimMethod_=optimMethod_,
       optimStart_=optimStart_,
       optimStartValues_=optimStartValues_,
       nre = nre,
       paramDeps = paramDeps,
       innerCalcNodes = innerCalcNodes,
       calcNodes = calcNodes,
       wrtNodes = wrtNodes,
       reTrans = reTrans,
       npar = npar,
       p_indices = p_indices,
       quadRule = quadRule
       )
}

## A single Laplace approximation for only one scalar random effect node
buildOneLaplace1D <- function(model, paramNodes, randomEffectsNodes, calcNodes, paramDeps, control = list()) {
  buildOneAGHQuad1D(model, nQuad = 1, paramNodes, randomEffectsNodes, calcNodes, paramDeps, control)
}

buildOneAGHQuad1D <- nimbleFunction(
  contains = AGHQuad_BASE,
  setup = function(model, nQuad, paramNodes, randomEffectsNodes, calcNodes, paramDeps, control = list()) {
    ## Check the number of random effects is 1
    ## optimControl_ <- extractControlElement(control, 'optimControl', nimOptimDefaultControl())
    ## optimMethod_ <- extractControlElement(control, 'optimMethod', 'BFGS')
    ## optimStart_ <- extractControlElement(control, 'optimStart', 'constant')
    ## optimStartValues_ <- extractControlElement(control, 'optimStartValues', 0)
    nQuad_ <- nQuad
    S <- setup_OneAGHQuad(model, paramNodes, randomEffectsNodes, calcNodes, paramDeps, control)
    optimControl_ <- S$optimControl_
    optimMethod_ <- S$optimMethod_
    optimStart_ <- S$optimStart_
    optimStartValues_ <- S$optimStartValues_
    nre  <-  S$nre
    paramDeps  <-  S$paramDeps
    innerCalcNodes  <-  S$innerCalcNodes
    calcNodes  <-  S$calcNodes
    wrtNodes  <-  S$wrtNodes
    reTrans  <-  S$reTrans
    npar  <-  S$npar
    p_indices  <-  S$p_indices
    quadRule_ <- S$quadRule

    nreTrans <- 1 # must be the case
    if(length(reTrans) != 1) stop("buildOneAGHQuad1D: The length of transformed random effects must be 1.")

    ## nre  <- length(model$expandNodeNames(randomEffectsNodes, returnScalarComponents = TRUE))
    if(length(nre) != 1) stop("buildOneAGHQuad1D: Number of random effects for buildOneAGHQuad1D or buildOneLaplace1D must be 1")
    reTrans_indices <- as.numeric(c(npar+1, -1))
    reTrans_indices_inner <- as.numeric(c(1, -1))
    p_reTrans_indices <- as.numeric(1:(npar + 1))

    ## Set up start values for the inner optimization of Laplace approximation
    if(!is.character(optimStart_) | length(optimStart_) != 1) stop("buildOneAGHQuad1D: There is a problem with `optimStart`: ", optimStart_)
    startID <- switch(optimStart_, last=1, last.best=2, constant=3, random=4, model=5)
    if(startID == 5) {
      constant_init_reTrans <- c(values(model, randomEffectsNodes), -1)
      startID <- 3  
    } else constant_init_reTrans <- c(optimStartValues_, -1)
    ## Update and constant nodes for obtaining derivatives using AD
    inner_derivsInfo    <- makeModelDerivsInfo(model = model, wrtNodes = randomEffectsNodes, calcNodes = innerCalcNodes)
    inner_updateNodes   <- inner_derivsInfo$updateNodes
    inner_constantNodes <- inner_derivsInfo$constantNodes
    joint_derivsInfo    <- makeModelDerivsInfo(model = model, wrtNodes = wrtNodes, calcNodes = calcNodes)
    joint_updateNodes   <- joint_derivsInfo$updateNodes
    joint_constantNodes <- joint_derivsInfo$constantNodes

    ## The following is used to ensure the one_time_fixes are run when needed.
    one_time_fixes_done <- FALSE

    ## Flags used for managing update and reset of groups of related derivs calls.
    ## Update means the "*_updateNodes" will be updated in the relevant tapes.
    ## Reset means the "*_constantNodes" will be updated, WHICH REQUIRES RE-TAPING and is thus costly.
    ## See multivariate version below for description of these flags
    ## and also the NOMENCLATURE such as "gr_RE".
    ##
    ## Flags for all gradients as a function of random effects only, i.e. inner gradients.
    gr_RE_update_once <- TRUE
    gr_RE_update_always <- FALSE
    gr_RE_reset_once <- TRUE

    ## Flags for all Hessians as a function of random effects only, i.e. inner Hessians.
    he_RE_update_once <- TRUE
    he_RE_update_always <- FALSE
    he_RE_reset_once <- TRUE

    ## Flags for all gradients as a function of parameters and random effects.
    gr_P_RE_update_once <- TRUE
    gr_P_RE_update_always <- FALSE
    gr_P_RE_reset_once <- TRUE

    ## Flags for all gradients as a function of parameters and random effects, but only wrt random effects.
    gr_P_RE_wrt_RE_update_once <- TRUE
    gr_P_RE_wrt_RE_update_always <- FALSE
    gr_P_RE_wrt_RE_reset_once <- TRUE

    ## Flags for all Hessians as a function of parameters and random effects wrt random effects, flattened upper triangular.
    he_P_RE_wrt_RE2_uptri_update_once <- TRUE
    he_P_RE_wrt_RE2_uptri_update_always <- FALSE
    he_P_RE_wrt_RE2_uptri_reset_once <- TRUE

    ## Caches for results of inner optimization:
    cache_inner_max <- TRUE
    saved_inner_argmax <- constant_init_reTrans
    saved_inner_max_value <- -Inf
    saved_inner_max_p <- if(npar > 1) rep(Inf, npar) else as.numeric(c(Inf, -1))
    saved_inner_negHess <- matrix(0, nrow = 1, ncol = 1)
    saved_inner_logdetNegHess <- 0

    ## Cache for set_P
    current_P_for_inner <- saved_inner_max_p

    ## Cache to ensure taping is done from init (RE)
    reInitTrans_for_taping <- constant_init_reTrans

    ## Caches to help with outer optimization:
    ## Record the maximum Laplace loglikelihood value for obtaining inner optimization start values
    max_margLogLik<- -Inf
    max_margLogLik_inner_argmax <- constant_init_reTrans
    margLogLik_saved_value <- -Inf

    ## Values to save when max inner log lik reached.
    max_outer_logLik <- -Inf
    outer_mode_inner_negHess <- matrix(0, nrow = 1, ncol = 1)
    outer_mode_inner_argmax <- as.numeric(c(1, -1))
    outer_param_max <- if(npar > 1) rep(Inf, npar) else as.numeric(c(Inf, -1))

    ## Cached gradients for AGHQ.
    gr_sigmahatwrtre <- numeric(1)
    gr_sigmahatwrtp <- if(npar > 1) numeric(npar) else as.numeric(c(1, -1))
    gr_rehatwrtp <- if(npar > 1) numeric(npar) else as.numeric(c(1, -1)) # double(1)
#    gr_QuadSum_value <- if(npar > 1) numeric(npar) else as.numeric(c(1, -1)) # not used
    AGHQuad_saved_gr <- if(npar > 1) numeric(npar) else as.numeric(c(1, -1))
    quadrature_previous_p <- if(npar > 1) rep(Inf, npar) else as.numeric(c(Inf, -1))

    ## Build AGHQ grid for 1D:
    ## This is set up to add other quad grids in the future. quadRule := "AGHQ" to start.
    quadGrid <- configureQuadGrid(d = 1, levels = nQuad_, quadRule = quadRule_)
    nodes <-  matrix(0, nrow = nQuad_, ncol = 1)
    wgts <- numeric(nQuad_)
    logDensity_quad <- numeric(nQuad_)
    if(nQuad_ == 1) {
      wgts <- c(0,-1)
      logDensity_quad <- c(0,-1)
    }
    quadTransform_ <- extractControlElement(control, "quadTransform", "cholesky")

    ## Convergence check for outer function.
    converged <- 0

    warn_optim <- extractControlElement(control, 'optimWarning', FALSE) ## Warn about inner optimization issues
  },
  run = function(){},
  methods = list(
    fix_one_vec = function(x = double(1)) {
      if(length(x) == 2) {
        if(x[2] == -1) {
          ans <- numeric(length = 1, value = x[1])
          return(ans)
        }
      }
      return(x)
      returnType(double(1))
    },
    one_time_fixes = function() {
      ## Run this once after compiling; remove extraneous -1 if necessary
      if(one_time_fixes_done) return()
      reTrans_indices <<- fix_one_vec(reTrans_indices)
      reTrans_indices_inner <<- fix_one_vec(reTrans_indices_inner)
      saved_inner_argmax <<- fix_one_vec(saved_inner_argmax)
      outer_mode_inner_argmax <<-  fix_one_vec(outer_mode_inner_argmax)
      max_margLogLik_inner_argmax <<- fix_one_vec(max_margLogLik_inner_argmax)
      constant_init_reTrans <<- fix_one_vec(constant_init_reTrans)
      reInitTrans_for_taping <<- fix_one_vec(reInitTrans_for_taping)
      #      if(startID == 3) optStart <<- fix_one_vec(optStart)
      if(npar == 1) {
        p_indices <<- fix_one_vec(p_indices)
        saved_inner_max_p <<- fix_one_vec(saved_inner_max_p)
        current_P_for_inner <<- fix_one_vec(current_P_for_inner)
        outer_param_max <<- fix_one_vec(outer_param_max)
        gr_sigmahatwrtp <<- fix_one_vec(gr_sigmahatwrtp)
        gr_rehatwrtp <<- fix_one_vec(gr_rehatwrtp)
#        gr_QuadSum_value <<- fix_one_vec(gr_QuadSum_value) # not used
        AGHQuad_saved_gr <<- fix_one_vec(AGHQuad_saved_gr)
        quadrature_previous_p <<- fix_one_vec(quadrature_previous_p)
      }
      if(nQuad_ == 1) {
        wgts <<- fix_one_vec(wgts)
        logDensity_quad <<- fix_one_vec(logDensity_quad)
      }
      reInit <- values(model, randomEffectsNodes)
      set_reInit(reInit)
      one_time_fixes_done <<- TRUE
    },
    updateSettings = function(optimMethod = character(0, default="NULL"),
                              optimStart = character(0, default="NULL"),
                              optimStartValues = double(1, default=Inf),
                              optimWarning = integer(0, default = -1),
                              useInnerCache = integer(0, default=-1),
                              nQuad = integer(0, default=-1),
                              quadTransform = character(0, default="NULL"),
                              optimControl = optimControlNimbleList(default=nimOptimDefaultControl()),
                              replace_optimControl = logical(0, default=FALSE)) {
      # Checking should have been done already. Or, if this is being called directly,
      # it will be for development or advanced uses and we can skip checking.
      if(optimMethod != "NULL") optimMethod_ <<- optimMethod
      if(optimStart != "NULL") {
        if(optimStart == "last") startID <<- 1 # last
        else if(optimStart == "last.best") startID <<- 2 # last.best
        else if(optimStart == "constant") startID <<- 3 # use fixed vector optimStart provided at setup time
        else if(optimStart == "random") startID <<- 4
        else if(optimStart == "model") {
          startID <<- 3
          constant_init_reTrans <<- reTrans$transform(values(model, randomEffectsNodes))
        }
      }
      if((length(optimStartValues) != 1) | (optimStartValues[1] != Inf) ) {
        if((length(optimStartValues) == 1) & (optimStartValues[1] == -Inf) ) { # numeric code for "model" setting
          constant_init_reTrans <<- reTrans$transform(values(model, randomEffectsNodes))
        } else {
          if(startID <= 3) {
            constant_init_reTrans <<- optimStartValues
            if(length(constant_init_reTrans) == 1)
              if(nre > 1)
                constant_init_reTrans <<- rep(constant_init_reTrans, nre)
          }
        }
      }
      if((!one_time_fixes_done) & (length(constant_init_reTrans) == 1)) {
         constant_init_reTrans <<- c(constant_init_reTrans, -1)
      }
      if(optimWarning != -1) {
        warn_optim <<- optimWarning != 0
      }
      if(useInnerCache != -1) {
        cache_inner_max <<- useInnerCache != 0
      }
      if(nQuad != -1) {
        quadGrid$buildGrid(method = quadRule_, nQuad = nQuad)
        nQuad_ <<- nQuad
      }
      if(quadTransform != "NULL") {
        quadTransform_ <<- quadTransform
      }
      if(replace_optimControl) {
        if(optimControl$fnscale == 1) optimControl$fnscale <- -1
        optimControl_ <<- optimControl
      }
    },
    set_reInit = function(re = double(1)) {
      reInitTrans <- reTrans$transform(re)
      saved_inner_argmax <<- reInitTrans
    },
    get_reInitTrans = function() {
      if(startID == 1) ans <- saved_inner_argmax                        ## last
      else if(startID == 2) ans <- max_margLogLik_inner_argmax          ## last.best
      else if(startID == 3) ans <- constant_init_reTrans                ## constant
      else if(startID == 4){                                            ## random (prior).
        model$simulate(randomEffectsNodes)
        ans <- reTrans$transform(values(model, randomEffectsNodes))     ## From prior:
      }
      return(ans)
      returnType(double(1))
    },
    get_reTransLength = function() {
      returnType(double(0))
      return(nreTrans) # must be 1 in this version
    },
    ## See comments in multivariate version below for NOMENCLATURE of method names.
    ##
    ## Joint log-likelihood with values of parameters fixed: used only for inner optimization
    set_P = function(p = double(1)) {
      values(model, paramNodes) <<- p
      model$calculate(paramDeps)
      gr_RE_update_once <<- TRUE
      he_RE_update_once <<- TRUE
      current_P_for_inner <<- p
    },
    reset = function(gr_RE = logical(0, default = TRUE),
                     he_RE = logical(0, default = TRUE),
                     gr_P_RE = logical(0, default = TRUE),
                     gr_P_RE_wrt_RE = logical(0, default = TRUE),
                     he_P_RE_wrt_RE2_uptri = logical(0, default = TRUE)){
      gr_RE_reset_once <<- gr_RE
      he_RE_reset_once <<- he_RE
      gr_P_RE_reset_once <<- gr_P_RE
      gr_P_RE_wrt_RE_reset_once <<- gr_P_RE_wrt_RE
      he_P_RE_wrt_RE2_uptri_reset_once <<- he_P_RE_wrt_RE2_uptri
      ## Reset the inner optimization cache.
    },
    ## The first set of functions are for "inner" steps, where
    ## P has already been set, so these are functions only of reTransform.
    logLik_RE = function(reTransform = double(1)) {
      re <- reTrans$inverseTransform(reTransform)
      values(model, randomEffectsNodes) <<- re
      ans <- model$calculate(innerCalcNodes) + reTrans$logDetJacobian(reTransform)
      return(ans)
      returnType(double())
    },
    # Gradient of the joint log-likelihood (p fixed) w.r.t. transformed random effects: used only for inner optimization
    gr_RE_a = function(reTransform = double(1),
                       forceUpdate = logical(0, default = FALSE),
                       forceReset = logical(0, default = FALSE)) {
      # renamed: previously had "internal" suffix
      #  previously gr_inner_logLik_internal
      do_reset <- forceReset | gr_RE_reset_once
      ans <- derivs(logLik_RE(reTransform), wrt = reTrans_indices_inner, order = 1, model = model,
                    updateNodes = inner_updateNodes, constantNodes = inner_constantNodes,
                    do_update = gr_RE_update_once | gr_RE_update_always | forceUpdate | do_reset,
                    reset=do_reset)
      gr_RE_update_once <<- FALSE
      gr_RE_reset_once <<- FALSE
      return(ans$jacobian[1,])
      returnType(double(1))
    },
    ## Double taping for efficiency
    gr_RE_b = function(reTransform = double(1),
                       forceUpdate = logical(0, default = FALSE),
                       forceReset = logical(0, default = FALSE)) {
      ## renamed: previously had no suffix
      ## previusly gr_inner_logLik
      do_reset <- forceReset | gr_RE_reset_once
      do_update <- gr_RE_update_once | gr_RE_update_always | forceUpdate | do_reset
      ans <- derivs(gr_RE_a(reTransform, forceUpdate=do_update, forceReset=do_reset),
                    wrt = reTrans_indices_inner, order = 0, model = model,
                    updateNodes = inner_updateNodes, constantNodes = inner_constantNodes,
                    do_update = do_update,
                    reset=do_reset)
      gr_RE_update_once <<- FALSE
      gr_RE_reset_once <<- FALSE
      return(ans$value)
      returnType(double(1))
    },
    gr_for_optim = function(reTransform = double(1)) {
      ## If the tape will be reset, we ensure we record it at the init params.
      ## I am not sure why except this came from experience.
      if(gr_RE_reset_once) {
        gr_RE_b(reInitTrans_for_taping)
      }
      return(gr_RE_b(reTransform))
      returnType(double(1))
    },
    # Hessian of the joint log-likelihood (p fixed) w.r.t. transformed random effects: used only for inner optimization
    he_RE_b = function(reTransform = double(1),
                       forceUpdate = logical(0, default = FALSE),
                       forceReset = logical(0, default = FALSE)) {
      ## renamed: previously had "internal" suffix
      ## previously he_inner_logLik_internal
      ## reimplemented: now uses order(1) from gr_inner_logLik
      do_reset <- forceReset | he_RE_reset_once
      do_update <- he_RE_update_once | he_RE_update_always | forceUpdate | do_reset
      ans <- derivs(gr_RE_a(reTransform, forceUpdate=do_update, forceReset=do_reset),
                    wrt = reTrans_indices_inner, order = 1, model = model,
                    updateNodes = inner_updateNodes, constantNodes = inner_constantNodes,
                    do_update = do_update,
                    reset=do_reset)
      he_RE_update_once <<- FALSE
      he_RE_reset_once <<- FALSE
      res <- ans$jacobian
      return(res)
      returnType(double(2))
    },
    he_RE_b_asvec = function(reTransform = double(1),
                             forceUpdate = logical(0, default = FALSE),
                             forceReset = logical(0, default = FALSE)) {
      ans <- he_RE_b(reTransform, forceUpdate=forceUpdate, forceReset=forceReset)
      res <- nimNumeric(value = ans, length = length(reTransform)*length(reTransform))
      return(res)
      returnType(double(1))
    },
    he_RE_c = function(reTransform = double(1),
                       forceUpdate = logical(0, default = FALSE),
                       forceReset = logical(0, default = FALSE)) {
      ## renamed: previously had no suffix
      # previously he_inner_logLik
      do_reset <- forceReset | he_RE_reset_once
      do_update <- he_RE_update_once | he_RE_update_always | forceUpdate | do_reset
      ans <- derivs(he_RE_b_asvec(reTransform, forceUpdate=do_update, forceReset=do_reset),
                    wrt = reTrans_indices_inner,
                    order = 0, model = model,
                    updateNodes = inner_updateNodes, constantNodes = inner_constantNodes,
                    do_update = do_update,
                    reset=do_reset)
      he_RE_update_once <<- FALSE
      he_RE_reset_once <<- FALSE
      res <- matrix(value = ans$value, nrow = nreTrans, ncol = nreTrans)
      return(res)
      returnType(double(2))
    },
    he_for_optim = function(reTransform = double(1)) {
      ## If the tape will be reset, we ensure we record it at the init params.
      ## I am not sure why except this came from experience.
      if(he_RE_reset_once) {
        he_RE_c(reInitTrans_for_taping)
      }
      return(he_RE_c(reTransform))
      returnType(double(2))
    },
    negHess = function(p = double(1),
                       reTransform = double(1),
                       forceReset = logical(0, default = FALSE)) {
      set_P(p) # This sets the update flag to TRUE.
      ans <- -he_RE_c(reTransform, forceUpdate=TRUE, forceReset=forceReset)
      return(ans)
      returnType(double(2))
    },
    ## The next set of functions are for "outer" steps,
    ## which are functions of both p and reTransform.
    ## Joint log-likelihood in terms of parameters and transformed random effects
    logLik_P_RE = function(p = double(1), reTransform = double(1)) {
        re <- reTrans$inverseTransform(reTransform)
        values(model, paramNodes) <<- p
        values(model, randomEffectsNodes) <<- re
        ans <- model$calculate(calcNodes) +  reTrans$logDetJacobian(reTransform)
        return(ans)
        returnType(double())
    },
    gr_P_RE_a = function(p = double(1), reTransform = double(1),
                         forceUpdate = logical(0, default = FALSE),
                         forceReset = logical(0, default = FALSE)) {
        # previously gr_joint_logLik_wrt_p_re_internal (?)
        do_reset <- forceReset | gr_P_RE_reset_once
        do_update <- gr_P_RE_update_once | gr_P_RE_update_always | forceUpdate | do_reset
        ans <- derivs(logLik_P_RE(p, reTransform), wrt = p_reTrans_indices, order = 1, model = model,
                      updateNodes = joint_updateNodes, constantNodes = joint_constantNodes,
                      do_update = do_update,
                      reset=do_reset)
        gr_P_RE_update_once <<- FALSE
        gr_P_RE_reset_once <<- FALSE
        return(ans$jacobian[1,])
        returnType(double(1))
    },
    gr_P_RE_b = function(p = double(1), reTransform = double(1),
                         forceUpdate = logical(0, default = FALSE),
                         forceReset = logical(0, default = FALSE)) {
      ## previously gr_joint_logLik_wrt_p_re
        do_reset <- forceReset | gr_P_RE_reset_once
        do_update <- gr_P_RE_update_once | gr_P_RE_update_always | forceUpdate | do_reset
        ans <- derivs(gr_P_RE_a(p, reTransform, forceUpdate=do_update, forceReset=do_reset), wrt = p_reTrans_indices, order = 0, model = model,
                      updateNodes = joint_updateNodes, constantNodes = joint_constantNodes,
                      do_update = do_update,
                      reset=do_reset)
        gr_P_RE_update_once <<- FALSE
        gr_P_RE_reset_once <<- FALSE
        return(ans$value)
        returnType(double(1))
    },
    gr_P_RE_wrt_RE_a = function(p = double(1), reTransform = double(1),
                                forceUpdate = logical(0, default = FALSE),
                                forceReset = logical(0, default = FALSE)) {
        do_reset <- forceReset | gr_P_RE_reset_once
        ans <- derivs(logLik_P_RE(p, reTransform), wrt = reTrans_indices, order = 1, model = model,
                      updateNodes = joint_updateNodes, constantNodes = joint_constantNodes,
                      do_update = gr_P_RE_wrt_RE_update_once | gr_P_RE_wrt_RE_update_always | forceUpdate | do_reset,
                      reset=do_reset)
        gr_P_RE_wrt_RE_update_once <<- FALSE
        gr_P_RE_wrt_RE_reset_once <<- FALSE
        return(ans$jacobian[1,])
        returnType(double(1))
    },
    jac_gr_P_RE_wrt_RE_outDir_b= function(p = double(1), reTransform = double(1),
                                                    outDir = double(1),
                                                    forceUpdate = logical(0, default = FALSE),
                                                    forceReset = logical(0, default = FALSE)) {
      do_reset <- forceReset | gr_P_RE_wrt_RE_reset_once
      do_update <- gr_P_RE_wrt_RE_update_once | gr_P_RE_wrt_RE_update_always | forceUpdate | do_reset
      ans <- derivs(gr_P_RE_wrt_RE_a(p, reTransform, forceUpdate=do_update, forceReset=do_reset),
                    wrt = p_reTrans_indices,
                    outDir = outDir,
                    order = 1, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes,
                    do_update = do_update,
                      reset=do_reset)
      gr_P_RE_wrt_RE_update_once <<- FALSE
      gr_P_RE_wrt_RE_reset_once <<- FALSE
      return(ans$jacobian)
      returnType(double(2))
    },
    he_P_RE_wrt_RE_wrt_P_b= function(p = double(1), reTransform = double(1),
                                                    forceUpdate = logical(0, default = FALSE),
                                                    forceReset = logical(0, default = FALSE)) {
      do_reset <- forceReset | gr_P_RE_wrt_RE_reset_once
      do_update <- gr_P_RE_wrt_RE_update_once | gr_P_RE_wrt_RE_update_always | forceUpdate | do_reset
      ans <- derivs(gr_P_RE_wrt_RE_a(p, reTransform, forceUpdate=do_update, forceReset=do_reset),
                    wrt = p_indices,
                    order = 1, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes,
                    do_update = do_update,
                      reset=do_reset)
      gr_P_RE_wrt_RE_update_once <<- FALSE
      gr_P_RE_wrt_RE_reset_once <<- FALSE
      return(ans$jacobian)
      returnType(double(2))
    },
    jac_gr_P_RE_wrt_RE_inDir= function(p = double(1), reTransform = double(1),
                                         inDir = double(1),
                                     forceUpdate = logical(0, default = FALSE),
                                     forceReset = logical(0, default = FALSE)) {
      do_reset <- forceReset | gr_P_RE_wrt_RE_reset_once
      do_update <- gr_P_RE_wrt_RE_update_once | gr_P_RE_wrt_RE_update_always | forceUpdate | do_reset
      ans <- derivs(gr_P_RE_wrt_RE_a(p, reTransform, forceUpdate=do_update, forceReset=do_reset),
                    inDir = inDir,
                    order = 1, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes,
                    do_update = do_update,
                    reset=do_reset)
      gr_P_RE_wrt_RE_update_once <<- FALSE
      gr_P_RE_wrt_RE_reset_once <<- FALSE
      return(ans$jacobian)
      returnType(double(2))
    },
    he_P_RE_wrt_RE2_uptri_b = function(p = double(1), reTransform = double(1),
                                       forceUpdate = logical(0, default = FALSE),
                                       forceReset = logical(0, default = FALSE)) {
      do_reset <- forceReset | he_P_RE_wrt_RE2_uptri_reset_once
      do_update <- he_P_RE_wrt_RE2_uptri_update_once | he_P_RE_wrt_RE2_uptri_update_always | forceUpdate | do_reset
      ans <- derivs(gr_P_RE_wrt_RE_a(p, reTransform, forceUpdate=do_update, forceReset=do_reset),
       wrt = reTrans_indices, order = 1, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes,
                    do_update = do_update,
                    reset=do_reset)
      he_P_RE_wrt_RE2_uptri_update_once <<- FALSE
      he_P_RE_wrt_RE2_uptri_reset_once <<- FALSE
      n <- 1L
      n <- nreTrans
      if(n != dim(ans$jacobian)[1]) stop("error (1) with dimensions in joint hessian")
      if(n != dim(ans$jacobian)[2]) stop("error (2) with dimensions in joint hessian")
      res <- nimNumeric(length = 0.5*n*(n+1), init=FALSE)
      ires <- 1L
      i <- 1L
      j <- 1L
      for(j in 1:n) {
        for(i in 1:j) {
          res[ires] <- ans$jacobian[i, j]
          ires <- ires+1
        }
      }
      return(res)
      returnType(double(1))
    },
    jac_he_P_RE_wrt_RE2_uptri_outDir_c = function(p = double(1), reTransform = double(1),
                                                  outDir = double(1),
                                                  forceUpdate = logical(0, default = FALSE),
                                                  forceReset = logical(0, default = FALSE)) {
      do_reset <- forceReset | he_P_RE_wrt_RE2_uptri_reset_once
      do_update <- he_P_RE_wrt_RE2_uptri_update_once | he_P_RE_wrt_RE2_uptri_update_always | forceUpdate | do_reset
      ans <- derivs(he_P_RE_wrt_RE2_uptri_b(p, reTransform, forceUpdate=do_update, forceReset=do_reset),
                    wrt = p_reTrans_indices,
                    outDir = outDir,
                    order = 1, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes,
                    do_update = do_update,
                    reset=do_reset)
      he_P_RE_wrt_RE2_uptri_update_once <<- FALSE
      he_P_RE_wrt_RE2_uptri_reset_once <<- FALSE
      return(ans$jacobian)
      returnType(double(2))
    },
    #################
    ## Solve the inner optimization for Laplace approximation
    max_logLik_RE = function(p = double(1)) {
      if(any(p != current_P_for_inner)) {
        set_P(p)
      }
      reInitTrans <- get_reInitTrans()
      fn_init <- logLik_RE(reInitTrans)
      if((fn_init == Inf) | (fn_init == -Inf) | (is.nan(fn_init)) | (is.na(fn_init))) {
        optRes <- optimResultNimbleList$new()
        optRes$par <- reInitTrans
        optRes$value <- -Inf
        optRes$convergence <- -1
        return(optRes)
      }
      reInitTrans_for_taping <<- reInitTrans
      optRes <- optim(reInitTrans, logLik_RE, gr = gr_for_optim, he = he_for_optim,
                      method = optimMethod_, control = optimControl_)
      if(optRes$convergence != 0 & warn_optim){
        print("  [Warning] `optim` did not converge for the inner optimization of AGHQ or Laplace approximation.")
      }
      converged <<- optRes$convergence
      return(optRes)
      returnType(optimResultNimbleList())
    },
    ## Outer check for inner convergence
    check_convergence = function(){
      returnType(double())
      return(converged)
    },
    ## These two update methods for max_logLik_RE use the same member data caches
    update_max_logLik_RE = function(p = double(1)) {
      optRes <- max_logLik_RE(p)
      saved_inner_argmax <<- optRes$par
      saved_inner_max_value <<- optRes$value
      saved_inner_max_p <<- p
      saved_inner_negHess <<- -he_RE_c(saved_inner_argmax)
      saved_inner_logdetNegHess <<- log(saved_inner_negHess[1,1])
      return(saved_inner_argmax)
      returnType(double(1))
    },
    ## Laplace approximation (version "2" for historical reasons)
    calcLogLik2 = function(p = double(1)){
      if(!one_time_fixes_done) one_time_fixes()
      if(any(p != saved_inner_max_p) | !cache_inner_max) {
        update_max_logLik_RE(p)
      }
      reTransform <- saved_inner_argmax
      maxValue <- saved_inner_max_value
      if(maxValue == -Inf) return(-Inf) # This would mean inner optimization failed

      if(nQuad_ == 1){
        ## Laplace approximation.
        margLogLik_saved_value <<- maxValue - 0.5 * saved_inner_logdetNegHess + 0.5 * 1 * log(2*pi)
      }else{
        ## Do Quadrature:
        margLogLik_saved_value <<- calcLogLik_AGHQuad(p)
      }

      if(margLogLik_saved_value > max_margLogLik) {
        max_margLogLik<<- margLogLik_saved_value
        max_margLogLik_inner_argmax <<- saved_inner_argmax
      }
      return(margLogLik_saved_value)
      returnType(double())
    },
    calcLogLik_AGHQuad = function(p = double(1)) {
      # This should be called ONLY from calcLogLik2. It is not for stand-along calls
      # because it assumes ths saved_inner_* values are already set.
      #
      ## AGHQ Approximation:  3 steps. build grid (happens once), transform z to re, do quad sum.
      quadGrid$buildGrid(method = quadRule_, nQuad = nQuad_)
      nQ <- quadGrid$gridSize()
      SD <- 1/sqrt(saved_inner_negHess[1,1])
      nodes <<- quadGrid$nodes(0)
      wgts <<- quadGrid$weights(0)
      logDensity_quad <<- numeric(value = 0, length = nQ)

      modeIndex <- quadGrid$modeIndex() ## if even, this is -1
      ans <- 0
      if(any(p != current_P_for_inner)) { # Needed for logLik_RE calls below.
        set_P(p)
      }
      for(i in 1:nQ) {
        if(i != modeIndex) {
          if(quadTransform_ != "identity")
            nodes[i,] <<- saved_inner_argmax + SD*nodes[i,]
          logDensity_quad[i] <<- logLik_RE(reTransform = nodes[i,])
          ans <- ans + exp(logDensity_quad[i] - saved_inner_max_value)*wgts[i]
        }else{
          if(quadTransform_ != "identity")
            nodes[i,] <<- saved_inner_argmax
          logDensity_quad[i] <<- saved_inner_max_value
          ans <- ans + wgts[i]
        }
      }
      ## Given all the saved values, weights and log density, do quadrature sum.
      res <- log(ans) + saved_inner_max_value - 0.5 * saved_inner_logdetNegHess
      quadrature_previous_p <<- p ## Cache this to make sure you have it for
      return(res)
      returnType(double())
    },
    ## Gradient of the Laplace approximation (version 2) w.r.t. parameters
    gr_logLik2 = function(p = double(1)){
      if(!one_time_fixes_done) one_time_fixes()
      if(any(p != saved_inner_max_p) | !cache_inner_max) {
        update_max_logLik_RE(p)
      }
      reTransform <- saved_inner_argmax
      negHessian <- saved_inner_negHess[1,1]
      invNegHessian <- 1/negHessian

      # invNegHessian <- inverse(negHessian)
      gr_logdetNegHess <- -jac_he_P_RE_wrt_RE2_uptri_outDir_c(p, reTransform, c(invNegHessian))
      grlogdetNegHesswrtp <- gr_logdetNegHess[1, p_indices]
      grlogdetNegHesswrtre <- gr_logdetNegHess[1, reTrans_indices][1]
      #hess_joint_logLik_wrt_p_wrt_re(p, reTransform)[,1]
      outDir <- numeric(length = 0) # outDir is not used in this case.
      hesslogLikwrtpre <- jac_gr_P_RE_wrt_RE_outDir_b(p, reTransform, outDir)[1, p_indices]
      if( nQuad_ == 1 ){
        ## Gradient of Laplace Approx
        p1 <- gr_P_RE_b(p, reTransform)[p_indices]
        AGHQuad_saved_gr <<- p1 - 0.5 * (grlogdetNegHesswrtp + hesslogLikwrtpre * (grlogdetNegHesswrtre / negHessian))
      }else{
        ## Gradient of AGHQ Approx.
        ## dre_hat/dp = d^2ll/drep / d^2ll/dre^2
        gr_rehatwrtp <<- hesslogLikwrtpre/negHessian
        ## dsigma_hat/dp (needed at real scale)
        sigma_hat <- 1/sqrt(negHessian)
        gr_sigmahatwrtp <<- -0.5*grlogdetNegHesswrtp*sigma_hat
        gr_sigmahatwrtre <<- -0.5*grlogdetNegHesswrtre*sigma_hat
        ## Sum gradient of each node.
        grp_AGHQuad_sum <- gr_AGHQuad_nodes(p = p, method = 2)
        AGHQuad_saved_gr <<- grp_AGHQuad_sum - 0.5 * (grlogdetNegHesswrtp + grlogdetNegHesswrtre * gr_rehatwrtp)
      }
      # N.B. An extra negation is built into gr_logdet because this is gradient of hessian, but the uptri_Omega_invNegHess is from the negative Hessian.

      return(AGHQuad_saved_gr)
      returnType(double(1))
    },
    ## Partial gradient of AGHQ nodes w respect to p.
    gr_AGHQuad_nodes = function(p = double(1), method = double()){

      ## Need to have quadrature sum for gradient:
      if(any(p != quadrature_previous_p)){
        calcLogLik_AGHQuad(p)
      }

      ## Method 2 implies double taping.
      modeIndex <- quadGrid$modeIndex()
      nQ <- quadGrid$gridSize()
      gr_wgted_wrt_p <- numeric(value = 0, length = dim(p)[1])
      wgts_lik <- numeric(value = 0, length = nQ)
      for(i in 1:nQ) {
        wgts_lik[i] <- exp(logDensity_quad[i] - saved_inner_max_value)*wgts[i]

        ## At the mode (z = 0, don't have additional z*sigma_hat gr complication).
	      if( modeIndex == i ){
          gr_jointlogLikwrtp <- gr_P_RE_b(p, nodes[i,])[p_indices]
          gr_wgted_wrt_p <- gr_wgted_wrt_p + wgts_lik[i]*gr_jointlogLikwrtp
        }else{
          ## Chain Rule: dll/dre * ( dre_hat/dp + dsigma_hat/dp*z_i )
          ## dll/dp
          gr_jointlogLik_i <- gr_P_RE_b(p, nodes[i,])
          gr_logLikwrtrewrtp_i <- gr_jointlogLik_i[reTrans_indices][1] *
                            ( (1 + gr_sigmahatwrtre*quadGrid$nodes(i)[1,1]) * gr_rehatwrtp  +  gr_sigmahatwrtp*quadGrid$nodes(i)[1,1] )
          ## The weighted gradient for the ith sum.
          gr_wgted_wrt_p <- gr_wgted_wrt_p + wgts_lik[i]*( gr_jointlogLik_i[p_indices] +  gr_logLikwrtrewrtp_i )
        }
      }
      return(gr_wgted_wrt_p / sum(wgts_lik[1:nQ]))
      returnType(double(1))
    },
    get_inner_mode = function(atOuterMode = integer(0, default = 0)){
      returnType(double(1))
      if(atOuterMode) return(outer_mode_inner_argmax)
      return(saved_inner_argmax)
    },
    get_inner_negHessian = function(atOuterMode = integer(0, default = 0)){
      returnType(double(2))
      if(atOuterMode) return(outer_mode_inner_negHess)
      return(saved_inner_negHess)
    },
    get_inner_negHessian_chol = function(atOuterMode = integer(0, default = 0)){
      returnType(double(2))
      if(atOuterMode) return(sqrt(outer_mode_inner_negHess))
      return(sqrt(saved_inner_negHess))
    },
    ## Update the maximum mode and neg hess based on the log likelihood passed via optim.
    ##  For efficient saving of values for calculating MLE values of random-effects and INLA simulation of them.
    save_outer_logLik = function(logLikVal = double()){
      if(logLikVal >= max_outer_logLik) {
        max_outer_logLik <<- logLikVal
        outer_mode_inner_negHess <<- saved_inner_negHess
        outer_mode_inner_argmax <<- saved_inner_argmax
        outer_param_max <<- saved_inner_max_p
      }
    },
    get_param_value = function(atOuterMode = integer(0, default = 0)){
      returnType(double(1))
      ## Ensures that the inner value will not match and cached values will not be used.
      if(!cache_inner_max) return(numeric(value = Inf, length = npar))
      if(atOuterMode) return(outer_param_max)
      return(saved_inner_max_p)
    },
    ## Need to reset every time optim is called to recache.
    reset_outer_logLik = function(){
      max_outer_logLik <<- -Inf
    },
    set_randomeffect_values = function(p = double(1)){
      foundIt <- FALSE
      ## Last value called:
      if(all(p == saved_inner_max_p)) {
        re <- reTrans$inverseTransform(saved_inner_argmax)
        foundIt <- TRUE
      }
      ## Best value called:
      if(all(p == outer_param_max)) {
        re <- reTrans$inverseTransform(outer_mode_inner_argmax)
        foundIt <- TRUE
      }
      if(foundIt){
        values(model, paramNodes) <<- p
        model$calculate(paramDeps)
      }else{
        # It would be nice to emit a message here, but different optimizers (e.g. BFGS vs nlminb)
        # behave differently as to whether the previous (last) parameters were always the MLE.
        # print("  [Warning] Have not cached the inner optimization. Running optimization now.")
        update_max_logLik_RE(p)
        re <- reTrans$inverseTransform(saved_inner_argmax)
      }
      ## Ensure the model is up to date for all nodes.
      values(model, randomEffectsNodes) <<- re
      model$calculate(innerCalcNodes)
    }
    ## set_inner_cache = function(cache = logical(0, default = TRUE)){
    ##   cache_inner_max <<- cache
    ## }
  ),
  buildDerivs = list(logLik_RE              = list(),
                     gr_RE_a                = list(),
                     he_RE_b                = list(),
                     he_RE_b_asvec          = list(),
                     logLik_P_RE            = list(),
                     gr_P_RE_a              = list(),
                     gr_P_RE_wrt_RE_a      = list(),
                     he_P_RE_wrt_RE2_uptri_b = list())
) ## End of buildOneAGHQuad1D


## A single Laplace approximation for models with more than one scalar random effect node
buildOneLaplace <- function(model, paramNodes, randomEffectsNodes, calcNodes, paramDeps, control = list()) {
  buildOneAGHQuad(model, nQuad = 1, paramNodes, randomEffectsNodes, calcNodes, paramDeps, control)
}

buildOneAGHQuad <- nimbleFunction(
  contains = AGHQuad_BASE,
  setup = function(model, nQuad = 1, paramNodes, randomEffectsNodes, calcNodes, paramDeps, control = list()) {
    ## Check and add necessary (upstream) deterministic nodes into calcNodes
    ## This ensures that deterministic nodes between paramNodes and calcNodes are used.
    ## optimControl_ <- extractControlElement(control, 'optimControl', nimOptimDefaultControl())
    ## optimMethod_ <- extractControlElement(control, 'optimMethod', 'BFGS')
    ## optimStart_ <- extractControlElement(control, 'optimStart', 'constant')
    ## optimStartValues_ <- extractControlElement(control, 'optimStartValues', 0)
    nQuad_ <- nQuad
    S <- setup_OneAGHQuad(model, paramNodes, randomEffectsNodes, calcNodes, paramDeps, control)
    optimControl_ <- S$optimControl_
    optimMethod_ <- S$optimMethod_
    optimStart_ <- S$optimStart_
    optimStartValues_ <- S$optimStartValues_
    nre  <-  S$nre
    paramDeps  <-  S$paramDeps
    innerCalcNodes  <-  S$innerCalcNodes
    calcNodes  <-  S$calcNodes
    wrtNodes  <-  S$wrtNodes
    reTrans  <-  S$reTrans
    npar  <-  S$npar
    p_indices  <-  S$p_indices
    quadRule_ <- S$quadRule

    ## OVERVIEW OF INDEXING SCHEME:
    ## We have two situations: Sometimes we need the logLik as a function of random effects only ("inner" or "_RE_"),
    ##        and sometimes as a function of parameters and random effects ("outer" or "_P_RE_").
    ## The random effects are always transformed, so we label them "reTrans".
    ## When using a function of params and reTrans, the order is always c(params, reTrans).
    ## SIZES
    ## nreTrans: length of reTrans
    ## nre: length of randomEffectsNodes (before transformation)
    ## npar: length of params
    ## reTrans_indices: indices of reTrans when used without params (1:nreTrans, possibly after one_time_fixes).
    ## reTrans_indices_inner: indices of reTrans when used with params (npar + (1:nreTrans), possibly after one_time_fixes).
    ## p_reTrans_indices: indices of params and reTrans when used together (1:(npar + nreTrans)).
    ## p_indices: indices of params when used with reTrans (1:npar). (set above)
    nreTrans <- reTrans$getTransformedLength()
    if(nreTrans > 1) reTrans_indices <- as.numeric((npar+1):(npar+nreTrans))
    else reTrans_indices <- as.numeric(c(npar+1, -1))
    ## Indices of randomEffectsNodes inside randomEffectsNodes for use in getting the derivative of
    ## the inner log-likelihood (paramNodes fixed) w.r.t. randomEffectsNodes.
    if(nreTrans > 1) reTrans_indices_inner <- as.numeric(1:nreTrans)
    else reTrans_indices_inner <- as.numeric(c(1, -1))
    p_reTrans_indices <- as.numeric(1:(npar + nreTrans))

    ## Set up start values for the inner optimization of Laplace approximation
    if(!is.character(optimStart_) | length(optimStart_) != 1) stop("problem with optimStart ", optimStart_)
    startID <- switch(optimStart_, last=1, last.best=2, constant=3, random=4, model=5)
    if(startID == 5) {
      constant_init_reTrans <- reTrans$transform(c(values(model, randomEffectsNodes)))
      startID <- 3  
    } else {
      if(length(optimStartValues_) == 1)
        constant_init_reTrans <- rep(optimStartValues_, nreTrans)
      else
        constant_init_reTrans <- optimStartValues_
    }
    if(length(constant_init_reTrans) != nreTrans)
      stop("buildOneAGHQuad: Found ", length(constant_init_reTrans), " initial values for inner optimization in Laplace or AGHQuad when expecting ", nreTrans)
    if(length(constant_init_reTrans) == 1) constant_init_reTrans <- c(constant_init_reTrans, -1)

    ## Update and constant nodes info for obtaining derivatives using AD
    inner_derivsInfo    <- makeModelDerivsInfo(model = model, wrtNodes = randomEffectsNodes, calcNodes = innerCalcNodes)
    inner_updateNodes   <- inner_derivsInfo$updateNodes
    inner_constantNodes <- inner_derivsInfo$constantNodes
    joint_derivsInfo    <- makeModelDerivsInfo(model = model, wrtNodes = wrtNodes, calcNodes = calcNodes)
    joint_updateNodes   <- joint_derivsInfo$updateNodes
    joint_constantNodes <- joint_derivsInfo$constantNodes

    ## Flags used to manage various steps
    ##
    ## one_time_fixes starts FALSE and then after being set TRUE never needs to be set FALSE.
    one_time_fixes_done <- FALSE
    ##
    ## For each group of related deriv functions, there are three flags:
    ## *_update_once: On the next call (only), do_update will be TRUE (updating the updateNodes).
    ##    This is used to ensure AD tapes are updated when needed but not otherwise.
    ## *_update_always: Over-ride *_update_once and instead update on every call.
    ##    This is not actively used, and can only be manually changed (not part of the API).
    ##    It is really for debugging purposes, in case there is a glitch with when updating is done.
    ## *_reset_once: On the next call (only), reset the AD tape, meaning tape it again from scratch.
    ##    This needs to be done if the constantNodes have changed, or if somehow the tapes
    ##    were made with bad inputs results in NaNs, which sometimes end up baked into a tape and make it useless.
    ##
    ## Note that each deriv function can also accept an argument for do_update and reset, which then over-ride
    ##   any of these flags.
    ##
    ## Within a group of deriv calls, flags are passed through nested calls nested for meta-taping.
    ##   In the case of reset, this really makes sense: when an outer tape is reset, the inner tapes should be too.
    ##   In the case of update, this is not really necessary when reset==FALSE, because only the outermost update is used
    ##     when simply playing a tape, but it doesn't matter because the inner update arguments are ignored if reset==FALSE.
    ##
    ## See NOMENCLATURE note below for labeling of deriv cases.
    ##
    ## Flags for all gradients as a function of random effects only, i.e. inner gradients.
    gr_RE_update_once <- TRUE
    gr_RE_update_always <- FALSE
    gr_RE_reset_once <- TRUE

    ## Flags for all Hessians as a function of random effects only, i.e. inner Hessians.
    he_RE_update_once <- TRUE
    he_RE_update_always <- FALSE
    he_RE_reset_once <- TRUE

    ## Flags for all gradients as a function of parameters and random effects.
    gr_P_RE_update_once <- TRUE
    gr_P_RE_update_always <- FALSE
    gr_P_RE_reset_once <- TRUE

    ## Flags for all gradients as a function of parameters and random effects, but only wrt random effects.
    gr_P_RE_wrt_RE_update_once <- TRUE
    gr_P_RE_wrt_RE_update_always <- FALSE
    gr_P_RE_wrt_RE_reset_once <- TRUE

    ## Flags for all Hessians as a function of parameters and random effects wrt random effects, flattened upper triangular.
    he_P_RE_wrt_RE2_uptri_update_once <- TRUE
    he_P_RE_wrt_RE2_uptri_update_always <- FALSE
    he_P_RE_wrt_RE2_uptri_reset_once <- TRUE

    ## Caches for results of inner optimization:
    cache_inner_max <- TRUE
    saved_inner_argmax <- constant_init_reTrans
    saved_inner_max_value <- -Inf #numeric(1)
    saved_inner_max_p <- if(npar > 1) rep(Inf, npar) else as.numeric(c(Inf, -1))
    saved_inner_negHess <- matrix(0, nrow = nre, ncol = nre)
    saved_inner_negHess_chol <- matrix(0, nrow = nre, ncol = nre)
    saved_inner_logdetNegHess <- 0

    ## Cache for set_P
    ## Any methods that will call logLik_RE or derivs that take argument RE only should
    ## check if all P match current_P_for_inner and call set_P if not.
    current_P_for_inner <- saved_inner_max_p

    ## Cache to ensure taping is done from init (RE)
    reInitTrans_for_taping <- constant_init_reTrans

    ## Caches to help with outer optimization:
    ## Record the maximum Laplace or AGHQ loglikelihood value for obtaining inner optimization start values
    ## These are labeled "margLogLik" for clarity.
    max_margLogLik<- -Inf
    max_margLogLik_inner_argmax <- constant_init_reTrans #if(nreTrans > 1) rep(Inf, nreTrans) else as.numeric(c(0, -1))
    margLogLik_saved_value <- -Inf


    ## Cache values for relevant to outer calls.
    max_outer_logLik <- -Inf
    outer_mode_inner_negHess <- matrix(0, nrow = nre, ncol = nre)
    outer_mode_inner_negHess_chol <- matrix(0, nrow = nre, ncol = nre)
    outer_mode_inner_argmax <- if(nreTrans > 1) numeric(nreTrans) else as.numeric(c(0, -1))
    outer_param_max <- if(npar > 1) rep(Inf, npar) else as.numeric(c(Inf, -1))

    ## Build Quadrature grid for any dimension:
    ## This is set up to add other quad grids in the future. quadRule := "AGHQ" to start.
    quadGrid <- configureQuadGrid(d = nreTrans, levels = nQuad_, quadRule = quadRule_)
    nodes <-  matrix(0, nrow = nQuad_, ncol = nreTrans)
    wgts <- numeric(nQuad_)
    logDensity_quad <- numeric(nQuad_)
    if(nQuad_ == 1) {
      wgts <- c(0,-1)
      logDensity_quad <- c(0,-1)
    }
    quadTransform_ <- extractControlElement(control, "quadTransform", "cholesky")

    converged <- 0
    warn_optim <- extractControlElement(control, 'optimWarning', FALSE) ## Warn about inner optimization issues
  },
  run = function(){},
  methods = list(
    fix_one_vec = function(x = double(1)) {
      if(length(x) == 2) {
        if(x[2] == -1) {
          ans <- numeric(length = 1, value = x[1])
          return(ans)
        }
      }
      return(x)
      returnType(double(1))
    },
    one_time_fixes = function() {
      if(one_time_fixes_done) return()
      if(nre == 1) {
        reTrans_indices <<- fix_one_vec(reTrans_indices)
        reTrans_indices_inner <<- fix_one_vec(reTrans_indices_inner)
        saved_inner_argmax <<- fix_one_vec(saved_inner_argmax)
        reInitTrans_for_taping <<- fix_one_vec(reInitTrans_for_taping)
        max_margLogLik_inner_argmax <<- fix_one_vec(max_margLogLik_inner_argmax)
        constant_init_reTrans <<- fix_one_vec(constant_init_reTrans)
        outer_mode_inner_argmax <<- fix_one_vec(outer_mode_inner_argmax)
      }
      if(npar == 1) {
        p_indices <<- fix_one_vec(p_indices)
        saved_inner_max_p <<- fix_one_vec(saved_inner_max_p)
        current_P_for_inner <<- fix_one_vec(current_P_for_inner)
        outer_param_max <<- fix_one_vec(outer_param_max)
      }
      if(nQuad_ == 1) {
        wgts <<- fix_one_vec(wgts)
        logDensity_quad <<- fix_one_vec(logDensity_quad)
      }
      reInit <- values(model, randomEffectsNodes)
      set_reInit(reInit)
      one_time_fixes_done <<- TRUE
    },
    updateSettings = function(optimMethod = character(0, default="NULL"),
                              optimStart = character(0, default="NULL"),
                              optimStartValues = double(1, default=Inf),
                              optimWarning = integer(0, default = -1),
                              useInnerCache = integer(0, default=-1),
                              nQuad = integer(0, default=-1),
                              quadTransform = character(0, default="NULL"),
                              optimControl = optimControlNimbleList(default=nimOptimDefaultControl()),
                              replace_optimControl = logical(0, default=FALSE)) {
      # Checking should have been done already. Or, if this is being called directly,
      # it will be for development or advanced uses and we can skip checking.
      if(optimMethod != "NULL") optimMethod_ <<- optimMethod
      if(optimStart != "NULL") {
        if(optimStart == "last") startID <<- 1 # last
        else if(optimStart == "last.best") startID <<- 2 # last.best
        else if(optimStart == "constant") startID <<- 3 # use fixed vector optimStart provided at setup time
        else if(optimStart == "random") startID <<- 4
        else if(optimStart == "model") {
          startID <<- 3
          constant_init_reTrans <<- reTrans$transform(values(model, randomEffectsNodes))
        }
      }
      if((length(optimStartValues) != 1) | (optimStartValues[1] != Inf) ) {
        if((length(optimStartValues) == 1) & (optimStartValues[1] == -Inf) ) { # numeric code for "model" setting
          constant_init_reTrans <<- reTrans$transform(values(model, randomEffectsNodes))
        } else {
          if(startID <= 3) {
            constant_init_reTrans <<- optimStartValues
            if(length(constant_init_reTrans) == 1)
              if(nreTrans > 1)
                constant_init_reTrans <<- rep(constant_init_reTrans, nreTrans)
          }
        }
      }
      if((!one_time_fixes_done) & (length(constant_init_reTrans) == 1)){
         constant_init_reTrans <<- c(constant_init_reTrans, -1)
      }
      if(optimWarning != -1) {
        warn_optim <<- optimWarning != 0
      }
      if(useInnerCache != -1) {
        cache_inner_max <<- useInnerCache != 0
      }
      if(nQuad != -1) {
        quadGrid$buildGrid(method = quadRule_, nQuad = nQuad)
        nQuad_ <<- nQuad
      }
      if(quadTransform != "NULL") {
        quadTransform_ <<- quadTransform
      }
      if(replace_optimControl) {
        if(optimControl$fnscale == 1) optimControl$fnscale <- -1
        optimControl_ <<- optimControl
      }
    },
    set_reInit = function(re = double(1)) {
      reInitTrans <- reTrans$transform(re)
      saved_inner_argmax <<- reInitTrans
    },
    get_reInitTrans = function() {
      if(startID == 1) ans <- saved_inner_argmax                ## last
      else if(startID == 2) ans <- max_margLogLik_inner_argmax            ## last best
      else if(startID == 3) ans <- constant_init_reTrans                      ## constant
      else if(startID == 4){                                              ## random
        model$simulate(randomEffectsNodes)
        ans <- reTrans$transform(values(model, randomEffectsNodes))
      }
      return(ans)
      returnType(double(1))
    },
    get_reTransLength = function(){
      returnType(double(0))
      return(nreTrans)
    },
    ## NOMENCLATURE:
    ## "P" = parameters, "RE" = random effects
    ## Together, "P", "RE", or "P_RE" indicate the arguments.
    ## "gr" indicates a gradient, "he" indicates a Hessian, and "jac" indicates a Jacobian.
    ## (Gradients and Jacobians are both first-order derivatives, but Gradients return a vector based on one output and Jacobians return a matrix based on multiple outputs.)
    ## LOG LIKELIHOOD FUNCTIONS:
    ## - logLik_RE is the log-likelihood as a function of parameters only.
    ##    N.B. logLik_RE must be preceded by calling set_P(p) to set the parameters.
    ##         Derivatives of logLik_RE below must be updated once after set_P(p) using the update flags.
    ## - logLik_P_RE is the log-likelihood as a function of parameters and random effects.
    ## - "logLik" is always the "inner" log-likelihood; what differs is whether viewed as a function of RE or P_RE.
    ## - In some flags and/or caches, the Laplace or AGQH result is called "margLogLik", i.e. marginal approximation.
    ## DERIVATIVE FUNCTIONS:
    ## All derivatives are ultimately of the logLik, so this is omitted from their names.
    ## "a" = first-level of taping, "b" = second-level of taping, "c" = third-level of taping
    ## "wrt" = with respect to, only included if not wrt all arguments.
    ## "outDir" = output direction (which will be used in reverse mode AD).
    ##.     If outDir is not included (or empty), then each output direction is used, resulting in a full Jacobian or Hessian.
    ##      There are only a couple of uses of outDir, for efficient calculation of the outer gradient (gradient of Laplace approx.)
    ## N.B. When an "outDir" is allowed, we indicate "gr_gr" or "gr_he".
    ## - gr_RE_a is the gradient of logLik_RE with respect to the random effects, first-level taping.
    ## - gr_RE_b is the gradient of logLik_RE with respect to the random effects, second-level taping (gradient of the gradient).
    ## - gr_for_optim checks for forced reset of taping, and then calls gr_RE_b.
    ## - he_RE_b is the Hessian of logLik_RE with respect to the random effects, second-level taping (gradient of the gradient).
    ##   (There is no he_RE_a, as we only do the Hessian by double-taping)
    ## - he_RE_b_asvec is he_RE_b returned as a flattened vector.
    ## - he_RE_c is the Hessian of logLik_RE with respect to the random effects, third-level taping (value (0) of the gradient (1) of the gradient (1)).
    ## - he_for_optim checks for forced reset of taping, and then calls he_RE_c.
    ## - negHess(P, reTransform) does set_P, and then returns -he_RE_c. The "_P_RE" is omitted
    ##      because it is called from the outer level.
    ## - gr_P_RE_a is the gradient of logLik_P_RE with respect to parameters and random effects, first-level taping.
    ## - gr_P_RE_b is the gradient of logLik_P_RE with respect to parameters and random effects, second-level taping (value of the gradient).
    ## - gr_P_RE_wrt_RE_a is the gradient of logLik_P_RE with respect to the random effects only, first-level taping.
    ## - jac_gr_P_RE_wrt_RE_outDir_b is the Jacobian of gr_P_RE_wrt_RE_a (i.e. 2nd derivs), from output direction outDir, second-level taping.
    ##   gr_P_RE_wrt_RE_a has length nreTrans, and
    ##  jac_gr_P_RE_wrt_RE_outDir_b is 1 x ntotal if outDir is provided and nreTrans x ntotal if not.
    ## - he_P_RE_wrt_RE2_uptri_b is the Hessian of logLik_P_RE with respect to the random effects, only the upper triangular part, flattened to a vector, second-level taping (gradient of the gradient).
    ## - jac_he_P_RE_wrt_RE2_uptri_outDir_c is the Jacobian of he_P_RE_wrt_RE2_uptri_b, from output direction outDir, third-level taping.
    ##   he_P_RE_wrt_RE2_uptri_b has length nreTrans*(nreTrans+1)/2, and
    ##     jac_he_P_RE_wrt_RE2_uptri_outDir_c is 1 x ntotal if outDir is provided and nreTrans*(nreTrans+1)/2 x ntotal if not.
    ## - he_P_RE_wrt_RE_wrt_P_b is the Hessian of logLik_P_RE with respect to the random effects x parameters, second-level taping (gradient of the gradient). This is called from the outer level.
    ## GUIDANCE ON WHAT TO CALL:
    ## For working with the inner optimization, use:
    ## - set_P(p) to set the parameters.
    ## - logLik_RE(reTransform) to get the log-likelihood.
    ## - gr_RE_b(reTransform) or gr_for_optim(reTransform) [if reset check is needed] get the gradient.
    ## - he_RE_c(reTransform) or he_for_optim(reTransform) [if reset check is needed] to get the Hessian.
    ## For working from the outer level, i.e. outer gradient stuff, use
    ## - logLik_P_RE(p, reTransform) to get the log-likelihood.
    ## - gr_P_RE_b(p, reTransform) to get the gradient.
    ## - negHess to get the negative Hessian. This is called from the outer level.

    set_P = function(p = double(1)) {
      values(model, paramNodes) <<- p
      model$calculate(paramDeps)
      gr_RE_update_once <<- TRUE
      he_RE_update_once <<- TRUE
      current_P_for_inner <<- p
    },
    reset = function(gr_RE = logical(0, default = TRUE),
                     he_RE = logical(0, default = TRUE),
                     gr_P_RE = logical(0, default = TRUE),
                     gr_P_RE_wrt_RE = logical(0, default = TRUE),
                     he_P_RE_wrt_RE2_uptri = logical(0, default = TRUE)){
      gr_RE_reset_once <<- gr_RE
      he_RE_reset_once <<- he_RE
      gr_P_RE_reset_once <<- gr_P_RE
      gr_P_RE_wrt_RE_reset_once <<- gr_P_RE_wrt_RE
      he_P_RE_wrt_RE2_uptri_reset_once <<- he_P_RE_wrt_RE2_uptri
      ## Reset the inner optimization cache.
    },
    ## Joint log-likelihood with values of parameters fixed: used only for inner optimization
    logLik_RE = function(reTransform = double(1)) {
      # previously inner_logLik
      re <- reTrans$inverseTransform(reTransform)
      values(model, randomEffectsNodes) <<- re
      ans <- model$calculate(innerCalcNodes) + reTrans$logDetJacobian(reTransform)
      return(ans)
      returnType(double())
    },
    gr_RE_a = function(reTransform = double(1),
                       forceUpdate = logical(0, default = FALSE),
                       forceReset = logical(0, default = FALSE)) {
      # renamed: previously had "internal" suffix
      #  previously gr_inner_logLik_internal
      do_reset <- forceReset | gr_RE_reset_once
      ans <- derivs(logLik_RE(reTransform), wrt = reTrans_indices_inner, order = 1, model = model,
                    updateNodes = inner_updateNodes, constantNodes = inner_constantNodes,
                    do_update = gr_RE_update_once | gr_RE_update_always | forceUpdate | do_reset,
                    reset=do_reset)
      gr_RE_update_once <<- FALSE
      gr_RE_reset_once <<- FALSE
      return(ans$jacobian[1,])
      returnType(double(1))
    },
    gr_RE_b = function(reTransform = double(1),
                       forceUpdate = logical(0, default = FALSE),
                       forceReset = logical(0, default = FALSE)) {
      ## renamed: previously had no suffix
      ## previusly gr_inner_logLik
      do_reset <- forceReset | gr_RE_reset_once
      do_update <- gr_RE_update_once | gr_RE_update_always | forceUpdate | do_reset
      ans <- derivs(gr_RE_a(reTransform, forceUpdate=do_update, forceReset=do_reset),
                    wrt = reTrans_indices_inner, order = 0, model = model,
                    updateNodes = inner_updateNodes, constantNodes = inner_constantNodes,
                    do_update = do_update,
                    reset=do_reset)
      gr_RE_update_once <<- FALSE
      gr_RE_reset_once <<- FALSE
      return(ans$value)
      returnType(double(1))
    },
    gr_for_optim = function(reTransform = double(1)) {
      ## If the tape will be reset, we ensure we record it at the init params.
      ## I am not sure why except this came from experience.
      if(gr_RE_reset_once) {
        gr_RE_b(reInitTrans_for_taping)
      }
      return(gr_RE_b(reTransform))
      returnType(double(1))
    },
    he_RE_b = function(reTransform = double(1),
                       forceUpdate = logical(0, default = FALSE),
                       forceReset = logical(0, default = FALSE)) {
      ## renamed: previously had "internal" suffix
      ## previously he_inner_logLik_internal
      ## reimplemented: now uses order(1) from gr_inner_logLik
      do_reset <- forceReset | he_RE_reset_once
      do_update <- he_RE_update_once | he_RE_update_always | forceUpdate | do_reset
      ans <- derivs(gr_RE_a(reTransform, forceUpdate=do_update, forceReset=do_reset),
                    wrt = reTrans_indices_inner, order = 1, model = model,
                    updateNodes = inner_updateNodes, constantNodes = inner_constantNodes,
                    do_update = do_update,
                    reset=do_reset)
      he_RE_update_once <<- FALSE
      he_RE_reset_once <<- FALSE
      res <- ans$jacobian
      return(res)
      returnType(double(2))
    },
    he_RE_b_asvec = function(reTransform = double(1),
                             forceUpdate = logical(0, default = FALSE),
                             forceReset = logical(0, default = FALSE)) {
      ans <- he_RE_b(reTransform, forceUpdate=forceUpdate, forceReset=forceReset)
      res <- nimNumeric(value = ans, length = nreTrans * nreTrans)
      return(res)
      returnType(double(1))
    },
    he_RE_c = function(reTransform = double(1),
                       forceUpdate = logical(0, default = FALSE),
                       forceReset = logical(0, default = FALSE)) {
      ## renamed: previously had no suffix
      # previously he_inner_logLik
      do_reset <- forceReset | he_RE_reset_once
      do_update <- he_RE_update_once | he_RE_update_always | forceUpdate | do_reset
      ans <- derivs(he_RE_b_asvec(reTransform, forceUpdate=do_update, forceReset=do_reset),
                    wrt = reTrans_indices_inner,
                    order = 0, model = model,
                    updateNodes = inner_updateNodes, constantNodes = inner_constantNodes,
                    do_update = do_update,
                    reset=do_reset)
      he_RE_update_once <<- FALSE
      he_RE_reset_once <<- FALSE
      res <- matrix(value = ans$value, nrow = nreTrans, ncol = nreTrans)
      return(res)
      returnType(double(2))
    },
    he_for_optim = function(reTransform = double(1)) {
      ## If the tape will be reset, we ensure we record it at the init params.
      ## I am not sure why except this came from experience.
      if(he_RE_reset_once) {
        he_RE_c(reInitTrans_for_taping)
      }
      return(he_RE_c(reTransform))
      returnType(double(2))
    },
    negHess = function(p = double(1),
                        reTransform = double(1),
                        forceReset = logical(0, default = FALSE)) {
      set_P(p) # This sets the update flag to TRUE.
      ans <- -he_RE_c(reTransform, forceUpdate=TRUE, forceReset=forceReset)
      return(ans)
      returnType(double(2))
    },
    logLik_P_RE = function(p = double(1), reTransform = double(1)) {
        re <- reTrans$inverseTransform(reTransform)
        values(model, paramNodes) <<- p
        values(model, randomEffectsNodes) <<- re
        ans <- model$calculate(calcNodes) +  reTrans$logDetJacobian(reTransform)
        return(ans)
        returnType(double())
    },
    gr_P_RE_a = function(p = double(1), reTransform = double(1),
                         forceUpdate = logical(0, default = FALSE),
                         forceReset = logical(0, default = FALSE)) {
        # previously gr_joint_logLik_wrt_p_re_internal (?)
        do_reset <- forceReset | gr_P_RE_reset_once
        do_update <- gr_P_RE_update_once | gr_P_RE_update_always | forceUpdate | do_reset
        ans <- derivs(logLik_P_RE(p, reTransform), wrt = p_reTrans_indices, order = 1, model = model,
                      updateNodes = joint_updateNodes, constantNodes = joint_constantNodes,
                      do_update = do_update,
                      reset=do_reset)
        gr_P_RE_update_once <<- FALSE
        gr_P_RE_reset_once <<- FALSE
        return(ans$jacobian[1,])
        returnType(double(1))
    },
    gr_P_RE_b = function(p = double(1), reTransform = double(1),
                         forceUpdate = logical(0, default = FALSE),
                         forceReset = logical(0, default = FALSE)) {
      ## previously gr_joint_logLik_wrt_p_re
        do_reset <- forceReset | gr_P_RE_reset_once
        do_update <- gr_P_RE_update_once | gr_P_RE_update_always | forceUpdate | do_reset
        ans <- derivs(gr_P_RE_a(p, reTransform, forceUpdate=do_update, forceReset=do_reset), wrt = p_reTrans_indices, order = 0, model = model,
                      updateNodes = joint_updateNodes, constantNodes = joint_updateNodes,
                      do_update = do_update,
                      reset=do_reset)
        gr_P_RE_update_once <<- FALSE
        gr_P_RE_reset_once <<- FALSE
        return(ans$value)
        returnType(double(1))
    },
    gr_P_RE_wrt_RE_a = function(p = double(1), reTransform = double(1),
                                forceUpdate = logical(0, default = FALSE),
                                forceReset = logical(0, default = FALSE)) {
        do_reset <- forceReset | gr_P_RE_reset_once
        ans <- derivs(logLik_P_RE(p, reTransform), wrt = reTrans_indices, order = 1, model = model,
                      updateNodes = joint_updateNodes, constantNodes = joint_constantNodes,
                      do_update = gr_P_RE_wrt_RE_update_once | gr_P_RE_wrt_RE_update_always | forceUpdate | do_reset,
                      reset=do_reset)
        gr_P_RE_wrt_RE_update_once <<- FALSE
        gr_P_RE_wrt_RE_reset_once <<- FALSE
        return(ans$jacobian[1,])
        returnType(double(1))
    },
    jac_gr_P_RE_wrt_RE_outDir_b= function(p = double(1), reTransform = double(1),
                                                    outDir = double(1),
                                                    forceUpdate = logical(0, default = FALSE),
                                                    forceReset = logical(0, default = FALSE)) {
      do_reset <- forceReset | gr_P_RE_wrt_RE_reset_once
      do_update <- gr_P_RE_wrt_RE_update_once | gr_P_RE_wrt_RE_update_always | forceUpdate | do_reset
      ans <- derivs(gr_P_RE_wrt_RE_a(p, reTransform, forceUpdate=do_update, forceReset=do_reset),
                    wrt = p_reTrans_indices,
                    outDir = outDir,
                    order = 1, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes,
                    do_update = do_update,
                      reset=do_reset)
      gr_P_RE_wrt_RE_update_once <<- FALSE
      gr_P_RE_wrt_RE_reset_once <<- FALSE
      return(ans$jacobian)
      returnType(double(2))
    },
    he_P_RE_wrt_RE_wrt_P_b= function(p = double(1), reTransform = double(1),
                                                    forceUpdate = logical(0, default = FALSE),
                                                    forceReset = logical(0, default = FALSE)) {
      do_reset <- forceReset | gr_P_RE_wrt_RE_reset_once
      do_update <- gr_P_RE_wrt_RE_update_once | gr_P_RE_wrt_RE_update_always | forceUpdate | do_reset
      ans <- derivs(gr_P_RE_wrt_RE_a(p, reTransform, forceUpdate=do_update, forceReset=do_reset),
                    wrt = p_indices,
                    order = 1, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes,
                    do_update = do_update,
                      reset=do_reset)
      gr_P_RE_wrt_RE_update_once <<- FALSE
      gr_P_RE_wrt_RE_reset_once <<- FALSE
      return(ans$jacobian)
      returnType(double(2))
    },
    jac_gr_P_RE_wrt_RE_inDir= function(p = double(1), reTransform = double(1),
                                         inDir = double(1),
                                     forceUpdate = logical(0, default = FALSE),
                                     forceReset = logical(0, default = FALSE)) {
      do_reset <- forceReset | gr_P_RE_wrt_RE_reset_once
      do_update <- gr_P_RE_wrt_RE_update_once | gr_P_RE_wrt_RE_update_always | forceUpdate | do_reset
      ans <- derivs(gr_P_RE_wrt_RE_a(p, reTransform, forceUpdate=do_update, forceReset=do_reset),
                    inDir = inDir,
                    order = 1, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes,
                    do_update = do_update,
                    reset=do_reset)
      gr_P_RE_wrt_RE_update_once <<- FALSE
      gr_P_RE_wrt_RE_reset_once <<- FALSE
      return(ans$jacobian)
      returnType(double(2))
    },
    he_P_RE_wrt_RE2_uptri_b = function(p = double(1), reTransform = double(1),
                                       forceUpdate = logical(0, default = FALSE),
                                       forceReset = logical(0, default = FALSE)) {
      do_reset <- forceReset | he_P_RE_wrt_RE2_uptri_reset_once
      do_update <- he_P_RE_wrt_RE2_uptri_update_once | he_P_RE_wrt_RE2_uptri_update_always | forceUpdate | do_reset
      ans <- derivs(gr_P_RE_wrt_RE_a(p, reTransform, forceUpdate=do_update, forceReset=do_reset),
       wrt = reTrans_indices, order = 1, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes,
                    do_update = do_update,
                    reset=do_reset)
      he_P_RE_wrt_RE2_uptri_update_once <<- FALSE
      he_P_RE_wrt_RE2_uptri_reset_once <<- FALSE
      n <- 1L
      n <- nreTrans
      if(n != dim(ans$jacobian)[1]) stop("error (1) with dimensions in joint hessian")
      if(n != dim(ans$jacobian)[2]) stop("error (2) with dimensions in joint hessian")
      res <- nimNumeric(length = 0.5*n*(n+1), init=FALSE)
      ires <- 1L
      i <- 1L
      j <- 1L
      for(j in 1:n) {
        for(i in 1:j) {
          res[ires] <- ans$jacobian[i, j]
          ires <- ires+1
        }
      }
      return(res)
      returnType(double(1))
    },
    jac_he_P_RE_wrt_RE2_uptri_outDir_c = function(p = double(1), reTransform = double(1),
                                                  outDir = double(1),
                                                  forceUpdate = logical(0, default = FALSE),
                                                  forceReset = logical(0, default = FALSE)) {
      do_reset <- forceReset | he_P_RE_wrt_RE2_uptri_reset_once
      do_update <- he_P_RE_wrt_RE2_uptri_update_once | he_P_RE_wrt_RE2_uptri_update_always | forceUpdate | do_reset
      ans <- derivs(he_P_RE_wrt_RE2_uptri_b(p, reTransform, forceUpdate=do_update, forceReset=do_reset),
                    wrt = p_reTrans_indices,
                    outDir = outDir,
                    order = 1, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes,
                    do_update = do_update,
                    reset=do_reset)
      he_P_RE_wrt_RE2_uptri_update_once <<- FALSE
      he_P_RE_wrt_RE2_uptri_reset_once <<- FALSE
      return(ans$jacobian)
      returnType(double(2))
    },
    ################
    ## Solve the inner optimization for Laplace approximation
    max_logLik_RE = function(p = double(1)) {
      if(any(p != current_P_for_inner)) {
        set_P(p)
      }
      reInitTrans <- get_reInitTrans()
      fn_init <- logLik_RE(reInitTrans)
      if((fn_init == Inf) | (fn_init == -Inf) | (is.nan(fn_init)) | (is.na(fn_init))) {
        optRes <- optimResultNimbleList$new()
        optRes$par <- reInitTrans
        optRes$value <- -Inf
        optRes$convergence <- -1
        return(optRes)
      }
      ## This ensure that on the very first calls, we record the AD tapes
      ## from the first reInitTrans, presumably safe.
      ## Otherwise, we don't actually know if an optim method will call
      ## the gradient or Hessian at the initial parameters.
      ## However, we defer these steps until they are actually called, because
      ## an optimizer might not even use gr or he.
      reInitTrans_for_taping <<- reInitTrans

      optRes <- optim(reInitTrans, logLik_RE, gr = gr_for_optim, he = he_for_optim,
                      method = optimMethod_, control = optimControl_)
      if(optRes$convergence != 0 & warn_optim){
        print("  [Warning] `optim` did not converge for the inner optimization of AGHQ or Laplace approximation.")
      }
      converged <<- optRes$convergence
      return(optRes)
      returnType(optimResultNimbleList())
    },
    ## Outer check on innner convergence.
    check_convergence = function(){
      returnType(double())
      return(converged)
    },
    update_max_logLik_RE = function(p = double(1)) {
      optRes <- max_logLik_RE(p)
      saved_inner_argmax <<- optRes$par
      saved_inner_max_value <<- optRes$value
      saved_inner_max_p <<- p
      saved_inner_negHess <<- -he_RE_c(saved_inner_argmax)
      saved_inner_negHess_chol <<- chol(saved_inner_negHess)
      saved_inner_logdetNegHess <<- 2 * sum(log(diag(saved_inner_negHess_chol)))
      return(saved_inner_argmax)
      returnType(double(1))
    },
    ## Laplace approximation (version "2" for historical reasons)
    calcLogLik2 = function(p = double(1)){
      if(!one_time_fixes_done) one_time_fixes()
      if(any(p != saved_inner_max_p) | !cache_inner_max) {
        update_max_logLik_RE(p)
      }
      maxValue <- saved_inner_max_value
      if(maxValue == -Inf) return(-Inf) # This would mean inner optimization failed
      if(nQuad_ == 1){
        margLogLik_saved_value <<- maxValue - 0.5 * saved_inner_logdetNegHess + 0.5 * nreTrans * log(2*pi)
      }else{
        margLogLik_saved_value <<- calcLogLik_AGHQuad(p)
      }
      if(margLogLik_saved_value > max_margLogLik) {
        max_margLogLik<<- margLogLik_saved_value
        max_margLogLik_inner_argmax <<- saved_inner_argmax
      }

      return(margLogLik_saved_value)
      returnType(double())
    },
    transformNode = function(z = double(1), eigenvec = double(2),
                             eigenval = double(1), method = character(0, "spectral")){
      if(method == "spectral"){
        theta <- numeric(value = 0, length = nreTrans)
        for( i in 1:nreTrans ){
          theta[i] <- saved_inner_argmax[i] + sum(eigenvec[i,] * z/sqrt(eigenval))
        }
      } else{
        if(method == "identity")
          theta <- z
        else ## Cholesky
          theta <- saved_inner_argmax + backsolve(saved_inner_negHess_chol, z)
      }
      returnType(double(1))
      return(theta)
    },
    calcLogLik_AGHQuad = function(p = double(1)) {
      ## AGHQ Approximation:  3 steps. build grid (happens once), transform z to re, save log density.
      quadGrid$buildGrid(method = quadRule_, nQuad = nQuad_)
      modeIndex <- quadGrid$modeIndex()

      nQ <- quadGrid$gridSize()
      nodes <<- quadGrid$nodes(0)  ## On standard scale but will be transformed.
      wgts <<- quadGrid$weights(0)
      logDensity_quad <<- numeric(value = 0, length = nQ)

      if(quadTransform_ == "spectral"){
        negH <- t(saved_inner_negHess_chol) %*% saved_inner_negHess_chol ## negHess isn't always computed.
				E <- eigen(negH, symmetric = TRUE) ## Should be symmetric...
				L <- E$values	# Always biggest to smallest.
				V <- E$vectors
      }else{
        L <- numeric(0, length = 1)
        V <- matrix(0, nrow = 1, ncol = 1)
      }
      ans <- 0
      if(any(p != current_P_for_inner)) {
        set_P(p)
      }
      for(i in 1:nQ) {
        if(i != modeIndex) {
          if(quadTransform_ != "identity")
            nodes[i,] <<- transformNode(z = nodes[i,], eigenvec = V, eigenval = L, method = quadTransform_)
          logDensity_quad[i] <<- logLik_RE(reTransform = nodes[i,]) ## Was joint_logLik, but p is constant?
          ans <- ans + exp(logDensity_quad[i] - saved_inner_max_value)*wgts[i]
        }else{
          if(quadTransform_ != "identity")
            nodes[i,] <<- saved_inner_argmax
          logDensity_quad[i] <<- saved_inner_max_value
          ans <- ans + wgts[i]
        }
      }
      ## Given all the saved values, weights and log density, do quadrature sum.
      res <- log(ans) + saved_inner_max_value - 0.5 * saved_inner_logdetNegHess
      return(res)
      returnType(double())
    },
    ## Gradient of the Laplace approximation w.r.t. parameters
    gr_logLik2 = function(p = double(1)){
      if(!one_time_fixes_done) one_time_fixes()
      if(any(p != saved_inner_max_p) | !cache_inner_max) {
        update_max_logLik_RE(p)
      }
      reTransform <- saved_inner_argmax
      invNegHess <- inverse(saved_inner_negHess) # want to change and do this by solves with cholesky
      gr_LL <- gr_P_RE_b(p, reTransform)
      uptri_Omega_invNegHess <- nimNumeric(length = nreTrans*(nreTrans+1)/2, fillZeros=FALSE, init=FALSE)
      i <- 1L; j <- 1L
      uptri_Omega_invNegHess[1] <- invNegHess[1, 1]
      k <- 2L
      for(j in 2:nreTrans) {
        for(i in 1:(j-1)) {
          uptri_Omega_invNegHess[k] <- 2*invNegHess[i, j]
          k <- k+1
        }
        uptri_Omega_invNegHess[k] <- invNegHess[j, j]
        k <- k+1
      }
      # N.B. An extra negation is built into gr_logdet because this is gradient of hessian, but the uptri_Omega_invNegHess is from the negative Hessian.
      gr_logdet <- jac_he_P_RE_wrt_RE2_uptri_outDir_c(p, reTransform, uptri_Omega_invNegHess)
      SigmaInv_gr_logdet <- (invNegHess %*% gr_logdet[1, reTrans_indices])[,1] # want to change to do this using solve of cholesky
      lastpiece <- jac_gr_P_RE_wrt_RE_outDir_b(p, reTransform,
                                              outDir=SigmaInv_gr_logdet) # return as matrix
      res <- gr_LL[p_indices] + 0.5*gr_logdet[1, p_indices] + 0.5*lastpiece[1, p_indices]
      return(res)
      returnType(double(1))
    },
    get_inner_mode = function(atOuterMode = integer(0, default = 0)){
      returnType(double(1))
      if(atOuterMode) return(outer_mode_inner_argmax)
      return(saved_inner_argmax)
    },
    get_inner_negHessian = function(atOuterMode = integer(0, default = 0)){
      returnType(double(2))
      if(atOuterMode) return(outer_mode_inner_negHess)
      return(saved_inner_negHess)
    },
    get_inner_negHessian_chol = function(atOuterMode = integer(0, default = 0)){
      returnType(double(2))
      if(atOuterMode) return(outer_mode_inner_negHess_chol)
      return(saved_inner_negHess_chol)
    },
    ## Update the maximum mode and neg hess based on the log likelihood passed via optim.
    ## For efficient saving of values for calculating MLE values of random-effects.
    save_outer_logLik = function(logLikVal = double()){
      if(logLikVal >= max_outer_logLik) {
        max_outer_logLik <<- logLikVal
        outer_mode_inner_negHess <<- saved_inner_negHess
        outer_mode_inner_argmax <<- saved_inner_argmax
        outer_mode_inner_negHess_chol <<- saved_inner_negHess_chol
        outer_param_max <<- saved_inner_max_p
      }
    },
    get_param_value = function(atOuterMode = integer(0, default = 0)){
      returnType(double(1))
      ## Ensures that the inner value will not match and cached values will not be used.
      if(!cache_inner_max) return(numeric(value = Inf, length = npar))
      if(atOuterMode) return(outer_param_max)
      return(saved_inner_max_p)
    },
    ## Need to reset every call optim to recache.
    reset_outer_logLik = function(){
      max_outer_logLik <<- -Inf
    },
    set_randomeffect_values = function(p = double(1)){
      foundIt <- FALSE
      ## Last value called:
      if(all(p == saved_inner_max_p)) {
        re <- reTrans$inverseTransform(saved_inner_argmax)
        foundIt <- TRUE
      }
      ## Best value called:
      if(all(p == outer_param_max)) {
        re <- reTrans$inverseTransform(outer_mode_inner_argmax)
        foundIt <- TRUE
      }
      if(foundIt){
        values(model, paramNodes) <<- p
        model$calculate(paramDeps)
      }else{
        # It would be nice to emit a message here, but different optimizers (e.g. BFGS vs nlminb)
        # behave differently as to whether the previous (last) parameters were always the MLE.
        # print("  [Warning] Have not cached the inner optimization. Running optimization now.")
        update_max_logLik_RE(p)
        re <- reTrans$inverseTransform(saved_inner_argmax)
      }
      ## Ensure the model is up to date for all nodes.
      values(model, randomEffectsNodes) <<- re
      model$calculate(innerCalcNodes)
    }
    ## set_inner_cache = function(cache = logical(0, default = TRUE)){
    ##   cache_inner_max <<- cache
    ## }
  ),
  buildDerivs = list(logLik_RE              = list(),
                     gr_RE_a                = list(),
                     he_RE_b                = list(),
                     he_RE_b_asvec          = list(),
                     logLik_P_RE            = list(),
                     gr_P_RE_a              = list(),
                     gr_P_RE_wrt_RE_a      = list(),
                     he_P_RE_wrt_RE2_uptri_b = list())
) ## End of buildOneAGHQuad


## Main function for Laplace approximation
#' @rdname laplace
#' @export
buildLaplace <- function(model, paramNodes, randomEffectsNodes, calcNodes, calcNodesOther, control = list()) {
  buildAGHQ(model, nQuad = 1, paramNodes, randomEffectsNodes, calcNodes, calcNodesOther, control)
}

## Main function for Adaptive Gauss-Hermite Quadrature
#' @rdname laplace
#' @export
buildAGHQ <- nimbleFunction(
  name = 'AGHQ',
  setup = function(model, nQuad = 1, paramNodes, randomEffectsNodes, calcNodes,
                   calcNodesOther, control = list()) {
    split <- extractControlElement(control, 'split', TRUE)
    check <- extractControlElement(control, 'check', TRUE)
    innerOptimWarning <- extractControlElement(control, 'innerOptimWarning', FALSE)

    if(!is.Rmodel(model))
        stop("`model` must be an R model, created by calling `nimbleModel`")

    if(nQuad %% 2 == 0)
      messageIfVerbose("  [Note] For computational efficiency, it is recommended to use an odd number of quadrature points.")
    if(nQuad > 35) {
      messageIfVerbose("  [Note] Currently only a maximum of 35 quadrature points are allowed; setting nQuad to 35.")
      nQuad <- 35
    }
    nQuad_ <- nQuad
    MargNodes <- NULL
    if(!missing(paramNodes)) {
      if(is.list(paramNodes)) {
        # The user called setupMargNodes and provided a list of that format to paramNodes.
        MargNodes <- paramNodes
      }
    }
    if(is.null(MargNodes)) {
      MargNodes <- setupMargNodes(model = model, paramNodes = paramNodes,
                                  randomEffectsNodes = randomEffectsNodes,
                                  calcNodes = calcNodes,
                                  calcNodesOther = calcNodesOther,
                                  split = split,
                                  check = check)
    }
    paramNodes <- MargNodes$paramNodes
    randomEffectsNodes <- MargNodes$randomEffectsNodes
    calcNodes <- MargNodes$calcNodes
    calcNodesOther <- MargNodes$calcNodesOther
    num_calcNodesOther <- length(calcNodesOther)
    # MargNodes$randomEffectsSets will be extracted below if needed

    if(length(calcNodesOther)) {
      otherLogLik_derivsInfo    <- makeModelDerivsInfo(model = model, wrtNodes = paramNodes, calcNodes = calcNodesOther)
      otherLogLik_updateNodes   <- otherLogLik_derivsInfo$updateNodes
      otherLogLik_constantNodes <- otherLogLik_derivsInfo$constantNodes
    }
    else { ## calcNodesOther is empty
      otherLogLik_updateNodes   <- character(0)
      otherLogLik_constantNodes <- character(0)
    }

    calcPrior_derivsInfo <- makeModelDerivsInfo(model, paramNodes, paramNodes)
    calcPrior_updateNodes   <- calcPrior_derivsInfo$updateNodes
    calcPrior_constantNodes <- calcPrior_derivsInfo$constantNodes

    ## Out and inner optimization settings
    outerOptimControl_   <- nimOptimDefaultControl()
    innerOptimControl_ <- nimOptimDefaultControl()
    optimControlArgNames <- c("trace", "fnscale", "parscale", "ndeps", "maxit", "abstol", "reltol", "alpha",
                              "beta", "gamma", "REPORT", "type", "lmm", "factr", "pgtol", "temp", "tmax")
    if(!is.null(control$outerOptimControl)){
      validNames <- intersect(names(control$outerOptimControl), optimControlArgNames)
      numValidNames <- length(validNames)
      if(numValidNames > 0){
        for(i in 1:numValidNames){
          outerOptimControl_[[validNames[i]]] <- control$outerOptimControl[[validNames[i]]]
        }
      }
    }
    if(!is.null(control$innerOptimControl)) {
      validNames_inner <- intersect(names(control$innerOptimControl), optimControlArgNames)
      numValidNames_inner <- length(validNames_inner)
      if(numValidNames_inner > 0){
        for(i in 1:numValidNames_inner)
          innerOptimControl_[[validNames_inner[i]]] <- control$innerOptimControl[[validNames_inner[i]]]
      }
    }
    outerOptimControl_$fnscale <- -1
    innerOptimControl_$fnscale <- -1
    if(!is.null(control$innerOptimMethod) &&
       ((control$innerOptimMethod %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B")) ||
        control$innerOptimMethod %in% ls(nimbleUserNamespace$.optimizers))){  # .optimizers by default contains 'nlminb'.
      innerOptimMethod <- control$innerOptimMethod
    } else innerOptimMethod <- "nlminb"

    if(!is.null(control$outerOptimMethod) &&
       ((control$outerOptimMethod %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B")) ||
        control$outerOptimMethod %in% ls(nimbleUserNamespace$.optimizers))){  # .optimizers by default contains 'nlminb'.
      outerOptimMethod_ <- control$outerOptimMethod
    } else outerOptimMethod_ <- "nlminb"

    innerOptimStart <- extractControlElement(control, "innerOptimStart", "last.best")
    if(!is.character(innerOptimStart) |
       length(innerOptimStart) != 1 |
       !(innerOptimStart %in% (validIOS <- c("last", "last.best", "constant", "random", "model", "zero"))))
      stop("buildAGHQ: `control$innerOptimStart` must be one of ", paste0('\'', validIOS, '\'', collapse=","))

    innerOptimStartValues <- NULL
    if(innerOptimStart == "model") {
      innerOptimStartValues <- 0 # will be ignored but is need to trigger next message
    }
    if(innerOptimStart == "zero") {
      innerOptimStart <- "constant"
      innerOptimStartValues <- 0
    }
    if(!is.null(innerOptimStartValues) &
       !is.null(control$innerOptimStartValues)) {
        messageIfVerbose("  [Note] Ignoring `control$innerOptimStartValues` because `control$innerOptimStart` is `", innerOptimStart, "`.")
    } else {
        innerOptimStartValues <- extractControlElement(control, "innerOptimStartValues", 0)
        if(is.character(innerOptimStartValues))
            if(length(innerOptimStartValues) != 1 |
               !(innerOptimStartValues == "model"))
                stop("buildAGHQ: The only valid character value for `control$innerOptimStartValues` is 'model'")
    }

    outerOptimUseAD <- extractControlElement(control, "outerOptimUseAD", TRUE)

    ## Create an AGHQuad (Adaptive Gauss-Hermite Quadrature) nimbleFunctionList
    AGHQuad_nfl <- nimbleFunctionList(AGHQuad_BASE)
    scalarRENodes <- model$expandNodeNames(randomEffectsNodes, returnScalarComponents = TRUE)
    nre <- length(scalarRENodes)
    multiSetsCheck <- FALSE ## AGHQ vs Laplace Check in findMLE.
    quadTransform <- extractControlElement(control, "quadTransform", "cholesky")
    innerControlList <- list(optimControl=innerOptimControl_,
                             optimMethod=innerOptimMethod,
                             optimStart=innerOptimStart,
                             optimStartValues=innerOptimStartValues,
                             optimWarning=innerOptimWarning,
                             quadTransform=quadTransform)
    if(nre > 0){
      ## Record the order of random effects processed internally
      internalRandomEffectsNodes <- NULL
      ## lenInternalRENodeSets <- NULL
      if(isFALSE(split)) { ## Do all randomEffectsNodes in one set
        internalRandomEffectsNodes <- randomEffectsNodes
        ## lenInternalRENodeSets <- nre
        reNodesAsScalars <- model$expandNodeNames(internalRandomEffectsNodes, returnScalarComponents = TRUE)
        # old default was "model". new default is "last.best" with values=0
        # Essentially the following steps will now b done in buildOneAGHQuad[1D]
        ## if(is.null(control$innerOptimStart)) innerOptimStart <- values(model, randomEffectsNodes)
        ## else {
        ##   providedStart <- control$innerOptimStart
        ##   if(any(providedStart %in% c("last", "last.best"))) innerOptimStart <- providedStart
        ##   else if(is.numeric(sum(providedStart)) && (length(providedStart) == nre)) innerOptimStart <- providedStart
        ##   else innerOptimStart <- values(model, randomEffectsNodes)
        ## }
        ## ## In case random effects are not properly initialized
        ## if(!any(innerOptimStart %in% c("last", "last.best")) & any(is.infinite(innerOptimStart) | is.na(innerOptimStart) | is.nan(innerOptimStart))){
        ##   all_reTransform <- parameterTransform(model, randomEffectsNodes)
        ##   all_nreTrans <- all_reTransform$getTransformedLength()
        ##   innerOptimStart <- all_reTransform$inverseTransform(rep(0, all_nreTrans))
        ## }
        ## Build AGHQuad
        if(nre > 1 | isTRUE(control[['force_nDim']])) {
          AGHQuad_nfl[[1]] <- buildOneAGHQuad(model, nQuad = nQuad_, paramNodes, randomEffectsNodes,
                                              calcNodes, control = innerControlList)
          multiSetsCheck <- TRUE
        } else {
          AGHQuad_nfl[[1]] <- buildOneAGHQuad1D(model, nQuad = nQuad_, paramNodes, randomEffectsNodes,
                                                calcNodes, control = innerControlList)
        }
        num_reSets <- 1
      }
      else {## Split randomEffectsNodes into conditionally independent sets
        reSets <- MargNodes$randomEffectsSets
        num_reSets <- length(reSets)
        reNodesAsScalars <- character()
        if(num_reSets == 0){
          stop("buildAGHQ: There was a problem determining conditionally independent random effects sets for this model")
        }
        if(nQuad_ == 1)
            msg <- "Laplace" else msg <- "AGHQ (extended Laplace)"
        if(length(reSets) > 1) {
            messageIfVerbose("Building ", num_reSets, " individual ", msg, " approximations (one dot for each): ", appendLF = FALSE)
        } else {
          messageIfVerbose("Building ", msg, " approximation.")
        }
        ## Do this once as shared across all individual AGHQs.
        paramDeps <- model$getDependencies(paramNodes, determOnly = TRUE, self=FALSE)
        for(i in seq_along(reSets)){
          ## Work with one conditionally independent set of latent states
          these_reNodes <- reSets[[i]]
          internalRandomEffectsNodes <- c(internalRandomEffectsNodes, these_reNodes)
          ## find paramNodes and calcNodes for this set of reNodes
          ## paramNodes are the same for all AGHQuad_nfl elements. In the future this could be customized.
          these_reDeps <- model$getDependencies(these_reNodes)  ## candidate calcNodes via reNodes
          these_calcNodes <- intersect(calcNodes, these_reDeps) ## definite calcNodes
          these_reNodesAsScalars <- model$expandNodeNames(these_reNodes, returnScalarComponents = TRUE)
          reNodesAsScalars <- c(reNodesAsScalars, these_reNodesAsScalars)
          nre_these <- length(these_reNodesAsScalars)
          ## lenInternalRENodeSets <- c(lenInternalRENodeSets, nre_these)
          ## Process start values for inner optimisation
          if(is.numeric(innerOptimStartValues) && (length(innerOptimStartValues) == nre)) {
            # input was a vector of all REs, so we must split it accordingly
            if(is.null(names(innerOptimStartValues))) {
              # split by order, because there are no names. use internalRandomEffectsNodes order.
              # the last portion will be for the present set of values
              these_reNodes_inds <- length(reNodesAsScalars) - nre_these + (1:nre_these)
            } else {
              # split by names
              these_reNodes_inds <- match(these_reNodesAsScalars, innerOptimStartValues, nomatch=0)
              if((length(these_reNodes_inds) != nre_these) |
                   (length(unique(these_reNodes_inds)) != length(these_reNodes_inds)) |
                   any(these_reNodes_inds==0))
                messageIfVerbose("  [Warning] There appears to be an incorrect name in `control$innerOptimStartValues`.")
            }
            these_innerOptimStartValues <- innerOptimStartValues[these_reNodes_inds]
          } else {
            these_innerOptimStartValues <- innerOptimStartValues
          }
          ## if(is.null(control$innerOptimStart)) innerOptimStart <- values(model, these_reNodes)
          ## else {
          ##   providedStart <- control$innerOptimStart
          ##   if(any(providedStart %in% c("last", "last.best"))) innerOptimStart <- providedStart
          ##   else if(is.numeric(sum(providedStart)) && (length(providedStart) == nre)){
          ##     these_reNodes_inds <- unlist(lapply(model$expandNodeNames(these_reNodes, returnScalarComponents = TRUE), function(x) {which(scalarRENodes == x)}))
          ##     innerOptimStart <- providedStart[these_reNodes_inds]
          ##   }
          ##   else innerOptimStart <- values(model, these_reNodes)
          ## }
          ## ## In case random effects are not properly initialized
          ## if(!any(innerOptimStart %in% c("last", "last.best")) & any(is.infinite(innerOptimStart) | is.na(innerOptimStart) | is.nan(innerOptimStart))){
          ##   these_reTransform <- parameterTransform(model, these_reNodes)
          ##   these_nreTrans <- these_reTransform$getTransformedLength()
          ##   innerOptimStart <- these_reTransform$inverseTransform(rep(0, these_nreTrans))
          ## }
          ## Build AGHQuad for each set
          ## **** WZ below the innerControlList is still the full one I do not know why the code above was commented out
          if(nre_these > 1 | isTRUE(control[['force_nDim']])){
            AGHQuad_nfl[[i]] <- buildOneAGHQuad(model, nQuad = nQuad_, paramNodes, these_reNodes, these_calcNodes,
                                                paramDeps, innerControlList)
            multiSetsCheck <- TRUE
          }
          else {
            AGHQuad_nfl[[i]] <- buildOneAGHQuad1D(model, nQuad = nQuad_, paramNodes, these_reNodes, these_calcNodes,
                                                  paramDeps, innerControlList)
          }
          if(length(reSets) > 1) messageIfVerbose(".", appendLF = FALSE)
        }
        if(length(reSets) > 1) messageIfVerbose("")
      }

      ## if(length(lenInternalRENodeSets) == 1) lenInternalRENodeSets <- c(lenInternalRENodeSets, -1)
      reTransform <- parameterTransform(model, internalRandomEffectsNodes)
      nreTrans <- reTransform$getTransformedLength()
      if(nreTrans > 1) reTransform_indices <- 1:nreTrans
      else reTransform_indices <- c(1, -1)

      reNodesAsScalars_vec <- reNodesAsScalars
      if(nre == 1) reNodesAsScalars_vec <- c(reNodesAsScalars, "_EXTRA_")
      reNodesAsScalars_first <- reNodesAsScalars[1]
      ## When transformed, put "_trans" to the end of the random effects names
      reTransNodeNames <- paste0("re_trans_", 1:nreTrans)
      if(nreTrans == 1) reTransNodeNames <- c(reTransNodeNames, "_EXTRA_")
    }
    else{
      ## No random effects
      ## lenInternalRENodeSets <- numeric(2)
      num_reSets <- 0
      nreTrans <- 0
      reTransform <- parameterTransform(model, paramNodes[1], control = list(allowDeterm = FALSE)) ## Won't be needed at all
      reTransform_indices <- numeric(2)
      reNodesAsScalars_vec <- character(0)
      reNodesAsScalars_first <- character(1)
      reTransNodeNames <- character(2)
      if(num_calcNodesOther == 0)
        stop("buildAGHQ: Both `calcNodesOther` and `randomEffectsNodes` are empty for Laplace or AGHQ for the given model")
    }

    paramNodesAsScalars <- model$expandNodeNames(paramNodes, returnScalarComponents = TRUE)
    npar <- length(paramNodesAsScalars)
    paramNodesAsScalars_vec <- paramNodesAsScalars
    if(npar == 1) paramNodesAsScalars_vec <- c(paramNodesAsScalars, "_EXTRA_")
    paramNodesAsScalars_first <- paramNodesAsScalars[1]
    if(npar == 1) p_indices <- c(1, -1)
    else p_indices <- 1:npar
    ## setupOutputs(reNodesAsScalars, paramNodesAsScalars)

    ## Automated transformation for parameters
    paramsTransform <- parameterTransform(model, paramNodes, control = list(allowDeterm = FALSE))
    nparTrans <- paramsTransform$getTransformedLength()
    if(nparTrans > 1) pTransform_indices <- 1:nparTrans
    else pTransform_indices <- c(1, -1)

    paramTransNodeNames <- paste0("param_trans_", seq_len(nparTrans))
    if(nparTrans == 1) paramTransNodeNames <- c(paramTransNodeNames, "_EXTRA_")

    ## Indicator for removing the redundant index -1 in pTransform_indices
    one_time_fixes_done <- FALSE
    ## Default calculation method for AGHQuad
    computeMethod_ <- extractControlElement(control, "computeMethod", 2)

    useInnerCache_ <- extractControlElement(control, "useInnerCache", TRUE)

    ## Set cached values for calculating prior and posterior in log density.
    includePrior_ <- TRUE
    includeJacobian_ <- TRUE

    ## Set up cached values for doing profile likelihood construction:
    pTransform_fixed <- 0
    pTransform_index_fixed <- 1
    pTransform_indices_other <- numeric(2)

    summaryCalcMethod = 2
    ## The nimbleList definitions AGHQuad_params and AGHQuad_summary
    ## have moved to predefined nimbleLists.
  },## End of setup
  run = function(){},
  methods = list(
    getREtransLength = function() {
      numre <- numeric(num_reSets)
      for(i in seq_along(AGHQuad_nfl))
        numre[i] <- AGHQuad_nfl[[i]]$get_reTransLength()
      returnType(double(1))
      return(numre)
    },
    getNodeNamesVec = function(returnParams = logical(0, default = TRUE)) {
      one_time_fixes()
      returnType(character(1))
      if(returnParams) return(paramNodesAsScalars_vec)
      else return(reNodesAsScalars_vec)
    },
    getNodeNameSingle = function(returnParams = logical(0, default = TRUE)) {
      returnType(character())
      if(returnParams) return(paramNodesAsScalars_first)
      else return(reNodesAsScalars_first)
    },
    updateSettings = function(innerOptimMethod = character(0, default="NULL"),
                              innerOptimStart = character(0, default="NULL"),
                              innerOptimStartValues = double(1, default=Inf),
                              innerOptimWarning = integer(0, default = -1),
                              useInnerCache = integer(0, default=-1),
                              nQuad = integer(0, default=-1),
                              quadTransform = character(0, default="NULL"),
                              innerOptimControl = optimControlNimbleList(default=nimOptimDefaultControl()),
                              outerOptimMethod = character(0, default="NULL"),
                              replace_innerOptimControl = logical(0, default=FALSE),
                              outerOptimControl = optimControlNimbleList(default=nimOptimDefaultControl()),
                              replace_outerOptimControl = logical(0, default=FALSE),
                              computeMethod = integer(0, default=-1)
                              ) {
      # checks
      if(innerOptimStart != "NULL") {
        if(innerOptimStart=="zero") {
          stop("updateSettings: `innerOptimStart` choice of 'zero' is not supported in `updateSettings`. Use `innerOptimStart='constant'` and `innerOptimStartValues = 0` to achieve 'zero' behavior")
        }
        if(innerOptimStart != "last" & innerOptimStart != "last.best" &
           innerOptimStart != "constant" & innerOptimStart != "random" &
           innerOptimStart != "model")
            stop("updateSettings: invalid value for `innerOptimStart`")
      }
      if(length(innerOptimStartValues) > 1) {
        if(length(innerOptimStartValues) != nre)
          stop("updateSettings: length of `innerOptimStartValues` must be 1 or total number of random effects")
      }
      if(nQuad != -1)  {
        if(nQuad < 1) stop("updateSettings: choose a positive number of grid points")
        if(nQuad > 35) stop("updateSettings: currently only a maximum of 35 quadrature points is allowed")
        threshold <- log(50000) # in text below too
        for(i in seq_along(AGHQuad_nfl)) {
          numre <- AGHQuad_nfl[[i]]$get_reTransLength()
          if(nQuad * log(numre) > threshold) {
              print("updateSettings: choice of `nQuad` would yield >50000 nodes for ", numre, " integration dimensions in conditionally independent set ", i, ".")
              stop("too many integration nodes")
          }
        }
      }
      if(computeMethod != -1) {
        if(!any(c(1, 2, 3) == computeMethod))  ## Cannot use `%in%` in nf code.
          stop("updateSettings: `computeMethod` must be 1, 2, or 3")
      }
      if(quadTransform != "NULL") {
        if(quadTransform != "spectral" & quadTransform != "cholesky")
          stop("`quadTransform` must be either cholesky or spectral.")
      }
      # actions
      one_time_fixes()
      if(nQuad != -1) nQuad_ <<- nQuad
      these_startValues <- innerOptimStartValues
      iStart <- 1
      for(i in seq_along(AGHQuad_nfl)) {
        if(length(innerOptimStartValues) > 1) {
          numre <- AGHQuad_nfl[[i]]$get_reTransLength()
          these_startValues <- innerOptimStartValues[iStart:(iStart + numre - 1)]
          iStart <- iStart + numre
        }
        AGHQuad_nfl[[i]]$updateSettings(optimMethod = innerOptimMethod,
                                        optimStart = innerOptimStart,
                                        optimStartValues = these_startValues,
                                        optimWarning = innerOptimWarning,
                                        useInnerCache = useInnerCache,
                                        nQuad = nQuad_,
                                        quadTransform = quadTransform,
                                        optimControl = innerOptimControl,
                                        replace_optimControl = replace_innerOptimControl)
      }
      # TO-DO: create useInnerCache_ and allow control arg.
      if(useInnerCache != -1) useInnerCache_ <<- useInnerCache != 0
      if(computeMethod != -1) computeMethod_ <<- computeMethod
      if(replace_outerOptimControl) {
        if(outerOptimControl$fnscale == 1) outerOptimControl$fnscale <- -1
        outerOptimControl_ <<- outerOptimControl
      }
      if(outerOptimMethod != "NULL")
        outerOptimMethod_ <<- outerOptimMethod
    },
    one_time_fixes = function() {
      if(one_time_fixes_done) return()
      if(nparTrans == 1){
        if(length(pTransform_indices) == 2){
          pTransform_indices <<- numeric(length = 1, value = 1)
        }
      }
      if(npar == 1){
        if(length(p_indices) == 2){
          p_indices <<- numeric(length = 1, value = 1)
        }
      }
      one_time_fixes_done <<- TRUE
    },
    reset = function(gr_RE = logical(0, default = TRUE),
                     he_RE = logical(0, default = TRUE),
                     gr_P_RE = logical(0, default = TRUE),
                     gr_P_RE_wrt_RE = logical(0, default = TRUE),
                     he_P_RE_wrt_RE2_uptri = logical(0, default = TRUE)){
      for(i in seq_along(AGHQuad_nfl)){
          AGHQuad_nfl[[i]]$reset(gr_RE = gr_RE, he_RE = he_RE,
                                gr_P_RE = gr_P_RE, gr_P_RE_wrt_RE = gr_P_RE_wrt_RE,
                                he_P_RE_wrt_RE2_uptri = he_P_RE_wrt_RE2_uptri)
        }
    },
    ## Check to see if the inner optimizations converged.
    checkInnerConvergence = function(message = logical(0, default = FALSE)){
      converged <- 0
      for(i in seq_along(AGHQuad_nfl)){
        conCheck <- AGHQuad_nfl[[i]]$check_convergence()
        if(conCheck != 0) {
          converged <- 1
          if(message) print("  [Warning] Inner optimization did not converge for conditionally independent set ", i, " with code ", conCheck, ".")
        }
      }
      returnType(double())
      return(converged)
    },
    ## Other log-likelihood (parts not involving random effects, i.e. simply
    ## additional calculations in the model) in terms of original parameters
    otherLogLik = function(p = double(1)) {
      if(num_calcNodesOther == 0) stop("`calcNodesOther` is empty: there is no exact likelihood component for the model")
      values(model, paramNodes) <<- p
      ans <- model$calculate(calcNodesOther)
      return(ans)
      returnType(double())
    },
    ## Gradient of the exact log-likelihood w.r.t parameters
    gr_otherLogLik_internal = function(p = double(1)) {
      if(num_calcNodesOther == 0) stop("`calcNodesOther` is empty: there is no exact likelihood component for the model")
      if(!one_time_fixes_done) one_time_fixes()
      ans <- derivs(otherLogLik(p), wrt = p_indices, order = 1, model = model,
                    updateNodes = otherLogLik_updateNodes, constantNodes = otherLogLik_constantNodes)
      return(ans$jacobian[1,])
      returnType(double(1))
    },
    ## Double taping for efficiency
    gr_otherLogLik = function(p = double(1)) {
      if(num_calcNodesOther == 0) stop("`calcNodesOther` is empty: there is no exact likelihood component for the model")
      if(!one_time_fixes_done) one_time_fixes()
      ans <- derivs(gr_otherLogLik_internal(p), wrt = p_indices, order = 0, model = model,
                    updateNodes = otherLogLik_updateNodes, constantNodes = otherLogLik_constantNodes)
      return(ans$value)
      returnType(double(1))
    },
    ## AGHQuad approximation in terms of original parameters
    calcLogLik = function(p = double(1), trans = logical(0, default = FALSE)) {
      if(!one_time_fixes_done) one_time_fixes()
      checkInterrupt()
      if(trans) {
          if(length(p) != nparTrans) {
              ## We cannot have variables in a nimStop.
              print("For `calcLogLik` (or `calcLaplace`) with `trans = TRUE`, `p` should be length ", nparTrans, " but was provided with length ", length(p), ".")
              stop("incorrect length for `p`")
          }
          p <- paramsTransform$inverseTransform(p)
      }
      if(length(p) != npar) {
          print("For `calcLogLik` (or `calcLaplace`), `p` should be length ", npar, " but is length ", length(p), ".")
          stop("incorrect length for `p`")
      }
      if(num_calcNodesOther > 0) ans <- otherLogLik(p)
      else ans <- 0
      if(nre > 0){
        for(i in seq_along(AGHQuad_nfl)){
          # Currently computeMethod_ is ignored. We retain it as an option for future use.
          ans <- ans + AGHQuad_nfl[[i]]$calcLogLik2(p) # Everything is method "2" at inner level.
        }
      }
      if(is.nan(ans) | is.na(ans)) ans <- -Inf
      return(ans)
      returnType(double())
    },
    calcLaplace = function(p = double(1), trans = logical(0, default = FALSE)) {
      if(nQuad_ > 1) {
        stop("`nQuad` must be equal to 1 to use `calcLaplace`. Either call `calcLogLik` or use `updateSettings()` to change `nQuad`")
      }
      ans <- calcLogLik(p, trans)
      return(ans)
      returnType(double())
    },
    ## Gradient of the AGHQuad approximation w.r.t. parameters
    gr_logLik = function(p = double(1), trans = logical(0, default=FALSE)) {
      if(!one_time_fixes_done) one_time_fixes()
      if(trans) {
        if(length(p) != nparTrans) {
            print("for `gr_logLik` (or `gr_Laplace`) with `trans = TRUE`, `p` should be length ", nparTrans, " but was provided with length ", length(p), ".")
            stop("incorrect length for `p`")
        }
        pDerivs <- derivs_pInverseTransform(p, c(0, 1))
        p <- pDerivs$value
      }
      if(length(p) != npar) {
          print("for `gr_logLik` (or `gr_Laplace`), `p` should be length ", npar, " but is length ", length(p), ".")
          stop("incorrect length for `p`")
      }
      if(num_calcNodesOther > 0) ans <- gr_otherLogLik(p) else ans <- numeric(length = npar)
      if(nre > 0){
        for(i in seq_along(AGHQuad_nfl)) {
          # Everything is computeMethod "2" and computeMethod_ is ignored.
          ans <- ans + AGHQuad_nfl[[i]]$gr_logLik2(p)
        }
      }
      if(trans) {
        ans <- (ans %*% pDerivs$jacobian)[1,]
      }
      return(ans)
      returnType(double(1))
    },
    gr_Laplace = function(p = double(1), trans = logical(0, default=FALSE)) {
      if(nQuad_ > 1)
        stop("`nQuad` must be equal to 1 to use `calcLaplace`. Either call `calcLogLik` or use `updateSettings()` to change `nQuad`")
      ans <- gr_logLik(p, trans)
      return(ans)
      returnType(double(1))
    },
    ## AGHQuad approximation in terms of transformed parameters
    calcLogLik_pTransformed = function(pTransform = double(1)) {
      ans <- calcLogLik(pTransform, trans = TRUE)
      ## if(!one_time_fixes_done) one_time_fixes()
      ## p <- paramsTransform$inverseTransform(pTransform)
      ## ans <- calcLogLik(p)
      ## if(is.nan(ans) | is.na(ans)) ans <- -Inf
      cache_outer_logLik(ans) ## Save outer in the inner to cache values at outer mode.
      return(ans)
      returnType(double())
    },
    pTransform_internal = function(p = double(1)){
      returnType(double(1))
      return(paramsTransform$transform(p))
    },
    ## Inverse transform parameters to original scale
    pInverseTransform = function(pTransform = double(1)) {
      p <- paramsTransform$inverseTransform(pTransform)
      return(p)
      returnType(double(1))
    },
    ## Jacobian of the inverse transformation for parameters
    derivs_pInverseTransform = function(pTransform = double(1), order = double(1)) {
      if(!one_time_fixes_done) one_time_fixes()
      ans <- derivs(pInverseTransform(pTransform), wrt = pTransform_indices, order = order)
      return(ans)
      returnType(ADNimbleList())
    },
    ## Inverse transform random effects to original scale
    reInverseTransform = function(reTrans = double(1)) {
      if(nre == 0) stop("no random effects in the model")
      re <- reTransform$inverseTransform(reTrans)
      return(re)
      returnType(double(1))
    },
    ## Jacobian of the inverse transformation
    derivs_reInverseTransform = function(reTrans = double(1), order = double(1)) {
      if(!one_time_fixes_done) one_time_fixes()
      if(nre == 0) stop("no random effects in the model")
      ans <- derivs(reInverseTransform(reTrans), wrt = reTransform_indices, order = order)
      return(ans)
      returnType(ADNimbleList())
    },
    ## Gradient of the AGHQuad approximation in terms of transformed parameters
    gr_logLik_pTransformed = function(pTransform = double(1)) {
      ans <- gr_logLik(pTransform, trans = TRUE)
      ## if(!one_time_fixes_done) one_time_fixes()
      ## pDerivs <- derivs_pInverseTransform(pTransform, c(0, 1))
      ## gr <- gr_logLik(pDerivs$value) ## pDerivs$value gives original param values
      ## ans <- (gr %*% pDerivs$jacobian)[1,]
      return(ans)
      returnType(double(1))
    },
    ## Prior contribution to the posterior
    calcPrior_p = function(p = double(1)){
      ## Prior log likelihood:
      values(model, paramNodes) <<- p
      ans <- model$calculate(paramNodes)
      return(ans)
      returnType(double())
    },
    ## Prior contribution to the posterior on the transformed scale.
    calcPrior_pTransformed = function(pTransform = double(1)) {
      p <- paramsTransform$inverseTransform(pTransform)
      ans <- calcPrior_p(p) + logDetJacobian(pTransform)
      return(ans)
      returnType(double())
    },
    ## Calculate posterior density at p log likelihood + log prior.
    calcLogDens = function(p = double(1), trans = logical(0, default = FALSE),
                           includeJacobian = logical(0, default = TRUE),
                           includePrior = logical(0, default = TRUE)) {
      ans <- 0
      if(trans) {
        pstar <- paramsTransform$inverseTransform(p)  ## Just want to do this once.
        if(includeJacobian) ans <- ans + logDetJacobian(p)  ## p is transformed, add Jacobian here.
      }else{
        pstar <- p
      }
      ans <- ans + calcLogLik(pstar, FALSE)
      if(includePrior) ans <- ans + calcPrior_p(pstar)

      returnType(double())
      return(ans)
    },
    ## Calculate posterior density at p transformed, log likelihood + log prior (transformed).
    calcLogDens_pTransformed = function(pTransform = double(1)) {
      ans <- calcLogDens(pTransform, trans = TRUE,
                         includeJacobian = includeJacobian_,
                         includePrior = includePrior_)
      cache_outer_logLik(ans) ## Update internal cache w/ prior.

      if(is.nan(ans) | is.na(ans)) ans <- -Inf
      returnType(double())
			return(ans)
    },
    calcLogDens_pTransformedFix1 = function(pTransform = double(1)){
      pTransform_star <- replaceOneVec(pTransform)
      ans <- calcLogDens(pTransform_star, trans = TRUE,
                         includeJacobian = includeJacobian_,
                         includePrior = includePrior_)
      cache_outer_logLik(ans) ## Update internal cache w/ prior.

      if(is.nan(ans) | is.na(ans)) ans <- -Inf
      returnType(double())
			return(ans)
    },
    ## Gradient of log det jacobian for parameter transformations.
    gr_logDetJacobian = function(pTransform = double(1)){
      ans <- derivs(logDetJacobian(pTransform), wrt = pTransform_indices, order = 1)
      return(ans$jacobian[1,])
      returnType(double(1))
    },
    ## Gradient of prior distribution.
    gr_prior = function(p = double(1)){
      ans <- derivs(calcPrior_p(p), wrt = p_indices, order = 1, model = model,
                    updateNodes = calcPrior_updateNodes, constantNodes = calcPrior_constantNodes)
      return(ans$jacobian[1,])
      returnType(double(1))
    },
    ## Gradient of posterior density on the transformed scale.
    gr_logDens = function(p = double(1), trans = logical(0, default = FALSE),
                          includeJacobian = logical(0, default = TRUE),
                          includePrior = logical(0, default = TRUE)){
      if(trans) {
        pDerivs <- derivs_pInverseTransform(p, c(0, 1))
        pstar <- pDerivs$value
      }else {
        pstar <- p
      }
      ## Gradient of log likelihood:
      ans <- gr_logLik(pstar, FALSE)
      if(includePrior) ans <- ans + gr_prior(pstar)

      if(trans){
        ans <- (ans %*% pDerivs$jacobian)[1,]
        if(includeJacobian) ans <- ans + gr_logDetJacobian(p)
      }

      return(ans)
      returnType(double(1))
    },
    gr_logDens_pTransformed = function(pTransform = double(1)){
      ans <- gr_logDens(pTransform, trans = TRUE,
                        includeJacobian = includeJacobian_,
                        includePrior = includePrior_)
      return(ans)
      returnType(double(1))
    },
    gr_logDens_pTransformedFix1 = function(pTransform = double(1)){
      pTransform_star <- replaceOneVec(pTransform)
      ans <- gr_logDens(pTransform_star, trans = TRUE,
                        includeJacobian = includeJacobian_,
                        includePrior = includePrior_)

      return(ans[pTransform_indices_other])
      returnType(double(1))
    },
    setLogDensType = function(includeJacobian = logical(0, default = TRUE),
                              includePrior = logical(0, default = TRUE)){
      includeJacobian_ <<- includeJacobian
      includePrior_ <<- includePrior
    },
    ## For internal purposes of building the gradient
    logDetJacobian = function(pTransform = double(1)){
      ans <- paramsTransform$logDetJacobian(pTransform)
      return(ans)
      returnType(double())
    },
    ## Calculate MLE of parameters
    findMLE = function(pStart  = double(1, default = Inf),
                       hessian = logical(0, default = TRUE) ){
      mleRes <- optimize(pStart  = pStart,
                         includePrior = FALSE,
                         includeJacobian = FALSE,
                         hessian = hessian,
                         parscale = "real")
      return(mleRes)
      returnType(optimResultNimbleList())
    },
    ## Calculate posterior mode of parameters
    findMAP = function(pStart  = double(1, default = Inf),
                       hessian = logical(0, default = TRUE) ){
      mapRes <- optimize(pStart  = pStart,
                         includePrior = TRUE,
                         includeJacobian = TRUE,
                         hessian = hessian,
                         parscale = "real")
      return(mapRes)
      returnType(optimResultNimbleList())
    },
    replaceOneVec = function(pTransform = double(1)){
      pTransform_star <- numeric(value = 0, length = nparTrans)
      pTransform_star[pTransform_index_fixed] <- pTransform_fixed
      pTransform_star[pTransform_indices_other] <- pTransform[1:(nparTrans-1)]
      returnType(double(1))
      return(pTransform_star)
    },
    findMax_fixedp = function(pStartTransform = double(1, default = Inf),
                              pTransformIndex = integer(),
                              pTransformValue = double(),
                              includePrior = logical(0, default = FALSE),
                              includeJacobian = logical(0, default = FALSE),
                              hessian = logical(0, default = TRUE)){
      pTransform_index_fixed <<- pTransformIndex
      pTransform_fixed <<- pTransformValue
      pTransform_indices_other <<- pTransform_indices[pTransform_indices != pTransform_index_fixed]

      if(length(pStartTransform) == (nparTrans-1)) {
        pStartTransform_star <- replaceOneVec(pStartTransform)
        pStart <- paramsTransform$inverseTransform(pStartTransform_star)
      }else{
        if(length(pStartTransform) == nparTrans) pStart <- paramsTransform$inverseTransform(pStartTransform)
        else pStart <- numeric(value = Inf, length = 1)
      }
      maxRes <- optimize(pStart = pStart,
                         includePrior = includePrior,
                         includeJacobian = includeJacobian,
                         hessian = hessian,
                         parscale = "transformed",
                         keepOneFixed = TRUE)
      return(maxRes)
      returnType(optimResultNimbleList())
    },
    ## General Maximization Function
    optimize = function(pStart = double(1, default = Inf),
                        includePrior = logical(0, default = FALSE),
                        includeJacobian = logical(0, default = FALSE),
                        hessian = logical(0, default = TRUE),
                        parscale = character(0, default = "transformed"),
                        keepOneFixed = logical(0, default = FALSE)) {
      if(!one_time_fixes_done) one_time_fixes() ## Otherwise summary will look bad.
      if(multiSetsCheck & nQuad_ > 1) stop("Currently only Laplace (`nQuad = 1`) is supported for maximization when integrations have more than one dimension at a time. Use `updateSettings(nQuad = 1)` to change.")
      if(any(abs(pStart) == Inf)) pStart <- values(model, paramNodes)
      ## Catch for a model that hasn't been initiated...
      if(any(is.na(pStart))) pStart <- numeric(value = 0, length = npar)
      if(length(pStart) != npar) {
        print("  [Warning] For maximization, `pStart` should be length ", npar, " but is length ", length(pStart), ".")
        ans <- optimResultNimbleList$new()
        return(ans)
      # stop("Wrong length for pStart in findMLE.")
      }
      ## Reset log likelihood internally for cache.
      reset_outer_inner_logLik()

      if(includeJacobian & !includePrior)
        stop("Should not include a Jacobian transformation when not including the prior distribution in the log density calculation.")

      pStartTransform <- paramsTransform$transform(pStart)
      ## In case parameter nodes are not properly initialized.
      ## We need to check on transformed scale as that is where 0 as initial value makes sense generally.
      invalidStart <- is.na(pStartTransform) | is.nan(pStartTransform) | abs(pStartTransform) == Inf
      pStartTransform[invalidStart] <- 0

      ## Choose the MLE, or the MAP, or a penalized MLE (:= no Jacobian MAP).
      # optRes <- optim(pStartTransform, calcLogLik_pTransformed, gr_logLik_pTransformed, method = outerOptimMethod_, control = outerOptimControl_, hessian = hessian)

      setLogDensType(includeJacobian = includeJacobian, includePrior = includePrior)
      if( !keepOneFixed ){
        if(outerOptimUseAD) {
            optRes <- optim(pStartTransform, calcLogDens_pTransformed, gr_logDens_pTransformed,
                            method = outerOptimMethod_, control = outerOptimControl_, hessian = hessian)
        } else optRes <- optim(pStartTransform, calcLogDens_pTransformed,
                               method = outerOptimMethod_, control = outerOptimControl_, hessian = hessian)
        p <- paramsTransform$inverseTransform(optRes$par)
        if(parscale == "real") optRes$par <- p
      } else {
        if(outerOptimUseAD) {
            optRes <- optim(pStartTransform[pTransform_indices_other], calcLogDens_pTransformedFix1, gr_logDens_pTransformedFix1,
                            method = outerOptimMethod_, control = outerOptimControl_, hessian = hessian)
        } else optRes <- optim(pStartTransform[pTransform_indices_other], calcLogDens_pTransformedFix1,
                               method = outerOptimMethod_, control = outerOptimControl_, hessian = hessian)
        fullpar <- replaceOneVec(optRes$par)
        p <- paramsTransform$inverseTransform(fullpar)
      }
      setLogDensType()  ## Reset it to default to posterior.

      if(optRes$convergence != 0)
          print("  [Warning] In maximizing the Laplace/AGHQ approximation,\n",
                "            `optim` has a non-zero convergence code: ", optRes$convergence, ",\n",
                "            with the message '", optRes$message, "'.\n",
                "            The control parameters of `optim` can be adjusted using the `outerOptimControl`\n",
                "            list component of the `control` list argument of `buildLaplace` or `buildAGHQ`.")

      ## Print out warning about inner convergence.
      if( checkInnerConvergence(FALSE) != 0 )
          print("  [Warning] Inner optimization had a non-zero convergence code.\n",
                "            Use the `checkInnerConvergence(TRUE)` method of the Laplace object to see details.")

      ## Back transform results to original scale if requested.

      setModelValues(p) ## Make sure the model object contains all the updated parameter values.

      ## Returns on transformed scale just like optim.
      return(optRes)
      returnType(optimResultNimbleList())
    },
    profileLogDens = function(pTransformValue = double(),
                              pTransformIndex = integer(),
                              pStart = double(1, default = Inf),
                              maxLogDens = double(),
                              limit = double(),
                              includePrior = logical(0, default = FALSE),
                              includeJacobian = logical(0, default = FALSE)){
      pTransform_index_fixed <<- pTransformIndex ## Assuming no dimension changes so that pIndex =: pTransformIndex!
      pTransform_fixed <<- pTransformValue
      pTransform_indices_other <<- pTransform_indices[pTransform_indices != pTransform_index_fixed]

      maxRes <- optimize(pStart = pStart,
                         includePrior = includePrior,
                         includeJacobian = includeJacobian,
                         hessian = FALSE,
                         parscale = "transformed",
                         keepOneFixed = TRUE)
      ans <- maxRes$value - maxLogDens + limit/2
      return(ans)
      returnType(double())
    },
    ## User can update whether or not a warning is set for inner optimization.
    ## setInnerOptimWarning = function(warn = logical(0, default = FALSE)){
    ##   for(i in seq_along(AGHQuad_nfl)){
    ##     AGHQuad_nfl[[i]]$set_warning(warn)
    ##   }
    ## },
    ## Grab the inner Cholesky from the cached last values.
    cache_outer_logLik = function(logLikVal = double()){
      for(i in seq_along(AGHQuad_nfl)){
        AGHQuad_nfl[[i]]$save_outer_logLik(logLikVal)
      }
    },
    ## Set cached log lik values to -Inf internally.
    reset_outer_inner_logLik = function(){
      for(i in seq_along(AGHQuad_nfl)){
        AGHQuad_nfl[[i]]$reset_outer_logLik()
      }
    },
    ## Grab the inner Cholesky from the cached last values.
    get_inner_cholesky = function(atOuterMode = integer(0, default = 0)){
      if(nre == 0) stop("no random effects in the model")
      cholesky <- matrix(value = 0, nrow = nreTrans, ncol = nreTrans)
      tot <- 0
      for(i in seq_along(AGHQuad_nfl)){
        ## numre <- lenInternalRENodeSets[i]
        numre <- AGHQuad_nfl[[i]]$get_reTransLength()
        cholesky[(tot+1):(tot+numre), (tot+1):(tot+numre)] <- AGHQuad_nfl[[i]]$get_inner_negHessian_chol(atOuterMode)
        tot <- tot + numre
      }
      return(cholesky)
      returnType(double(2))
    },
    ## Grab the inner mode from the cached last values.
    get_inner_mode = function(atOuterMode = integer(0, default = 0)){
      if(nre == 0) stop("no random effects in the model")
      raneff <- numeric(nreTrans)
      tot <- 0
      for(i in seq_along(AGHQuad_nfl)){
        ## numre <- lenInternalRENodeSets[i]
        numre <- AGHQuad_nfl[[i]]$get_reTransLength()
        raneff[(tot+1):(tot+numre)] <- AGHQuad_nfl[[i]]$get_inner_mode(atOuterMode)
        tot <- tot + numre
      }
      return(raneff)
      returnType(double(1))
    },
    ## Optimized random effects given transformed parameter values
    optimRandomEffects = function(pTransform = double(1)){
      if(nre == 0) stop("No random effects in the model")
      p <- paramsTransform$inverseTransform(pTransform)
      raneff <- numeric(nreTrans)
      tmp <- numeric(nreTrans) ## Not sure this is needed.
      tot <- 0

      computeMethod <- -1
      if(useInnerCache_){
        pMLE <- AGHQuad_nfl[[1]]$get_param_value(atOuterMode = 1)
        pLast <- AGHQuad_nfl[[1]]$get_param_value(atOuterMode = 0)
        ## Cache check for either last value or MLE
        if(all(p == pMLE)) computeMethod <- 1
        else if(all(p == pLast)) computeMethod <- 0
      }

      for(i in seq_along(AGHQuad_nfl)){
        if(computeMethod == -1 ){ ## Is this valid?
          tmp <- AGHQuad_nfl[[i]]$update_max_logLik_RE(p)
        }else{
          tmp <- AGHQuad_nfl[[i]]$get_inner_mode(atOuterMode = computeMethod)
        }
        numre <- dim(tmp)[1]
        raneff[(tot+1):(tot+numre)] <- tmp
        tot <- tot + numre
      }
      return(raneff)
      returnType(double(1))
    },
    ## Inverse of the negative Hessian of log-likelihood wrt transformed random effects
    inverse_negHess = function(p = double(1), reTransform = double(1)){
      if(nre == 0) stop("no random effects in the model")
      invHess <- matrix(value = 0, nrow = nreTrans, ncol = nreTrans)
      tot <- 0

      outer_mode_case <- -1
      if(useInnerCache_){
        pMLE <- AGHQuad_nfl[[1]]$get_param_value(atOuterMode = 1)
        pLast <- AGHQuad_nfl[[1]]$get_param_value(atOuterMode = 0)
        ## Cache check for either last value or MLE
        if(all(p == pMLE)) outer_mode_case <- 1
        else if(all(p == pLast)) outer_mode_case <- 0
      }

      for(i in seq_along(AGHQuad_nfl)){
        ## numre <- lenInternalRENodeSets[i]
        numre <- AGHQuad_nfl[[i]]$get_reTransLength()
        if(outer_mode_case == -1){
          tmp <- AGHQuad_nfl[[i]]$negHess(p, reTransform[(tot+1):(tot+numre)])
        }else{
          U <- AGHQuad_nfl[[i]]$get_inner_negHessian_chol(atOuterMode = outer_mode_case)
          tmp <- t(U) %*% U
        }
        invHess[(tot+1):(tot+numre), (tot+1):(tot+numre)] <- inverse(tmp)
        tot <- tot + numre
      }
      return(invHess)
      returnType(double(2))
    },
    ## Hessian of joint log-likelihood wrt parameters and (transformed) random effects
    hess_logLik_wrt_p_wrt_re = function(p = double(1), reTransform = double(1)){
      if(nre == 0) stop("no random effects in the model")
      ans <- matrix(value = 0, nrow = npar, ncol = nreTrans)
      tot <- 0
      for(i in seq_along(AGHQuad_nfl)){
        ## numre <- lenInternalRENodeSets[i]
        numre <- AGHQuad_nfl[[i]]$get_reTransLength()
        # Alternative future computeMethod_ settings could be used here.
#        tmp <- AGHQuad_nfl[[i]]$hess_joint_logLik_wrt_p_wrt_re(p, reTransform[(tot+1):(tot+numre)])
        tmp <- AGHQuad_nfl[[i]]$he_P_RE_wrt_RE_wrt_P_b(p, reTransform[(tot+1):(tot+numre)])
        ans[1:npar, (tot+1):(tot+numre)] <- t(tmp)
        tot <- tot + numre
      }
      return(ans)
      returnType(double(2))
    },
    jac_gr_logLik_wrt_re_inDir_p = function(p = double(1), reTransform = double(1),
                                              inDir2D = double(2)) {
      if(nre == 0) stop("no random effects in the model")
      num_inDirs <- dim(inDir2D)[2] # should be nparTrans when called from summary
      ans <- matrix(value = 0, nrow = num_inDirs, ncol = nreTrans)
      tot <- 0
      for(i in seq_along(AGHQuad_nfl)){
        ## numre <- lenInternalRENodeSets[i]
        numre <- AGHQuad_nfl[[i]]$get_reTransLength()
        # Alternative future computeMethod_ settings could be used here.
        #        tmp <- AGHQuad_nfl[[i]]$hess_joint_logLik_wrt_p_wrt_re(p, reTransform[(tot+1):(tot+numre)])
        length_each_inDir <- npar + numre
        inDir1D <- nimNumeric(value = 0, length = length_each_inDir * num_inDirs)
        for(j in 1:num_inDirs) {
          iStart <- (j-1)*length_each_inDir
          inDir1D[(iStart+1):(iStart + npar)] <- inDir2D[1:npar, j]
        }
        tmp <- AGHQuad_nfl[[i]]$jac_gr_P_RE_wrt_RE_inDir(p,
                                                          reTransform[(tot+1):(tot+numre)],
                                                          inDir1D)
        ans[1:num_inDirs, (tot+1):(tot+numre)] <- t(tmp)
        tot <- tot + numre
      }
      return(ans)
      returnType(double(2))

    },
    ## Gives the user control to start fresh by removing internally saved values.
    ## setInnerCache = function(useCache = logical(0, default = TRUE)){
    ##   innerCache <<- useCache
    ##   for(i in seq_along(AGHQuad_nfl)) AGHQuad_nfl[[i]]$set_inner_cache(useCache)
    ## },
    ## Set all model values after finding the MLE. Function will repeat inner optimization if the inner cached values
    ## the inner cached values don't match p.
    setModelValues = function(p = double(1)){
      for(i in seq_along(AGHQuad_nfl))
         AGHQuad_nfl[[i]]$set_randomeffect_values(p)
    },
    ## Summarise AGHQuad MLE results
    summary = function(MLEoutput             = optimResultNimbleList(),
                       originalScale         = logical(0, default = TRUE),
                       randomEffectsStdError = logical(0, default = TRUE),
                       jointCovariance       = logical(0, default = FALSE)){
      if(dim(MLEoutput$hessian)[1] == 0) stop("Hessian matrix was not calculated for Laplace or AGHQ MLE")
      ## Output lists
      ans <- AGHQuad_summary$new()
      pres <- AGHQuad_params$new()
      ranres <- AGHQuad_params$new()
      ## Parameters
      p <- MLEoutput$par
      pTransform <- paramsTransform$transform(p)
      vcov_pTransform <- -inverse(MLEoutput$hessian)
      stdErr_pTransform <- sqrt(diag(vcov_pTransform))
      if(nre == 0) { ## No random effects
        ranres$estimate <- numeric(0)
        ranres$stdError <- numeric(0)
        if(originalScale){
          derivspInvTransform  <- derivs_pInverseTransform(pTransform, c(0, 1))
          JacobpInvTransform   <- derivspInvTransform$jacobian
          stdErr_p <- numeric(npar)
          if(jointCovariance) {
            vcov <- JacobpInvTransform %*% vcov_pTransform %*% t(JacobpInvTransform)
            stdErr_p <- sqrt(diag(vcov))
            ans$vcov <- vcov
          }
          else{
            for(i in 1:npar){
              var_p_i <- (JacobpInvTransform[i,,drop=FALSE] %*% vcov_pTransform %*% t(JacobpInvTransform[i,,drop=FALSE]))[1,1]
              stdErr_p[i] <- sqrt(var_p_i)
            }
            ans$vcov <- matrix(nrow = 0, ncol = 0)
          }
          pres$estimate <- p
          pres$stdError <- stdErr_p
        }
        else {
          pres$estimate <- pTransform
          pres$stdError <- stdErr_pTransform
          if(jointCovariance) ans$vcov <- vcov_pTransform
          else ans$vcov <- matrix(0, nrow = 0, ncol = 0)
        }
      }
      else{
        ## Random effects
        optreTransform <- optimRandomEffects(pTransform)  
        optre <- reInverseTransform(optreTransform)
        ntot <- npar + nreTrans
        if(jointCovariance) {
          ## Inverse of the negative Hessian of log-likelihood wrt transformed random effects at MLEs
          inv_negHess <- inverse_negHess(p, optreTransform)   
          jointInvNegHessZero <- matrix(0, nrow = ntot, ncol = ntot)
          jointInvNegHessZero[(npar+1):ntot, (npar+1):ntot] <- inv_negHess
          ## Derivative of inverse transformation for params
          derivspInvTransform  <- derivs_pInverseTransform(pTransform, c(0, 1))
          JacobpInvTransform   <- derivspInvTransform$jacobian
          ## Jacobian of optimized random effects wrt transformed parameters
          ## Hessian of log-likelihood wrt to params and transformed random effects
          if(summaryCalcMethod == 1) {
            hessLoglikwrtpre <- hess_logLik_wrt_p_wrt_re(p, optreTransform)
            JacobOptreWrtParams <- inv_negHess %*% t(hessLoglikwrtpre) %*% JacobpInvTransform
          } else {
            hessLLcrossterms_by_JpIT <- jac_gr_logLik_wrt_re_inDir_p(p, optreTransform, JacobpInvTransform)
            JacobOptreWrtParams <- inv_negHess %*% t(hessLLcrossterms_by_JpIT)
          }
          jointJacob <- matrix(init = FALSE, nrow = ntot, ncol = npar)
          jointJacob[(npar+1):ntot, 1:npar] <- JacobOptreWrtParams
          jointJacob[1:npar, 1:npar] <- diag(npar)
          ## Joint covariance matrix on transformed scale
          vcov_Transform <- jointInvNegHessZero + jointJacob %*% vcov_pTransform %*% t(jointJacob)
          if(originalScale){
            derivs_reInvTransform <- derivs_reInverseTransform(optreTransform, c(0, 1))
            Jacob_reInvTransform  <- derivs_reInvTransform$jacobian
            ## In case the number of random effects differs after transformation
            ntot2 <- npar + nre
            Jacob_JointInvTransform <- matrix(0, nrow = ntot2, ncol = ntot)
            Jacob_JointInvTransform[(npar+1):ntot2, (npar+1):ntot] <- Jacob_reInvTransform
            Jacob_JointInvTransform[1:npar, 1:npar] <- JacobpInvTransform
            vcov <- Jacob_JointInvTransform %*% vcov_Transform %*% t(Jacob_JointInvTransform)

            stdErr_p_re <- sqrt(diag(vcov))
            stdErr_p <- stdErr_p_re[1:npar]
            if(randomEffectsStdError){
              ranres$stdError <- stdErr_p_re[(npar+1):ntot2]
            }
            else{
              ranres$stdError <- numeric(0)
            }
            ans$vcov <- vcov
            pres$estimate <- p
            pres$stdError <- stdErr_p
            ranres$estimate <- optre
          }## End of if(originalScale)
          else { ## On transformed scale
            if(randomEffectsStdError){
              stdErr_reTransform <- sqrt(diag(vcov_Transform)[(npar+1):ntot])
              ranres$stdError <- stdErr_reTransform
            }
            else{
              ranres$stdError <- numeric(0)
            }
            ans$vcov <- vcov_Transform
            pres$estimate <- pTransform
            pres$stdError <- sqrt(diag(vcov_Transform)[1:npar])
            ranres$estimate <- optreTransform
          }
        }## End of if(jointCovariance)
        else { ## Do not return joint covariance matrix
          if(originalScale){## On original scale
            pres$estimate <- p
            ranres$estimate <- optre
            if(randomEffectsStdError){
              ## Joint covariance matrix on transform scale
              inv_negHess <- inverse_negHess(p, optreTransform)
              # jointInvNegHessZero <- matrix(0, nrow = ntot, ncol = ntot)
              # jointInvNegHessZero[1:nre, 1:nre] <- inv_negHess
              ## Derivative of inverse transformation for params
              derivspInvTransform  <- derivs_pInverseTransform(pTransform, c(0, 1))
              JacobpInvTransform   <- derivspInvTransform$jacobian
              ## Covariance matrix for params on the original scale
              vcov_p <- JacobpInvTransform %*% vcov_pTransform %*% t(JacobpInvTransform)
              ## Jacobian of optimized random effects wrt transformed parameters
              if(summaryCalcMethod == 1) {
              ## Hessian of log-likelihood wrt to params and transformed random effects
                hessLoglikwrtpre <- hess_logLik_wrt_p_wrt_re(p, optreTransform)
                JacobOptreWrtParams <- inv_negHess %*% t(hessLoglikwrtpre) %*% JacobpInvTransform
              } else {
                hessLLcrossterms_by_JpIT <- jac_gr_logLik_wrt_re_inDir_p(p, optreTransform, JacobpInvTransform)
                JacobOptreWrtParams <- inv_negHess %*% t(hessLLcrossterms_by_JpIT)
              }
              # jointJacob <- matrix(NA, nrow = ntot, ncol = npar)
              # jointJacob[1:nre, 1:npar] <- JacobOptreWrtParams
              # jointJacob[(nre+1):ntot, 1:npar] <- diag(npar)
              ## Join covariance matrix on transformed scale
              # vcov_Transform <- jointInvNegHessZero + jointJacob %*% vcov_pTransform %*% t(jointJacob)
              ## Covariance matrix for random effects (transformed)
              vcov_reTransform <- inv_negHess + JacobOptreWrtParams %*% vcov_pTransform %*% t(JacobOptreWrtParams)
              ## Derivatives information
              derivs_reInvTransform <- derivs_reInverseTransform(optreTransform, c(0, 1))
              Jacob_reInvTransform  <- derivs_reInvTransform$jacobian
              # Jacob_JointInvTransform <- matrix(0, nrow = ntot, ncol = ntot)
              # Jacob_JointInvTransform[1:nre, 1:nre] <- Jacob_reInvTransform
              # Jacob_JointInvTransform[(nre+1):ntot, (nre+1):ntot] <- JacobpInvTransform
              stdErr_re <- numeric(nre)
              for(i in 1:nre){
                var_i <- (Jacob_reInvTransform[i,,drop=FALSE] %*% vcov_reTransform %*% t(Jacob_reInvTransform[i,,drop=FALSE]))[1,1]
                stdErr_re[i] <- sqrt(var_i)
              }
              stdErr_p <- sqrt(diag(vcov_p))
              pres$stdError   <- stdErr_p
              ranres$stdError <- stdErr_re
              ans$vcov <- vcov_p
            }## End of if(randomEffectsStdError)
            else { ## Do not calculate standard errors of random effects estimates
              derivspInvTransform  <- derivs_pInverseTransform(pTransform, c(0, 1))
              JacobpInvTransform   <- derivspInvTransform$jacobian
              ## Covariance matrix for params on the original scale
              vcov_p <- JacobpInvTransform %*% vcov_pTransform %*% t(JacobpInvTransform)
              # stdErr_p <- numeric(npar)
              # for(i in 1:npar){
              #   var_p_i <- (JacobpInvTransform[i,,drop=FALSE] %*% vcov_pTransform %*% t(JacobpInvTransform[i,,drop=FALSE]))[1,1]
              #   stdErr_p[i] <- sqrt(var_p_i)
              # }
              stdErr_p <- sqrt(diag(vcov_p))
              pres$stdError <- stdErr_p
              ranres$stdError <- numeric(0)
              ans$vcov <- vcov_p
            }
          }## End of if(originalScale)
          else {## On transformed scale
            pres$estimate <- pTransform
            pres$stdError <- stdErr_pTransform
            ranres$estimate <- optreTransform
            ans$vcov <- vcov_pTransform
            if(randomEffectsStdError){
              inv_negHess <- inverse_negHess(p, optreTransform)
              ## jointInvNegHessZero <- matrix(0, nrow = ntot, ncol = ntot)
              ## jointInvNegHessZero[1:nreTrans, 1:nreTrans] <- inv_negHess
              ## Derivative of inverse transformation for params
              derivspInvTransform  <- derivs_pInverseTransform(pTransform, c(0, 1))
              JacobpInvTransform   <- derivspInvTransform$jacobian
              ## Jacobian of optimized random effects wrt transformed parameters
              ## Hessian of log-likelihood wrt to params and transformed random effects
              if(summaryCalcMethod == 1) {
                hessLoglikwrtpre <- hess_logLik_wrt_p_wrt_re(p, optreTransform)
                JacobOptreWrtParams <- inv_negHess %*% t(hessLoglikwrtpre) %*% JacobpInvTransform
              } else {
                hessLLcrossterms_by_JpIT <- jac_gr_logLik_wrt_re_inDir_p(p, optreTransform, JacobpInvTransform)
                JacobOptreWrtParams <- inv_negHess %*% t(hessLLcrossterms_by_JpIT)
              }
              stdErr_reTransform <- numeric(nreTrans)
              for(i in 1:nreTrans){
                var_reTransform_i <- inv_negHess[i, i] + (JacobOptreWrtParams[i,,drop=FALSE] %*% vcov_pTransform %*% t(JacobOptreWrtParams[i,,drop=FALSE]))[1,1]
                stdErr_reTransform[i] <- sqrt(var_reTransform_i)
              }
              ranres$stdError <- stdErr_reTransform
            }
            else{
              ranres$stdError <- numeric(0)
            }
          }
        }
      }
      if(originalScale) {
          pres$names <- paramNodesAsScalars_vec
          ranres$names <- reNodesAsScalars_vec
      } else {
        pres$names <- paramTransNodeNames
        ranres$names <- reTransNodeNames
      }
      ans$params <- pres
      ans$randomEffects <- ranres
      ans$originalScale <- originalScale
      return(ans)
      returnType(AGHQuad_summary())
    }
  ),
  buildDerivs = list(pInverseTransform  = list(),
                     reInverseTransform = list(),
                     otherLogLik = list(),
                     gr_otherLogLik_internal = list(),
                     logDetJacobian = list(),
                     calcPrior_p = list()
                     )
)

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
#'   space there are no constraints. (default = TRUE)
#'
#' @param randomEffectsStdError If TRUE, calculate the standard error of the
#'   estimates of random effects values. (default = TRUE)
#'
#' @param jointCovariance If TRUE, calculate the joint covariance matrix of
#'   the parameters and random effects together. If FALSE, calculate the
#'   covariance matrix of the parameters. (default = FALSE)
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
                           originalScale = TRUE,
                           randomEffectsStdError = TRUE,
                           jointCovariance = FALSE) {
  summary <- laplace$summary(MLEoutput, originalScale = originalScale,
                             randomEffectsStdError = randomEffectsStdError,
                             jointCovariance = jointCovariance)
  paramNames <- summary$params$names
  paramEsts <- summary$params$estimate
  if(length(paramEsts) < length(paramNames)) paramNames <- paramNames[1:(length(paramNames)-1)]
  names(paramEsts) <- paramNames
  stdErrParams <- summary$params$stdError
  paramsDF <- data.frame(estimate = paramEsts, stdError = stdErrParams, row.names = paramNames)

  REnames <- summary$randomEffects$names
  REests <- summary$randomEffects$estimate
  if(length(REests) < length(REnames)) REnames <- REnames[1:(length(REnames)-1)]
  REstdErrs <- summary$randomEffects$stdError
  if(length(REstdErrs))
    REDF <- data.frame(estimate = REests, stdError = REstdErrs, row.names = REnames)
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
       logLik = MLEoutput$value,
       df = length(paramEsts),
       originalScale = originalScale)
}

#' @rdname summaryLaplace
#' @export
summaryAGHQ <- function(AGHQ, MLEoutput,
                        originalScale =TRUE,
                        randomEffectsStdError = TRUE,
                        jointCovariance = FALSE) {
  summaryLaplace(AGHQ, MLEoutput, originalScale, randomEffectsStdError, jointCovariance)
}

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
#' @param pStart Initial values for parameters to begin optimization search for
#'   the maximum likelihood estimates. If omitted, the values currently in the
#'   (compiled or uncompiled) model object will be used.
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
runLaplace <- function(laplace, pStart,
                       originalScale = TRUE,
                       randomEffectsStdError = TRUE,
                       jointCovariance = FALSE) {
  if(missing(pStart)) pStart <- Inf # code to use values in model
  runAGHQ(AGHQ = laplace, pStart, originalScale, randomEffectsStdError, jointCovariance)
}

#' @rdname runLaplace
#' @export
runAGHQ <- function(AGHQ, pStart,
                    originalScale = TRUE,
                    randomEffectsStdError = TRUE,
                    jointCovariance = FALSE) {
  if(missing(AGHQ)) stop('runAGHQ: must provide a NIMBLE Laplace or AGHQ algorithm')
  if(!inherits(AGHQ, c("AGHQ","AGHQ_refClass")))
    stop('runAGHQ: AGHQ or Laplace argument must be a NIMBLE Laplace or AGHQ algorithm (compiled or uncompiled) from `buildLaplace` or `buildAGHQ`.')
  if(!is.Cnf(AGHQ)) {
    messageIfVerbose('  [Warning] Running an uncompiled Laplace or AGHQ algorithm.\n',
                     '            Use `compileNimble()` for faster execution.')
    tmp <- AGHQ$gr_logDens_pTransformed
    tmp <- AGHQ$calcLogDens_pTransformed
    for(i in seq_along(AGHQ$AGHQuad_nfl)) {
        tmp <- AGHQ$AGHQuad_nfl[[i]]$gr_inner_logLik
        tmp <- AGHQ$AGHQuad_nfl[[i]]$he_inner_logLik
    }
  }

  if(missing(pStart)) pStart <- Inf # code to use values in the model

  opt <- try(AGHQ$findMLE(pStart = pStart, hessian = TRUE))
  if(inherits(opt, "try-error"))
    stop("method findMLE had an error.")

  summary  <- try(summaryLaplace(laplace=AGHQ, MLEoutput=opt,
                              originalScale=originalScale,
                              randomEffectsStdError=randomEffectsStdError,
                              jointCovariance=jointCovariance))
  if(inherits(summary, "try-error")) {
    messageIfVerbose("  [Warning] `summaryLaplace` had an error. Only the MLE result will be returned.")
    summary <- NULL
  }
  list(MLE = opt, summary=summary)
}

#' Laplace approximation and adaptive Gauss-Hermite quadrature
#'
#' Build a Laplace or AGHQ approximation algorithm for a given NIMBLE model.
#'
#' @param model a NIMBLE model object, such as returned by \code{nimbleModel}.
#'   The model must have automatic derivatives (AD) turned on, e.g. by using
#'   \code{buildDerivs=TRUE} in \code{nimbleModel}.
#' @param nQuad number of quadrature points for AGHQ (in one dimension). Laplace approximation is
#'   AGHQ with `nQuad=1`. Only odd numbers of nodes really
#'   make sense. Often only one or a few nodes can achieve high accuracy. A maximum of
#'   35 nodes is supported. Note that for multivariate quadratures, the number
#'   of nodes will be (number of dimensions)^nQuad.
#' @param paramNodes a character vector of names of parameter nodes in the
#'   model; defaults are provided by \code{\link{setupMargNodes}}.
#'   Alternatively, \code{paramNodes} can be a list in the format returned by
#'   \code{setupMargNodes}, in which case \code{randomEffectsNodes},
#'   \code{calcNodes}, and \code{calcNodesOther} are not needed (and will be
#'   ignored).
#' @param randomEffectsNodes a character vector of names of continuous
#'   unobserved (latent) nodes to marginalize (integrate) over using Laplace/AGHQ
#'   approximation; defaults are provided by \code{\link{setupMargNodes}}.
#' @param calcNodes a character vector of names of nodes for calculating the
#'   integrand for Laplace/AGHQ approximation; defaults are provided by
#'   \code{\link{setupMargNodes}}. There may be deterministic nodes between
#'   \code{paramNodes} and \code{calcNodes}. These will be included in
#'   calculations automatically and thus do not need to be included in
#'   \code{calcNodes} (but there is no problem if they are).
#' @param calcNodesOther a character vector of names of nodes for calculating
#'   terms in the log-likelihood that do not depend on any
#'   \code{randomEffectsNodes}, and thus are not part of the marginalization,
#'   but should be included for purposes of finding the MLE. This defaults to
#'   stochastic nodes that depend on \code{paramNodes} but are not part of and
#'   do not depend on \code{randomEffectsNodes}. There may be deterministic
#'   nodes between \code{paramNodes} and \code{calcNodesOther}. These will be
#'   included in calculations automatically and thus do not need to be included
#'   in \code{calcNodesOther} (but there is no problem if they are).
#' @param control a named list for providing additional settings used in Laplace/AGHQ
#'   approximation. See \code{control} section below. Most of these can be
#'   updated later with the `updateSettings` method.
#'
#' @section \code{buildLaplace} and \code{buildAGHQ}:
#'
#' \code{buildLaplace} creates an object that can run Laplace approximation
#'   for a given model or part of a model. \code{buildAGHQ} creates an object
#'   that can run adaptive Gauss-Hermite quadrature (AGHQ, sometimes called
#'   "adaptive Gaussian quadrature") for a given model or part of a model.
#'   Laplace approximation is AGHQ with one quadrature point, hence
#'   `buildLaplace` simply calls `buildAGHQ` with `nQuad=1`. These methods
#'   approximate the integration over continuous random effects in a
#'   hierarchical model to calculate the (marginal) likelihood.
#'
#' \code{buildAGHQ} and \code{buildLaplace} will by default (unless changed
#' manually via `control$split`) determine from the model which random effects
#' can be integrated over (marginalized) independently. For example, in a GLMM
#' with a grouping factor and an independent random effect intercept for each
#' group, the random effects can be marginalized as a set of univariate
#' approximations rather than one multivariate approximation. On the other hand,
#' correlated or nested random effects would require multivariate marginalization.
#'
#' Maximum likelihood estimation is available for Laplace approximation
#' (`nQuad=1`) with univariate or multivariate integrations. With `nQuad > 1`,
#' maximum likelihood estimation is available only if all integrations are
#' univariate (e.g., a set of univariate random effects). If there are
#' multivariate integrations, these can be calculated at chosen input parameters
#' but not maximized over parameters. For example, one can find the MLE based on
#' Laplace approximation and then increase `nQuad` (using the `updateSettings`
#' method below) to check on accuracy of the marginal log likelihood at the MLE.
#'
#' Beware that quadrature will use `nQuad^k` quadrature points, where `k` is the
#' dimension of each integration. Therefore quadrature for `k` greater that 2 or
#' 3 can be slow. As just noted, `buildAGHQ` will determine independent
#' dimensions of quadrature, so it is fine to have a set of univariate random
#' effects, as these will each have k=1. Multivariate quadrature (k>1) is only
#' necessary for nested, correlated, or otherwise dependent random effects.
#'
#' The recommended way to find the maximum likelihood estimate and associated
#' outputs is by calling \code{\link{runLaplace}} or \code{\link{runAGHQ}}. The
#' input should be the compiled Laplace or AGHQ algorithm object. This would be
#' produced by running \code{\link{compileNimble}} with input that is the result
#' of \code{buildLaplace} or \code{buildAGHQ}.
#'
#' For more granular control, see below for methods \code{findMLE} and
#'   \code{summary}. See function \code{\link{summaryLaplace}} for an easier way
#'   to call the \code{summary} method and obtain results that include node
#'   names. These steps are all done within \code{runLaplace} and
#'   \code{runAGHQ}.
#'
#' The NIMBLE User Manual at r-nimble.org also contains an example of Laplace
#' approximation.
#'
#' @section How input nodes are processed:
#'
#' \code{buildLaplace} and \code{buildAGHQ} make good tries at deciding what
#' to do with the input model and any (optional) of the node arguments. However,
#' random effects (over which approximate integration will be done) can be
#' written in models in multiple equivalent ways, and customized use cases may
#' call for integrating over chosen parts of a model. Hence, one can take full
#' charge of how different parts of the model will be used.
#'
#' Any of the input node vectors, when provided, will be processed using
#'   \code{nodes <- model$expandNodeNames(nodes)}, where \code{nodes} may be
#'   \code{paramNodes}, \code{randomEffectsNodes}, and so on. This step allows
#'   any of the inputs to include node-name-like syntax that might contain
#'   multiple nodes. For example, \code{paramNodes = 'beta[1:10]'} can be
#'   provided if there are actually 10 scalar parameters, 'beta[1]' through
#'   'beta[10]'. The actual node names in the model will be determined by the
#'   \code{exapndNodeNames} step.
#'
#' In many (but not all) cases, one only needs to provide a NIMBLE model object
#'   and then the function will construct reasonable defaults necessary for
#'   Laplace approximation to marginalize over all continuous latent states
#'   (aka random effects) in a model. The default values for the four groups of
#'   nodes are obtained by calling \code{\link{setupMargNodes}}, whose arguments
#'   match those here (except for a few arguments which are taken from control
#'   list elements here).
#'
#' \code{setupMargNodes} tries to give sensible defaults from
#'   any combination of \code{paramNodes}, \code{randomEffectsNodes},
#'   \code{calcNodes}, and \code{calcNodesOther} that are provided. For example,
#'   if you provide only \code{randomEffectsNodes} (perhaps you want to
#'   marginalize over only some of the random effects in your model),
#'   \code{setupMargNodes} will try to determine appropriate choices for the
#'   others.
#'
#' \code{setupMargNodes} also determines which integration dimensions are
#' conditionally independent, i.e., which can be done separately from each
#' other. For example, when possible, 10 univariate random effects will be split
#' into 10 univariate integration problems rather than one 10-dimensional
#' integration problem.
#'
#' The defaults make general assumptions such as that
#'   \code{randomEffectsNodes} have \code{paramNodes} as parents. However, The
#'   steps for determining defaults are not simple, and it is possible that they
#'   will be refined in the future. It is also possible that they simply don't
#'   give what you want for a particular model. One example where they will not
#'   give desired results can occur when random effects have no prior
#'   parameters, such as `N(0,1)` nodes that will be multiplied by a scale
#'   factor (e.g. sigma) and added to other explanatory terms in a model. Such
#'   nodes look like top-level parameters in terms of model structure, so
#'   you must provide a \code{randomEffectsNodes} argument to indicate which
#'   they are.
#'
#' It can be helpful to call \code{setupMargNodes} directly to see exactly how
#'   nodes will be arranged for Laplace approximation. For example, you may want
#'   to verify the choice of \code{randomEffectsNodes} or get the order of
#'   parameters it has established to use for making sense of the MLE and
#'   results from the \code{summary} method. One can also call
#'   \code{setupMargNodes}, customize the returned list, and then provide that
#'   to \code{buildLaplace} as \code{paramNodes}. In that case,
#'   \code{setupMargNodes} will not be called (again) by \code{buildLaplace}.
#'
#' If \code{setupMargNodes} is emitting an unnecessary warning, simply use
#'   \code{control=list(check=FALSE)}.
#'
#' @section Managing parameter transformations that may be used internally:
#'
#' If any \code{paramNodes} (parameters) or \code{randomEffectsNodes} (random
#'   effects / latent states) have constraints on the range of valid values
#'   (because of the distribution they follow), they will be used on a
#'   transformed scale determined by \code{parameterTransform}. This means the
#'   Laplace approximation itself will be done on the transformed scale for
#'   random effects and finding the MLE will be done on the transformed scale
#'   for parameters. For parameters, prior distributions are not included in
#'   calculations, but they are used to determine valid parameter ranges and
#'   hence to set up any transformations. For example, if \code{sigma} is a
#'   standard deviation, you can declare it with a prior such as \code{sigma ~
#'   dhalfflat()} to indicate that it must be greater than 0.
#'
#' For default determination of when transformations are needed, all parameters
#'   must have a prior distribution simply to indicate the range of valid
#'   values. For a param \code{p} that has no constraint, a simple choice is
#'   \code{p ~ dflat()}.
#'
#' @section Understanding inner and outer optimizations:
#'
#' Note that there are two numerical optimizations when finding maximum
#' likelihood estimates with a Laplace or (1D) AGHQ algorithm: (1) maximizing
#' the joint log-likelihood of random effects and data given a parameter value
#' to construct the approximation to the marginal log-likelihood at the given
#' parameter value; (2) maximizing the approximation to the marginal
#' log-likelihood over the parameters. In what follows, the prefix 'inner'
#' refers to optimization (1) and 'outer' refers to optimization (2). Currently
#' both optimizations default to using method \code{"nlminb"}. However, one can
#' use other optimizers or simply run optimization (2) manually from R; see the
#' example below. In some problems, choice of inner and/or outer optimizer can
#' make a big difference for obtaining accurate results, especially for standard
#' errors. Hence it is worth experimenting if one is in doubt.
#'
#' @section \code{control} list arguments:
#'
#' The \code{control} list allows additional settings to be made using named
#' elements of the list. Most (or all) of these can be updated later using the
#' `updateSettings` method. Supported elements include:
#'
#' \itemize{
#'
#'   \item \code{split}. If TRUE (default), \code{randomEffectsNodes} will be
#'         split into conditionally independent sets if possible. This
#'         facilitates more efficient Laplace or AGHQ approximation because each
#'         conditionally independent set can be marginalized independently. If
#'         FALSE, \code{randomEffectsNodes} will be handled as one multivariate
#'         block, with one multivariate approximation. If \code{split} is a
#'         numeric vector, \code{randomEffectsNodes} will be split by calling
#'         \code{split}(\code{randomEffectsNodes}, \code{control$split}). The
#'         last option allows arbitrary control over how
#'         \code{randomEffectsNodes} are blocked.
#'
#'   \item \code{check}. If TRUE (default), a warning is issued if
#'         \code{paramNodes}, \code{randomEffectsNodes} and/or \code{calcNodes}
#'         are provided but seem to have missing or unnecessary
#'         elements based on some default inspections of the model. If
#'         unnecessary warnings are emitted, simply set \code{check=FALSE}.
#'
#'   \item \code{innerOptimControl}. An `optimControlNimbleList` list
#'         of control parameters (an R list is sufficient for uncompiled operation) for the inner
#'         optimization of Laplace approximation using \code{nimOptim}. See
#'         'Details' of \code{\link{nimOptim}} for further information. Default
#'         is `nimOptimDefaultControl()`.
#'
#'   \item \code{innerOptimMethod}. Optimization method to be used in
#'         \code{nimOptim} for the inner optimization. See 'Details' of
#'         \code{\link{nimOptim}}. Currently \code{nimOptim} in NIMBLE supports:
#'         \code{"Nelder-Mead"}", \code{"BFGS"}, \code{"CG"}, \code{"L-BFGS-B"},
#'         \code{"nlminb"}, \code{"bobyqa"}, and user-provided optimizers. By default, method
#'         \code{"nlminb"} is used for both univariate and multivariate cases. For
#'         \code{"nlminb"}, \code{"bobyqa"}, or user-provided optimizers, only a subset of
#'         elements of the \code{innerOptimControlList} are supported. (Note
#'         that control over the outer optimization method is available as an
#'         argument to `findMLE`). Choice of optimizers can be important and so
#'         can be worth exploring.
#'
#'   \item \code{innerOptimStart}. Method for determining starting values for
#'         the inner optimization. Options are:
#'
#' \itemize{
#'
#' \item \code{"last.best"} (default): use optimized random effects values corresponding to
#'         the best outer optimization (i.e. the largest marginal log likelihood value) so far
#'         for each conditionally independent part of the approximation;
#'
#' \item \code{"last"}: use the result of the last inner optimization;
#'
#' \item \code{"zero"}: use all zeros;
#'
#' \item \code{"constant"}: always use the same values, determined by
#'         \code{innerOptimStartValues};
#'
#' \item \code{"random"}: randomly draw new starting values from the
#'       model (i.e., from the prior);
#'
#' \item \code{"model"}: use values for random effects stored in the
#'          model, which are determined from the first call.
#'
#' }
#'
#'       Note that \code{"model"} and \code{"zero"} are shorthand for
#'         \code{"constant"} with particular choices of
#'         \code{innerOptimStartValues}. Note that \code{"last"} and
#'         \code{"last.best"} require a choice for the very first values, which will
#'         come from \code{innerOptimStartValues}. The default is
#'         \code{innerOptimStart="zero"} and may change in the future.
#'
#'   \item \code{innerOptimStartValues}. Values for some of
#'         \code{innerOptimStart} approaches. If a scalar is provided, that
#'         value is used for all elements of random effects for each
#'         conditionally independent set. If a vector is provided, it must be
#'         the length of *all* random effects. If these are named (by node
#'         names), the names will be used to split them correctly among each
#'         conditionally independent set of random effects. If they are not
#'         named, it is not always obvious what the order should be because it
#'         may depend on the conditionally independent sets of random
#'         effects. It should match the order of names returned as part of
#'         `summaryLaplace`.
#'
#'   \item \code{innerOptimWarning}. If FALSE (default), do not emit warnings
#'   from the inner optimization. Optimization methods may sometimes emit a
#'   warning such as for bad parameter values encountered during the
#'   optimization search. Often, a method can recover and still find the
#'   optimum. In the approximations here, sometimes the inner optimization
#'   search can fail entirely, yet the outer optimization see this as one failed
#'   parameter value and can recover. Hence, it is often desirable to silence
#'   warnings from the inner optimizer, and this is done by default. Set
#'   \code{innerOptimWarning=TRUE} to see all warnings.
#'
#'   \item \code{useInnerCache}. If TRUE (default), use caching system for
#'     efficiency of inner optimizations. The caching system records one set of
#'     previous parameters and uses the corresponding results if those parameters
#'    are used again (e.g., in a gradient call). This should generally not be
#'    modified.
#'
#'   \item \code{outerOptimMethod}. Optimization method to be used in
#'         \code{nimOptim} for the outer optimization. See 'Details' of
#'         \code{\link{nimOptim}}. Currently \code{nimOptim} in NIMBLE supports:
#'         \code{"Nelder-Mead"}", \code{"BFGS"}, \code{"CG"}, \code{"L-BFGS-B"},
#'         \code{"nlminb"}, \code{"bobyqa"}, and user-provided optimizers. By default, method
#'         \code{"nlminb"} is used for both univariate and multivariate cases,
#'         although some problems may benefit from other choices. For
#'         \code{"nlminb"}, \code{"bobyqa"}, or user-provided optimizers, only a subset of
#'         elements of the \code{innerOptimControlList} are supported. (Note
#'         that control over the outer optimization method is available as an
#'         argument to `findMLE`). Choice of optimizers can be important and so
#'         can be worth exploring.
#'
#' \item \code{outerOptimControl}. An `optimControlNimbleList` of control parameters
#'         (an R list is sufficient for uncompiled operation) for maximizing
#'         the Laplace log-likelihood using \code{nimOptim}. See 'Details' of
#'         \code{\link{nimOptim}} for further information.
#'
#'   \item \code{computeMethod}. There are three approaches available for
#'   internal details of how the approximations, and specifically derivatives
#'   involved in their calculation, are handled. These are labeled simply 1, 2,
#'   and 3, and the default is 2. The relative performance of the methods will
#'   depend on the specific model. Users wanting to explore efficiency can try
#'   switching from method 2 (default) to methods 1 or 3 and comparing
#'   performance. The first Laplace approximation with each method will be
#'   (much) slower than subsequent Laplace approximations. Further details are
#'   not provided at this time.
#'
#'  \item \code{quadTransform} (relevant only \code{nQuad>1}). For multivariate AGHQ,
#'  a grid must be constructed based on the Hessian at the inner mode. Options
#'  include "cholesky" (default) and "spectral" (i.e., eigenvectors and
#'  eigenvalues) for the corresponding matrix decompositions on which the grid
#'  can be based.
#'
#' \item \code{outerOptimUseAD}. The optimization of the (hyper)parameters (the
#' "outer" optimization can provide an AD-based gradient to the chosen outer
#' optimization function or can omit this, causing any derivative-based
#' optimization method to use finite differences. Turning this off allows one
#' to avoid any complexity associatend with use of AD applied to the inner
#' Laplace/AGHQ approximation.
#'
#' } # end itemize
#'
#' @section Available methods:
#'
#' The object returned by \code{buildLaplace} or \code{buildAGHQ} is a nimbleFunction object with
#' numerous methods (functions). Here these are described in three tiers of user
#' relevance.
#'
#' @section Most useful methods:
#'
#' The most relevant methods to a user are:
#'
#' \itemize{
#'
#' \item \code{calcLogLik(p, trans=FALSE)}. Calculate the approximation to the
#'       marginal log-likelihood function at parameter value \code{p}, which (if
#'       \code{trans} is FALSE) should match the order of \code{paramNodes}. For
#'       any non-scalar nodes in \code{paramNodes}, the order within the node is
#'       column-major. The order of names can be obtained from method
#'       \code{getNodeNamesVec(TRUE)}. Return value is the scalar (approximate,
#'       marginal) log likelihood.
#'
#'       If \code{trans} is TRUE, then \code{p} is the vector of parameters on
#'       the transformed scale, if any, described above. In this case, the
#'       parameters on the original scale (as the model was written) will be
#'       determined by calling the method \code{pInverseTransform(p)}. Note that
#'       the length of the parameter vector on the transformed scale might not
#'       be the same as on the original scale (because some constraints of
#'       non-scalar parameters result in fewer free transformed parameters than
#'       original parameters).
#'
#' \item \code{calcLaplace(p, trans)}. This is the same as \code{calcLogLik} but
#'        requires that the approximation be Laplace (i.e \code{nQuad} is 1),
#'        and results in an error otherwise.
#'
#' \item \code{findMLE(pStart, hessian)}. Find the maximum likelihood
#'         estimates of parameters using the approximated marginal likelihood.
#'         This can be used if \code{nQuad} is 1 (Laplace case) or if
#'         \code{nQuad>1} and all marginalizations involve only univariate
#'         random effects. Arguments are \code{pStart}: initial parameter
#'         values (defaults to parameter values currently in the model);
#'          and \code{hessian}: whether to calculate and return the
#'         Hessian matrix (defaults to \code{TRUE}, which is required for
#'         subsequent use of \code{summary} method). Second derivatives in the
#'         Hessian are determined by finite differences of the gradients
#'         obtained by automatic differentiation (AD). Return value is a
#'         nimbleList of type \code{optimResultNimbleList}, similar to what is
#'         returned by R's optim. See \code{help(nimOptim)}. Note that
#'         parameters (\code{par}) are returned for the natural parameters, i.e. how
#'         they are defined in the model. But the \code{hessian}, if requested, is
#'         computed for the parameters as transformed for optimization if
#'         necessary. Hence one must be careful interpreting `hessian` if any
#'         parameters have constraints, and the safest next step is to use the
#'         \code{summary} method or \code{summaryLaplace} function.
#'
#' \item \code{summary(MLEoutput, originalScale, randomEffectsStdError,
#'        jointCovariance)}. Summarize the maximum likelihood estimation
#'        results, given object \code{MLEoutput} that was returned by
#'        \code{findMLE}. The summary can include a covariance matrix for the
#'        parameters, the random effects, or both), and these can be returned on
#'        the original parameter scale or on the (potentially) transformed
#'        scale(s) used in estimation. It is often preferred instead to call
#'        function (not method) `summaryLaplace` because this will attach
#'        parameter and random effects names (i.e., node names) to the results.
#'
#' In more detail, \code{summary} accepts the following optional arguments:
#'
#'        \itemize{
#'
#'           \item \code{originalScale}. Logical. If TRUE, the function returns
#'           results on the original scale(s) of parameters and random effects;
#'           otherwise, it returns results on the transformed scale(s). If there
#'           are no constraints, the two scales are identical. Defaults to TRUE.
#'
#'           \item \code{randomEffectsStdError}. Logical. If TRUE, standard
#'           errors of random effects will be calculated.
#'           Defaults to TRUE.
#'
#'           \item \code{jointCovariance}. Logical. If TRUE, the joint
#'           variance-covariance matrix of the parameters and the random effects
#'           will be returned. If FALSE, the variance-covariance matrix of the
#'           parameters will be returned. Defaults to FALSE.
#'
#'        }
#'
#'        The object returned by \code{summary} is an \code{AGHQuad_summary}
#'        nimbleList with elements:
#'
#'        \itemize{
#'
#'           \item \code{params}. A nimbleList that contains estimates and
#'           standard errors of parameters (on the original or transformed
#'           scale, as chosen by \code{originalScale}).
#'
#'           \item \code{randomEffects}. A nimbleList that contains estimates of
#'           random effects and, if requested
#'           (\code{randomEffectsStdError=TRUE}) their standard errors, on
#'           original or transformed scale. Standard errors are calculated
#'           following the generalized delta method of Kass and Steffey (1989).
#'
#'           \item \code{vcov}. If requested (i.e.
#'           \code{jointCovariance=TRUE}), the joint variance-covariance
#'           matrix of the parameters and random effects, on original or
#'           transformed scale. If \code{jointCovariance=FALSE}, the
#'           covariance matrix of the parameters, on original or transformed
#'           scale.
#'
#'           \item \code{scale}. \code{"original"} or \code{"transformed"}, the
#'           scale on which results were requested.
#'        }
#'     }
#'
#'
#' @section Methods for more advanced uses:
#'
#' Additional methods to access or control more details of the Laplace/AGHQ
#' approximation include:
#'
#' \itemize{
#'
#'   \item \code{updateSettings}. This provides a single function through which
#'   many of the settings described above (mostly for the \code{control} list)
#'   can be later changed. Options that can be changed include:
#'   \code{innerOptimMethod}, \code{innerOptimStart},
#'   \code{innerOptimStartValues}, \code{useInnerCache}, \code{nQuad},
#'   \code{quadTransform}, \code{innerOptimControl}, \code{outerOptimControl}, and
#'   \code{computeMethod}. For \code{innerOptimStart}, method "zero" cannot be
#'   specified but can be achieved by choosing method "constant" with
#'   \code{innerOptimStartValues=0}. Only provided options will be modified. The
#'   exceptions are \code{innerOptimControl}, \code{outerOptimControl}, which
#'   are replaced only when \code{replace_innerOptimControl=TRUE} or
#'   \code{replace_outerOptimControl=TRUE}, respectively.
#'
#'   \item \code{getNodeNamesVec(returnParams)}. Return a vector (>1) of names
#'   of parameters/random effects nodes, according to \code{returnParams =
#'   TRUE/FALSE}. Use this if there is more than one node.
#'
#'   \item \code{getNodeNameSingle(returnParams)}. Return the name of a
#'   single parameter/random effect node, according to \code{returnParams =
#'   TRUE/FALSE}. Use this if there is only one node.
#'
#'   \item \code{checkInnerConvergence(message)}. Checks whether all internal
#'   optimizers converged. Returns a zero if everything converged and one
#'   otherwise. If \code{message = TRUE}, it will print more details about
#'   convergence for each conditionally independent set.
#'
#'   \item \code{gr_logLik(p, trans)}. Gradient of the (approximated)
#'   marginal log-likelihood at parameter value \code{p}. Argument \code{trans}
#'   is similar to that in \code{calcLaplace}. If there are multiple parameters,
#'   the vector \code{p} is given in the order of parameter names returned by
#'   \code{getNodeNamesVec(returnParams=TRUE)}.
#'
#'   \item \code{gr_Laplace(p, trans)}. This is the same as \code{gr_logLik}.
#'
#'   \item \code{otherLogLik(p)}. Calculate the \code{calcNodesOther}
#'   nodes, which returns the log-likelihood of the parts of the model that are
#'   not included in the Laplace or AGHQ approximation.
#'
#'   \item \code{gr_otherLogLik(p)}. Gradient (vector of derivatives with
#'   respect to each parameter) of \code{otherLogLik(p)}. Results should
#'   match \code{gr_otherLogLik_internal(p)} but may be more efficient after
#'   the first call.
#'
#' }
#'
#' @section Internal or development methods:
#'
#' Some methods are included for calculating the (approximate) marginal log
#' posterior density by including the prior distribution of the parameters. This
#' is useful for finding the maximum a posteriori probability (MAP) estimate.
#' Currently these are provided for point calculations without estimation methods.
#'
#' \itemize{
#'
#'   \item \code{calcPrior_p(p)}. Log density of prior distribution.
#'
#'   \item \code{calcPrior_pTransformed(pTransform)}. Log density of prior distribution on transformed scale, includes the Jacobian.
#'
#'   \item \code{calcPostLogDens(p)}. Marginal log posterior density in terms of the parameter p.
#'
#'   \item \code{calcPostLogDens_pTransformed (pTransform)}. Marginal log posterior density in terms of the transformed
#'   parameter, which includes the Jacobian transformation.
#'
#'   \item \code{gr_postLogDens_pTransformed(pTransform)}. Graident of marginal log posterior density on the transformed scale.
#'   Other available options that are used in the derivative for more flexible include \code{logDetJacobian(pTransform)} and
#'   \code{gr_logDeJacobian(pTransform)}, as well as \code{gr_prior(p)}.
#' }
#'
#' Finally, methods that are primarily for internal use by other methods include:
#'
#' \itemize{
#'
#'    \item \code{gr_logLik_pTransformed}. Gradient of the Laplace
#'     approximation (\code{calcLogLik_pTransformed(pTransform)}) at transformed
#'     (unconstrained) parameter value \code{pTransform}.
#'
#'    \item \code{pInverseTransform(pTransform)}. Back-transform the transformed
#'    parameter value \code{pTransform} to original scale.
#'
#'    \item \code{derivs_pInverseTransform(pTransform, order)}. Derivatives of
#'    the back-transformation (i.e. inverse of parameter transformation) with
#'    respect to transformed parameters at \code{pTransform}. Derivative order
#'    is given by \code{order} (any of 0, 1, and/or 2).
#'
#'    \item \code{reInverseTransform(reTrans)}. Back-transform the transformed
#'    random effects value \code{reTrans} to original scale.
#'
#'    \item \code{derivs_reInverseTransform(reTrans, order)}. Derivatives of the
#'    back-transformation (i.e. inverse of random effects transformation) with
#'    respect to transformed random effects at \code{reTrans}. Derivative order
#'    is given by \code{order} (any of 0, 1, and/or 2).
#'
#'    \item \code{optimRandomEffects(pTransform)}. Calculate the optimized
#'    random effects given transformed parameter value \code{pTransform}. The
#'    optimized random effects are the mode of the conditional distribution of
#'    random effects given data at parameters \code{pTransform}, i.e. the
#'    calculation of \code{calcNodes}.
#'
#'    \item \code{inverse_negHess(p, reTransform)}. Calculate the inverse of the
#'    negative Hessian matrix of the joint (parameters and random effects)
#'    log-likelihood with respect to transformed random effects, evaluated at
#'    parameter value \code{p} and transformed random effects
#'    \code{reTransform}.
#'
#'    \item \code{hess_logLik_wrt_p_wrt_re(p, reTransform)}. Calculate the
#'    Hessian matrix of the joint log-likelihood with respect to parameters and
#'    transformed random effects, evaluated at parameter value \code{p} and
#'    transformed random effects \code{reTransform}.
#'
#'   \item \code{one_time_fixes()}. Users never need to run this. Is is called
#'   when necessary internally to fix dimensionality issues if there is only
#'   one parameter in the model.
#'
#'   \item \code{calcLogLik_pTransformed(pTransform)}. Laplace approximation at
#'         transformed (unconstrained) parameter value \code{pTransform}. To
#'         make maximizing the Laplace likelihood unconstrained, an automated
#'         transformation via \code{\link{parameterTransform}} is performed on
#'         any parameters with constraints indicated by their priors (even
#'         though the prior probabilities are not used).
#'
#'   \item \code{gr_otherLogLik_internal(p)}. Gradient (vector of
#'   derivatives with respect to each parameter) of \code{otherLogLik(p)}.
#'   This is obtained using automatic differentiation (AD) with single-taping.
#'   First call will always be slower than later calls.
#'
#'   \item \code{cache_outer_logLik(logLikVal)}. Save the marginal log likelihood value
#'   to the inner Laplace mariginlization functions to track the outer maximum internally.
#'
#'   \item \code{reset_outer_inner_logLik()}. Reset the internal saved maximum marginal log likelihood.
#'
#'   \item \code{get_inner_cholesky(atOuterMode = integer(0, default = 0))}. Returns the cholesky
#'   of the negative Hessian with respect to the random effects. If \code{atOuterMode = 1} then returns
#'   the value at the overall best marginal likelihood value, otherwise \code{atOuterMode = 0} returns the last.
#'
#'   \item \code{get_inner_mode(atOuterMode = integer(0, default = 0))}. Returns the mode of the random effects
#'   for either the last call to the innner quadrature functions (\code{atOuterMode = 0} ), or the last best
#'   value for the marginal log likelihood, \code{atOuterMode = 1}.
#'
#' }
#'
#' @author Wei Zhang, Perry de Valpine, Paul van Dam-Bates
#'
#' @name laplace
#'
#' @aliases Laplace buildLaplace AGHQuad buildAGHQ AGHQ
#'
#' @examples
#' pumpCode <- nimbleCode({
#'   for (i in 1:N){
#'     theta[i] ~ dgamma(alpha, beta)
#'     lambda[i] <- theta[i] * t[i]
#'     x[i] ~ dpois(lambda[i])
#'   }
#'   alpha ~ dexp(1.0)
#'   beta ~ dgamma(0.1, 1.0)
#' })
#' pumpConsts <- list(N = 10, t = c(94.3, 15.7, 62.9, 126, 5.24, 31.4, 1.05, 1.05, 2.1, 10.5))
#' pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22))
#' pumpInits <- list(alpha = 0.1, beta = 0.1, theta = rep(0.1, pumpConsts$N))
#' pump <- nimbleModel(code = pumpCode, name = "pump", constants = pumpConsts,
#'                     data = pumpData, inits = pumpInits, buildDerivs = TRUE)
#'
#' # Build Laplace approximation
#' pumpLaplace <- buildLaplace(pump)
#'
#' \dontrun{
#' # Compile the model
#' Cpump <- compileNimble(pump)
#' CpumpLaplace <- compileNimble(pumpLaplace, project = pump)
#' # Calculate MLEs of parameters
#' MLEres <- CpumpLaplace$findMLE()
#' # Calculate estimates and standard errors for parameters and random effects on original scale
#' allres <- CpumpLaplace$summary(MLEres, randomEffectsStdError = TRUE)
#'
#' # Change the settings and also illustrate runLaplace
#' innerControl <- nimOptimDefaultControl()
#' innerControl$maxit <- 1000
#' CpumpLaplace$updateSettings(innerOptimControl = innerControl,
#'                             replace_innerOptimControl = TRUE)
#' newres <- runLaplace(CpumpLaplace)
#'
#' # Illustrate use of the component log likelihood and gradient functions to
#' # run an optimizer manually from R.
#' # Use nlminb to find MLEs
#' MLEres.manual <- nlminb(c(0.1, 0.1),
#'                         function(x) -CpumpLaplace$calcLogLik(x),
#'                         function(x) -CpumpLaplace$gr_Laplace(x))
#' }
#'
#' @references
#'
#' Kass, R. and Steffey, D. (1989). Approximate Bayesian inference in
#' conditionally independent hierarchical models (parametric empirical Bayes
#' models). \emph{Journal of the American Statistical Association}, 84(407),
#' 717-726.
#'
#' Liu, Q. and Pierce, D. A. (1994). A Note on Gauss-Hermite Quadrature. \emph{Biometrika}, 81(3) 624-629.
#'
#' Jackel, P. (2005). A note on multivariate Gauss-Hermite quadrature. London: \emph{ABN-Amro. Re.}
#'
#' Skaug, H. and Fournier, D. (2006). Automatic approximation of the marginal
#' likelihood in non-Gaussian hierarchical models. \emph{Computational
#' Statistics & Data Analysis}, 56, 699-709.
#'
NULL
