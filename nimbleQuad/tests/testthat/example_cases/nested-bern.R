library(nimbleQuad)
crossed <- FALSE # crossed or nested REs?
## most (initially reported) results below for Laplace are for crossed

## Consider introducing imbalance

qpts <- c(.025,.25,.5,.75,.975)

set.seed(1)
n <- 1e5 # 1000, 10000, 100000, 1e6, 7283575

ntowns <- 262
nstates <- 47

gender <- sample(c("M","F"), n, replace = TRUE)
livingarrangement <- sample(c("H","D","I"), n, prob = c(.5,.2,.3), replace = TRUE)
race <- sample(c("AI","AS","B","H", "W"), n, prob = c(.03,.05,.1,.32,.5), replace = TRUE)

sex <- as.numeric(gender == "F")
live1 <- as.numeric(livingarrangement == "D")
live2 <- as.numeric(livingarrangement == "I")

race1 <- as.numeric(race == "AI")
race2 <- as.numeric(race == "AS")
race3 <- as.numeric(race == "B")
race4 <- as.numeric(race == "H")

if(crossed) {
    town <- rcat(n, prob = rep(1/ntowns, ntowns))
    if(n == 1000)
        town <- (rep(1:ntowns, 5))[1:n]
    state <- rcat(n, prob = rep(1/nstates, nstates))
} else {
    town <- rcat(n, prob = rep(1/ntowns, ntowns))
    town2state <- sample(seq_len(nstates), ntowns, replace = TRUE)
    state <- town2state[town]
}

sigma_state <- 0.5
sigma_town <- 0.25

beta0 <- 0.1
beta_sex <- 0.2
beta_live <- c(-.1, .05)
beta_race <- c(-.2, .1, -.05, 0)

u1 <- rnorm(nstates, 0, sigma_state)
u2 <- rnorm(ntowns, 0, sigma_town)

eta <- beta0 + beta_sex*sex + beta_live[1]*live1 + beta_live[2]*live2 +
            beta_race[1]*race1 + beta_race[2]*race2 + beta_race[3]*race3 + beta_race[4]*race4 + u1[state] + u2[town]

p <- expit(eta)

y <- rbinom(n, 1, p)


code <- nimbleCode({
    for(j in 1:nstates)
        u1[j] ~ dnorm(0, sd = sigma_state)
    for(j in 1:ntowns)
        u2[j] ~ dnorm(0, sd = sigma_town)
    sigma_state ~ dhalfflat()
    sigma_town ~ dhalfflat()

    for(i in 1:n) {
        eta[i] <- beta0 + beta_sex*sex[i] + beta_live[1]*live1[i] + beta_live[2]*live2[i] +
            beta_race[1]*race1[i] + beta_race[2]*race2[i] + beta_race[3]*race3[i] + beta_race[4]*race4[i] +
            u1[state[i]] + u2[town[i]]
        y[i] ~ dbern(expit(eta[i]))
    }
    beta0 ~ dflat()
    beta_sex ~ dflat()
    for(j in 1:2)
        beta_live[j] ~ dflat()
    for(j in 1:4)
        beta_race[j] ~ dflat()
})

outerOptimUseAD <- TRUE  # To isolate behavior to inner optimization.

set.seed(1)
m <- nimbleModel(code, data = list(y=y), constants = list(ntowns = ntowns, nstates = nstates, n=n, state = state, town = town, sex = sex, race1=race1,race2=race2,race3=race3,race4=race4,live1=live1,live2=live2), inits = list(beta0 = 0, beta_sex = 0, beta_live = rep(0,2), beta_race = rep(0,4), sigma_state = 1, sigma_town = 1, u1 = rnorm(nstates), u2 = rnorm(ntowns)), calculate = FALSE, buildDerivs = TRUE)  # slow with n=1e6

cm <- compileNimble(m)

## This won't be necessary once we automatically exclude latents with no deps.
if(!crossed && n==10000) {
    u1nodes <- paste0('u1[', (1:nstates)[-c(9,28)], ']')
} else u1nodes <- 'u1'

## NIMBLE Laplace

system.time(mLaplace <- buildLaplace(m, randomEffectsNodes = c('u1','u2','beta0','beta_sex','beta_race','beta_live'), paramNodes = c('sigma_state','sigma_town'), control = list(outerOptimUseAD = TRUE))) # n=1e5: 29 sec. (32 sec. for nested); n=1e6: 350 sec.
system.time(cLaplace <- compileNimble(mLaplace, project = m))  # n=1e5: 170 sec.; 14 sec. second time?; 151 sec. for nested; n=1e6: 373 sec.
## system.time(cLaplace$findMLE(c(0,0)))
system.time(result <- runLaplace(cLaplace))
## n=1e5: 41 sec; 5 sec if done a second time; 28 sec. for nested; n=1e6: 448 sec.
## comparable to glmmTMB

system.time(mLaplace <- buildLaplace(m, randomEffectsNodes = c('u1','u2','beta0','beta_sex','beta_race','beta_live'), paramNodes = c('sigma_state','sigma_town'), control = list(outerOptimUseAD = FALSE)))
## TODO: for nested case (perhaps others) now getting compilation failure related to `outDir` arg in laplace.
system.time(cLaplace <- compileNimble(mLaplace, project = m))
## system.time(cLaplace$findMLE(c(0,0)))
system.time(result <- runLaplace(cLaplace))
## n=1e5: 45 sec; 11 sec if done a second time; 34 sec. for nested

system.time(mLaplace <- buildLaplace(m, randomEffectsNodes = c('u1','u2'), paramNodes = c('sigma_state','sigma_town','beta0','beta_sex','beta_race','beta_live'), control = list(outerOptimUseAD = TRUE)))
system.time(cLaplace <- compileNimble(mLaplace, project = m))  # 153 sec.
system.time(result <- runLaplace(cLaplace))  
## n=1e5: 68 sec; 24 sec second time
## very similar to glmmTMB; presumably this is the approach used in glmmTMB

## for nested, have individual Laplace approxs; buildLaplace is slow 5300 sec.
## results are quite similar to glmmTMB but I don't think as similar as with crossed


## glmmTMB Laplace

dat <- data.frame(y = y, gender = factor(gender, levels = c('M','F')), race = factor(race, levels=c("W","AI","AS","B","H")), livingarrangement=factor(livingarrangement, levels=c("H","D","I")), state = as.character(state), town = as.character(town))
system.time(result <- glmmTMB::glmmTMB(
    y ~ gender + livingarrangement + race + (1|state) + (1|town), # same as (1|state)+(1|state:town) for nested
    data = dat,
    family = binomial ,
    verbose = FALSE #,
    # sparseX = c('cond' = TRUE)
    ))
## glmmTMB::ranef(result)
## (for n=1e5) 26 s. w/ or w/o sparseX (much more than 9 s. for full stringer in paper...); little memory use
## additional time may be because of serial computation, but actually top shows ~300% cpu use
## `openmp()` not available (not sure what package it is in) so can't seem to do in parallel

## Crossed:
## n=1e5: 25 s. (elapsed; 41 s. user)
## numerical results for var comps and fixed effects comparable to NIMBLE with fixed in randomEffectsNodes;
## very similar to NIMBLE with fixed in params

## n=1e6: 226 s. (elapsed; 241 s. user)

## Nested
## n=1e5: 5 sec.

## nested approx

system.time(approx <- buildNestedApprox(m, latentNodes = c(u1nodes,'u2','beta0','beta_sex','beta_race','beta_live'), paramNodes = c('sigma_state','sigma_town'), control = list(outerOptimUseAD = outerOptimUseAD)))  
system.time(capprox <- compileNimble(approx, project = m)) 
system.time(result <- runNestedApprox(capprox)) 

system.time(latent_sample <- result$sampleLatents(1000)) 
system.time(result$improveParamMarginals(c("sigma_state","sigma_town"), nMarginalGrid = 5, nQuad = 3)) 

t(apply(latent_sample[,1:12], 2, quantile, qpts))

## Timing (build, compile run, sample, improve):
## 10000 crossed: 4, 154, 4, 2, 13 
## 10000 nested: 4, 199, 5,  2, 14
## 100000 crossed: 40, 163, 45, 9, 109 (redo?)
## 100000 nested: 42, 172, 33, 9, 94

## For crossed, n=10000: inference good for VCs (marginal improvement has limited effect) and FEs; pretty good for REs compared to HMC for n=10000

## HMC

library(nimbleHMC)

## Use noncentered for better HMC mixing.
code <- nimbleCode({
    for(j in 1:nstates)
        u1[j] ~ dnorm(0, sd = 1)
    for(j in 1:ntowns)
        u2[j] ~ dnorm(0, sd = 1)
    sigma_state ~ dhalfflat()
    sigma_town ~ dhalfflat()

    for(i in 1:n) {
        eta[i] <- beta0 + beta_sex*sex[i] + beta_live[1]*live1[i] + beta_live[2]*live2[i] +
            beta_race[1]*race1[i] + beta_race[2]*race2[i] + beta_race[3]*race3[i] + beta_race[4]*race4[i] +
            sigma_state*u1[state[i]] + sigma_town*u2[town[i]]
        y[i] ~ dbern(expit(eta[i]))
    }
    beta0 ~ dflat()
    beta_sex ~ dflat()
    for(j in 1:2)
        beta_live[j] ~ dflat()
    for(j in 1:4)
        beta_race[j] ~ dflat()
})

set.seed(1)
m <- nimbleModel(code, data = list(y=y), constants = list(ntowns = ntowns, nstates = nstates, n=n, state = state, town = town, sex = sex, race1=race1,race2=race2,race3=race3,race4=race4,live1=live1,live2=live2), inits = list(beta0 = 0, beta_sex = 0, beta_live = rep(0,2), beta_race = rep(0,4), sigma_state = 1, sigma_town = 1, u1 = rnorm(nstates), u2 = rnorm(ntowns)), calculate = TRUE, buildDerivs = TRUE)

cm <- compileNimble(m)
mcmc <- buildHMC(m)
cmcmc <- compileNimble(mcmc, project = m)

system.time(out <- runMCMC(cmcmc, niter = 21000, nburnin = 1000)) #  n=1e4: 1367 sec.

t(apply(out[,1:10], 2, quantile, qpts))

out[,11:57] <- out[,11:57]*out[,'sigma_state']
out[,58:319] <- out[,58:319]*out[,'sigma_town']
t(apply(out[,11:14], 2, quantile, qpts))
t(apply(out[,58:62], 2, quantile, qpts))

## INLA

library(INLA)
dat <- data.frame(y = y, gender = gender, race = race, livingarrangement=livingarrangement, state = as.character(state), town = as.character(town))

## Need to make priors consistent for comparison of statistical results.
system.time(inlamod <- tryCatch(inla(
      y ~ gender + race + livingarrangement + 
        f(state,model = "iid",hyper = list(prec = list(prior = "pc.prec",param = c(.5,.5)))) +
        f(town,model = "iid",hyper = list(prec = list(prior = "pc.prec",param = c(.5,.5)))),
      data = dat,
      family = 'binomial')))  #  7 sec (14 sec user) for n=10000,crossed; 4 sec (11 for user) for n=10000 nested
## That seems to run in parallel with inla.mkl executable.


## Try in parallel and time
system.time(inlamod <- tryCatch(inla(
      y ~ gender + race + livingarrangement + 
        f(state,model = "iid",hyper = list(prec = list(prior = "pc.prec",param = c(.5,.5)))) +
        f(town,model = "iid",hyper = list(prec = list(prior = "pc.prec",param = c(.5,.5)))),
      data = dat,
      family = 'binomial',
      control.compute = list(
        openmp.strategy = 'pardiso',
        smtp = 'pardiso'
      ),
      control.inla = list(
        strategy = 'gaussian',
        int.strategy = 'ccd'
      ))))
