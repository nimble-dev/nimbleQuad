#!/usr/bin/env Rscript

requirements <- c(
    'igraph',
    'coda',
    'R6',
    'testthat',
    'pracma',     ## for AD
    'numDeriv',    ## for AD
    'mvQuad',
    'RTMB',
    'polynom'
    ## 'lme4'    
    )

for(package in requirements) {
    install.packages(package)
}

## We need NIMBLE's tests installed for `{AD_,}test_utils.R`.
install.packages('nimble', INSTALL_opts = '--install-tests')

## Apparently a bug in Matrix (as of early 2024) is causing an issue (https://bioconductor.org/packages/devel/bioc/vignettes/dreamlet/inst/doc/errors.html) that is causing Laplace test failures when fitting a model with lmer.
install.packages('lme4', type = 'source')


