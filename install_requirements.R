#!/usr/bin/env Rscript

requirements <- c(
    'igraph',
    'coda',
    'R6',
    ##'nimble',      ## add back in later
    'testthat',
    'pracma',     ## for AD
    'numDeriv',    ## for AD
    'mvQuad',
    'RTMB'
    ## 'lme4'     ## for test-ADlaplace.R
    )

for(package in requirements) {
    install.packages(package)
}

## Apparently a bug in Matrix (as of early 2024) is causing an issue (https://bioconductor.org/packages/devel/bioc/vignettes/dreamlet/inst/doc/errors.html) that is causing test-ADlaplace.R failures when fitting a model with lmer.
install.packages('lme4', type = 'source')


## later, remove below:
install.packages('devtools')
library(devtools)
devtools::install_github('nimble-dev/nimble', ref = 'devel', subdir = 'packages/nimble', INSTALL_opts = c("--install-tests"))
## remove until here

