#!/usr/bin/env Rscript

library(methods)
library(testthat)
library(nimble)
library(nimbleQuad)

testFiles <-
    grep('test-.+\\.R$',
         list.files('nimbleQuad/tests/testthat', full.names = TRUE),
         value = TRUE)

## Put test-AGHQ last to see if that resolves failure on Windows.
testFiles <- testFiles[c(2:length(testFiles),1)]

## See issues 65 and 66 for strange error preventing running some
## Laplace tests after earlier tests.
#if(Sys.info()['sysname'] == "Windows")
#    testFiles <- testFiles[!grepl("laplace[2-9]", testFiles)]

for(test in testFiles) {
    cat('===========================================================\n')
    cat(paste0('Running test-', gsub('.*test-', '', test), '\n'))
    cat('===========================================================\n')
    source(test)
}

