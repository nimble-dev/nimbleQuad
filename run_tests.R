#!/usr/bin/env Rscript

library(methods)
library(testthat)
library(nimble)
library(nimbleQuad)

testFiles <- "nimbleQuad/tests/testthat/test-nested-interface-partial.R"

for(test in testFiles) {
    cat('===========================================================\n')
    cat(paste0('Running test-', gsub('.*test-', '', test), '\n'))
    cat('===========================================================\n')
    source(test)
}

