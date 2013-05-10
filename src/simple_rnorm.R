#!/usr/bin/env Rscript
args <- commandArgs(TRUE)

stopifnot(length(args)==3)
args = as.numeric(args)
n = args[1]; mean = args[2]; sd = args[3]
rnorm(n, mean, sd)

