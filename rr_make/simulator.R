## Generate some samples from a distribution
args <- commandArgs(TRUE)

source('params.R')                      #should define n
dat <- do.call(args[1], list(n=n))
write(dat, file='')                     #go to stdout


