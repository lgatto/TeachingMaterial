## simulator.R --- Generate some samples from a distribution
## The first argument is the name of the function to generate random samples
## (e.g. "rnorm").
## The number of samples to generate is taken from the file params.R

args <- commandArgs(TRUE)

source('params.R')                      #should define n
dat <- do.call(args[1], list(n=n))
write(dat, file='')                     #go to stdout


