## Draw a histogram of data stored in file.
## Assumes that the file ends in .dat
args <- commandArgs(TRUE)

infile <- args[1]
outfile <- sub('.dat$', '.png', infile)

data <- scan(infile, quiet=TRUE)
main <- sprintf('Histogram of %s (%d samples)',
                infile, length(data))
png(file=outfile)
par(las=1, bty='n')
hist(data, main=main)
invisible(dev.off())



