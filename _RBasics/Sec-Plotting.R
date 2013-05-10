### R code from vignette source 'Sec-Plotting.Rnw'

###################################################
### code chunk number 1: env
###################################################
library(Biobase)
data(sample.ExpressionSet)


###################################################
### code chunk number 2: plotcode
###################################################
plot(exprs(sample.ExpressionSet[, 1]), 
     exprs(sample.ExpressionSet[, 2]), 
     log = "xy", 
     xlab = sampleNames(sample.ExpressionSet)[1],
     ylab = sampleNames(sample.ExpressionSet)[2])
abline(0, 1)
grid()


###################################################
### code chunk number 3: plotfig1
###################################################
plot(exprs(sample.ExpressionSet[, 1]), 
     exprs(sample.ExpressionSet[, 2]), 
     log = "xy", 
     xlab = sampleNames(sample.ExpressionSet)[1],
     ylab = sampleNames(sample.ExpressionSet)[2])
abline(0, 1)
grid()


###################################################
### code chunk number 4: plotcode
###################################################
pairs(log2(exprs(sample.ExpressionSet)[, 1:4]),
      pch = 19,
      col = "#0000FF20")


###################################################
### code chunk number 5: plotfig2
###################################################
pairs(log2(exprs(sample.ExpressionSet)[, 1:4]),
      pch = 19,
      col = "#0000FF20")


###################################################
### code chunk number 6: plotcode
###################################################
boxplot(log2(exprs(sample.ExpressionSet)))


###################################################
### code chunk number 7: plotfig3
###################################################
boxplot(log2(exprs(sample.ExpressionSet)))


###################################################
### code chunk number 8: plotcode
###################################################
smoothScatter(log2(exprs(sample.ExpressionSet)[, 1:2]))


###################################################
### code chunk number 9: plotfig4
###################################################
smoothScatter(log2(exprs(sample.ExpressionSet)[, 1:2]))


