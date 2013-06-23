
## @knitr env, echo=FALSE, massage = FALSE
suppressPackageStartupMessages(library(Biobase))
data(sample.ExpressionSet)


## @knitr plotcode, eval = FALSE, tidy = FALSE
## plot(exprs(sample.ExpressionSet[, 1]),
##      exprs(sample.ExpressionSet[, 2]),
##      log = "xy",
##      xlab = sampleNames(sample.ExpressionSet)[1],
##      ylab = sampleNames(sample.ExpressionSet)[2])
## abline(0, 1)
## grid()


## @knitr label=plotfig1, echo=FALSE, fig.width=5, fig.height=4, tidy=FALSE, warning = FALSE
plot(exprs(sample.ExpressionSet[, 1]), 
     exprs(sample.ExpressionSet[, 2]), 
     log = "xy", 
     xlab = sampleNames(sample.ExpressionSet)[1],
     ylab = sampleNames(sample.ExpressionSet)[2])
abline(0, 1)
grid()


## @knitr plotcode0, eval = FALSE, tidy = FALSE
## pairs(log2(exprs(sample.ExpressionSet)[, 1:4]),
##       pch = 19,
##       col = "#0000FF20")


## @knitr label=plotfig2,echo=FALSE,fig.width=3.5,fig.height=3.5,tidy=FALSE, warning = FALSE
pairs(log2(exprs(sample.ExpressionSet)[, 1:3]),
      pch = 19,
      col = "#0000FF20")


## @knitr plotcode1, eval = FALSE
## boxplot(log2(exprs(sample.ExpressionSet)))


## @knitr label=plotfig3,echo=FALSE,fig.width=5,fig.height=4,tidy=FALSE
boxplot(log2(exprs(sample.ExpressionSet)))


## @knitr label=plotfig4,echo=TRUE,fig.width=5,fig.height=4,tidy=FALSE, warning = FALSE, message = FALSE
smoothScatter(log2(exprs(sample.ExpressionSet)[, 1:2]))


