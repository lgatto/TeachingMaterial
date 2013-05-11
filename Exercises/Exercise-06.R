

library("Biobase")
library("affy")
library("AnnotationDbi")



dir()
flnms <- dir(pattern = "*CEL")
flnms



rawdata <- ReadAffy(filenames = flnms)
class(rawdata)



par(mfrow = c(2, 2))
image(rawdata[, 1:4])



boxplot(log2(exprs(rawdata)))



deg <- AffyRNAdeg(rawdata)
plotAffyRNAdeg(deg)



eset <- rma(rawdata)
class(eset)
eset
head(exprs(eset))



boxplot(exprs(eset))



ctrl <- grep("AFFX", featureNames(eset))
heatmap(exprs(eset[ctrl, ]))



pca <- prcomp(t(exprs(eset)))
plot(pca$x[, 1:2], pch = 19, col = "#AB000030", cex = 3)
text(pca$x[, 1:2], labels = 1:5)
grid()



## save(eset, file = "eset.rda")



## library("arrayQualityMetrics")
## arrayQualityMetrics(eset, outdir = "aqmReport")



## vignette(package = "arrayQualityMetrics")
## vignette("arrayQualityMetrics", package = "arrayQualityMetrics")


