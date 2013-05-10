library("Biobase")
library("affy")


dir()
flnms <- dir(pattern = "*CEL")

rawdata <- ReadAffy(filenames = flnms)
class(rawdata)
rawdata

boxplot(rawdata)
hist(rawdata)
par(mfrow = c(1, 2))
image(rawdata[, 1:2])
dev.off()

deg <- AffyRNAdeg(rawdata)
class(deg)
names(deg)
plotAffyRNAdeg(deg)

eset <- rma(rawdata)
class(eset)
eset
head(exprs(eset))

boxplot(exprs(eset))

ctrl <- grep("AFFX", featureNames(eset))
heatmap(exprs(eset[ctrl, ]))

pca <- prcomp(t(exprs(eset)))
biplot(pca)
plot(pca$x[, 1:2], pch = 19, col = "#AB000030", cex = 3)
text(pca$x[, 1:2], labels = 1:5)
grid()

## save(eset, file = "eset.rda")

library("arrayQualityMetrics")
?arrayQualityMetrics
vignette(package = "arrayQualityMetrics")
vignette("arrayQualityMetrics", package = "arrayQualityMetrics")
arrayQualityMetrics(eset, outdir = "aqmReport")

