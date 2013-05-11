

exp <- read.csv("MAdata1.csv", row.names = 1)
head(exp)



class(exp)
exp <- as.matrix(exp)
class(exp)



smeta <- read.csv("smeta1.csv", row.names = 1)
fmeta <- read.csv("fmeta1.csv", row.names = 1)
class(smeta)



marray <- list(expression = exp,
               featuremeta = fmeta,
               samplemeta = smeta)
str(marray)



save(marray, file = "marray.rda")
rm(list = ls()) ## clean working environment
load("marray.rda")
str(marray)



names(marray$featuremeta)
de <- marray$featuremeta[, "bh"] < 0.05
class(de)
table(de)



boxplot(marray$expression)



plot(marray$expression[, 1],
     marray$expression[, 4])



exp <- marray$expression
plot(exp[, 1], exp[, 4],
     xlab = colnames(exp)[1],
     ylab = colnames(exp)[4])
title(main = "Comparing expression levels")
grid()
points(exp[de, 1], exp[de, 4], col = "red", pch = 19)
abline(0, 1)



## identify(exp[, 1], exp[, 4], labels = fmeta$genes)



pairs(exp)



smoothScatter(exp[, 1], exp[, 4])



hist(exp[, 1])



par(mfrow = c(3, 2))
for (i in 1:6) {
  hist(exp[, i],
       xlab = marray$samplemeta[i, "sample"],
       main = "Histrogram of intensities")
  rug(exp[, i])
}



colnames(exp) <- marray$samplemeta$sample
rownames(exp) <- marray$featuremeta$genes



heatmap(exp)
heatmap(exp[de, ])



dir()
dir(pattern = "fmeta")



for (fmetafile in dir(pattern = "fmeta")) {
  fmeta <- read.csv(fmetafile, row.names = 1)
  print(dim(fmeta))
}



for (i in 1:3) {
  fmetafile <- paste0("fmeta", i, ".csv")
  cat(fmetafile, ":\n")
  fmeta <- read.csv(fmetafile, row.names = 1)
  de <- which(fmeta$bh < 0.05)
  if (length(de) > 0) {
    expfile <- paste0("MAdata", i, ".csv")
    exp <- read.csv(expfile, row.names = 1)
    exp <- as.matrix(exp)
    pdffile <- paste0("heatmap", i, ".pdf")
	cat(length(de), "DE genes.\n")
	cat("Saved", pdffile, ".\n")
    pdf(pdffile)
    heatmap(exp[de, ])
    dev.off()
  } else {
    cat("No DE found.\n")
  }
  cat("\n")
}


