

fls <- dir(pattern = "Exercise-05-DE")
length(fls)
fls



l <- lapply(fls, read.delim, row.names = 1)
length(l)
class(l)



allgenes <- lapply(l, rownames)
allgenes <- unique(unlist(allgenes))
length(allgenes)
head(allgenes)



eset <- matrix(0, nrow = length(allgenes), ncol = length(l))
rownames(eset) <- allgenes



for (i in 1:length(l)) {
  ## find the row indices
  idx <- match(rownames(l[[i]]), rownames(eset))
  eset[idx, i] <- l[[i]][, 1]
}



heatmap(eset)



tab <- read.delim("Exercise-05-table.tsv", stringsAsFactors = FALSE)
tab



genes <- unlist( strsplit(tab$genes, ",") )
genes



all(genes %in% rownames(eset))



heatmap(eset[genes, ])


