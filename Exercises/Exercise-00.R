fls <- dir(pattern = "Exercise-05-DE")

l <- lapply(fls, read.delim, row.names = 1)
class(l)
length(l)
class(l[[1]])
l[[1]]

allgenes <- lapply(l, rownames)
allgenes
allgenes <- unique(unlist(allgenes))
allgenes

eset <- matrix(0, nrow = length(allgenes), ncol = length(l))
rownames(eset) <- allgenes

for (i in 1:length(l)) {
  idx <- match(rownames(l[[i]]), rownames(eset))
  eset[idx, i] <- l[[i]][, 1]
}

heatmap(eset)

tab <- read.delim("Exercise-05-table.tsv",
                  stringsAsFactors = FALSE)

genes <- unlist( strsplit(tab$genes, ",") )

all(genes %in% rownames(eset))

heatmap(eset[genes, ])
