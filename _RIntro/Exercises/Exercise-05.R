#################################
## (1) Combine expression data ##
#################################

## a vector for the file names of interest
fls <- dir(pattern = "Exercise-05-DE")

## iterate over file names, import as data.frame
## and return as a list (of data.frames)
l <- lapply(fls, read.delim, row.names = 1)

## obtain a vector of unique gene names
allgenes <- lapply(l, rownames)
allgenes <- unique(unlist(allgenes))

## initialise matrix
eset <- matrix(0, nrow = length(allgenes), ncol = length(l))
rownames(eset) <- allgenes

## populate matrix
for (i in 1:length(l)) {
  ## find the row indices
  idx <- match(rownames(l[[i]]), rownames(eset))
  eset[idx, i] <- l[[i]][, 1]
}

heatmap(eset)


##################################################
## (2) Extract and visualise genes of interest  ##
##################################################

tab <- read.delim("Exercise-05-table.tsv", stringsAsFactors = FALSE)

## split the comma separated genes names 
## store as a character vector 
genes <- unlist( strsplit(tab$genes, ",") )

## check that all gene names are in the matrix
all(genes %in% rownames(eset))

heatmap(eset[genes, ])
