#################################################################
## 11_CFAcountData.R                                           ##
## Exercise in writing scripts with multiple steps             ##
## Script assigns several factor objects.  Reads in 3 data     ##
## files.  Undertakes some analysis and generates results      ##
##                                                             ##
## Robert Stojnic. rs550@cam.ac.uk                             ## 
#################################################################
##
##
### START ###
###############################

########################################
## Part 1.  Build a list of all genes ##
########################################

# start out with an empty collection of genes
genes <- c() 
for(fileNum in 1:5){
	# load in files 13_DiffGenes1.tsv ...
	t <- read.delim(paste("13_DiffGenes", fileNum, ".tsv", sep=""), as.is=T, header=F)
	# label the input columns to help code readability
	names(t) <- c("gene", "expression")
	# 
	genes <- union(genes, t$gene)
}

# for tidiness order our genes by name
genes <- sort(genes)


###############################################
## Part 2.  Combine data into a single table ##
###############################################

# make the destination table
values <- matrix(0, nrow=length(genes), ncol=5)
# row are going to the complete set of genes
rownames(values) <- genes

for(fileNum in 1:5){
	# read in the file again
	t <- read.delim(paste("13_DiffGenes", fileNum, ".tsv", sep=""), as.is=T, header=F)
	names(t) <- c("gene", "expression")
	
	# match the names of the genes to the rows in our big table
	index <- match(t$gene, rownames(values))
	# use the matched indices to copy the expression levels
	values[index,fileNum] <- t$expression
}

# make the heat map with hierarchical clustering
#png("geneClusters.png")
heatmap(values, scale="none", col = cm.colors(256))
#dev.off()

##################################################
## Part 3.  Redo the heatmap on subset of genes ##
##################################################

# load in the experimentally verified genes
t.exp <- read.delim("13_ExperimentalGenes.tsv", as.is=T)
# split all gene names by "," and then flatten it out into a single vector
experim.genes <- unlist( strsplit(t.exp$genes, ",") )

# redo the heatmap by using just the genes in the experimentally verified set
is.experimental <- rownames(values) %in% experim.genes
#png("geneClusterExp.png")
x11() # make a new window
heatmap(values[ is.experimental, ], scale="none", col = cm.colors(256))
#dev.off()


