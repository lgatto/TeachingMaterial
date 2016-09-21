
source("02_GenericSeq.R")


#### S3 methods

# generics
alphabet <- function(x) UseMethod("alphabet")

# methods
alphabet.GenericSeq <- function(x) x$alphabet
length.GenericSeq <- function(x) nchar(seq(x))
rev.GenericSeq <- function(x) paste(rev(unlist(strsplit(x$seq, ""))), collapse="")

##### Test code
print(alphabet(s))
print(length(s))
print(rev(s))
