##' Reads sequences data in fasta and create \code{DnaSeq}
##' and \code{RnaSeq} instances. 
##'
##' This funtion reads DNA and RNA fasta files and generates
##' valid \code{"DnaSeq"} and \code{"RnaSeq"} instances.
##' 
##' @title Read fasta files.
##' @param infile  the name of the fasta file which the data are to be read from.
##' @return an instance of \code{DnaSeq} or \code{RnaSeq}.
##' @seealso \code{\linkS4class{GenericSeq}}, \code{\linkS4class{DnaSeq}} and \code{\linkS4class{RnaSeq}}.
##' @examples
##' f <- dir(system.file("extdata",package="sequences"),pattern="fasta",full.names=TRUE)
##' f
##' aa <- readFasta(f[1])
##' aa
##' @author Laurent Gatto \email{lg390@@cam.ac.uk}
##' @keywords IO, file
readFasta <- function(infile){
  lines <- readLines(infile)
  header <- grep("^>", lines)
  if (length(header)>1) {
    warning("Reading first sequence only.")
    lines <- lines[header[1]:(header[2]-1)]
    header <- header[1]
  }
  .id <- sub("^> *","",lines[header],perl=TRUE)
  .sequence <- toupper(paste(lines[(header+1):length(lines)],collapse=""))
  .alphabet <- toupper(unique(strsplit(.sequence,"")[[1]]))
  if (all(.alphabet %in% c("A","C","G","T"))) {
    newseq <- new("DnaSeq",
                  id=.id,
                  sequence=.sequence)
  } else if (all(.alphabet %in% c("A","C","G","U"))) {
    newseq <- new("RnaSeq",
                  id=.id,
                  sequence=.sequence)
  } else {
    stop("Alphabet ",.alphabet," is unknown.")
  }
  if (validObject(newseq))
    return(newseq)
}


gccount <- function(inseq) {
  .Call("gccount",
        inseq,
        PACKAGE="sequences")
}

gccount2 <- function(inseq) {
  .Call("gccount2",
        inseq,
        PACKAGE="sequences")
}


debugme <- function() {
    f <- list.files(system.file("scripts", package = "sequences"),
                    full.names = TRUE)
    source(f)
}
