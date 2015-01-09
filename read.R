#' Reads sequences data in fasta and create \code{GenericSeq}
#'
#' This funtion reads DNA and RNA fasta files and generates
#' valid \code{"GenericSeq"} instances.
#'
#' @title Read fasta files.
#' @param infile the name of the fasta file which the data are to be read from.
#' @return an instance of \code{GenericSeq}.
#' @seealso \code{\linkS4class{GenericSeq}}
#' @examples
#' f <- dir(system.file("extdata",package="sequences"), pattern="fasta",full.names=TRUE)
#' f
#' aa <- readFasta(f)
#' aa
#' @author Laurent Gatto \email{lg390@@cam.ac.uk}
#' @export
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
  
  newseq <- new("GenericSeq", id=.id, alphabet=.alphabet,
                  sequence=.sequence)

  return(newseq)
}
