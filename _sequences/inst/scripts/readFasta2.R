##' Read a FASTA file, possibly with multiple sequences
##'
##' @param infile the input filename
##' @return a \code{list} of sequences read from the filename, each
##' sequence is an instance of class \code{DnaSeq}.
##' @examples
##' readFasta2("moreDnaSeqs.fasta")
readFasta2 <- function(infile) {
  raw <- readLines(infile)
  header.inx <- grep("^>", raw)
  ## find out how to group lines together
  if (length(header.inx) == 1) {
    seq.groups = rep(1, length(raw))
  } else {
    ## number of lines each distinct fasta entry takes up
    seq.groups.rep <- c(header.inx[2:length(header.inx)], length(raw)+1) - header.inx
    seq.groups <- rep(1:length(header.inx), times=seq.groups.rep)
  }
  
  fasta.headers <- raw[header.inx]
  ## concatenete sequence lines
  fasta.seqs <- as.vector(by(raw, seq.groups , function(x) paste(x, collapse="" )))
  ## extract only the id part
  fasta.ids <- sub("^> *", "", fasta.headers)
	
  ## only one type of sequence supported in single file
  alphabet <- paste(sort(unique(unlist(strsplit(fasta.seqs, "")))), collapse="")
  
  out <- vector("list", length=length(fasta.ids))
  
  for (i in 1:length(fasta.ids))
      out[[i]] <- new("DnaSeq",
                      id = fasta.ids[i],
                      alphabet = alphabet,
                      seq = fasta.seqs[i])  
  return(out)
}
