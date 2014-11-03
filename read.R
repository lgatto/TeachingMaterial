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
