# Load in a single sequence from a fasta file
readFasta <- function(infile){
  lines <- readLines(infile)
  header <- grep("^>", lines)
 
  if ( length(header) > 1 ) {
    warning("Reading first sequence only.")
    lines <- lines[header[1]:(header[2]-1)]
    header <- header[1]
  }
  
  # extract the first sequence
  id <- sub("^> *","",lines[header],perl=TRUE)
  sequence <- toupper(paste(lines[(header+1):length(lines)],collapse=""))
  alphabet <- unique(strsplit(sequence,"")[[1]])

  # make the return object and assign class  
  return.value <- list(id=id, sequence=sequence, alphabet=alphabet)
  class(return.value) <- "GenericSeq"
  
  return.value
}

# based on alphabet set the appropriate classes for a sequence object
setSeqSubtype <- function(s){
  if (all( alphabet(s) %in% c("A","C","G","T") )) {
    class(s) <- c("DnaSeq", "GenericSeq")  
  } else if (all( alphabet(s) %in% c("A","C","G","U") )) {
    class(s) <- c("RnaSeq", "GenericSeq")
  } else {
    stop("Alphabet ", alphabet(s) ," is unknown.")
  }
  
  return(s)
}


#### define S3 methods

# generics
id <- function(x){ UseMethod("id") }
alphabet <- function(x) UseMethod("alphabet")
comp <- function(x){ UseMethod("comp") }

# methods
id.GenericSeq <- function(x){ x$id } 
seq.GenericSeq = function(x){ x$seq }
alphabet.GenericSeq <- function(x) x$alphabet
length.GenericSeq <- function(x) nchar(seq(x))
rev.GenericSeq <- function(x) paste(rev(unlist(strsplit(x$seq, ""))), collapse="")
comp.DnaSeq = function(x) chartr("ACGT", "TGCA", seq(x))

##### Test code
s <- readFasta("aDnaSeq.fasta")
s.dna <- setSeqSubtype(s)

print( comp(s.dna) )

