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

#### define S3 methods

# generics
id <- function(x){ UseMethod("id") }

# methods
id.GenericSeq <- function(x){ x$id } 
seq.GenericSeq = function(x){ x$seq }

##### Test code
s <- readFasta("aDnaSeq.fasta")
print(id(s))
print(seq(s))
