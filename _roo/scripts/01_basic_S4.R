
# Function to read in data
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

#### S4 class definition
setClass("GenericSeq",
         representation = representation(
           id = "character",
           alphabet = "character",
           sequence =  "character"
          ))
          

### S4 methods

# generics

setGeneric("rev", function(x) standardGeneric("rev"))
setGeneric("id", function(object, ...) standardGeneric("id"))
setGeneric("id<-", function(object,value) standardGeneric("id<-"))


# methods
setMethod("rev", "GenericSeq",
          function(x) paste(rev(unlist(strsplit(x@sequence, ""))), collapse=""))          
setMethod("id", "GenericSeq", function(object, ...) object@id)
setReplaceMethod("id", signature(object="GenericSeq",
                           value="character"),
                 function(object, value) {
                   object@id <- value
                   return(object)
                 })

## test code
genseq <- new("GenericSeq", id="sequence name", 
       alphabet=c("A", "C", "G", "T"), sequence="AGATACCCCGAAACGA")

id(genseq) <- "new sequence name"
id(genseq)          
genseq
genseq@id
slot(genseq, "id")
showMethods("rev")
