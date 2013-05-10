#### functions

# Read in a single sequence from a FASTA file in an S4 object
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
  return(newseq)
}

#### S4 class definition
setClass("GenericSeq",
         representation = representation(
           id = "character",
           alphabet = "character",
           sequence =  "character",
           "VIRTUAL"),
         validity = function(object) {
           isValid <- TRUE
           if (nchar(object@sequence)>0) {
             chars <- casefold(unique(unlist(strsplit(object@sequence,""))))
             isValid <- all(chars %in% casefold(object@alphabet))
           }
           if (!isValid)
             cat("Some characters are not defined in the alphabet.\n")
           return(isValid)
         })

setClass("DnaSeq",
         contains="GenericSeq",
         prototype = prototype(
           id = paste("my DNA sequence",date()),
           alphabet = c("A","C","G","T"),
           sequence = character())
         )

setClass("RnaSeq",
         contains="GenericSeq",
         prototype = prototype(
           id = paste("my RNA sequence",date()),
           alphabet = c("A","C","G","U"),
           sequence = character())
         )          

### S4 methods

# generics

setGeneric("rev", function(x) standardGeneric("rev"))
setGeneric("id", function(object) standardGeneric("id"))
setGeneric("id<-", function(object,value) standardGeneric("id<-"))
setGeneric("alphabet", function(object) standardGeneric("alphabet"))

## There is already a S3 'seq' method (see ?seq),
## although not a generic one (see isGeneric(seq))
setGeneric("seq", function(...) standardGeneric("seq"))
setGeneric("seq<-", function(object,value) standardGeneric("seq<-"))
setGeneric("print", function(x,...) standardGeneric("print"))
setGeneric("comp",function(object,...) standardGeneric("comp"))

# methods
setMethod("rev", "GenericSeq",
          function(x) paste(rev(unlist(strsplit(x@sequence, ""))), collapse=""))          
setMethod("id", "GenericSeq", function(object) object@id)
                 
setMethod("alphabet", "GenericSeq", function(object) object@alphabet)
setMethod("length", "GenericSeq", function(x) nchar(x@sequence))
setMethod("seq", "GenericSeq", function(object,...) object@sequence)

setReplaceMethod("id",
                 signature(object="GenericSeq",
                           value="character"),
                 function(object, value) {
                   object@id <- value
                   if (validObject(object))
                     return(object)
                 })
                 
setReplaceMethod("seq",
                 signature(object="GenericSeq",
                           value="character"),
                 function(object, value) {
                   object@sequence <- value
                   if (validObject(object))
                     return(object)
                 }) 
                                                
setMethod("show",
          "GenericSeq",
          function(object) {
            cat("Object of class",class(object),"\n")
            cat(" Id:",id(object),"\n")
            cat(" Length:",length(object),"\n")
            cat(" Alphabet:",alphabet(object),"\n")
            cat(" Sequence:",seq(object), "\n")
          })
          

setMethod("print", "GenericSeq",
          function(x) {
            sq <- strsplit(seq(x),"")[[1]]
            cat(">",id(x),"\n"," 1\t")
            for (i in 1:length(x)) {
              if ((i %% 10)==0) {
                cat("\n",i,"\t")
              }
              cat(sq[i])
            }
            cat("\n")
          })                           

setMethod("initialize", "GenericSeq",
    function(.Object, ..., id="", sequence=""){
        .Object@id <- id
        .Object@sequence <- toupper(sequence)
        .Object@alphabet <- unique(strsplit(.Object@sequence, "")[[1]])
##        callNextMethod(.Object, ...) # call parent class initialize()
        .Object
    })	
    
    
setMethod("comp","DnaSeq",
          function(object, ...) {
            chartr("ACGT","TGCA",seq(object))
          })
          
setMethod("comp","RnaSeq",
          function(object, ....) {
            chartr("ACGU","UGCA",seq(object))
          })    

## test code

