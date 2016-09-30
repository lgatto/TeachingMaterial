
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
setGeneric("id", function(object) standardGeneric("id"))
setGeneric("id<-", function(object,value) standardGeneric("id<-"))


# methods
setMethod("rev", "GenericSeq",
          function(x) paste(rev(unlist(strsplit(x@sequence, ""))), collapse=""))          
setMethod("id", "GenericSeq", function(object) object@id)
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
