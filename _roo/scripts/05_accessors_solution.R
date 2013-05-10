source("04_basic_S4.R")

#### S4 methods

# generics
setGeneric("alphabet", function(object) standardGeneric("alphabet"))

## There is already a 'seq' method (see ?seq),
## although not a generic one (see isGeneric(seq))
setGeneric("seq", function(...) standardGeneric("seq"))
setGeneric("seq<-", function(object,value) standardGeneric("seq<-"))


# methods
setMethod("alphabet", "GenericSeq", function(object) object@alphabet)
setMethod("length", "GenericSeq", function(x) nchar(x@sequence))
setMethod("seq", "GenericSeq", function(object,...) object@sequence)

setReplaceMethod("seq",
                 signature(object="GenericSeq",
                           value="character"),
                 function(object, value) {
                   object@sequence <- value
                   return(object)
                 })

### test code

print( seq(genseq) )
print( alphabet(genseq) )
print( length(genseq) )

seq(genseq) <- "AAAATTT"
print( seq(genseq) )



