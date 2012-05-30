setMethod("show",
          "GenericSeq",
          function(object) {
            cat("Object of class",class(object),"\n")
            cat(" Id:",id(object),"\n")
            cat(" Length:",length(object),"\n")
            cat(" Alphabet:",alphabet(object),"\n")
            cat(" Sequence:",seq(object), "\n")
          })


setMethod("print",
          "GenericSeq",
          function(x) {
            sq <- strsplit(seq(x),"")[[1]]
            cat(">",id(x),"\n")
            cat(" 1   ")
            for (i in 1:length(x)) {
              if ((i %% 10)==0) {
                cat("\n")
                cat(i,"  ")
              }
              cat(sq[i])
            }
            cat("\n")
          })


setMethod("id", "GenericSeq", function(object) object@id)
setMethod("id<-", "GenericSeq",
          function(object,value) object@id <- value)
setReplaceMethod("id",
                 signature(object="GenericSeq",
                           value="character"),
                 function(object, value) {
                   object@id <- value
                   if (validObject(object))
                     return(object)
                 })


setMethod("alphabet", "GenericSeq", function(object) object@alphabet)
setMethod("length", "GenericSeq", function(x) nchar(x@sequence))
setMethod("seq", "GenericSeq", function(object,...) object@sequence)

setMethod("seq<-", "GenericSeq",
          function(object,value) object@sequence <- toupper(value))

setReplaceMethod("seq",
                 signature(object="GenericSeq",
                           value="character"),
                 function(object, value) {
                   object@sequence <- value
                   if (validObject(object))
                     return(object)
                 })

setMethod("rev","GenericSeq",
          function(x) paste(rev(strsplit(seq(x),"")[[1]]),collapse=""))
          

## this is only an example of initialize function,
## note Generic is in fact virtual          
setMethod("initialize", "GenericSeq",
          function(.Object, ..., id="", sequence=""){
            .Object@id <- id
            .Object@sequence <- toupper(sequence)
            callNextMethod(.Object, ...)
          })
	

setMethod("[","GenericSeq",
          function(x,i,j="missing",drop="missing") {
            if (any(i > length(x)))
              stop("subscript out of bounds")
            s <- seq(x)
            s <- paste(strsplit(s,"")[[1]][i], collapse="")
            x@sequence <- s
            if (validObject(x))
              return(x)
          })
