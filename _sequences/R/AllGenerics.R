setGeneric("id", function(object) standardGeneric("id"))
setGeneric("id<-", function(object,value) standardGeneric("id<-"))
setGeneric("alphabet", function(object) standardGeneric("alphabet"))

## There is already a 'seq' method (see ?seq),
## although not a generic one (see isGeneric(seq))
setGeneric("seq", function(...) standardGeneric("seq"))
setGeneric("seq<-", function(object,value) standardGeneric("seq<-"))
## Same note that above for print
setGeneric("print", function(x,...) standardGeneric("print"))

setGeneric("rev",function(x) standardGeneric("rev"))
setGeneric("comp",function(object) standardGeneric("comp"))

setGeneric("transcribe", function(object) standardGeneric("transcribe"))
