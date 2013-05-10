source("06_S4_complete.R")

setGeneric("transcribe", function(object, ...) standardGeneric("transcribe"))
setMethod("transcribe","DnaSeq",
          function(object, ...) {
            .sequence <- chartr("T","U",toupper(seq(object)))
            .id <-  paste(id(object),"-- transcribed")
            rna <- new("RnaSeq",
                       id=.id,
                       alphabet=c("A","C","G","U"),
                       sequence=.sequence)
            if (validObject(rna))
              return(rna)
          })
          
## test code

x <- readFasta("aDnaSeq.fasta")

print(transcribe(x))
     
