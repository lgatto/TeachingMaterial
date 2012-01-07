setMethod("comp","DnaSeq",
          function(object) {
            chartr("ACGT","TGCA",seq(object))
          })


setMethod("transcribe","DnaSeq",
          function(object) {
            .sequence <- chartr("T","U",toupper(seq(object)))
            .id <-  paste(id(object),"-- transcribed")
            rna <- new("RnaSeq",
                       id=.id,
                       alphabet=c("A","C","G","U"),
                       sequence=.sequence)
            return(rna)
          })

