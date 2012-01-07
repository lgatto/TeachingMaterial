setClass("GenericSeq",
         representation = representation(
           id = "character",
           alphabet = "character",
           sequence =  "character",
           "VIRTUAL"),
         validity = function(object) {
           isValid <- TRUE
           if (length(object@sequence)>0) {
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
