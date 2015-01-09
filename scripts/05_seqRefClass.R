Seq <- setRefClass("Seq",
                   fields = list(
                     id = "character",
                     alphabet = "character",
                     sequence = "character"),
                   methods = list(
                     generateAlphabet = function() {
                       'Generate the alphabet.'
                       alphabet <<- unique(strsplit(.self$sequence, "")[[1]])
                       if (!.self$valid())
                         stop("Object is not valid")
                     } ,
                     valid = function() {
                       'Checks validity of object.'
                       chars <- casefold(unique(unlist(strsplit(.self$sequence,""))))
                       isValid <- all(chars %in% casefold(.self$alphabet))
                       return(isValid)
                     },
                     show = function() {
                       'Method for automatically printing matrix editors'
                       cat("Reference class",
                           class(.self), "\n")
                       cat("Id: ", .self$id, "\n", sep="")
                       cat("Length: ", nchar(.self$sequence), "\n", sep="")
                       cat("Alphabet: ", .self$alphabet, "\n", sep="")
                       cat("Sequence: ", .self$sequence, "\n", sep="")
                     },
                     rev = function() {
                       'Reverses the sequence.'
                       sequence <<- paste(base::rev(strsplit(.self$sequence,"")[[1]]), collapse="")
                       id <<- paste(.self$id, "-- reversed")
                     },
                     comp = function() {
                       'Complements the (DNA) sequence'
                       sequence <<- chartr("ACGT","TGCA",.self$sequence)
                       id <<- paste(.self$id, "-- complemented")                     
                     },
                     seq = function() {
                       'Returns the object\'s sequence'
                       return(.self$sequence)
                     },
                     transcribe = function() {
                       'Transcribes the DNA sequence into RNA.'
                       sequence <<- chartr("T", "U", toupper(.self$sequence))
                       id <<- paste(.self$id, "-- transcribed")
                       .self$generateAlphabet()
                     }
                     ))
                     
#### Test code

## Seq is a generator object for the "Seq" class
Seq$fields()
Seq$methods()
## Inline help
Seq$help("rev")
## Making new objects
s <- Seq$new(id="foo", sequence="GATCATCA")
s
s$valid()
s$generateAlphabet()
s
s$valid()

## Changing object and accessing fields
s$rev()
s
s$seq()
s$sequence
s$field("sequence")

## Assignment
s2 <- s
s$seq()
s2$seq()
s$sequence <- "AAAAAAAAAAA"
s$seq()
s2$seq()
s2 <- s$copy()
s$sequence <- "GCAGCGCGCATC"
s$seq()
s2$seq()

s$sequence
s$field("sequence", "TTT")
s$sequence
s$sequence <- "GACTCATCA"
s$sequence

## More complex functions
s$comp()
s
s$transcribe()
s                

