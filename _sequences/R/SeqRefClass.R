Seq <- setRefClass("Seq",
                   fields = list(
                     id = "character",
                     alphabet = "character",
                     sequence = "character"),
                   methods = list(
                     ## initialize = function(id , sequence, alphabet=character()) {
                     ##   'Initialises the \'Seq\' instance.'
                     ##   sequence <<- sequence
                     ##   id <<- id
                     ##   alphabet <<- alphabet
                     ##   if (length(.self$alphabet) == 0)
                     ##     alphabet <<- unique(strsplit(.self$sequence, "")[[1]])
                     ##   if (!.self$valid())
                     ##     stop("Object is not valid!")
                     ## },
                     setAlphabet = function() {
                       'Sets the alphabet.'
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
                       cat("Reference class ",
                           classLabel(class(.self)), "\n")
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
                       .self$setAlphabet()
                     }
                     ))

