setMethod("comp","RnaSeq",
          function(object) {
            chartr("ACGU","UGCA",seq(object))
          })
