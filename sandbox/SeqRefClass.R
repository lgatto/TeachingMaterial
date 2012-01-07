library("sequences")

Seq$fields()
Seq$methods()
Seq$help("rev")
s <- Seq$new(id="foo", sequence="GATCATCA")
s
s$valid()
s$setAlphabet()
s
s$valid()
##s <- Seq$new(id="foo", sequence="GATCATCA")
##s
s$rev()
s
s$seq()
s$sequence
s$field("sequence")
## !!!
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

s$comp()
s
s$transcribe()
s
