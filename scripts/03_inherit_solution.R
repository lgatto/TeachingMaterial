
source("03_inherit.R")

#### S3 methods
# method definition
comp.RnaSeq <- function(x) chartr("ACGU","UGCA", seq(x))

#### Test code

s.rna <- structure(list(id="RNA sequence", sequence="AUGCUAGUAC", alphabet=c("A","C","G","U")), class="GenericSeq")
s.rna <- setSeqSubtype(s.rna)

print(s.rna)
print(comp(s.rna))

# assigning more classes
class(s.rna) = c("lm", class(s.rna))
print(class(s.rna))
print(s.rna) # produces strange output
print(seq(s.rna)) # but seq() function still works..
