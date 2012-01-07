\name{gccount}
\alias{gccount}

\title{Count bases in sequence}

\usage{gccount(inseq)}

\description{Returns the number of 'A', 'C', 'G' and 'T' bases in the
  'inseq' sequence string. }

\details{This function calls a C primitive}

\value{A numeric of length 4}

\author{Laurent Gatto \email{lg390@cam.ac.uk}}

\arguments{\item{inseq}{ a DNA sequence string. }}

\examples{
s <- "AAAACCCGGT"
cnt <- gccount(s)
cnt
stopifnot(cnt==table(strsplit(s,"")))
}
