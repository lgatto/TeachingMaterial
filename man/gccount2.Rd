\name{gccount2}
\alias{gccount2}

\title{Count bases in sequence (using Rcpp)}

\usage{gccount2(inseq)}

\description{Returns the number of 'A', 'C', 'G' and 'T' bases in the
  'inseq' sequence string. }

\details{This function calls a C primitive}

\value{A numeric of length 4}

\author{Robert Stojnic \email{rs550@cam.ac.uk}}

\arguments{\item{inseq}{ a DNA sequence string. }}

\examples{
s <- "AAAACCCGGT"
cnt <- gccount2(s)
cnt
stopifnot(cnt==table(strsplit(s,"")))
}
