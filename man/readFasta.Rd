\name{readFasta}
\alias{readFasta}
\title{Read fasta files.}
\usage{readFasta(infile)}
\description{Reads sequences data in fasta and create \code{DnaSeq}
and \code{RnaSeq} instances.}
\details{This funtion reads DNA and RNA fasta files and generates
valid \code{"DnaSeq"} and \code{"RnaSeq"} instances.}
\value{an instance of \code{DnaSeq} or \code{RnaSeq}.}
\seealso{\code{\linkS4class{GenericSeq}}, \code{\linkS4class{DnaSeq}} and \code{\linkS4class{RnaSeq}}.}
\author{Laurent Gatto \email{lg390@cam.ac.uk}}
\keyword{IO,}
\keyword{file}
\arguments{\item{infile}{the name of the fasta file which the data are to be read from.}}
\examples{
f <- dir(system.file("extdata",package="sequences"),pattern="fasta",full.names=TRUE)
f
aa <- readFasta(f)
aa}
