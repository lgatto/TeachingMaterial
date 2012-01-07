\name{Seq-class}
\Rdversion{1.1}
\docType{class}
\alias{Seq-class}
\alias{class:Seq}
\alias{Seq}

\title{Class \code{"Seq"}}

\description{
  The \code{Seq} class implements a nucleic acid sequence as well as
  various sequence transformations. It's aim is to illustrate the usage
  of reference classes.
  
}
\section{Extends}{
  All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.
}


\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

\seealso{
  \code{\linkS4class{DnaSeq}} for a S4 implementation of the same subject.
}

\examples{
s <- Seq$new(id="My sequence",
             sequence="AGCAGCTACGACT",
             alphabet=c("A","C","G","T"))
s
s$rev()
s
}

\keyword{classes}
\section{Fields}{
  \describe{
    \item{\code{id}:}{Object of class \code{character} defining the
      sequence identifier. }
    \item{\code{alphabet}:}{Object of class \code{character} defining
      the sequence alphabet. }
    \item{\code{sequence}:}{Object of class \code{character} with the
      actual sequence string. }
  }
}
\section{Methods}{
  \describe{
    \item{\code{transcribe()}:}{ Transcribes a DNA sequence into RNA. }
    \item{\code{seq()}:}{ Retuns the sequence string. }
    \item{\code{comp()}:}{ Complements the sequence. }
    \item{\code{rev()}:}{ Reverses the sequence. }
    \item{\code{valid()}:}{ Checks sequence validity. }
    \item{\code{setAlphabet()}:}{ Sets the sequence alphabet. }
    \item{\code{show()}:}{ Shows the 'Seq' instance. }
  }
}
