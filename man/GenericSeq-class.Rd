\name{GenericSeq-class}
\Rdversion{1.1}
\docType{class}
\alias{GenericSeq-class}
\alias{class:GenericSeq}

\alias{alphabet}
\alias{alphabet-method}
\alias{alphabet,GenericSeq-method}
\alias{alphabet,DnaSeq-method}
\alias{alphabet,RnaSeq-method}
\alias{id<-}
\alias{id<-,GenericSeq,ANY-method}
\alias{id<-,GenericSeq,character-method}
\alias{id}
\alias{id-method}
\alias{id,GenericSeq-method}
\alias{id,DnaSeq-method}
\alias{id,RnaSeq-method}
\alias{length,GenericSeq-method}
\alias{length,DnaSeq-method}
\alias{length,RnaSeq-method}
\alias{print,GenericSeq-method}
\alias{print,DnaSeq-method}
\alias{print,RnaSeq-method}
\alias{rev}
\alias{rev-method}
\alias{rev,GenericSeq-method}
\alias{rev,RnaSeq-method}
\alias{rev,DnaSeq-method}
\alias{seq}
\alias{seq-method}
\alias{seq,GenericSeq-method}
\alias{seq,DnaSeq-method}
\alias{seq,RnaSeq-method}
\alias{seq<-}
\alias{seq<-,GenericSeq,ANY-method}
\alias{seq<-,GenericSeq,character-method}
\alias{show,GenericSeq-method}
\alias{show,DnaSeq-method}
\alias{show,RnaSeq-method}
\alias{[,GenericSeq-method}
\alias{[,DnaSeq-method}
\alias{[,RnaSeq-method}



\title{Class "GenericSeq"}

\description{
  This virtual class defines a generic sequence object composed
  by an identifier, a sequence and an alphabet that defines the
  letters of the sequence. A set of accessor methods and setters
  are described below.
}

\section{Objects from the Class}{A virtual Class: No objects may be
  created from it.} 
\section{Slots}{
  \describe{
    \item{\code{id}:}{Object of class \code{"character"} that identifies
      the \code{sequence} instance.}
    \item{\code{alphabet}:}{Object of class \code{"character"} that
      defines all the letters that a valid \code{sequence} can be made
      of. }
    \item{\code{sequence}:}{Object of class \code{"character"} storing
      the actual sequence string of the \code{sequence} instance. }
  }
}

\section{Methods}{
  \describe{
    \item{alphabet}{\code{signature(object = "GenericSeq")}: retrieves
      the alphabet of valid sequence letters, returned as
      \code{characher}. }
    \item{id}{\code{signature(object = "GenericSeq")}: retrieves the
      identifier of the \code{sequence} instance. }
    \item{id<-}{\code{signature(object = "GenericSeq", value =
	"character")}: sets the identifier of the \code{sequence}
      instance to \code{value}. }
    \item{length}{\code{signature(x = "GenericSeq")}: returns the length
      the \code{sequence} instance, i.e. the number of characters of the
      \code{sequence} string.}
    \item{print}{\code{signature(x = "GenericSeq")}: prints the
      \code{sequence} of the instance.}
    \item{rev}{\code{signature(x = "GenericSeq")}: return the reverse of
      the \code{sequence} string of the \code{sequence} instance. }
    \item{seq}{\code{signature(... = "GenericSeq")}: returns the
      \code{sequence} character string of the instance. }
    \item{seq<-}{\code{signature(object = "GenericSeq", value =
	"character")}: changes the \code{sequence} string of the
      \code{sequence} instance to \code{value}. The new string must be
      compatbile with the \code{sequence} instance alphabet.}
    \item{show}{\code{signature(object = "GenericSeq")}: displays the
      object content. }
    \item{[}{\code{signature(x = "GenSeq")}: Subsets the sequence of
      \code{x} and returns the correspong subsetted object of same class
      as \code{x}. }
  }
}


\author{
  Laurent Gatto <lg390@cam.ac.uk>
}



\seealso{
  See \code{\linkS4class{DnaSeq}} and \code{\linkS4class{RnaSeq}}
  classes for examples of subclasses and \code{\link{readFasta}} 
  on how to create instances.
}

\examples{
showClass("GenericSeq") 
}

\keyword{classes}
