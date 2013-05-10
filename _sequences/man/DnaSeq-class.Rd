\name{DnaSeq-class}
\Rdversion{1.1}
\docType{class}
\alias{DnaSeq-class}
\alias{class:DnaSeq}
\alias{RnaSeq-class}
\alias{class:RnaSeq}

\alias{comp}
\alias{comp-method}
\alias{comp,DnaSeq-method}
\alias{comp,RnaSeq-method}
\alias{transcribe}
\alias{transcribe-method}
\alias{transcribe,DnaSeq-method}


\title{Class "DnaSeq" and "RnaSeq"}

\description{
  The \code{DnaSeq} and \code{RnaSeq} are instances of the virtual
  \code{\linkS4class{GenericSeq}} class for "DNA" and "RNA"
  respectively. Their alphabets are defined accordingly. 
  
}

\section{Objects from the Class}{
%% Objects can be created by calls of the form \code{new("DnaSeq",
%% ...)}.
Instance of \code{DnaSeq} and \code{RnaSeq} can be created by reading
\code{fasta} files using the \code{\link{readFasta}} function. 
}

\section{Slots}{
  See the \code{\linkS4class{GenericSeq}} class for details about the
  generic slots and methods.
}
\section{Extends}{
Class \code{\linkS4class{GenericSeq}}, directly.
}

\section{Methods}{
  See \code{\linkS4class{GenericSeq}} for inherited methods.
  \describe{
    \item{comp}{\code{signature(object = "DnaSeq")}: return the complement of
      the \code{sequence} string of the \code{sequence} instance. }
    \item{transcribe}{\code{signature(object = "DnaSeq")}: transcribes
      the \code{DnaSeq} instance into a \code{RnaSeq} instance. }
  }
}

\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

% \note{
% %%  ~~further notes~~
% }

% %% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  See also the \code{\linkS4class{GenericSeq}} virutal class.
}

\examples{
data("dnaseq")
dnaseq
rnaseq <- transcribe(dnaseq)
rnaseq
}

\keyword{classes}
