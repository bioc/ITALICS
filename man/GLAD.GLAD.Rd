\name{analyseCGH}
\alias{analyseCGH}
\title{GLAD analysis} 
\description{
 Glad Analysis ot the genomic profile
}
\usage{
analyseCGH(data, amplicon, deletion, deltaN, forceGL, param, nbsigma, ...)
}
\arguments{
  \item{data}{A data frame containing SNP's intensity, chromosome and
    position on the genome. data must have a Chr, X and LogRatio
    columns} 
  \item{amplicon}{see the amplicon parameter in the daglad function}
  \item{deletion}{see the deletion parameter in the daglad function}
  \item{deltaN}{see the deltaN parameter in the daglad function}
  \item{forceGL}{see the forceGL parameter in the daglad function}
  \item{param}{see the param parameter in the daglad function}
  \item{nbsigma}{see the nbsigma parameter in the daglad function}
  \item{...}{Other daglad parameters.}
}
\value{
  An object of class profileCGH
}

\note{People interested in tools dealing with array CGH analysis and DNA
  copy number analysis can 
  visit our web-page \url{http://bioinfo.curie.fr}.}

\author{Guillem Rigaill, \email{italics@curie.fr}.}
\source{Institut Curie, \email{italics@curie.fr}.}
\keyword{}

