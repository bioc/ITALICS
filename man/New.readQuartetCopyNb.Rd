\name{readQuartetCopyNb}
\alias{readQuartetCopyNb}
\alias{readQuartetCopyNb}
\alias{readQuartetCopyNb}
\title{Read PM probes of selected quartets and compute the quartet intensity} 
\description{
This function read the cel files and return the raw-value of each quartet = mean of allele A and B}
\usage{
  readQuartetCopyNb(tmpExprs)
}
\arguments{
  \item{tmpExprs}{A vector of the perfect match intensity of allele A and B of the quartets. This vector should be sorted in a specific order. See the example given in the help of the ITALICS function.}
  
}
\value{
  return a vector with the raw-value of each quartet
}
\note{People interested in tools dealing with array CGH analysis and DNA copy number analysis can
  visit our web-page \url{http://bioinfo.curie.fr}.}

\author{Guillem Rigaill, \email{italics@curie.fr}.}
\source{Institut Curie, \email{italics@curie.fr}.}
\keyword{}
