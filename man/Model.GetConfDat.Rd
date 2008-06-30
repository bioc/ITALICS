\name{getConfDat}
\alias{getConfDat}
\title{Elimination of badly predicted probes} 
\description{
 This function eliminate badly predicted probes using a regression table and an estimated model given by the function \code{getModel} or \code{getBestBICModelLight}. Then it computes the corrected intensity.
}
\usage{
getConfDat(confidence, quartetInfo, model)
}
\arguments{
  \item{confidence}{The confidence interval : 0.95}
  \item{quartetInfo}{A Regression table containing the variables in the model}
  \item{model}{The class lm object given by the function \code{getModel} }
}
\value{
  \item{A data frame with the corrected intensity. Only goodly predicted probes are taken into account. SNP's with more than 8 badly predicted probes get a NA.} 
}
\note{People interested in tools dealing with array CGH analysis and DNA copy number analysis can
  visit our web-page \url{http://bioinfo.curie.fr}.}

\author{Guillem Rigaill, \email{italics@curie.fr}.}
\source{Institut Curie, \email{italics@curie.fr}.}
\keyword{}


