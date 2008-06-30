\name{addInfo}
\alias{addInfo}
\alias{addInfo}
\alias{addInfo}
\title{add info to quartet annotation} 
\description{
This function merge information obtain from the getQuartet function and a given table
}
\usage{
  addInfo(quartet, dat)
}
\arguments{
  \item{quartet}{list obtain through the getQuartet Function}
  \item{dat}{a data.frame with additionnal information it must contain a fsetid and fid column}
}
\value{
  a data.frame similar to the quartetInfo item of quartet plus additionnal column
}
\note{People interested in tools dealing with array CGH analysis and DNA copy number analysis can
  visit our web-page \url{http://bioinfo.curie.fr}.}

\author{Guillem Rigaill, \email{italics@curie.fr}.}
\source{Institut Curie, \email{italics@curie.fr}.}
\keyword{}
