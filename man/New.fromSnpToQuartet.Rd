\name{fromSnpToQuartet}
\alias{fromSnpToQuartet}
\alias{fromSnpToQuartet}
\alias{fromSnpToQuartet}
\title{Function to get from snp to quartet} 
\description{
This function put the smoothing value of each SNP in front of its corresponding quartet
in the quartetInfo data.frame.

}
\usage{
 fromSnpToQuartet(quartetInfo, profilSNP)
}
\arguments{
  \item{quartetInfo}{a data frame containing all the quartet values plus there GC content, fragment length and GC content and Quartet effect}
  \item{profilSNP}{a data frame, corresponding to the profileValues argument of a profilCGH object (see GLAD)}
}
\value{return the data.frame quartetInfo with an additionnal column: "Smoothing" corresponding to the estimated smoothing value.}

\note{People interested in tools dealing with array CGH analysis and DNA copy number analysis can
  visit our web-page \url{http://bioinfo.curie.fr}.}

\author{Guillem Rigaill, \email{italics@curie.fr}.}
\source{Institut Curie, \email{italics@curie.fr}.}
\keyword{}
