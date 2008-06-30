\name{fromQuartetToSnp}
\alias{fromQuartetToSnp}
\alias{fromQuartetToSnp}
\alias{fromQuartetToSnp}
\title{Compute the copy number of each SNP from its quartets intensities} 
\description{
This function removes the LogRatio column of the snpInfo data.frame.
Then compute the copy number of each SNP having its quartet intensities.
And return the snpInfo data.frame with the newly computed LogRatio.
}
\usage{
  fromQuartetToSnp(quartetInfo, snpInfo, cIntensity="quartetLogRatio", nLog=1)
}
\arguments{
 
  \item{quartetInfo}{A table containing the quartet intensities and other quartet information. 
  It must have a column called : fsetid.}
  \item{snpInfo}{A table containing snp information.}
  \item{cIntensity}{A vector containing the names of the quartet information to be aggregate. For example quartetLogRatio. }
  \item{nLog}{The position of the field which will be named LogRatio in the snpInfo data.frame. For example if cIntensity = c("a", "b") and you want b to be considered as the LogRatio you should set nLog=2}
}
\value{
  \item{return the data.frame snpInfo with additionnal columns. } 
}
\examples{}
\note{People interested in tools dealing with array CGH analysis and DNA copy number analysis can
  visit our web-page \url{http://bioinfo.curie.fr}.}

\author{Guillem Rigaill, \email{italics@curie.fr}.}
\source{Institut Curie, \email{italics@curie.fr}.}
\keyword{}
