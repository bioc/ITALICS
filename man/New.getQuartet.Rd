\name{getQuartet}
\alias{getQuartet}
\alias{getQuartet}
\alias{getQuartet}
\title{Function to retrieve the information of each quartet } 
\description{
This function retrieve information of each quartet. 
This function use the pd.mapping50k.xba240, pd.mapping50k.hind240, pd.mapping250k.sty and pd.mapping250k.nsp package.
}
\usage{
  getQuartet(pkgname, snpInfo)
}
\arguments{
  \item{pkgname}{the chip type pd.mapping50k.xba240, pd.mapping50k.hind240, pd.mapping250k.sty or pd.mapping250k.nsp}
  \item{snpInfo}{a data frame containing SNPs position along the genome}
}
\value{
  \item{return a list with two fields. 
   fid : containing the position of each quartet on the CEL file. 
   quartetInfo : a data fame containing the columns : fsetid, fid, FL (fragment length) and GC (content of the quartet)} 
}
\note{People interested in tools dealing with array CGH analysis and DNA copy number analysis can
  visit our web-page \url{http://bioinfo.curie.fr}.}

\author{Guillem Rigaill, \email{italics@curie.fr}.}
\source{Institut Curie, \email{italics@curie.fr}.}
\keyword{}





 
