\name{getSnpInfo}
\alias{getSnpInfo}
\alias{getSnpInfo}
\alias{getSnpInfo}
\title{Function to retrieve the chromosome and the position of each SNP on a given Affymetrix SNP array} 
\description{
This function retrieve the chromosome and position in bp of each SNP of a given Affymetrix SNP array.
This function use the pd.mapping50k.xba240, pd.mapping50k.hind240, pd.mapping250k.sty and pd.mapping250k.nsp package.
}
\usage{
  getSnpInfo(pkgname)
}
\arguments{
  \item{pkgname}{the chip type  pd.mapping50k.xba240, pd.mapping50k.hind240, pd.mapping250k.sty or pd.mapping250k.nsp}
}
\value{
  \item{Return a data.frame with five columns : fsetid, dbsnp_rs_id, Chr, X and fragment_length corresponding to the fsetid, the rs_id, the chromosome,
the position on the chromosome and the PCR amplified fragment length respectively.} 
}
\note{People interested in tools dealing with array CGH analysis and DNA copy number analysis can
  visit our web-page \url{http://bioinfo.curie.fr}.}

\author{Guillem Rigaill, \email{italics@curie.fr}.}
\source{Institut Curie, \email{italics@curie.fr}.}
\keyword{}


 
