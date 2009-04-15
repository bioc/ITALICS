\name{trainITALICS}
\alias{trainITALICS}
\alias{trainItalics}
\alias{trainitalics}
\title{ITALICS training} 
\description{
 Estimation of the quartet effect based on several normal sample chips
}
\usage{
trainITALICS (dir,  amplicon=2.1, deletion=-3.5, deltaN=0.15, forceGL=c(-0.2,0.2), param=c(d=2), nbsigma=1, ...)
}
\arguments{
  \item{dir}{The directory containing the normal sample chips. All theses chips should be of the same type hind, xba, nsp or sty. Only .CEL files be considered}
  \item{amplicon}{see the amplicon parameter in the daglad function}
  \item{deletion}{see the deletion parameter in the daglad function}
  \item{deltaN}{see the deltaN parameter in the daglad function}
  \item{forceGL}{see the forceGL parameter in the daglad function}
  \item{param}{see the param parameter in the daglad function}
  \item{nbsigma}{see the nbsigma parameter in the daglad function}
  \item{...}{Other daglad parameters.}
}
\value{
  a data.frame with two column fsetid and quartetEffect
}

\details{
The ITALICS function take into account a quartet effect which is computed
on a reference data set of normal women samples. The ITALICSData provide 
quartetEffect for the Xba, Hind, Sty and Nsp chip computed on our own reference data
set. 

We recommand that you use your own reference data set to compute the quartet Effect 
by using the trainITALICS function. ITALICS reference data should contain only woman normal samples.
Furthermore we recommand that you check that none of these chip have obvious spatial artifact. To so
read the cel files using the read.affybatch (form the affy package). Then use the image function on the obtain
affybatch object.
}
\note{People interested in tools dealing with array CGH analysis and DNA copy number analysis can
  visit our web-page \url{http://bioinfo.curie.fr}.}

\author{Guillem Rigaill, \email{italics@curie.fr}.}
\source{Institut Curie, \email{italics@curie.fr}.}
\keyword{}


