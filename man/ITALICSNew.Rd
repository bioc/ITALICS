\name{ITALICS}
\alias{ITALICS}
\alias{Italics}
\alias{italics}
\title{Affymetrix SNP chip normalization} 
\description{
 Normalize and analyse Affymetrix SNP array 100K and 500K set (see the vignette)
}
\usage{
ITALICS(quartetInfo, snpInfo, confidence=0.95, iteration=2, 
    formule="Smoothing+QuartetEffect+FL+I(FL^2)+I(FL^3)+GC+I(GC^2)+I(GC^3)", prc=0.3,
    amplicon=2.1, deletion=-3.5, deltaN=0.15, forceGL=c(-0.2,0.2), param=c(d=2), nbsigma=1, ... )
}
\arguments{
  \item{quartetInfo}{a data frame containing all the raw quartet intensities plus there GC content, fragment length, and Quartet effect}
  \item{snpInfo}{a data frame containing SNPs position along the genome and raw copy number}
  \item{confidence}{The confidence interval. After the last bias estimation step, quartets outside this confidence interval are flagged. The lower confidence is, the more quartets will be flagged. See also the parameter prc.}
  \item{iteration}{The number of iteration you d'like to do}
  \item{formule}{A symbolic description of the term of the model. The default value of formule means that we want correct the observed quartetLogRatio using the estimated copy number (Smoothing), the Quartet Effect, the quartet Fragment Length (FL) and the quartet GC content.}
  \item{prc}{prc is a frequence (between 0 and 1). After the final iteration of ITALICS, badly predicted probes are flagged (see also the parameter confidence). Only SNPs having more than prc of their probes non-flagged are kept for the final GLAD analysis. The higher prc is, the more SNPs are removed before the final GLAD analysis.}
  \item{amplicon}{see the amplicon parameter in the daglad function}
  \item{deletion}{see the deletion parameter in the daglad function}
  \item{deltaN}{see the deltaN parameter in the daglad function}
  \item{forceGL}{see the forceGL parameter in the daglad function}
  \item{param}{see the param parameter in the daglad function}
  \item{nbsigma}{see the nbsigma parameter in the daglad function}
  \item{...}{Other daglad parameters.}
}
\value{
  Return an object of class profileCGH
}
\examples{
\dontrun{
## step to get the path of the HF0844_Hind.CEL file
ITALICSDataPATH <- attr(as.environment(match("package:ITALICSData",search())),"path")
filename <- paste(ITALICSDataPATH,"/extdata/HF0844_Hind.CEL", sep="")
quartetEffectFile <- paste(ITALICSDataPATH,"/data/Hind.QuartetEffect.csv", sep="")

## load quartet effect
quartetEffect <- read.table(quartetEffectFile, sep=";", header=TRUE)

## load annotation using the pd.mapping50k.xba24 or pd.mapping50k.hind240 or  pd.mapping250k.sty or pd.mapping250k.nsp package
headdetails <- readCelHeader(filename[1])
pkgname <- cleanPlatformName(headdetails[["chiptype"]])
snpInfo <- getSnpInfo(pkgname)
quartet <- getQuartet(pkgname, snpInfo)

## read cel files and format data
tmpExprs <- readCelIntensities(filename, indices=quartet$fid)
quartet$quartetInfo$quartetLogRatio <- readQuartetCopyNb(tmpExprs)
quartet$quartetInfo <- addInfo(quartet, quartetEffect)
snpInfo <- fromQuartetToSnp(cIntensity="quartetLogRatio", quartetInfo=quartet$quartetInfo, snpInfo=snpInfo)


## ITALICS normalization
profilSNPHind <- ITALICS(quartet$quartetInfo, snpInfo,
    formule="Smoothing+QuartetEffect+FL+I(FL^2)+I(FL^3)+GC+I(GC^2)+I(GC^3)")

## plot the profile
data(cytoband)
plotProfile(profilSNPHind, Smoothing="Smoothing", Bkp=TRUE, cytoband = cytoband)
}
}
\details{

  The function \code{ITALICS} implements the methodology which
  is described in the article : ITALICS: an algorithm for normalization and DNA copy number calling for Affymetrix SNP arrays   (Rigaill et al., Bioinformatics Advance Access published on February 5, 2008).

  The principle of the ITALICS algorithm:
  ITALICS, is a normalization method that estimates both the biological and the non-relevant effects
  in an alternate and iterative way to accurately remove the non-relevant effects.

  ITALICS deals with known systematic sources of variation such as the GC-content of the quartets, the PCR amplified fragment length and the GC-content of the PCR amplified fragment . It also takes into account the quartet effect which corresponds to the fact that some quartets systematically have a small intensity while others tend to have a high intensity. 
  ITALICS is also able to correct spatial artifacts which sometimes arise on Affymetrix SNP arrays 100K and 500K set.

}
\note{People interested in tools dealing with array CGH analysis and DNA copy number analysis can
  visit our web-page \url{http://bioinfo.curie.fr}.}

\author{Guillem Rigaill, \email{italics@curie.fr}.}
\source{Institut Curie, \email{italics@curie.fr}.}
\keyword{}


 
