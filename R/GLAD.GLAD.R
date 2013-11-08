analyseCGH <- function(data, amplicon, deletion, deltaN, forceGL, param, nbsigma, ...){
	ITALICSPATH <- attr(as.environment(match("package:ITALICS",search())),"path")
	source(paste(ITALICSPATH,"/scripts/tmp_filterBkp2.R", sep=""))
	source(paste(ITALICSPATH,"/scripts/tmp_filterBkpStep2.R", sep=""))
    ## format data in order to be used by as.profileCGH
    InputFields <- names(data)	
    data$Chromosome <- data$Chr
    data$PosBase <- data$X
    data <- data[order(data$Chromosome, data$PosBase),]
    data$PosOrder <- 1:length(data[, 1])

    ## use daglad on data
    profileCGH <- as.profileCGH(data)
    t1 <- system.time(profileCGH <- daglad(profileCGH, param=param,
                                       amplicon=amplicon, deletion=deletion, deltaN=deltaN,  forceGL=forceGL, 
                                       nbsigma=nbsigma, ...))
    
    return(profileCGH)
}

