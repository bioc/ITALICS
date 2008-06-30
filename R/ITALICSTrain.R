

trainITALICS <- function(dir,  amplicon=2.1, deletion=-3.5, deltaN=0.15, forceGL=c(-0.2,0.2), param=c(d=2), nbsigma=1, ...){
    filenames <- list.files(dir, full.name=TRUE)
    filenames <- filenames[grep("CEL", filenames)]
    i <- TRUE
    for(filename in filenames){
        print(filename)
        #### first chip ( get info )
        headdetails <- readCelHeader(filename)
        pkgname <- cleanPlatformName(headdetails[["chiptype"]])
        if( i == TRUE){
            print("First Time")
            print(pkgname)
            snpInfo <- getSnpInfo(pkgname)
            quartet <- getQuartet(pkgname, snpInfo)
            
       }
       #### LOAD DATA
       tmpExprs <- readCelIntensities(filename, indices=quartet$fid)
       quartet$quartetInfo$quartetLogRatio <- readQuartetCopyNb(tmpExprs)
       snpInfo <- fromQuartetToSnp(cIntensity="quartetLogRatio", quartetInfo=quartet$quartetInfo, snpInfo=snpInfo)
       profilSNP <- analyseCGH(snpInfo, amplicon, deletion, deltaN, forceGL, param, nbsigma)


       profilSNP <- profilSNP$profileValues

       ## data scaling
       ec <- sd(profilSNP$LogRatio- profilSNP$Smoothing)
       mu <- mean(profilSNP$LogRatio)
       if(i == TRUE){
           quartetLogRatio  <- (quartet$quartetInfo$quartetLogRatio - mu)/ec
           i <- FALSE
       } else {
           quartetLogRatio  <- quartetLogRatio + (quartet$quartetInfo$quartetLogRatio - mu)/ec
       }
       
   } 
   
   res <- data.frame(quartet$quartetInfo$fsetid, quartet$quartetInfo$fid, quartetLogRatio/length(filenames))
   colnames(res) <- c("fsetid", "fid", "QuartetEffect")
   return(res)
     

}
