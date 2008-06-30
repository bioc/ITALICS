ITALICS <- function(quartetInfo, snpInfo, confidence=0.95, iteration=2, 
    formule="Smoothing+QuartetEffect+FL+I(FL^2)+I(FL^3)+GC+I(GC^2)+I(GC^3)", prc=0.3,
    amplicon=2.1, deletion=-3.5, deltaN=0.15, forceGL=c(-0.2,0.2), param=c(d=2), nbsigma=1, ... ){
    if(iteration >=2){

    cat("####### FIRST ROUND #######\n")
    profilSNP <- analyseCGH(snpInfo, amplicon, deletion, deltaN, forceGL, param, nbsigma)
    profilSNP <- profilSNP$profileValues

    quartetInfo <- fromSnpToQuartet(quartetInfo, profilSNP)

    ## data scaling
    ec <- sd(profilSNP$LogRatio- profilSNP$Smoothing)
    mu <- mean(profilSNP$LogRatio)
    quartetInfo$quartetLogRatio  <- (quartetInfo$quartetLogRatio - mu)/ec
    quartetInfo$Smoothing <- (quartetInfo$Smoothing - mu)/ec
    rm(profilSNP)

    ## bias estimation 1
    cat("Bias Estimation\n")
    model <- getModel(formule, "quartetLogRatio", quartetInfo)

    quartetInfo$quartetCorrection <- getCorrection("Smoothing", model, quartetInfo)
    snpInfo   <- fromQuartetToSnp(cIntensity="quartetCorrection", quartetInfo=quartetInfo, snpInfo=snpInfo)

    ## glad 2
    i <- 2
    while(i <= iteration-1){
        cat(paste("####### ROUND ",i," #######\n",sep=""))
        profilSNP <- analyseCGH(snpInfo, amplicon, deletion, deltaN, forceGL, param, nbsigma)
        profilSNP <- profilSNP$profileValues

        quartetInfo <- fromSnpToQuartet(quartetInfo, profilSNP)
        rm(profilSNP)

        ## bias estimation 2
        cat("Bias Estimation\n")
        model <- getModel(formule, "quartetLogRatio", quartetInfo)

        quartetInfo$quartetCorrection <- getCorrection("Smoothing", model, quartetInfo)
        snpInfo   <- fromQuartetToSnp(cIntensity="quartetCorrection", quartetInfo=quartetInfo, snpInfo=snpInfo)
        i <- i+1
    }
    cat("####### FINAL ROUND  #######\n")
    profilSNP <- analyseCGH(snpInfo, amplicon, deletion, deltaN, forceGL, param, nbsigma)
    profilSNP <- profilSNP$profileValues

    quartetInfo <- fromSnpToQuartet(quartetInfo, profilSNP)
    rm(profilSNP)

    ## bias estimation 2
    cat("Bias Estimation\n")
    model <- getModel(formule, "quartetLogRatio", quartetInfo)

    quartetInfo$quartetCorrection <- getCorrection("Smoothing", model, quartetInfo)

    ## Removing badly predicted probes
    cat("Elimination of badly predicted probes\n")
    quartetInfo$quartetCorrectionFlagged <- getConfDat(confidence, quartetInfo, model)
    ## Count percentage of eliminated quartet per SNP
    quartetInfo$quartetNbKept <- 1
    quartetInfo$quartetNbKept[is.na(quartetInfo$quartetCorrectionFlagged)] <- 0

    snpInfo   <- fromQuartetToSnp(cIntensity=c("quartetCorrection", "quartetCorrectionFlagged", "quartetNbKept"), nLog=2, quartetInfo=quartetInfo, snpInfo=snpInfo)
    
    ## removing SNP with too many eliminated quartets
    snpInfo$LogRatio[snpInfo$quartetNbKept < prc] <- NA

    ## glad 3 analyse
    cat("####### ANALYSIS #######\n")
    profilSNP <- analyseCGH(snpInfo, amplicon, deletion, deltaN, forceGL, param, nbsigma)

    return(profilSNP)
    }else{
    cat("You must choose a number of iteration equal or greater than 2\n")
    return(NA)
    }
}

