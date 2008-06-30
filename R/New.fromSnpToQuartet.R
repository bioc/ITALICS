#############################################
### get mean snp intensity for each quartet
fromSnpToQuartet <- function(quartetInfo, profilSNP){
    #which(colnames(profilSNP) == "fsetid" | colnames(profilSNP) == "Smoothing")
    e <- which(colnames(quartetInfo) == "Smoothing")
    if(length(e) > 0) quartetInfo <- quartetInfo[, -e]
    profilSNP <- profilSNP[,  which(colnames(profilSNP) == "fsetid" | colnames(profilSNP) == "Smoothing")]
    dat <- merge(profilSNP, quartetInfo)
    return(dat)
}

