############################################
### get copy nb per snp from quartetInfo
fromQuartetToSnp <- function(quartetInfo, snpInfo, cIntensity="quartetLogRatio", nLog=1){
    e <- which(colnames(snpInfo) == "LogRatio")
    if(length(e) > 0) snpInfo <- snpInfo[, -e]
    e <-  unlist(lapply(paste("^", cIntensity, "$", sep=""), grep, colnames(quartetInfo)))
    dat <- aggregate(data.frame(quartetInfo[, e]), by= list(quartetInfo$fsetid), FUN=mean, na.rm=TRUE)
    colnames(dat)[c(1, nLog+1)] <- c("fsetid", "LogRatio")
    dat <- merge(snpInfo, dat)
    dat <- dat[order(dat$Chr, dat$X), ]
    return(dat)
}

