addInfo <- function(quartet, dat){
    quartet$quartetInfo$ORDER <- 1: nrow(quartet$quartetInfo)
    dat <- merge(dat, quartet$quartetInfo)
    dat <- dat[order(dat$ORDER), ]
    dat <- dat[, -which(colnames(dat) ==  "ORDER")]
    return(dat)
}
