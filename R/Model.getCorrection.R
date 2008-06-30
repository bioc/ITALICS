getResidu <- function(model){
    return(as.vector(model$residuals))
}


getEffet <- function(effet, model, regTab){
    colonReg <- grep(effet,colnames(regTab))
    placeCoef <- grep(effet,names(model$coefficients))
    coef <- model$coefficient[placeCoef]
    effet <- coef * regTab[,colonReg]
    return(effet)
}

## fonction to eliminate the bias effect
getCorrection <- function(effet, model, regTab){
    res <- getResidu(model)
    effet <- getEffet(effet, model, regTab)
    return(effet+res)
}

## fonction which eliminate badly predicted probes and eliminate the bias effect
getConfDat <- function(confidence, quartetInfo, model){
    residu <- getResidu(model)
    logRatio <- getCorrection("Smoothing", model, quartetInfo)
    foo <- predict.lm(model, quartetInfo, level=confidence, interval="predict")
    marge <- (foo[, 1]-foo[, 2]) 
    rm(foo)

    e <- which(residu < -marge | residu > marge)
    rm(marge)
    rm(residu)
    logRatio[e] <- NA
    return(logRatio)
}

######################
