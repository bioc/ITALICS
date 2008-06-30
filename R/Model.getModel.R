## fonction to computes the lm object and get the BIC criteria of a model
getModel <- function(formule, response, regTab){
    model <- lm(formula = as.formula(paste(response, "~", formule,  sep="")), data=regTab, na.action=na.omit)
    return(model)
}
####################