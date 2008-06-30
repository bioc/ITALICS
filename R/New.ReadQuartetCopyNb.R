################# READ
### read only pm intensities given their position on the cel file
readQuartetCopyNb <-  function(tmpExprs){
   b <- 2 * 1:(length(tmpExprs)/2) 
   return(log2(tmpExprs[b, ]+ tmpExprs[b-1, ]))
}

