####################################
####################################
### Get snp information of the specified affymetrix snp array
getSnpInfo <- function(pkgname){
    ### open connection to database
    require(pkgname, character.only=TRUE, quietly=TRUE)
    conn <- db(get(pkgname))

    ### SNP informationp
    #snps <- paste("(", paste(unique(snpInfo$fsetid), collapse=","), ")", sep="")
    cat("Loading SNP annotation form the pd.mapping50k.xba240, pd.mapping50k.hind240, pd.mapping250k.sty or pd.mapping250k.nsp packages\n")
    
    sql <- paste("SELECT fsetid, dbsnp_rs_id, chrom, physical_pos, fragment_length",  
        "FROM featureSet ")
    #    "WHERE fsetid IN", snps)   
            
    system.time(snpInfo <- dbGetQuery(conn, sql))
    snpInfo <- snpInfo[!is.na(snpInfo$dbsnp_rs_id), ]
    snpInfo <- snpInfo[!is.na(snpInfo$fragment_length), ]
    colnames(snpInfo)[3:4] <- c("Chr", "X")
    return(snpInfo)
}



