############################ INFO
############################
####### Get position on the cel, and GC probe content, fragment length and quartetEffect
getQuartet <- function(pkgname, snpInfo){
    #### find the  cel type hind, xba... and load the corresponding annotation package
 
    require(pkgname, character.only=TRUE, quietly=TRUE)

    ### open connection to database
    conn <- db(get(pkgname))

    ### retrive information for each probes: allele, strand, fid, fsetid, sequence, fragment length...
    #snps <- paste("(", paste(unique(snpInfo$fsetid), collapse=","), ")", sep="")
    cat("Loading probes annotation from the pd.mapping50k.xba240, pd.mapping50k.hind240, pd.mapping250k.sty or pd.mapping250k.nsp packages\n")

    #sql <- paste("SELECT pmfeature.fsetid, pmfeature.fid, pmfeature.allele, pmfeature.strand, offset, fragment_length, seq", 
    #      " FROM sequence, pmfeature, featureSet ",
    #      "WHERE pmfeature.fid=sequence.fid AND pmfeature.fsetid=featureSet.fsetid AND pmfeature.fsetid IN", snps) 
    #AND featureSet.affy_snp_id""
    sql <- paste("SELECT pmfeature.fsetid, pmfeature.fid, pmfeature.allele, pmfeature.strand, offset, seq", 
          " FROM sequence, pmfeature", 
          "WHERE pmfeature.fid=sequence.fid ")
#AND pmfeature.fsetid IN", snps)
    system.time(quartetInfo <- dbGetQuery(conn, sql))

    quartetInfo <- merge(quartetInfo, snpInfo)
    ### order the probe information according to the SNP, strand, offset and allele
    ### this way all probes of a SNP are regrouped
    ### probes of a quartet are regrouped a line for allele A and the following for allele B
    quartetInfo <- quartetInfo[order(quartetInfo$fsetid, quartetInfo$strand, quartetInfo$offset, quartetInfo$allele), ]
 
    ### compute GC
    cat("Computing probes GC content. Please wait. This can be quite long especially for Sty and Nsp chips.\n")
    quartetInfo$GC <- apply(basecontent(quartetInfo$seq)[, 3:4], 1, sum)
 

    ###
    fid <- quartetInfo$fid
    b <- 2 * 1:(nrow(quartetInfo)/2) 
    GC <- quartetInfo$GC[b-1]+ quartetInfo$GC[b]
    quartetInfo <- quartetInfo[b, -c(3:9, 11)]
    quartetInfo$GC <- GC
    
    colnames(quartetInfo)[3] <- "FL"
    ### Result
    quartet <- list() 
    quartet$fid <- fid
    quartet$quartetInfo <- quartetInfo
    return(quartet)
}
