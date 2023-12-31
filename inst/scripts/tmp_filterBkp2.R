filterBkp <- function(...)
  {
    UseMethod("filterBkp")
  }

filterBkp.profileCGH <- function(profileCGH, MinBkpWeight=0.25, assignGNLOut=TRUE, ...)
  {


    if (is.data.frame(profileCGH$BkpInfo))
      {
        profileCGH$profileValues <- profileCGH$profileValues[order(profileCGH$profileValues$PosOrder),]
        RecomputeGNL <- FALSE


##################################################################################
###
### On supprime les Breakpoints qui sont situ�s au sein des r�gions amplifi�es
###
##################################################################################
        
        
        indexBkpToDel <- which(profileCGH$BkpInfo$GNLchange==0 & profileCGH$BkpInfo$ZoneGNL==2)
        if (length(indexBkpToDel)>0)
          {
            RecomputeGNL <- TRUE
            profileCGH$profileValues$Breakpoints[profileCGH$BkpInfo$PosOrder[indexBkpToDel]] <- -1
            profileCGH$BkpInfo <- profileCGH$BkpInfo[-indexBkpToDel,]
          }


        
##################################################################################
###
### On D�place les Bkp qui sont aussi Outliers et dont
### le GNL correspond � celui du BAC d'apr�s
### On d�place �galement les Bkp apr�s lequel il y a un outlier
### correspondant au statut du Bkp        
###
##################################################################################
        

        profileCGH$profileValues <- profileCGH$profileValues[order(profileCGH$profileValues$PosOrder),]

        nb <- length(profileCGH$profileValues[,1])-1
                                       


        moveBkp <- .C("moveBkp",
                      as.integer(profileCGH$profileValues$ZoneGNL),
                      Level=as.integer(profileCGH$profileValues$Level),
                      Breakpoints=as.integer(profileCGH$profileValues$Breakpoints),
                      OutliersTot=as.integer(profileCGH$profileValues$OutliersTot),
                      OutliersAws=as.integer(profileCGH$profileValues$OutliersAws),
                      as.integer(profileCGH$profileValues$Chromosome),
                      RecomputeSmt=as.integer(0),
                      as.integer(nb),
                      PACKAGE="GLAD")


        if (moveBkp$RecomputeSmt==1)
          {            
            rownames(profileCGH$profileValues) <- 0:(length(profileCGH$profileValues[,1])-1)
            RecomputeGNL <- TRUE
            profileCGH$profileValues$Level <- moveBkp$Level
            profileCGH$profileValues$Breakpoints <- moveBkp$Breakpoints
            profileCGH$profileValues$OutliersTot <- moveBkp$OutliersTot
            profileCGH$profileValues$OutliersAws <- moveBkp$OutliersAws
          }
        

##################################################################################
###        
### Suppression des Bkp dont est poids est inf�rieur � un seuil
### et qui ne correspondent pas � un changement de GNL
###
##################################################################################        


        indexWeightToSmall <- which(profileCGH$BkpInfo$Weight<MinBkpWeight & profileCGH$BkpInfo$GNLchange==0 & profileCGH$BkpInfo$ZoneGNL!=2)
        if (length(indexWeightToSmall)>0)
          {
            RecomputeGNL <- TRUE
            for (PosBkp in profileCGH$BkpInfo$PosOrder[indexWeightToSmall])
              {
                #indexPos <- which(profileCGH$profileValues$PosOrder==PosBkp)
                #profileCGH$profileValues$Breakpoints[indexPos] <- -1
                profileCGH$profileValues$Breakpoints[PosBkp] <- -1
              }
            if (length(indexWeightToSmall)==length(profileCGH$BkpInfo[,1]))
              {
                profileCGH$BkpInfo <- NA
              }
            else
              {
                profileCGH$BkpInfo <- profileCGH$BkpInfo[-indexWeightToSmall,]
              }

          }

##################################################################################
###        
### Suppression des Bkp dont le poids vaut 0
### et qui correspondent � un changement de GNL
### cette situation peut arriver apr�s �limination des Bkp
### dont le poids est inf�rieur � un seuil
###
##################################################################################


        indexWeightZero <- which(profileCGH$BkpInfo$Weight==0 & profileCGH$BkpInfo$GNLchange==1)
        if (length(indexWeightZero)>0)
          {
            RecomputeGNL <- TRUE
            for (PosBkp in profileCGH$BkpInfo$PosOrder[indexWeightZero])
              {
                #indexPos <- which(profileCGH$profileValues$PosOrder==PosBkp)
                #profileCGH$profileValues$Breakpoints[indexPos] <- -1
                profileCGH$profileValues$Breakpoints[PosBkp] <- -1
              }
            if (length(indexWeightZero)==length(profileCGH$BkpInfo[,1]))
              {
                profileCGH$BkpInfo <- NA
              }
            else
              {
                profileCGH$BkpInfo <- profileCGH$BkpInfo[-indexWeightZero,]
              }

          }        


### Quand je vais recalculer les Outliers, il faut le NormalRef
### Attention � ce que celui-ci soit bien transmis
### Normalement NormalRef vaut 0 puisqu'en sortie de gladLA
### les log-ratios sont centr�s sur NormalRef
        

        if (RecomputeGNL)
          {


            profileCGH$profileValues <- profileCGH$profileValues[order(profileCGH$profileValues$PosOrder),]
            
            l <- length(profileCGH$profileValues[,1])
            updateLevel <- .C("updateLevel",
                              as.integer(profileCGH$profileValues$Chromosome),
                              Breakpoints=as.integer(profileCGH$profileValues$Breakpoints),
                              Level=as.integer(profileCGH$profileValues$Level),
                              as.integer(profileCGH$profileValues$PosOrder),
                              NextLogRatio=as.double(profileCGH$profileValues$NextLogRatio),
                              as.double(profileCGH$profileValues$LogRatio),
                              as.integer(max(profileCGH$profileValues$Level)),
                              as.integer(l),
                              PACKAGE="GLAD")

            profileCGH$profileValues$Level <- updateLevel$Level
            profileCGH$profileValues$Breakpoints <- updateLevel$Breakpoints
            profileCGH$profileValues$NextLogRatio <- updateLevel$NextLogRatio

            updateOutliers <- .C("updateOutliers",
                                 OutliersAws=as.integer(profileCGH$profileValues$OutliersAws),
                                 Level=as.integer(profileCGH$profileValues$Level),
                                 Breakpoints=as.integer(profileCGH$profileValues$Breakpoints),
                                 Smoothing=as.double(profileCGH$profileValues$Smoothing),
                                 as.integer(l),
                                 PACKAGE="GLAD")
            profileCGH$profileValues$Level <- updateOutliers$Level
            profileCGH$profileValues$Breakpoints <- updateOutliers$Breakpoints
            profileCGH$profileValues$OutliersAws <- updateOutliers$OutliersAws
            profileCGH$profileValues$Smoothing <- updateOutliers$Smoothing
            
            

### Recalcul des Outliers
            class(profileCGH) <- "profileChr"
            profileCGH <- detectOutliers(profileCGH, region="Level", alpha=profileCGH$alpha, msize=profileCGH$msize)

            
### recalcul de la smoothing line
            agg <- aggregate(profileCGH$profileValues$LogRatio, list(Level=profileCGH$profileValues$Level), median)
            agg$Level <- as.numeric(as.character(agg$Level))
            names(agg) <- c("Level","Smoothing")
            profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"Smoothing"))
            profileCGH$profileValues <- merge(profileCGH$profileValues, agg, by="Level", all=TRUE)



####################
            profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"ZoneGNL"))
            
            indexNormalLevel <- which(abs(profileCGH$profileValues$Smoothing-profileCGH$NormalRef)<=profileCGH$deltaN)
            profileCGH$profileValues$NormalRange <- profileCGH$profileValues$Level
            profileCGH$profileValues$NormalRange[indexNormalLevel] <- 0

            
### le clustering est fait sur les niveaux NormalRange
            profileCGH <- findCluster(profileCGH, region="NormalRange", method=profileCGH$method, genome=TRUE,
                                      lambda=profileCGH$lambdaclusterGen,
                                      nmin=profileCGH$NbClusterOpt, nmax=profileCGH$NbClusterOpt)
            

### le cluster correspondant au normal est celui qui comprend
### le NormalRange 0
            indexNormalRange <- which(profileCGH$profileValues$NormalRange==0)
            NormalCluster <- unique(profileCGH$profileValues$ZoneGen[indexNormalRange])
            MedianCluster <- aggregate(profileCGH$profileValues$LogRatio, list(ZoneGen=profileCGH$profileValues$ZoneGen),median,na.rm=TRUE)
            MedianCluster$ZoneGen <- as.numeric(as.character(MedianCluster$ZoneGen))
            names(MedianCluster) <- c("ZoneGen","Median")
            RefNorm <- MedianCluster$Median[which(MedianCluster$ZoneGen==NormalCluster)]
            MedianCluster$ZoneGNL <- rep(0,length(MedianCluster[,1]))
            indexClusterGain <- which(MedianCluster$Median>RefNorm)
            MedianCluster$ZoneGNL[indexClusterGain] <- 1
            indexClusterLost <- which(MedianCluster$Median<RefNorm)
            MedianCluster$ZoneGNL[indexClusterLost] <- -1                                        
            profileCGH$profileValues <- merge(profileCGH$profileValues, MedianCluster[,c("ZoneGen","ZoneGNL")], all=TRUE, by="ZoneGen")


### on force les gains et les pertes pour certaines valeur de smoothing
            indexForceGain <- which((profileCGH$profileValues$Smoothing-profileCGH$NormalRef) >= profileCGH$forceGL[2])
            profileCGH$profileValues$ZoneGNL[indexForceGain] <- 1
            indexForceLost <- which((profileCGH$profileValues$Smoothing-profileCGH$NormalRef) <= profileCGH$forceGL[1])
            profileCGH$profileValues$ZoneGNL[indexForceLost] <- -1

### Amplicon et deletion
            indexAmp <- which((profileCGH$profileValues$Smoothing-profileCGH$NormalRef) >= profileCGH$amplicon)
            profileCGH$profileValues$ZoneGNL[indexAmp] <- 2
            indexDel <- which((profileCGH$profileValues$Smoothing-profileCGH$NormalRef) <= profileCGH$deletion)
            profileCGH$profileValues$ZoneGNL[indexDel] <- -10
            
            profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"ZoneGen"))
            profileCGH$profileValues <- subset(profileCGH$profileValues, select=setdiff(names(profileCGH$profileValues),"NormalRange"))

                        

            class(profileCGH) <- "profileCGH"

            profileCGH$BkpInfo <- BkpInfo(profileCGH)
            
### Mise � jour du GNL des Outliers
            class(profileCGH) <- "profileCGH"

            if(assignGNLOut)
              {
                profileCGH <- OutliersGNL(profileCGH, alpha=profileCGH$alpha,
                                          sigma=profileCGH$SigmaG$Value, NormalRef=profileCGH$NormalRef,
                                          amplicon=profileCGH$amplicon, deletion=profileCGH$deletion)
              }

          }
        
      }

    return(profileCGH)

    
  }


