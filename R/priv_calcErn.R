
# #############################################################################
# Calculates gene enrichment: i-CisTarget version
# Arguments should be kept in the same order as .calcEnr_Aprox
.calcEnr_iCisTarget <- function(gsRankings, maxRank,
                                signifRankingNames, plotCurve, nCores, nMean)
{
  # nMean: ignored

  # Calculate RCC for all motifs
  rccs <- .calcRCC(gsRankings, maxRank, nCores)

  # Estimate mean and mean+2sd at each rank
  rccMean <- apply(rccs, 1, mean)
  rccM2sd <- rccMean + (2*apply(rccs, 1, sd))
  # TO DO: Save/return somehow?

  # Max enrichment for the selected rankings[,srName]
  maxEnr <- sapply(signifRankingNames, function(sr) {
    x <- min(which.max(rccs[,sr]-rccM2sd))
    c(x=x, y=unname(rccs[x,sr]))
  })

  # Plot
  if(plotCurve)
  {
    na <- sapply(colnames(maxEnr), function(srName) {
      .plotRCC(rccMean, rccM2sd, srName, maxEnr[,srName], rccs[,srName])
    })
  }
  return(maxEnr)
}
# 50-110 secs

##############################################################################
# Calculates gene enrichment: Aproximated/faster version
# Arguments should be kept in the same order as .calcEnr_iCisTarget
.calcEnr_Aprox <- function(gsRankings, maxRank,
                           signifRankingNames, plotCurve, nCores, nMean)
{
  # Calculate aproximated-RCC across all motifs at each rank position
  maxRankExtra <- maxRank+nMean
  gsRankings.asMat <- as.matrix(gsRankings) # Much faster!
  globalMat <- matrix(0, nrow=nrow(gsRankings), ncol=maxRankExtra)
  for(i in 1:nrow(gsRankings)) # (TO DO: Paralellize?)
  {
    x <- gsRankings.asMat[i,]
    x <- sort(x[x<maxRankExtra])
    
    if(length(x) > 0){
      coords <- cbind(y=seq_along(x), x)
      globalMat[coords] <- globalMat[coords]+1
    }
  }
  

  # Estimate mean and mean+2sd at each rank
  rccStatsRaw <- apply(globalMat, 2, function(x){
    tmp <- x
    if(sum(x)>0) tmp <- rep(seq_along(x), x)
    rccMean <- mean(tmp)
    rccSd <- sd(tmp)
    c(mean=rccMean, sd=rccSd)
  })
  # Remove NAs (assign left value)
  nas <- which(is.na(rccStatsRaw), arr.ind=TRUE)
  if(any(nas[,2] == 1)) {  # First value cannot be assigned to left
    rccStatsRaw[nas[which(nas[,2]==1), ]] <- 0
    nas <- nas[which(nas[,2]!=1), , drop=FALSE]
  }
  if(nrow(nas)>0){
    for(i in seq_len(nrow(nas)))
    {
      x <- nas[i,]
      rccStatsRaw[x[1], x[2]] <- rccStatsRaw[x[1], x[2]-1]
    }
  }


  # Reduce noise in the stats with the rolling mean
  rccStats <- t(apply(rccStatsRaw, 1,
                      function(x)
                        c(x[1:5],   # Correct??
                          zoo::rollmean(x, nMean, align="center",
                                        fill="extend"))))[,1:(maxRank-1)]
  rccM2sd <- rccStats["mean",] + (2*rccStats["sd",])
  rccMean <- rccStats["mean",]
  rm(rccStats); rm(globalMat)


  # Calculate real RCC & max enrichment for selected rankings
  rccs <- .calcRCC(gsRankings[signifRankingNames,,drop=FALSE], maxRank, nCores)
  maxEnr <- sapply(signifRankingNames, function(sr) {
    x <- min(which.max(rccs[,sr]-rccM2sd))
    c(x=x, y=unname(rccs[x,sr]))
  })

  # Plot
  if(plotCurve)
  {
    # Global estimation plot
    plot(rccStatsRaw["mean",]+2*rccStatsRaw["sd",],
         type="l", col="lightgreen",
         xlab="Rank", ylab="#genes recovered",
         main="Global mean and SD estimation", 
         xlim=c(0,maxRank))
    lines(rccStatsRaw["mean",],type="l", col="pink")
    lines(rccMean, col="red")
    lines(rccM2sd, col="darkgreen")

    # RCC for each significant ranking
    na <- sapply(colnames(maxEnr), function(srName) {
      .plotRCC(rccMean, rccM2sd, srName, maxEnr[,srName], rccs[,srName])
    })
  }
  return(maxEnr)
}

###############################################################################
# Aux functions

# Calculates RCC (of the gene-set) ONE RANKING
.calcRCC.oneRanking <- function(x, maxRank)
{
  x <- unlist(x)
  x <- sort(x[x<maxRank])

  curranking <- c(x, maxRank)
  unlist(mapply(rep, seq_along(curranking)-1,
                c(curranking[1], diff(curranking))))[-1]
}

# Apply .calcRCC.oneRanking on all rankings (= each column)
.calcRCC <- function(gsRankings, maxRank, nCores)
{
  if(nCores==1)
  {
    rccs <- apply(gsRankings, 1, .calcRCC.oneRanking, maxRank)
  }else
  {
    # Split rankings into 10 groups and run in parallel
    doParallel::registerDoParallel()
    options(cores=nCores)

    rowsNam <- rownames(gsRankings)
    suppressWarnings(
      rowNamsGroups <- split(rowsNam, (seq_along(rowsNam)) %% nCores)) 
        # Expected warning: Not multiple
    
    # rccs <- foreach(colsGr=colsNamsGroups, .combine="cbind") %do%
    rowsGr <- NULL
    rccs <- foreach::"%do%"(foreach::foreach(rowsGr=rowNamsGroups,
                                             .combine="cbind"),
    {
      apply(gsRankings[rowsGr,, drop=FALSE],1, .calcRCC.oneRanking, maxRank)
    })
  }
  return(rccs)
}
