# Max ranking to calculate the enrichment. By default it is the same as to calculate AUC... but...

# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)
# @import data.table
#'
#' @title to do
#' @description to do
#' @param to do
#' @return to do
#' @details # required?
# @example # To do
#' @export
addSignificantGenes <- function(aucTable, geneSets, rankings, maxRank=5000, plotCurve=FALSE, genesFormat="geneList", method="aprox", nCores=1)
{
  # suppressPackageStartupMessages(library(data.table))
  method <- tolower(method[1])
  if(!method %in% c("icistarget", "icistargetaprox", "aprox")) stop("'method' should be either 'iCisTarget' or 'iCisTargetAprox'.")

  geneSetNames <- unique(aucTable$geneSet)

  rnkType <- c("ranking", "motif")
  rnkType <- rnkType[which(rnkType %in% colnames(aucTable))]

  signifMotifsAsList <- lapply(geneSetNames, function(gsn) # (Paralelized inside enrichment function)
  {
    enrRnkT_ByGs <- aucTable[geneSet==gsn]

    signifGenes <- getSignificantGenes(geneSet=geneSets[[unique(enrRnkT_ByGs$geneSet)]],
                                       rankings=rankings,
                                       signifRankingNames=unname(unlist(subset(enrRnkT_ByGs, select=rnkType))),
                                       method=method,
                                       maxRank=maxRank,
                                       plotCurve=plotCurve,
                                       genesFormat=genesFormat,
                                       nCores=nCores, digits=3)

    enrRnkT_ByGs <- cbind(enrRnkT_ByGs, signifGenes$enrStats)
    if("geneList" %in% genesFormat) enrRnkT_ByGs <- cbind(enrRnkT_ByGs, enrichedGenes=sapply(signifGenes$enrichedGenes, function(x) paste(unlist(x), collapse=";")))
    if("incidMatrix" %in% genesFormat) enrRnkT_ByGs <- cbind(enrRnkT_ByGs, signifGenes$incidMatrix)
    enrRnkT_ByGs
  })

  rbindlist(signifMotifsAsList)
}

# geneSet <- geneLists$hypoxia
# rankings <- motifRankings
# signifRankingNames <- subset(motifEnrichmentTable, geneSet=="hypoxia")$motif[1:4]
# signifRankingNames="homer__GCACGTACCC_HIF2a"
# method="aprox"
# maxRank= 5000 #0.05*nrow(rankings)
# plotCurve=FALSE
# genesFormat="geneList"
# digits=3

# signifGenes <- getSignificantGenes(geneLists$hypoxia, motifRankings,
#                                    signifRankingNames=signifRankingNames, genesFormat="geneList", method="aprox")
# TO DO: Faster?
# method: "aprox", is a quick scan function, faster if only for a few rankings. If for many, maybe better iCisTarget



#' @import data.table
#'
#' @title to do
#' @description to do
#' @param to do
#' @return to do
#' @examples # to do
#' @export
getSignificantGenes <- function(geneSet, rankings, signifRankingNames=NULL, method="iCisTarget", maxRank=5000, plotCurve=FALSE, genesFormat=c("geneList", "incidMatrix"), nCores=1, digits=3)
{
  ################################################################################
  # Argument checks & init. vars
  # aucThreshold <- 0.05*nrow(rankings)
  maxRank <- round(maxRank)
  if(is.null(signifRankingNames)) {
    signifRankingNames <- colnames(rankings)
    warning("'signifRankingNames' has not been provided. The significant genes will be calculated for all rankings.")
  }
  signifRankingNames <- unname(signifRankingNames)
  if(!all(genesFormat %in% c("geneList", "incidMatrix", "none"))) stop('"genesFormat" should be ither "geneList" and/or "incidMatrix".')

  method <- tolower(method[1])
  if(!method %in% c("icistarget", "icistargetaprox", "aprox")) stop("'method' should be either 'iCisTarget' or 'iCisTargetAprox'.")
  if(method == "icistarget") {
    calcEnrFunct <- .calcErn_iCisTarget
  } else {
    calcEnrFunct <- .calcErn_Aprox
    if(!"zoo" %in% rownames(installed.packages())) stop("Package 'zoo' is required to calculate the aproximate RCC distributions. To install it, run:\t install.packages('zoo')")
  }

  # Remove missing genes from geneSet...
  geneSet <- unique(geneSet)
  geneSet <- geneSet[which(geneSet %in% rankings$rn)]

  gSetRanks <- subset(rankings, rn %in% geneSet)
  # gSetRanks <- rankings[geneSet]
  geneSetIndex <- match(geneSet, gSetRanks$rn)
  rm(rankings)

  ################################################################################
  # Calculate enrichment
  enrStats <- t(calcEnrFunct(gSetRanks[,-"rn", with=FALSE], maxRank, signifRankingNames, plotCurve, nCores))
  enrStats <- enrStats[,c("y", "x"), drop=FALSE]
  colnames(enrStats) <- c("nErnGenes", "rankAtMax")

  enrichedGenes <- list()
  if("incidMatrix" %in% genesFormat) incidMatrix <- matrix(0, nrow=length(signifRankingNames), ncol=length(geneSet),
                                                              dimnames=list(signifRankingNames, sort(geneSet)))
  for(rankingName in signifRankingNames)
  {
    enrichedGenes[[rankingName]] <- sort(gSetRanks$rn[which(gSetRanks[,rankingName,with=F] <= enrStats[rankingName, "rankAtMax"])])
    if("incidMatrix" %in% genesFormat) incidMatrix[rankingName, enrichedGenes[[rankingName]]] <- 1
  }

  ################################################################################
  # Return
  ret <- list(enrStats=enrStats) #data.table(enrStats, keep.rownames=TRUE))
  if("geneList" %in% genesFormat) ret <- c(ret, enrichedGenes=list(enrichedGenes))
  if("incidMatrix" %in% genesFormat) ret <- c(ret, incidMatrix=list(incidMatrix))
  if(tolower(genesFormat) != "none") return(ret)
}
