# Max ranking to calculate the enrichment

# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)

#' @import data.table
##' @import GSEABase
#' @importFrom methods new
#' @importFrom graphics lines plot points polygon
#' @importFrom utils installed.packages
#'
#' @title Add significant genes
#' @description Identify which genes (of the gene-set) are highly ranked
#' for each motif.
#'
#' \itemize{
#'   \item addSignificantGenes(): adds them to the results table.
#'   \item getSignificantGenes():
#'   Calculates the significant genes for ONE gene set.
#'   It provides the plot and the gene list (it is used by addSignificantGenes).
#' }
#' @param resultsTable [addSignificantGenes]
#' Output table from \code{\link{addMotifAnnotation}}
#' @param geneSets [addSignificantGenes] List of gene-sets which was analyzed.
#' @param rankings Motif rankings used to analyze the gene list
#' (They should be the same as used for calcAUC in this same analysis).
#' @param maxRank Maximum rank to take into account for the recovery curve
#' (Default: 5000).
#' @param plotCurve Logical. Wether to plot the recovery curve (Default: FALSE).
#' @param genesFormat "geneList" or "incidMatrix". Format to return the genes
#' (Default: "geneList").
#' @param method "iCisTarget" or "aprox". There are two methods to identify the
#'  highly ranked genes:
#'   (1) equivalent to the ones used in iRegulon and i-cisTarget
#'   (method="iCisTarget", recommended if running time is not an issue),
#'   and (2) a faster implementation based on an approximate distribution
#'   using the average at each rank (method="aprox",
#'   useful to scan multiple gene sets). (Default: "aprox")
#' @param nMean Only used for "aprox" method: Interval to calculate the running
#' mean and sd. Default: 20 (aprox. nGenesInRanking/1000).
#' @param nCores Number of cores to use for parallelization (Default: 1).
#' @param geneSet [getSignificantGenes] Gene-set to analyze (Only one).
#' @param signifRankingNames [getSignificantGenes] Motif ranking name.
#' @param digits  [getSignificantGenes]
#' Number of digits to include in the output.
#' @return Output from \code{\link{addMotifAnnotation}}
#' adding the folowing columns:
#' \itemize{
#'   \item nEnrGenes: Number of genes highly ranked
#'   \item rankAtMax: Ranking at the maximum enrichment,
#'   used to determine the number of enriched genes.
#'   \item enrichedGenes: Genes that are highly ranked for the given motif.
#'   If genesFormat="geneList", the gene names are collapsed into a comma
#'   separated text field (alphabetical order). If genesFormat="incidMatrix",
#'   they are formatted as an indicence matrix, i.e. indicanting with 1 the
#'   genes present, and 0 absent.
#' }
#' @return If plotCurve=TRUE, the recovery curve is plotted.
#' @details
#' The highly ranked genes are selected based on the distribution of the
#' recovery curves of the gene set across all the motifs in the database.
#' In the plot, the red line indicates the average of the recovery curves of
#' all the motifs, the green line the average + standard deviation, and the
#' blue line the recovery curve of the current motif.
#' The point of maximum distance between the current motif and the green curve
#' (mean+sd), is the rank selected as maximum enrichment.
#' All the genes with lower rank will be considered enriched.
#'
#' Depending on whether the method is "iCisTarget" or "aprox", the mean and
#' SD at each rank are calculated slightly different.
#' "iCisTarget" method calculates the recovery curves for all the motifs, and
#'  then calculates the average and SD at each rank.
#' Due to the implementation of the function in R, this method is slower than
#' just subsetting the ranks of the genes in for each motif,
#' and calculating the average of the available ones at each position with a
#' sliding window.
#' Since there are over 18k motifs, the chances of getting several measures at
#' each rank are very high and highly resemble the results calculated
#' by iCisTarget, though they are often not exactly the same
#' (hence the name: "aprox" method).

#' @seealso Previous step in the workflow: \code{\link{addMotifAnnotation}}.
#'
#' See the package vignette for examples and more details:
#' \code{vignette("RcisTarget")}
#' @example inst/examples/example_addSignificantGenes.R
#' @export
#'
#' @export
setGeneric("addSignificantGenes", signature="geneSets",
  function(resultsTable, geneSets, rankings, maxRank=5000, plotCurve=FALSE,
    genesFormat="geneList", method="aprox", nMean=20, nCores=1)
  {
    standardGeneric("addSignificantGenes")
  })

#' @rdname addSignificantGenes
#' @aliases addSignificantGenes,list-method
setMethod("addSignificantGenes", "list",
  function(resultsTable, geneSets, rankings, maxRank=5000, plotCurve=FALSE,
    genesFormat="geneList", method="aprox", nMean=20, nCores=1)
  {
    .addSignificantGenes(resultsTable=resultsTable,
                         geneSets=geneSets,
                         rankings=rankings,
                         maxRank=maxRank,
                         plotCurve=plotCurve,
                         genesFormat=genesFormat,
                         method=method,
                         nMean=nMean,
                         nCores=nCores)
  })

#' @rdname addSignificantGenes
#' @aliases addSignificantGenes,character-method
setMethod("addSignificantGenes", "character",
  function(resultsTable, geneSets, rankings, maxRank=5000, plotCurve=FALSE,
    genesFormat="geneList", method="aprox", nMean=20, nCores=1)
  {
    geneSets <- list(geneSet=geneSets)

    .addSignificantGenes(resultsTable=resultsTable,
                         geneSets=geneSets,
                         rankings=rankings,
                         maxRank=maxRank,
                         plotCurve=plotCurve,
                         genesFormat=genesFormat,
                         method=method,
                         nMean=nMean,
                         nCores=nCores)
  })

# #' @rdname addSignificantGenes
# #' @aliases addSignificantGenes,GeneSet-method
# setMethod("addSignificantGenes", "GeneSet",
#   function(resultsTable, geneSets, rankings, maxRank=5000, plotCurve=FALSE,
#     genesFormat="geneList", method="aprox", nMean=20, nCores=1)
#   {
#     geneSets <- setNames(list(GSEABase::geneIds(geneSets)),
#                          GSEABase::setName(geneSets))
#
#     .addSignificantGenes(resultsTable=resultsTable,
#                          geneSets=geneSets,
#                          rankings=rankings,
#                          maxRank=maxRank,
#                          plotCurve=plotCurve,
#                          genesFormat=genesFormat,
#                          method=method,
#                          nMean=nMean,
#                          nCores=nCores)
#   })

# #' @rdname addSignificantGenes
# #' @aliases addSignificantGenes,GeneSetCollection-method
# setMethod("addSignificantGenes", "GeneSetCollection",
#   function(resultsTable, geneSets, rankings, maxRank=5000,
#     plotCurve=FALSE, genesFormat="geneList", method="aprox", nMean=20, nCores=1)
#   {
#     geneSets <- GSEABase::geneIds(geneSets)
#
#     .addSignificantGenes(resultsTable=resultsTable,
#                          geneSets=geneSets,
#                          rankings=rankings,
#                          maxRank=maxRank,
#                          plotCurve=plotCurve,
#                          genesFormat=genesFormat,
#                          method=method,
#                          nMean=nMean,
#                          nCores=nCores)
#   })

.addSignificantGenes <- function(resultsTable, geneSets, rankings,
                         maxRank=5000, plotCurve=FALSE, genesFormat="geneList",
                         method="aprox", nMean=20, nCores=1)
{
  if(isS4(rankings)) {
    if(getMaxRank(rankings) < Inf)
    {
      if(maxRank > getMaxRank(rankings))
        stop("maxRank (", maxRank, ") should not be bigger ",
             "than the maximum ranking available in the database (",
             getMaxRank(rankings),")")
    }

    rankings <- rankings@rankings
  }

  # suppressPackageStartupMessages(library(data.table))
  method <- tolower(method[1])
  if(!method %in% c("icistarget", "icistargetaprox", "aprox"))
    stop("'method' should be either 'iCisTarget' or 'iCisTargetAprox'.")

  geneSetNames <- unique(resultsTable$geneSet)

  rnkType <- c("ranking", "motif")
  rnkType <- rnkType[which(rnkType %in% colnames(resultsTable))]

  # (Paralelized inside enrichment function)
  signifMotifsAsList <- lapply(geneSetNames, function(gsn)
  {
    enrRnkT_ByGs <- resultsTable[resultsTable$geneSet==gsn]

    signifGenes <- .getSignificantGenes(
                       geneSet=as.character(
                         geneSets[[unique(enrRnkT_ByGs$geneSet)]]),
                       rankings=rankings,
                       signifRankingNames=unname(unlist(subset(enrRnkT_ByGs,
                                                            select=rnkType))),
                       method=method,
                       maxRank=maxRank,
                       plotCurve=plotCurve,
                       genesFormat=genesFormat,
                       nCores=nCores, digits=3, nMean=nMean)

    enrRnkT_ByGs <- cbind(enrRnkT_ByGs, signifGenes$enrStats)
    if("geneList" %in% genesFormat)
      enrRnkT_ByGs <- cbind(enrRnkT_ByGs,
                            enrichedGenes=sapply(signifGenes$enrichedGenes,
                                  function(x) paste(unlist(x), collapse=";")))
    if("incidMatrix" %in% genesFormat)
      enrRnkT_ByGs <- cbind(enrRnkT_ByGs, signifGenes$incidMatrix)
    enrRnkT_ByGs
  })

  rbindlist(signifMotifsAsList)
}

#' @import data.table
#' @rdname addSignificantGenes
#' @export
setGeneric("getSignificantGenes", signature="geneSet",
  function(geneSet, rankings, signifRankingNames=NULL, method="iCisTarget",
    maxRank=5000, plotCurve=FALSE, genesFormat=c("geneList", "incidMatrix"),
    nCores=1, digits=3, nMean=10)
  {
    standardGeneric("getSignificantGenes")
  })

#' @rdname addSignificantGenes
#' @aliases getSignificantGenes,list-method
setMethod("getSignificantGenes", "list",
  function(geneSet, rankings, signifRankingNames=NULL, method="iCisTarget",
    maxRank=5000, plotCurve=FALSE, genesFormat=c("geneList", "incidMatrix"),
    nCores=1, digits=3, nMean=10)
  {
    if(length(geneSet)>1) stop("Provide only one gene set.")
    geneSet <- as.character(unname(unlist(geneSet)))
    .getSignificantGenes(geneSet=geneSet,
                         rankings=rankings,
                         signifRankingNames=signifRankingNames,
                         method=method,
                         maxRank=maxRank,
                         plotCurve=plotCurve,
                         genesFormat=genesFormat,
                         nCores=nCores,
                         digits=digits,
                         nMean=nMean)
  })

#' @rdname addSignificantGenes
#' @aliases getSignificantGenes,character-method
setMethod("getSignificantGenes", "character",
  function(geneSet, rankings, signifRankingNames=NULL, method="iCisTarget",
    maxRank=5000, plotCurve=FALSE, genesFormat=c("geneList", "incidMatrix"),
    nCores=1, digits=3, nMean=10)
  {
    .getSignificantGenes(geneSet=geneSet,
                         rankings=rankings,
                         signifRankingNames=signifRankingNames,
                         method=method,
                         maxRank=maxRank,
                         plotCurve=plotCurve,
                         genesFormat=genesFormat,
                         nCores=nCores,
                         digits=digits,
                         nMean=nMean)
  })

#' @rdname addSignificantGenes
#' @aliases getSignificantGenes,factor-method
setMethod("getSignificantGenes", "factor",
  function(geneSet, rankings, signifRankingNames=NULL, method="iCisTarget",
    maxRank=5000, plotCurve=FALSE, genesFormat=c("geneList", "incidMatrix"),
    nCores=1, digits=3, nMean=10)
  {
    geneSet <- as.character(geneSet)
    .getSignificantGenes(geneSet=geneSet,
                         rankings=rankings,
                         signifRankingNames=signifRankingNames,
                         method=method,
                         maxRank=maxRank,
                         plotCurve=plotCurve,
                         genesFormat=genesFormat,
                         nCores=nCores,
                         digits=digits,
                         nMean=nMean)
  })

# #' @rdname addSignificantGenes
# #' @aliases getSignificantGenes,GeneSet-method
# setMethod("getSignificantGenes", "GeneSet",
#   function(geneSet, rankings, signifRankingNames=NULL, method="iCisTarget",
#     maxRank=5000, plotCurve=FALSE, genesFormat=c("geneList", "incidMatrix"),
#     nCores=1, digits=3, nMean=10)
#   {
#     geneSet <- GSEABase::geneIds(geneSet)
#
#     .getSignificantGenes(geneSet=geneSet,
#                          rankings=rankings,
#                          signifRankingNames=signifRankingNames,
#                          method=method,
#                          maxRank=maxRank,
#                          plotCurve=plotCurve,
#                          genesFormat=genesFormat,
#                          nCores=nCores,
#                          digits=digits,
#                          nMean=nMean)
#   })

# #' @rdname addSignificantGenes
# #' @aliases getSignificantGenes,GeneSetCollection-method
# setMethod("getSignificantGenes", "GeneSetCollection",
#   function(geneSet, rankings, signifRankingNames=NULL, method="iCisTarget",
#     maxRank=5000, plotCurve=FALSE, genesFormat=c("geneList", "incidMatrix"),
#     nCores=1, digits=3, nMean=10)
#   {
#     if(length(geneSet)>1) stop("Provide only one gene set.")
#     geneSet <- unlist(GSEABase::geneIds(geneSet))
#
#     .getSignificantGenes(geneSet=geneSet,
#                          rankings=rankings,
#                          signifRankingNames=signifRankingNames,
#                          method=method,
#                          maxRank=maxRank,
#                          plotCurve=plotCurve,
#                          genesFormat=genesFormat,
#                          nCores=nCores,
#                          digits=digits,
#                          nMean=nMean)
#   })
.getSignificantGenes <- function(geneSet,
                                 rankings,
                                 signifRankingNames=NULL,
                                 method="iCisTarget",
                                 maxRank=5000,
                                 plotCurve=FALSE,
                                 genesFormat=c("geneList", "incidMatrix"),
                                 nCores=1,
                                 digits=3,
                                 nMean=10)
{
  ############################################################################
  # Argument checks & init. vars
  # aucThreshold <- 0.05*nrow(rankings)

  maxRank <- round(maxRank)
  if(isS4(rankings)) {
    if(getMaxRank(rankings) < Inf)
    {
      if(maxRank > getMaxRank(rankings))
        stop("maxRank (", maxRank, ") should not be bigger ",
             "than the maximum ranking available in the database (",
             getMaxRank(rankings),")")
    }
  }
  if(isS4(rankings)) rankings <- rankings@rankings


  if(is.null(signifRankingNames)) {
    signifRankingNames <- colnames(rankings)
    warning("'signifRankingNames' has not been provided.",
            "The significant genes will be calculated for all rankings.")
  }

  signifRankingNames <- unname(signifRankingNames)
  if(!all(genesFormat %in% c("geneList", "incidMatrix", "none")))
    stop('"genesFormat" should be ither "geneList" and/or "incidMatrix".')

  method <- tolower(method[1])
  if(!method %in% c("icistarget", "icistargetaprox", "aprox"))
    stop("'method' should be either 'iCisTarget' or 'iCisTargetAprox'.")
  if(method == "icistarget") {
    calcEnrFunct <- .calcEnr_iCisTarget
  } else {
    calcEnrFunct <- .calcEnr_Aprox
    if(!"zoo" %in% rownames(utils::installed.packages()))
      stop("Package 'zoo' is required ",
           "to calculate the aproximate RCC distributions.",
           "To install it, run:\t install.packages('zoo')")
  }

  # Remove missing genes from geneSet...
  geneSet <- unique(geneSet)
  geneSet <- geneSet[which(geneSet %in% rankings$rn)]

  gSetRanks <- subset(rankings, rankings$rn %in% geneSet)
  # gSetRanks <- rankings[geneSet]
  geneSetIndex <- match(geneSet, gSetRanks$rn)
  rm(rankings)

  #############################################################################
  # Calculate enrichment
  enrStats <- t(calcEnrFunct(gSetRanks[,-"rn", with=FALSE],
                             maxRank,
                             signifRankingNames,
                             plotCurve,
                             nCores,
                             nMean))
  enrStats <- enrStats[,c("y", "x"), drop=FALSE]
  colnames(enrStats) <- c("nEnrGenes", "rankAtMax")

  enrichedGenes <- list()
  if("incidMatrix" %in% genesFormat)
    incidMatrix <- matrix(0, nrow=length(signifRankingNames),
                          ncol=length(geneSet),
                          dimnames=list(signifRankingNames, sort(geneSet)))
  for(rankingName in signifRankingNames)
  {
    enrichedGenes[[rankingName]] <- sort(gSetRanks$rn[which(
      gSetRanks[,rankingName,with=FALSE] <=
        enrStats[rankingName, "rankAtMax"])])
    if("incidMatrix" %in% genesFormat)
      incidMatrix[rankingName, enrichedGenes[[rankingName]]] <- 1
  }

  #############################################################################
  # Return
  ret <- list(enrStats=enrStats) #data.table(enrStats, keep.rownames=TRUE))
  if("geneList" %in% genesFormat)
    ret <- c(ret, enrichedGenes=list(enrichedGenes))
  if("incidMatrix" %in% genesFormat)
    ret <- c(ret, incidMatrix=list(incidMatrix))
  if(any(tolower(genesFormat) != "none"))
    return(ret)
}

