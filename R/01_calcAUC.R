#' @import GSEABase
#' @importFrom methods new
#'
#' @title Calculate AUC
#' @description Calculates the Area Under the Curve (AUC) of each gene-set
#' for each motif ranking.
#' This measure is used in the following steps to identify the DNA motifs
#' that are significantly over-represented in the gene-set.
#' @param geneSets List of gene-sets to analyze.
#' The gene-sets should be provided as \code{\link{GeneSet}},
#' \code{\link{GeneSetCollection}} or character list (see examples).
#' @param rankings 'Motif rankings' database for the required organism and
#' search-space (i.e. 10kbp around- or 500bp upstream the TSS).
#' These objects are provided in separate files, 
#' which can be imported with \code{importRankings()}:
#' \itemize{
#' \item \url{http://pyscenic.aertslab.org/databases/mm9-500bp-upstream-7species.mc9nr.feather}[mm9-500bp-upstream-7species.mc9nr] (Mouse, 500bp)
#' \item \url{http://pyscenic.aertslab.org/databases/mm9-tss-centered-10kb-7species.mc9nr.feather}[mm9-tss-centered-10kb-7species.mc9nr] (Mouse, 10kbp)
#' \item \url{http://pyscenic.aertslab.org/databases/hg19-500bp-upstream-7species.mc9nr.feather}[hg19-500bp-upstream-7species.mc9nr] (Human, 500bp)
#' \item \url{http://pyscenic.aertslab.org/databases/hg19-tss-centered-10kb-7species.mc9nr.feather}[hg19-tss-centered-10kb-7species.mc9nr] (Human, 10kbp)
#' \item -Coming soon- (Fly)
#' }
#' See \code{vignette("RcisTarget")} for an exhaustive list of databases.
#' 
#' Since the normalized enrichment score (NES) of the motif
#' depends on the total number of motifs in the database,
#' we highly recommend to use the full version of the databases (20k motifs).
#' A smaller version of the human databases,
#' containing only the 4.6k motifs from cisbp,
#' are available in Bioconductor:
#' \itemize{
#' \item RcisTarget.hg19.motifDBs.cisbpOnly.500bp (Human)
#' }
#' @param nCores Number of cores to use for computation.
#' Note: In general, using a higher number of cores (e.g. processes)
#' decreases overall running time.
#' However, it also deppends on the available memory and overall system load.
#' Setting nCores too high might also decrease performance.
#' @param aucMaxRank Threshold to calculate the AUC.
#' In a simplified way, the AUC value represents the fraction of genes,
#' within the top X genes in the ranking, that are included in the signature.
#' The parameter 'aucMaxRank' allows to modify the number of genes
#' (maximum ranking) that is used to perform this computation.
#' By default it is set to 5\% of the total number of genes in the rankings.
#' Common values range from 1 to 10\%.
#' See \code{vignette("RcisTarget")} for examples and more details.
#' @param verbose Should the function show progress messages? (TRUE / FALSE)
#' @return \code{\link{aucScores}} of gene-sets (columns) by motifs (rows)
#' with the value of AUC for each pair as content.
#' @seealso Next step in the workflow: \code{\link{addMotifAnnotation}}.
#'
#' See the package vignette for examples and more details:
#' \code{vignette("RcisTarget")}
#' @example inst/examples/example_workflow.R
#' @rdname calcAUC
#' @export
setGeneric("calcAUC", signature="geneSets",
  function(geneSets, rankings, nCores=1,
      aucMaxRank=0.03*getNumColsInDB(rankings), verbose=TRUE)
  {
    standardGeneric("calcAUC")
  })

#' @rdname calcAUC
#' @aliases calcAUC,list-method
setMethod("calcAUC", "list",
  function(geneSets, rankings, nCores=1,
           aucMaxRank=0.03*getNumColsInDB(rankings), verbose=TRUE)
  {
    .RcisTarget_calcAUC(geneSets=geneSets,
                        rankings=rankings,
                        nCores=nCores,
                        aucMaxRank=aucMaxRank,
                        verbose=verbose)
  })

#' @rdname calcAUC
#' @aliases calcAUC,character-method
setMethod("calcAUC", "character",
  function(geneSets, rankings, nCores=1,
           aucMaxRank=0.03*getNumColsInDB(rankings), verbose=TRUE)
  {
    geneSets <- list(geneSet=geneSets)

    .RcisTarget_calcAUC(geneSets=geneSets,
                        rankings=rankings,
                        nCores=nCores,
                        aucMaxRank=aucMaxRank,
                        verbose=verbose)
  })

#' @rdname calcAUC
#' @aliases calcAUC,GeneSet-method
setMethod("calcAUC", "GeneSet",
  function(geneSets, rankings, nCores=1,
           aucMaxRank=0.03*getNumColsInDB(rankings), verbose=TRUE)
  {
    geneSets <- setNames(list(GSEABase::geneIds(geneSets)),
                         GSEABase::setName(geneSets))

    .RcisTarget_calcAUC(geneSets=geneSets,
                        rankings=rankings,
                        nCores=nCores,
                        aucMaxRank=aucMaxRank,
                        verbose=verbose)
  })

#' @rdname calcAUC
#' @aliases calcAUC,GeneSetCollection-method
setMethod("calcAUC", "GeneSetCollection",
  function(geneSets, rankings, nCores=1,
           aucMaxRank=0.03*getNumColsInDB(rankings), verbose=TRUE)
  {
    geneSets <- GSEABase::geneIds(geneSets)

    .RcisTarget_calcAUC(geneSets=geneSets,
                        rankings=rankings,
                        nCores=nCores,
                        aucMaxRank=aucMaxRank,
                        verbose=verbose)
  })

.RcisTarget_calcAUC <- function(geneSets, rankings, nCores=1,
                                aucMaxRank=0.03*getNumColsInDB(rankings), verbose=TRUE)
{
  # Check the gene sets
  if(!is.list(geneSets))
    stop("geneSets should be a named list.")
  if(is.null(names(geneSets)))
    stop("geneSets should be a named list.")
  if(nCores > length(geneSets))
    nCores <- length(geneSets) # No point in using more...
  
  # Check the rankings
  if(! "rankingRcisTarget" %in% class(rankings))
    stop("The rankings should be of class rankingRcisTarget.")
  
  if(getMaxRank(rankings) < Inf)
  {
    if(aucMaxRank > getMaxRank(rankings))
      stop("aucMaxRank (", aucMaxRank,
        ") should not be bigger than the maximum ranking available ",
        "in the database (", getMaxRank(rankings),")")
  }

  rankingsInfo <- c(org="", genome="", description="")
  if(isS4(rankings)) {
    rankingsInfo <- c(org=rankings@org,
                      genome=rankings@genome,
                      description=rankings@description)
    rankings <- getRanking(rankings)
    
    # could be added as argument...
    if(!"features" %in% colnames(rankings)) 
      stop("No feature names in the ranking (column 'features')")
  }
  if(!is.data.frame(rankings))
    stop("Rankings does not have the right format.")

  # Check if the genes are included in the rankings
  allGenes <- unique(unlist(geneSets))
  if(sum(allGenes %in% colnames(rankings))/length(allGenes) < .80)
    stop("Fewer than 80% of the genes/features in the gene-sets ",
     "are included in the rankings.",
     "Check wether the IDs in the 'rankings' (columns) and 'geneSets' match.")


  ######################################################################
  #### 1. Calculate the AUC for each gene set
  if(nCores==1)
  {
    gSetName <- NULL
    aucMatrix <- vapply(names(geneSets), function(gSetName)
      .AUC.geneSet(geneSet=geneSets[[gSetName]],
                   rankings=rankings[,-1], # featureName
                   aucMaxRank=aucMaxRank,
                   gSetName=gSetName),
      FUN.VALUE=numeric(nrow(rankings)+2))
    aucMatrix <- t(aucMatrix)
    colnames(aucMatrix)[1:(ncol(aucMatrix)-2)]<-as.character(rankings$features)
  }else
  {
    # Run each geneSet in parallel
    doParallel::registerDoParallel()
    options(cores=nCores)

    if(verbose)
      message("Using ", foreach::getDoParWorkers(), " cores.")

    # aucMatrix <- foreach(gSetName=names(geneSets)) %dopar%
    gSetName <- NULL
    "%dopar%"<- foreach::"%dopar%"
    aucMatrix <- foreach::"%dopar%"(foreach::foreach(gSetName=names(geneSets)),
    {
      setNames(list(.AUC.geneSet(geneSet=geneSets[[gSetName]],
                                 rankings=rankings[,-1],  # featureName
                                 aucMaxRank=aucMaxRank,
                                 gSetName=gSetName)), gSetName)
    })
    aucMatrix <- do.call(rbind,
            unlist(aucMatrix, recursive = FALSE)[names(geneSets)])
    colnames(aucMatrix)[1:(ncol(aucMatrix)-2)]<-as.character(rankings$features)
  }

  ######################################################################
  ##### Messages for missing genes
  missingGenes <- as.matrix(aucMatrix[,c("missing", "nGenes") , drop=FALSE])
  missingPercent <- as.numeric(missingGenes[,"missing"])/
    as.numeric(missingGenes[,"nGenes"])
  missingPercent <- setNames(missingPercent, rownames(missingGenes))

  if(all(missingPercent>=.80))
    stop("Fewer than 20% of the genes in the gene sets",
         " are included in the rankings.",
         "Check wether the gene IDs in the 'rankings' and 'geneSets' match.")

  if(any(missingPercent>.80))
  {
    warning("The following gene sets will be excluded from the analysis",
            "(less than 20% of their genes are available):\n",
            paste(names(missingPercent)[which(missingPercent >= .80)],
                  collapse=", "),
            immediate.=TRUE)
    aucMatrix <- aucMatrix[which(missingPercent < .80),,drop=FALSE]
  }

  if(sum(missingGenes[,"missing"])>0)
  {
    msg1 <- "Genes in the gene sets NOT available in the dataset: \n"
    msg2 <-  vapply(rownames(missingGenes)[which(missingGenes[,"missing"]>0)],
              function(gSetName)
                paste("\t", gSetName, ": \t", missingGenes[gSetName,"missing"],
                " (",round(missingPercent[gSetName]*100),"% of ",
                missingGenes[gSetName,"nGenes"],")",sep=""),
              FUN.VALUE="")
    if(verbose)
      message(paste(msg1, paste(msg2, collapse="\n"), sep=""))
  }
  # (remove missing genes info from AUC matrix)
  aucMatrix <- aucMatrix[,1:(ncol(aucMatrix)-2), drop=FALSE]


  ######################################################################
  #### End: Return
  names(dimnames(aucMatrix)) <- c("gene-set", "motif")
  new("aucScores",
      SummarizedExperiment::SummarizedExperiment(assays=list(AUC=aucMatrix)),
      org = rankingsInfo["org"],
      genome = rankingsInfo["genome"],
      description = rankingsInfo["description"])
}

.AUC.geneSet <- function(geneSet, rankings, aucMaxRank, gSetName="")
{
  geneSet <- unique(geneSet)
  nGenes <- length(geneSet)
  geneSet <- geneSet[which(geneSet %in% colnames(rankings))]
  missing <- nGenes-length(geneSet)

  # gene names are no longer needed
  gSetRanks <- rankings[,geneSet]
  rm(rankings)

  aucThreshold <- round(aucMaxRank)
  maxAUC <- aucThreshold * ncol(gSetRanks)  

  auc <- apply(gSetRanks, 1, .auc, aucThreshold, maxAUC)

  c(auc, missing=missing, nGenes=nGenes)
}

# oneRanking <- gSetRanks[,3, with=FALSE]
.auc <- function(oneRanking, aucThreshold, maxAUC)
{
  x <- as.numeric(oneRanking)
  x <- sort(x[x<aucThreshold])

  y <- seq_along(x)
  sum(diff(c(x, aucThreshold)) * y)/maxAUC
}
