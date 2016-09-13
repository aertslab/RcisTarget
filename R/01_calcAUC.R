
# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)
#' @import data.table
#'
#' @title Calculate AUC
#' @description Calculates the Area Under the Curve (AUC) of each gene-set for each motif ranking. This measure is used in the following steps to identify the DNA motifs that are significantly over-represented in the gene-set.
#' @param geneSets List of gene-sets to analyze. The gene-sets should be provided as a 'named list' in which each element is a gene-set (i.e. \code{list(geneSet1=c("gene1", "gene2"))})
#' @param rankings 'Motif rankings' database for the required organism and search-space (i.e. 10kbp around- or 500bp upstream the TSS).
#' These objects are provided in separate packages: \link[http://bioconductor.org/packages/RcisTarget.hg19.motifDatabases]{RcisTarget.hg19.motifDatabases} (Human), \link[http://bioconductor.org/packages/RcisTarget.mm9.motifDatabases]{RcisTarget.mm9.motifDatabases} (Mouse).
#' See the help files for more information: i.e. \code{help(RcisTarget.hg19.motifDatabases)}.
#' @param nCores Number of cores to use for computation.
#' Note: In general, using a higher number of cores (e.g. processes) decreases overall running time. However, it also deppends on the available memory and overall system load. Setting nCores too high might also decrease performance.
#' @param aucMaxRank Threshold to calculate the AUC.
#' In a simplified way, the AUC value represents the fraction of genes, within the top X genes in the ranking, that are included in the signature.
#' The parameter 'aucMaxRank' allows to modify the number of genes (maximum ranking) that is used to perform this computation.
#' By default it is set to 5\% of the total number of genes in the rankings. Common values range from 1 to 10\%.
#' See \code{vignette("RcisTarget")} for examples and more details.
#' @param verbose Should the function show progress messages? (TRUE / FALSE)
#' @return Matrix of motifs (row) by gene-sets (columns) with the value of AUC for each pair as content.
#' @seealso \code{vignette("RcisTarget")}. Next step in the workflow: \link[addMotifAnnotation]{addMotifAnnotation}
#' @example inst/examples/workflow_example.R
#' @export
calcAUC <- function(geneSets, rankings, nCores=1, aucMaxRank=0.05*nrow(rankings), verbose=TRUE) #, seed=123, plotHist=TRUE
{
  if(!is.list(geneSets)) stop("geneSets should be a named list.")
  if(is.null(names(geneSets))) stop("geneSets should be a named list.")
  if(nCores > length(geneSets)) nCores <- length(geneSets) # No point in using more...

  if(!is.data.table(rankings)) stop("Rankings should be a data.table (i.e. genes x [cells or motifs])")
  # if(!key(rankings) == "rn") stop("The rankings key should be 'rn'.")

  ######################################################################
  #### 1. Calculate the AUC for each gene set
  if(nCores==1)
  {
    aucMatrix <- sapply(names(geneSets), function(gSetName)  .AUC.geneSet(geneSet=geneSets[[gSetName]], rankings=rankings, aucMaxRank=aucMaxRank, gSetName=gSetName))
  }else
  {
    # Run each geneSet in parallel
    suppressMessages(require("doMC", quietly=TRUE))
    doMC::registerDoMC(nCores)
    if(verbose) message(paste("Using", getDoParWorkers(), "cores."))

    aucMatrix <- foreach(gSetName=names(geneSets)) %dopar%
    {
      setNames(list(.AUC.geneSet(geneSet=geneSets[[gSetName]], rankings=rankings, aucMaxRank=aucMaxRank, gSetName=gSetName)), gSetName)
    }
    aucMatrix <- do.call(cbind, unlist(aucMatrix, recursive = FALSE)[names(geneSets)])
  }

  ######################################################################
  ##### Messages for missing genes
  missinGenes <- aucMatrix[c("missing", "nGenes"), , drop=FALSE]
  missingPercent <- missinGenes["missing",, drop=FALSE]/missinGenes["nGenes",, drop=FALSE]
  if(all(missingPercent>=.80)) stop("Fewer than 20% of the genes in the gene sets are included in the rankings. Check wether the gene IDs in the 'rankings' and 'geneSets' match.")

  if(any(missingPercent>.80))
  {
    warning(paste("The following gene sets will be excluded from the analysis (less than 20% of their genes are available):\n",
     paste(names(missingPercent)[which(missingPercent >= .80)], collapse=", "), sep=""), immediate.=TRUE)
    aucMatrix <- aucMatrix[,which(missingPercent < .80),drop=FALSE]
  }

  if(sum(missinGenes["missing",])>0)
  {
    msg1 <- "Genes in the gene sets NOT available in the dataset: \n"
    msg2 <-  sapply(colnames(missinGenes)[which(missinGenes["missing",]>0)], function(gSetName)
        paste("\t", gSetName, ": \t", missinGenes["missing",gSetName],
          " (",round(missingPercent[,gSetName]*100),"% of ", missinGenes["nGenes",gSetName],")",sep=""))
    if(verbose) message(paste(msg1, paste(msg2, collapse="\n"), sep=""))
  }
  # (remove missing genes info from AUC matrix)
  aucMatrix <- aucMatrix[1:(nrow(aucMatrix)-2),, drop=FALSE]


  ######################################################################
  #### End: Return
  aucMatrix
}

.AUC.geneSet <- function(geneSet, rankings, aucMaxRank, gSetName="")  # add?: AUCThreshold
{
  geneSet <- unique(geneSet)
  nGenes <- length(geneSet)
  geneSet <- geneSet[which(geneSet %in% rankings$rn)]
  missing <- nGenes-length(geneSet)

  # stop(paste(c(sum(!geneSet %in% rankings$rn), class(geneSet)), collapse=", "))
  # gSetRanks <- rankings[geneSet][,-"rn", with=FALSE] # gene names are no longer needed
  gSetRanks <- subset(rankings, rn %in% geneSet)[,-"rn", with=FALSE] # gene names are no longer needed
  rm(rankings)

  aucThreshold <- round(aucMaxRank)            # 5%: same as .ini, but Stein said 3%
  maxAUC <- aucThreshold * nrow(gSetRanks)     # database.gene_count  ->  IS THIS CORRECT?

  # Apply by columns (i.e. to each ranking)
  auc <- sapply(gSetRanks, .calcAUC, aucThreshold, maxAUC)

  c(auc, missing=missing, nGenes=nGenes)
}


# oneRanking <- gSetRanks[,3, with=F]
.calcAUC <- function(oneRanking, aucThreshold, maxAUC)
{
  x <- unlist(oneRanking)

  x <- sort(x[x<aucThreshold])
  y <- 1:length(x)
  sum(diff(c(x, aucThreshold)) * y)/maxAUC
}
