# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)
#' @import data.table
#'
#' @title Add motif annotation
#' @description Select significant motifs and/or annotate motifs to transcription factors.
#' The motifs are considered significantly enriched if they pass the the Normalized Enrichment Score (NES) threshold.
#' @param auc Output from calcAUC.
#' @param nesThreshold Numeric. NES threshold to calculate the motif significant (3.0 by default). The NES is calculated -for each motif- based on the AUC distribution of all the motifs for the gene-set [(x-mean)/sd].
#' @param digits Integer. Number of digits for the AUC and NES in the output table.
#' @param motifAnnot_direct Motif annotation database containing DIRECT annotations of the motif to transcription factors. The names should match the ranking column names.
#' @param motifAnnot_indirect Motif annotation database containing the expanded annotations of the motif to transcription factors based on motif similarity.
#' @param highlightTFs Character. If a list of transcription factors is provided, the column TFinDB in the otuput table will indicate whether any of those TFs are included within the direct annotation (two asterisks, **) or indirect annotation (one asterisk, *) of the motif.
#' The vector can be named to indicate which TF to highlight for each gene-set. Otherwise, all TFs will be used for all geneSets.
#' @return \code{\link[data.table]{data.table}} with the folowing columns:
#' \itemize{
#' \item geneSet: Name of the gene set
#' \item motif: ID of the motif (colnames of the ranking, it might be other kind of feature)
#' \item NES: Normalized enrichment score of the motif in the gene-set
#' \item AUC: Area Under the Curve (used to calculate the NES)
#' \item TFinDB: Indicates whether the highlightedTFs are included within the direct annotation (two asterisks, **) or indirect annotation (one asterisk, *)
#' \item TF_direct: Transcription factors annotated to the motif according to 'direct annotation'.
#' \item TF_indirect: Transcription factors annotated to the motif according to 'indirect annotation'.
#' }
#' @seealso Next step in the workflow: \code{\link{addSignificantGenes}}.
#'
#' Previous step in the workflow: \code{\link{calcAUC}}.
#'
#' See the package vignette for examples and more details: \code{vignette("RcisTarget")}
#' @example inst/examples/example_addMotifAnnotation.R
#' @export
addMotifAnnotation <- function(auc, nesThreshold=3.0, digits=3, motifAnnot_direct=NULL, motifAnnot_indirect=NULL, highlightTFs=NULL)
{
  auc <- getAuc(auc)
  #### Check inputs
  if(!is.null(highlightTFs))
  {
    if(is.null(motifAnnot_direct) && is.null(motifAnnot_indirect)) stop("To hightlight TFs, please provide a motif-TF annotation.")
    if(is.null(names(highlightTFs))) {
      warning("The input TFs are not named, all TFs will be used with all Gene Sets.")
      highlightTFs <- setNames(rep(list(highlightTFs), nrow(auc)), rownames(auc))
    }

    if(!all(names(highlightTFs) %in% rownames(auc))) warning("TFs 1")
    if(!all(rownames(auc) %in% names(highlightTFs))) warning("TFs 2")
  }

  #### Runs "auc.asTable" on each AUC columns i.e. signatures/cells
  ret <- lapply(rownames(auc), function(geneSet) {
      tfs <- highlightTFs[[geneSet]]
      aucTable <- .auc.asTable(auc[geneSet,], nesThreshold=nesThreshold, digits=digits)
      if(nrow(aucTable)>0)
      {
        aucTable <- .addTfs(aucTable, motifAnnot_direct=motifAnnot_direct, motifAnnot_indirect=motifAnnot_indirect, highlightTFs=tfs)
        aucTable <- data.table(geneSet=geneSet, aucTable)
      }else{
        aucTable <- NULL
      }
      aucTable
    })

  #### Merge the results from each signature/geneSet/regionSet into a single data.table
  # ret <- do.call(rbind, unname(ret))  # Slower?
  # library(data.table)
  ret <- data.table::rbindlist(ret)

  colnames(ret)[which(colnames(ret) == "ranking")] <- "motif"
  return(ret)
}


############ PRIVATE
.calcNES <- function(AUC)
{
  meanAUC <- mean(AUC)
  sdAUC <- sd(AUC)

  # NES = (AUC-mean)/sd
  NES <- sapply(AUC, function(x) (x-meanAUC)/sdAUC)
  return(NES)
}

#' @import data.table
.auc.asTable <- function(auc, nesThreshold=3.0, digits=3)
{
  nes <- .calcNES(auc)
  nes <- sort(nes, decreasing=TRUE)

  signifRankings <- names(nes)[which(nes >= nesThreshold)]
  aucTable <- data.table(motif=signifRankings,
                         NES=signif(nes[signifRankings], digits=digits),
                         AUC=signif(auc[signifRankings],digits=digits))
  aucTable
}

#' @import data.table
.addTfs <- function(aucTable, motifAnnot_direct=NULL, motifAnnot_indirect=NULL, highlightTFs=NULL)
{
  if(!is.null(highlightTFs))
  {
    aucTable <- data.table(aucTable, highlightedTFs=paste(highlightTFs, collapse=", ") , TFinDB="")
    tmp <- .tfInAnnot(aucTable$motif, inputTFs=highlightTFs, motifAnnot_direct=motifAnnot_direct, motifAnnot_indirect=motifAnnot_indirect)
    if(!is.null(motifAnnot_indirect)){
      wMotifs <- names(tmp$indirectAnnot)[which(tmp$indirectAnnot)]
      aucTable[which(aucTable$motif %in% wMotifs),"TFinDB"] <- "*"
    }
    if(!is.null(motifAnnot_direct)) {
      wMotifs <- names(tmp$directAnnot)[which(tmp$directAnnot)]
      aucTable[which(aucTable$motif %in% wMotifs),"TFinDB"] <- "**" # (Overrides indirect)
    }
  }

  if(!is.null(motifAnnot_direct))
  {
    TF_direct <- sapply(aucTable$motif, function(x) {
      paste(motifAnnot_direct[[x]][,1], collapse="; ")
    })
    aucTable <- data.table(aucTable, TF_direct=TF_direct)
  }

  if(!is.null(motifAnnot_indirect))
  {
    TF_indirect <- sapply(aucTable$motif, function(x) {
      paste(motifAnnot_indirect[[x]][,1], collapse="; ")
    })
    aucTable <- data.table(aucTable, TF_indirect=TF_indirect)
  }
  aucTable
}


# Not exclusive!
# TO DO: Optimize??
#' @import data.table
.tfInAnnot <- function(motifList, inputTFs, motifAnnot_direct=NULL, motifAnnot_indirect=NULL)
{
  in000 <- NULL
  in001 <- NULL
  # if(is.null(motifAnnot_direct) && is.null(motifAnnot_direct)) stop("Please provide the annotation.")
  if(!is.null(motifAnnot_direct))
  {
    motifs_00 <- motifList[which(motifList %in% names(motifAnnot_direct))]
    in000 <- sapply(motifAnnot_direct[motifs_00], function(x) any(inputTFs %in% x[,1]))
    if(length(in000)==0) in000 <- setNames(rep(FALSE, length(motifList)), motifList)
  }
  if(!is.null(motifAnnot_indirect))
  {
    motifs_001 <- motifList[which(motifList %in% names(motifAnnot_indirect))]
    in001 <- sapply(motifAnnot_indirect[motifs_001], function(x) any(inputTFs %in% x[,1]))
    if(length(in001)==0) in001 <- setNames(rep(FALSE, length(motifList)), motifList)
  }
  list(directAnnot=in000, indirectAnnot=in001)
}


