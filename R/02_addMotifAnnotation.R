

# Intput TFs: If not named, all TFs will be used with all geneSets


# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)
#' @import data.table
#'
#' @title to do
#' @description to do
#' @param to do
#' @return to do Columna TF
# @example to do
#' @export
addMotifAnnotation <- AUC.asDataTable <- function(AUCellOutput, nesThreshold=3.0, digits=3, motifAnnot_direct=NULL, motifAnnot_indirect=NULL, highlightTFs=NULL)
{
  #### Check inputs
  if(!is.null(highlightTFs))
  {
    if(is.null(motifAnnot_direct) && is.null(motifAnnot_indirect)) stop("To hightlight TFs, please provide a motif-TF annotation.")
    if(is.null(names(highlightTFs))) {
      warning("The input TFs are not named, all TFs will be used with all Gene Sets.")
      highlightTFs <- setNames(rep(list(highlightTFs), ncol(AUCellOutput)), colnames(AUCellOutput))
    }

    if(!all(names(highlightTFs) %in% colnames(AUCellOutput))) warning("TFs 1")
    if(!all(colnames(AUCellOutput) %in% names(highlightTFs))) warning("TFs 2")
  }

  #### Runs "auc.asTable" on each AUC columns i.e. signatures/cells
  ret <- lapply(colnames(AUCellOutput), function(geneSet) {
      tfs <- highlightTFs[[geneSet]]
      aucTable <- .auc.asTable(AUCellOutput[,geneSet], nesThreshold=nesThreshold, digits=digits)
      if(nrow(aucTable)>0)
      {
        aucTable <- .addTfs(aucTable, motifAnnot_direct=motifAnnot_direct, motifAnnot_indirect=motifAnnot_indirect, highlightTFs=tfs)
        aucTable <- data.table(geneSet=geneSet, aucTable)
      }else{
        aucTable <- NULL
      }
      aucTable
    })

  #### Merge the results from each signature/cell into a single data.table
  # ret <- do.call(rbind, unname(ret))  # Slower?
  # library(data.table)
  ret <- rbindlist(ret)
  return(ret)
}

#### Output:
#             geneSet                            motif  NES    AUC TF_direct
# 1:  Astrocyte_Cahoy                   stark__AAASTTT 5.70 0.0608
# 2:  Astrocyte_Cahoy                     cisbp__M0574 5.58 0.0603
# 3:  Astrocyte_Cahoy                     cisbp__M0648 4.29 0.0550
# 4:  Astrocyte_Cahoy                     cisbp__M0659 4.25 0.0549
# 5:  Astrocyte_Cahoy                     cisbp__M0664 4.19 0.0546
# ---
#   1075: RandGeneSet_100g                     cisbp__M6289 3.03 0.0823     Hoxa9
# 1076: RandGeneSet_100g homer__ATAGTGCCACCTGGTGGCCA_CTCF 3.02 0.0821      Ctcf
# 1077: RandGeneSet_100g                     cisbp__M0828 3.02 0.0821
# 1078: RandGeneSet_100g                elemento__TGGGGCC 3.02 0.0821
# 1079: RandGeneSet_100g  taipale__HNF4A_NRGTCCAAAGTCCANY 3.01 0.0819     Hnf4a


############ PRIVATE
.calcNES <- function(AUC)
{
  meanAUC <- mean(AUC)
  sdAUC <- sd(AUC)

  # NES = (AUC-mean)/sd
  NES <- sapply(AUC, function(x) (x-meanAUC)/sdAUC)
  return(NES)
}


.auc.asTable <- function(auc, nesThreshold=3.0, digits=3)
{
  nes <- .calcNES(auc)
  nes <- sort(nes, decreasing=TRUE)

  signifRankings <- names(nes)[which(nes >= nesThreshold)]
  aucTable <- data.table(ranking=signifRankings,
                         NES=signif(nes[signifRankings], digits=digits),
                         AUC=signif(auc[signifRankings],digits=digits))
  aucTable
}

.addTfs <- function(aucTable, motifAnnot_direct=NULL, motifAnnot_indirect=NULL, highlightTFs=NULL)
{
  if((!is.null(motifAnnot_direct)) || (!is.null(motifAnnot_indirect)) || (!is.null(highlightTFs))) colnames(aucTable)[which(colnames(aucTable) == "ranking")] <- "motif"

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


