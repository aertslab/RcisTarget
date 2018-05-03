# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)

#' @importFrom stats sd setNames
#'
#' @title Add motif annotation
#' @description Select significant motifs and/or annotate motifs to
#' transcription factors.
#' The motifs are considered significantly enriched if they pass the the
#' Normalized Enrichment Score (NES) threshold.
#' @param auc Output from calcAUC.
#' @param nesThreshold Numeric. NES threshold to calculate the motif significant
#' (3.0 by default). The NES is calculated -for each motif- based on the AUC
#' distribution of all the motifs for the gene-set [(x-mean)/sd].
#' @param digits Integer. Number of digits for the AUC and NES in the
#' output table.
#' @param motifAnnot Motif annotation database containing the 
#' annotations of the motif to transcription factors.
#' The names should match the ranking column names.
#' @param motifAnnot_highConfCat Categories considered as source for 
#' 'high confidence' annotations. By default, 
#' "directAnnotation" (annotated in the source database), and 
#' "inferredBy_Orthology" (the motif is annotated to an homologous/ortologous 
#' gene).
#' @param motifAnnot_lowConfCat Categories considered 
#' 'lower confidence' source for annotations. By default, the annotations 
#' inferred based on motif similarity ("inferredBy_MotifSimilarity", 
#' "inferredBy_MotifSimilarity_n_Orthology").
#' @param highlightTFs Character. If a list of transcription factors is
#' provided, the column TFinDB in the otuput table will indicate whether any
#' of those TFs are included within the 'high-confidence' annotation 
#' (two asterisks, **)
#' or 'low-confidence' annotation (one asterisk, *) of the motif.
#' The vector can be named to indicate which TF to highlight for each gene-set.
#' Otherwise, all TFs will be used for all geneSets.
#' 
#' @return \code{\link[data.table]{data.table}} with the folowing columns:
#' \itemize{
#' \item geneSet: Name of the gene set
#' \item motif: ID of the motif
#' (colnames of the ranking, it might be other kind of feature)
#' \item NES: Normalized enrichment score of the motif in the gene-set
#' \item AUC: Area Under the Curve (used to calculate the NES)
#' \item TFinDB: Indicates whether the highlightedTFs are included within the
#' high-confidence annotation (two asterisks, **)
#' or lower-confidence annotation (one asterisk, *)
#' \item TF_highConf: Transcription factors annotated to the motif 
#' based on high-confidence annotations.
#' \item TF_lowConf: Transcription factors annotated to the motif according to
#' based on lower-confidence annotations.
#' }
#' @seealso Next step in the workflow: \code{\link{addSignificantGenes}}.
#'
#' Previous step in the workflow: \code{\link{calcAUC}}.
#'
#' See the package vignette for examples and more details:
#' \code{vignette("RcisTarget")}
#' @example inst/examples/example_addMotifAnnotation.R
#' @export
addMotifAnnotation <- function(auc, nesThreshold=3.0, digits=3,
       motifAnnot=NULL, 
       motifAnnot_highConfCat=c("directAnnotation", "inferredBy_Orthology"), 
       motifAnnot_lowConfCat=c("inferredBy_MotifSimilarity", 
                               "inferredBy_MotifSimilarity_n_Orthology"), 
       highlightTFs=NULL)
{
  auc <- getAUC(auc)
  #### Check inputs
  if(!is.null(highlightTFs))
  {
    if(is.null(motifAnnot))
      stop("To hightlight TFs, please provide a motif-TF annotation.")
    if(is.null(names(highlightTFs))) {
      warning("The input TFs are not named, ",
              "all TFs will be used with all Gene Sets.")
      highlightTFs <- setNames(rep(list(highlightTFs), nrow(auc)),
                               rownames(auc))
    }
    
    if(!all(names(highlightTFs) %in% rownames(auc))) warning("TFs 1")
    if(!all(rownames(auc) %in% names(highlightTFs))) warning("TFs 2")
  }
  
  if(!is.null(motifAnnot))
  {
    if(!is.data.table(motifAnnot))
      stop("motifAnnot should be a data.table")
    if(!is.null(motifAnnot_highConfCat) && 
       any(!motifAnnot_highConfCat %in% levels(motifAnnot$annotationSource)))
      stop("'motifAnnot_highConfCat' 
           should be a value in the column 'annotationSource'.") 
    
    if(!is.null(motifAnnot_lowConfCat) && 
       any(!motifAnnot_lowConfCat %in% levels(motifAnnot$annotationSource)))
      stop("'motifAnnot_lowConfCat' 
           should be a value in the column 'annotationSource'.") 
    
    commonCat <- intersect(motifAnnot_highConfCat, motifAnnot_lowConfCat)
    if(length(commonCat)>0)
      warning("The following annotation types are both in", 
              "'motifAnnot_highConfCat' and 'motifAnnot_lowConfCat': ", commonCat)
  }
  
  #### Runs "auc.asTable" on each signature/geneset
  applyFun <- lapply
  if((nrow(auc)>4) && ("BiocParallel" %in% installed.packages())) 
    applyFun <- BiocParallel::bplapply
  
  ret <- applyFun(rownames(auc), function(geneSet) {
    tfs <- highlightTFs[[geneSet]]
    aucTable <- .auc.asTable(auc[geneSet,],
                             nesThreshold=nesThreshold, digits=digits)
    if(nrow(aucTable)>0)
    {
      aucTable <- .addTfs(aucTable,
                          motifAnnot=motifAnnot,
                          TFs=tfs, 
                          motifAnnot_highConfCat=motifAnnot_highConfCat,
                          motifAnnot_lowConfCat=motifAnnot_lowConfCat)
      aucTable <- data.table::data.table(geneSet=geneSet, aucTable)
    }else{
      aucTable <- NULL
    }
    aucTable
  })
  
  ## Merge the results from each signature/geneSet/regionSet into a single dt
  # ret <- do.call(rbind, unname(ret))  # Slower?
  # library(data.table)
  ret <- data.table::rbindlist(ret)
  
  if(nrow(ret)>0)
    colnames(ret)[which(colnames(ret) == "ranking")] <- "motif"
  return(ret)
}


############ PRIVATE
.calcNES <- function(AUC)
{
  meanAUC <- mean(AUC)
  sdAUC <- sd(AUC)
  
  # NES = (AUC-mean)/sd
  NES <- vapply(AUC, function(x) (x-meanAUC)/sdAUC,
                FUN.VALUE=numeric(1))
  return(NES)
}

.auc.asTable <- function(auc, nesThreshold=3.0, digits=3)
{
  nes <- .calcNES(auc)
  nes <- sort(nes, decreasing=TRUE)
  
  signifRankings <- names(nes)[which(nes >= nesThreshold)]
  aucTable <- data.table::data.table(motif=signifRankings,
                                     NES=signif(nes[signifRankings], digits=digits),
                                     AUC=signif(auc[signifRankings],digits=digits))
  aucTable
}

.addTfs <- function(aucTable,
                    motifAnnot=NULL,
                    TFs=NULL,
                    motifAnnot_highConfCat=NULL,
                    motifAnnot_lowConfCat=NULL)
{
  if(!is.null(TFs))
  {
    aucTable <- data.table::data.table(aucTable,
                                 highlightedTFs=paste(TFs, collapse=", ") ,
                                 TFinDB="")
    
    if(!is.null(motifAnnot)) {
      motifAnnot_subset <- motifAnnot[(motifAnnot$motif %in% aucTable$motif) 
                                      & (motifAnnot$TF %in% TFs), 
                                      c("motif", "TF", "annotationSource")]
      motifAnnot_subset <- split(motifAnnot_subset, motifAnnot_subset$motif)
      for(motifName in names(motifAnnot_subset))
      {
        if(any(as.character(motifAnnot_subset[[motifName]]$annotationSource) 
               %in% motifAnnot_lowConfCat))
          aucTable[aucTable$motif==motifName,"TFinDB"] <- "*"
        
        # overrides lowConf
        if(any(as.character(motifAnnot_subset[[motifName]]$annotationSource) 
               %in% motifAnnot_highConfCat))
          aucTable[aucTable$motif==motifName,"TFinDB"] <- "**"
      }
    }
  }
  
  if(!is.null(motifAnnot))
  {
    if(!is.null(motifAnnot_highConfCat))
    {
      TF_highConf <- .formatTfs(aucTable=aucTable, 
                                motifAnnot=motifAnnot,
                                annotCats=motifAnnot_highConfCat)
      
      aucTable <- data.table::data.table(aucTable, TF_highConf=TF_highConf)
    }
    
    if(!is.null(motifAnnot_lowConfCat))
    {
      TF_lowConf <- .formatTfs(aucTable=aucTable, 
                               motifAnnot=motifAnnot,
                               annotCats=motifAnnot_lowConfCat)
      
      aucTable <- data.table::data.table(aucTable, TF_lowConf=TF_lowConf)
    }
  }
  
  aucTable
}


.formatTfs <- function(aucTable, motifAnnot, annotCats)
{
  motifAnnot_subset <- motifAnnot[motifAnnot$annotationSource %in% annotCats, ] 
  motifAnnot_subset <- motifAnnot_subset[motifAnnot_subset$motif %in% aucTable$motif, ] 
  motifAnnot_Cats <- vapply(split(motifAnnot_subset, motifAnnot_subset$motif), 
              function(mat){
                mat <- split(mat$TF, factor(mat$annotationSource))
                tfsByCat <- vapply(names(mat),
                                   function(x) paste(paste(unlist(mat[[x]]),
                                                           collapse="; "),
                                                     " (",x,").",
                                                     sep=""), "")
                paste(tfsByCat, collapse="")
              }, FUN.VALUE="")
  
  ret <- setNames(rep("", nrow(aucTable)), aucTable$motif)
  ret[names(motifAnnot_Cats)] <- motifAnnot_Cats
  return(ret)
}
## Previous version:
# .formatTfs <- function(aucTable, motifAnnot, annotCats)
# {
#   motifAnnot_Cats <- motifAnnot[motifAnnot$annotationSource %in% annotCats, ] 
#   motifAnnot_Cats <- motifAnnot_Cats[motifAnnot_Cats$motif %in% aucTable$motif, ] 
#   
#   vapply(aucTable$motif, function(mot) {
#     motifAnnot_selected <- motifAnnot_Cats[motifAnnot_Cats$motif==mot, ]
#     motifAnnot_selected <- split(motifAnnot_selected$TF,
#                                  motifAnnot_selected$annotationSource)
#     motifAnnot_selected <- motifAnnot_selected[
#       which(lengths(motifAnnot_selected)>0)]
# 
#     if(length(motifAnnot_selected) > 0){
#       tfsByCat <- vapply(names(motifAnnot_selected),
#                          function(x) paste(paste(unlist(motifAnnot_selected[[x]]),
#                                                  collapse="; "),
#                                            " (",x,"). ",
#                                            sep=""), "")
#       paste(tfsByCat, collapse="")
#     }else
#     {
#       ""
#     }
#   }, FUN.VALUE="")
# }
