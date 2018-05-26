#' @title Imports the annotations of motifs to transcription factors 
#' @aliases importAnnotations
#' @param annotFile File containing the motif annotations corresponding to the
#' rankings. They should match in organism and version. 
#' @param motifsInRanking Subset of motifs to keep.
#' e.g. 
#' \code{motifsInRanking=getRanking(motifRankings)$features}
#' @param columnsToKeep Other colums from the file to keep
#' 
#' @description RcisTarget package includes the motif annotations for the rankings 
#' using motif collection version 9 ('mc9nr', 24453 motifs): 
#' 
#' \itemize{
#' \item{Human: }{ \code{data(motifAnnotations_hgnc)}}
#' \item{Mouse: }{ \code{data(motifAnnotations_mgi)}}
#' \item{Fly: }{ \code{data(motifAnnotations_dmel)}}
#' \item{Previous versions (Annotations for motif collection version 8:}{
#' motifAnnotations_hgnc_v8 (Human), motifAnnotations_mgi_v8 (Mouse), motifAnnotations_dmel_v8 (Fly)
#' }
#' }
#' 
#' This function (importAnnotations) allows to import annotations 
#' for other versions of the rankings, or to keep extra data columns.
#' 
#' e.g. Source of the annotations (motif collection 9 'mc9nr'):
#' \itemize{
#' \item{ Human: }{ https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl}
#' \item{ Mouse: }{ https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.mgi-m0.001-o0.0.tbl}
#' }
#' @seealso See \emph{iRegulon} paper and documentation for \bold{details} on how the
#' rankings and annotations were built.
#' @return
#' data.table with the annotations of the motifs to transcription factors
#' 
#' Columns:
#'  \itemize{
#'    \item{\bold{motif}: }{Motif ID.}
#'    \item{\bold{TF}: }{Transcription factor (or inferred gene).}
#'    \item{\bold{directAnnotation}, \bold{inferred_Orthology}, 
#'      \bold{inferred_MotifSimil}: }{
#'        Boolean values indicating whether the motif is 
#'        annotated to the TF in the source database ("directAnnotation"), 
#'        or whether it was inferred by orthology ("inferred_Orthology") 
#'        or motif similarity ("inferred_MotifSimil").}
#'    \item{\bold{Description}: }{Description of the source of the annotation.}
#'    \item{\bold{annotationSource}: }{Source of the annotation 
#'      formatted as factor (e.g. for subsetting). 
#'      Levels: directAnnotation, inferredBy_Orthology, 
#'      inferredBy_MotifSimilarity, inferredBy_MotifSimilarity_n_Orthology.}
#'  }
#' 
#' @examples
#' # motifAnnotations <- importAnnotations("motifs-v9-nr.hgnc-m0.001-o0.0.tbl")
#' 
#' ## To save (decrease from ~100MB to 1MB): 
#' # attr(motifAnnotations, "version") <- "motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
#' # save(motifAnnotations, file="hgnc_motifAnnotations.RData", compress='xz')
#' 
#' # This code would generate the equivalent to
#' data(motifAnnotations_hgnc)
##### Load/import the ranking from a feather file:
#' @rdname importAnnotations
#' @rawNamespace import(data.table, except = shift)
#' @export
importAnnotations <- function(annotFile, 
                              motifsInRanking=NULL, 
                              columnsToKeep=NULL)
{
  mAnnot <- data.table::fread(annotFile, 
                          colClasses = setNames("character","source_version"))
  colnames(mAnnot)[which(colnames(mAnnot)=="#motif_id")] <- "motif"
  colnames(mAnnot)[which(colnames(mAnnot)=="gene_name")] <- "TF"
  
  if(!is.null(motifsInRanking))
  {
    nMissingMotifs <- sum(unique(motifsInRanking) %in% mAnnot$motif)
    if(nMissingMotifs>0) 
      warning(nMissingMotifs,
              " of the requested motifs are not in the annotation file.")
    mAnnot <- mAnnot[mAnnot$motif %in% motifsInRanking]
  }
  
  # Re-format annotations
  mAnnot$directAnnotation <- unlist(
    mAnnot[,"description"]=="gene is directly annotated")
  mAnnot$inferred_Orthology <- mAnnot[,"orthologous_gene_name"]!="None"
  mAnnot$inferred_MotifSimil <- mAnnot[,"similar_motif_id"]!="None"
  
  annotationSource <- rep("", nrow(mAnnot))
  annotationSource[mAnnot$directAnnotation] <- "directAnnotation"
  annotationSource[mAnnot$inferred_Orthology] <- "inferredBy_Orthology"
  annotationSource[mAnnot$inferred_MotifSimil] <- "inferredBy_MotifSimilarity"
  annotationSource[mAnnot$inferred_Orthology 
      & mAnnot$inferred_MotifSimil] <- "inferredBy_MotifSimilarity_n_Orthology"
  mAnnot$annotationSource=factor(annotationSource); rm(annotationSource)
  #table(mAnnot$directAnnotation, mAnnot$annotationSource)
  
  selectedColumns <- c("motif", "TF", 
             "directAnnotation", "inferred_Orthology", "inferred_MotifSimil",
             "annotationSource", "description", columnsToKeep)
  mAnnot <- mAnnot[,selectedColumns,with=FALSE]
  data.table::setkeyv(mAnnot, c("motif", "TF"))
  
  return(mAnnot)
}

