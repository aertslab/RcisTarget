#' Documentation for the data

#' @rdname motifAnnotations
#' @name motifAnnotations_hgnc
#' @title Motifs to TF annotations (human)
#' @description Contains the annotations to HUMAN transcription factors for 
#' the rankings using motif collection version 9 ('mc9nr', 24453 motifs).
#' [Source: [motifs-v9-nr.hgnc-m0.001-o0.0](http://pyscenic.aertslab.org/resources/motifs-v9-nr.hgnc-m0.001-o0.0.tbl)]
#' 
#' This object is meant to be provided to RcisTarget without modification, 
#' but it can also be explored by the user 
#' to obtain further information about the motifs.
#' Columns:
#' \itemize{
#'  \item{\bold{motif}: }{Motif ID.}
#'  \item{\bold{TF}: }{Transcription factor (or inferred gene).}
#'  \item{\bold{directAnnotation}, \bold{inferred_Orthology}, 
#'        \bold{inferred_MotifSimil}: }
#'     {Boolean values indicating whether the motif is 
#'        annotated to the TF in the source database ("directAnnotation"), 
#'        or whether it was inferred by orthology ("inferred_Orthology") 
#'        or motif similarity ("inferred_MotifSimil").}
#'  \item{\bold{Description}: }{Description of the source of the annotation.}
#'    \item{\bold{annotationSource}: }{
#'      Source of the annotation formatted as factor (e.g. for subsetting). 
#'      Levels: directAnnotation, inferredBy_Orthology, 
#'      inferredBy_MotifSimilarity, inferredBy_MotifSimilarity_n_Orthology.}
#' }
#' @docType data
#' @seealso importAnnotations, RcisTarget
#' @keywords datasets
NULL

#' @rdname motifAnnotations
#' @name motifAnnotations_mgi
#' @title Motifs to TF annotations (mouse)
#' @description Contains the annotations to MOUSE transcription factors for 
#' the rankings using motif collection version 9 ('mc9nr', 24453 motifs).
#' [Source: [motifs-v9-nr.mgi-m0.001-o0.0](http://pyscenic.aertslab.org/resources/motifs-v9-nr.mgi-m0.001-o0.0.tbl)]
#' 
#' This object is meant to be provided to RcisTarget without modification, 
#' but it can also be explored by the user 
#' to obtain further information about the motifs.
#' Columns:
#' \itemize{
#'  \item{\bold{motif}: }{Motif ID.}
#'  \item{\bold{TF}: }{Transcription factor (or inferred gene).}
#'  \item{\bold{directAnnotation}, \bold{inferred_Orthology}, 
#'        \bold{inferred_MotifSimil}: }
#'     {Boolean values indicating whether the motif is 
#'        annotated to the TF in the source database ("directAnnotation"), 
#'        or whether it was inferred by orthology ("inferred_Orthology") 
#'        or motif similarity ("inferred_MotifSimil").}
#'  \item{\bold{Description}: }{Description of the source of the annotation.}
#'    \item{\bold{annotationSource}: }{
#'      Source of the annotation formatted as factor (e.g. for subsetting). 
#'      Levels: directAnnotation, inferredBy_Orthology, 
#'      inferredBy_MotifSimilarity, inferredBy_MotifSimilarity_n_Orthology.}
#' }
#' @docType data
#' @seealso importAnnotations, RcisTarget
#' @keywords datasets
NULL


