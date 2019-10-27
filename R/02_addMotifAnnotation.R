# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)

#' @importFrom stats sd setNames
#'
#' @title Add motif annotation
#' @description Select significant motifs and/or annotate motifs to
#' genes or transcription factors.
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
#' @param idColumn Annotation column containing the ID (e.g. motif, accession)
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
       idColumn="motif",
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
    # if(!is.null(motifAnnot_highConfCat) && 
    #    any(!motifAnnot_highConfCat %in% levels(motifAnnot$annotationSource)))
    #   warning("'motifAnnot_highConfCat' 
    #        should be a value in the column 'annotationSource'.") 
    # 
    # if(!is.null(motifAnnot_lowConfCat) && 
    #    any(!motifAnnot_lowConfCat %in% levels(motifAnnot$annotationSource)))
    #   warning("'motifAnnot_lowConfCat' 
    #        should be a value in the column 'annotationSource'.") 
    
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
                             nesThreshold=nesThreshold, digits=digits, idColumn=idColumn)
    if(nrow(aucTable)>0)
    {
      aucTable <- .addTfs(aucTable,
                          motifAnnot=motifAnnot,
                          TFs=tfs, 
                          motifAnnot_highConfCat=motifAnnot_highConfCat,
                          motifAnnot_lowConfCat=motifAnnot_lowConfCat,
                          idColumn=idColumn)
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

.auc.asTable <- function(auc, nesThreshold=3.0, digits=3, idColumn="motif")
{
  nes <- .calcNES(auc)
  nes <- sort(nes, decreasing=TRUE)
  
  signifRankings <- names(nes)[which(nes >= nesThreshold)]
  aucTable <- data.table::data.table(motif=signifRankings,
                                     NES=signif(nes[signifRankings], digits=digits),
                                     AUC=signif(auc[signifRankings],digits=digits))
  colnames(aucTable)[1] <- idColumn
  aucTable
}

.addTfs <- function(aucTable,
                    motifAnnot=NULL,
                    TFs=NULL,
                    motifAnnot_highConfCat=NULL,
                    motifAnnot_lowConfCat=NULL,
                    idColumn="motif")
{
  if(!is.null(TFs))
  {
    aucTable <- data.table::data.table(aucTable,
                                 highlightedTFs=paste(TFs, collapse=", ") ,
                                 TFinDB="")
    
    if(!is.null(motifAnnot)) {
      motifAnnot_subset <- motifAnnot[(motifAnnot[[idColumn]] %in% aucTable[[idColumn]]) 
                                      & (motifAnnot$TF %in% TFs), ][,c(idColumn, "TF", "annotationSource"),with=F]
      motifAnnot_subset <- split(motifAnnot_subset, motifAnnot_subset[[idColumn]])
      for(motifName in names(motifAnnot_subset))
      {
        if(any(as.character(motifAnnot_subset[[motifName]]$annotationSource) 
               %in% motifAnnot_lowConfCat))
          aucTable[aucTable[[idColumn]]==motifName,"TFinDB"] <- "*"
        
        # overrides lowConf
        if(any(as.character(motifAnnot_subset[[motifName]]$annotationSource) 
               %in% motifAnnot_highConfCat))
          aucTable[aucTable[[idColumn]]==motifName,"TFinDB"] <- "**"
      }
    }
  }
  
  if(!is.null(motifAnnot))
  {
    if(!is.null(motifAnnot_highConfCat))
    {
      TF_highConf <- .formatTfs(motifs=aucTable[[idColumn]], 
                                motifAnnot=motifAnnot,
                                annotCats=motifAnnot_highConfCat,
                                idColumn=idColumn)
      
      aucTable <- data.table::data.table(aucTable, TF_highConf=TF_highConf)
    }
    
    if(!is.null(motifAnnot_lowConfCat))
    {
      TF_lowConf <- .formatTfs(motifs=aucTable[[idColumn]], 
                               motifAnnot=motifAnnot,
                               annotCats=motifAnnot_lowConfCat,
                               idColumn=idColumn)
      
      aucTable <- data.table::data.table(aucTable, TF_lowConf=TF_lowConf)
    }
  }
  
  aucTable
}


## 26 apr 2019
# Replaced input: aucTable by motifs. In calls:  .formatTfs(motifs=aucTable[[idColumn]]
# aucTable$motif --> motifs
# nrow(aucTable) --> length(motifs)
.formatTfs <- function(motifs, motifAnnot, annotCats, idColumn)
{
  motifAnnot_subset <- motifAnnot[motifAnnot$annotationSource %in% annotCats, ] 
  motifAnnot_subset <- motifAnnot_subset[motifAnnot_subset[[idColumn]] %in% motifs, ] 
  motifAnnot_Cats <- vapply(split(motifAnnot_subset, motifAnnot_subset[[idColumn]]), 
              function(mat){
                mat <- split(mat$TF, factor(mat$annotationSource))
                tfsByCat <- vapply(names(mat),
                                   function(x) paste(paste(unlist(mat[[x]]),
                                                           collapse="; "),
                                                     " (",x,"). ",
                                                     sep=""), "")
                paste(tfsByCat, collapse="")
              }, FUN.VALUE="")
  
  ret <- setNames(rep("", length(motifs)), motifs)
  ret[names(motifAnnot_Cats)] <- motifAnnot_Cats
  return(ret)
}

#' @title Get motif annotation
#' @description Get the genes/transcription factors annotated to the given motifs

#' @param motifs Motif IDs 
#' @param motifAnnot Motif annotation database containing the 
#' annotations of the motif to genes or transcription factors.
#' @param annotCats Annotation categories to be considered:  
#' "directAnnotation" (annotated in the source database),  
#' "inferredBy_Orthology" (the motif is annotated to an homologous/ortologous 
#' gene), or inferred based on motif similarity ("inferredBy_MotifSimilarity", 
#' "inferredBy_MotifSimilarity_n_Orthology").
#' @param returnFormat Determines the output format. Choose one of the following values:
#' \itemize{
#' \item \code{asCharacter}: Named vector with the genes or TFs annotated to the given motifs (in the same order, including empty and duplicated values).
#' \item \code{subset}: Subset of the annotation table (list split by motif)
#' \item \code{list}: List of TF names (unique values), duplicated motifs or motifs without annotation are not returned.
#' }
#' @return See argument \code{returnFormat}
#' @seealso \code{\link{addMotifAnnotation}} add the annotation directly to the motif enrichment results.
#' 
#' See the package vignette for examples and more details:
#' \code{vignette("RcisTarget")}
#' @example inst/examples/example_addMotifAnnotation.R
#' @export
getMotifAnnotation <- function(motifs, 
                               motifAnnot, 
                               annotCats=c("directAnnotation",
                                           "inferredBy_MotifSimilarity",
                                           "inferredBy_Orthology",
                                           "inferredBy_MotifSimilarity_n_Orthology"),
                               idColumn="motif",
                               returnFormat=c("asCharacter","subset","list")[1])
{
  ## Check inputs:
  returnFormat <- tolower(returnFormat)
  if(!returnFormat %in% c("ascharacter","subset","list")) stop("returnFormat should be eihter 'asCharacter', 'subset', or 'list'.")
  if(length(returnFormat)>1) stop("Please, choose ONE returnFormat.")
  
  ## Run:
  if(returnFormat=="ascharacter"){
    ret <- .formatTfs(motifs=motifs, 
                      motifAnnot=motifAnnot,
                      annotCats=annotCats,
                      idColumn=idColumn) 
  }else{
    ret <- .getTfs(motifs=motifs, 
                   motifAnnot=motifAnnot,
                   annotCats=annotCats,
                   idColumn=idColumn,
                   returnFormat=returnFormat) 
  }
  return(ret)
}

.getTfs <- function(motifs, motifAnnot, annotCats, idColumn, returnFormat)
{
  motifAnnot_subset <- motifAnnot[motifAnnot$annotationSource %in% annotCats, ] 
  motifAnnot_subset <- motifAnnot_subset[motifAnnot_subset[[idColumn]] %in% motifs, ] 
  # motifAnnot_subset <- split(motifAnnot_subset[,c("TF", "directAnnotation", "inferred_Orthology", "inferred_MotifSimil","annotationSource")], motifAnnot_subset[[idColumn]])
  motifAnnot_subset <- split(motifAnnot_subset, motifAnnot_subset[[idColumn]])
  
  ret <- motifAnnot_subset # returnFormat=="subset"
  
  if(returnFormat=="list")
    ret <- lapply(ret, function(x) sort(unique(x[["TF"]])))
  
  return(ret)
}

