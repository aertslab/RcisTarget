# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)
# @import data.table
#' @title cisTarget
#' @description Identifies DNA motifs significantly over-represented in a gene-set.
#'
#' This is the main function to run RcisTarget (motif enrichment analysis):
#' It is a wrapper that calls \link{calcAUC}, \link{addMotifAnnotation}
#' and \link{addSignificantGenes}.
#' @param geneSets List of gene-sets to analyze. The gene-sets should be provided as a 'named list' in which each element is a gene-set (i.e. \code{list(geneSet1=c("gene1", "gene2"))})
#' @param motifRankings Database of the appropiate organism and search-space (i.e. 10kbp around- or 500bp upstream the TSS).
#' These objects are provided in separate packages: \link[http://bioconductor.org/packages/RcisTarget.hg19.motifDatabases]{RcisTarget.hg19.motifDatabases} (Human), \link[http://bioconductor.org/packages/RcisTarget.mm9.motifDatabases]{RcisTarget.mm9.motifDatabases} (Mouse).
#' See the help files for more information: i.e. \code{help(RcisTarget.hg19.motifDatabases)}.
#' @param motifAnnot_direct TO DO
#' @param motifAnnot_indirect TO DO
#' @param highlightTFs TO DO
#' @param nesThreshold TO DO
#' @param aucThresholdPERC Threshold to calculate the AUC.
#' In a simplified way, the AUC value represents the fraction of genes -within the top X genes in the ranking- that are included in the signature.
#' The parameter 'aucThresholdPERC' allows to modify the percentage of genes (of the top of the ranking) that is used to perform this computation.
#' By default it is set to 5\% of the total number of genes in the rankings. Common values range from 1 to 10\%.
#' @param geneErnMethod TO DO
#' @param geneErnMmaxRank TO DO
#' @param nCores Number of cores to use for computation.
#' Note: In general, using a higher number of cores (e.g. processes) decreases overall running time. However, it also deppends on the available memory and overall system load. Setting nCores too high might also decrease performance.
#' @param verbose Should the function show progress messages? (TRUE / FALSE)
#' See \code{vignette("RcisTarget")} for examples and more details.
#' @return Data table containing the over-represented motifs (according to the selected NES threshold),
#' their statistics, annotation to transcription factors and the genes with high enrichment of the motif.
#' @example inst/examples/cisTarget_example.R
#' @export
cisTarget <- function(geneSets, motifRankings,
            motifAnnot_direct=NULL, motifAnnot_indirect=NULL, highlightTFs=NULL,
            nesThreshold=3.0, aucThresholdPERC=0.05,
            geneErnMethod="aprox", geneErnMmaxRank=5000,
            nCores=1, verbose=TRUE)
{
  suppressPackageStartupMessages(library(data.table))
  # Calculate AUC
  motifs_AUC <- calcAUC(geneSets=geneSets, rankings=motifRankings, nCores=nCores, aucMaxRank=aucThresholdPERC*nrow(motifRankings), verbose=verbose)

  # Select significant motifs, add TF annotation & format as table
  motifEnrichmentTable <- addMotifAnnotation(AUCellOutput=motifs_AUC, nesThreshold=nesThreshold, motifAnnot_direct=motifAnnot_direct, motifAnnot_indirect=motifAnnot_indirect, highlightTFs=highlightTFs)

  # Identify significant genes for each motif (i.e. genes from the gene set in the top of the ranking)
  motifEnrichmentTable_wGenes <- addSignificantGenes(aucTable=motifEnrichmentTable,
                                                     geneSets=geneSets,
                                                     rankings=motifRankings,
                                                     maxRank=geneErnMmaxRank,
                                                     plotCurve=FALSE, genesFormat="geneList",
                                                     method=geneErnMethod, nCores=nCores)
}
