# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)
#' @import data.table
#' @title cisTarget
#' @rdname RcisTarget
#' @description Identifies DNA motifs significantly over-represented in a gene-set.
#'
#' This is the main function to run RcisTarget. It includes on the following steps:
#' \itemize{
#' \item 1. Motif enrichment analysis (\link{calcAUC})
#' \item 2. Motif-TF annotation (\link{addMotifAnnotation})
#' \item 3. Selection of significant genes (\link{addSignificantGenes})
#' }
#'
#' @param geneSets List of gene-sets to analyze. The gene-sets should be provided as a 'named list' in which each element is a gene-set (i.e. \code{list(geneSet1=c("gene1", "gene2"))})
#' @param motifRankings Database of the appropiate organism and search-space (i.e. 10kbp around- or 500bp upstream the TSS).
#' These objects are provided in separate packages:
#' \itemize{
#' \item \url{http://bioconductor.org/packages/RcisTarget.hg19.motifDatabases} (Human)
#' \item \url{http://bioconductor.org/packages/RcisTarget.mm9.motifDatabases} (Mouse)
#' }
#' See the help files for more information: i.e. \code{help(RcisTarget.hg19.motifDatabases)}
#' @param motifAnnot_direct Motif annotation database containing DIRECT annotations of the motif to transcription factors.
#' @param motifAnnot_indirect Motif annotation database containing the expanded annotations of the motif to transcription factors based on motif similarity.
#' @param highlightTFs Character. If a list of transcription factors is provided, the column TFinDB in the otuput table will indicate whether any of those TFs are included within the direct annotation (two asterisks, **) or indirect annotation (one asterisk, *) of the motif. The vector can be named to indicate which TF to highlight for each gene-set. Otherwise, all TFs will be used for all geneSets.
#' @param nesThreshold Numeric. NES threshold to calculate the motif significant (3.0 by default). The NES is calculated -for each motif- based on the AUC distribution of all the motifs for the gene-set [(x-mean)/sd]. The motifs are considered significantly enriched if they pass the the Normalized Enrichment Score (NES) threshold.
#' @param aucThresholdPERC Threshold to calculate the AUC.
#' In a simplified way, the AUC value represents the fraction of genes -within the top X genes in the ranking- that are included in the signature.
#' The parameter 'aucThresholdPERC' allows to modify the percentage of genes (of the top of the ranking) that is used to perform this computation.
#' By default it is set to 5\% of the total number of genes in the rankings. Common values range from 1 to 10\%.
#' @param geneErnMethod "iCisTarget" or "aprox". Method to identify the highly ranked genes (see \link{addSignificantGenes} for details).
#' @param geneErnMmaxRank Maximum rank to take into account for the gene enrichment recovery curve (see \link{addSignificantGenes} for details).
#' @param nCores Number of cores to use for computation.
#' Note: In general, using a higher number of cores (e.g. processes) decreases overall running time. However, it also deppends on the available memory and overall system load. Setting nCores too high might also decrease performance.
#' @param verbose Should the function show progress messages? (TRUE / FALSE)
#' @return \code{\link[data.table]{data.table}} containing the over-represented motifs (according to the selected NES threshold),
#' their statistics, annotation to transcription factors and the genes with high enrichment of the motif.
#' @seealso See the package vignette for examples and more details: \code{vignette("RcisTarget")}
#' @example inst/examples/example_cisTarget.R
#' @export
cisTarget <- function(geneSets, motifRankings,
            motifAnnot_direct=NULL, motifAnnot_indirect=NULL, highlightTFs=NULL,
            nesThreshold=3.0, aucThresholdPERC=0.05,
            geneErnMethod="aprox", geneErnMmaxRank=5000,
            nCores=1, verbose=TRUE)
{
  # suppressPackageStartupMessages(library(data.table))
  # Calculate AUC
  motifs_AUC <- calcAUC(geneSets=geneSets, rankings=motifRankings, nCores=nCores, aucMaxRank=aucThresholdPERC*nrow(motifRankings), verbose=verbose)

  # Select significant motifs, add TF annotation & format as table
  motifEnrichmentTable <- addMotifAnnotation(auc=motifs_AUC, nesThreshold=nesThreshold, motifAnnot_direct=motifAnnot_direct, motifAnnot_indirect=motifAnnot_indirect, highlightTFs=highlightTFs)

  # Identify significant genes for each motif (i.e. genes from the gene set in the top of the ranking)
  motifEnrichmentTable_wGenes <- addSignificantGenes(resultsTable=motifEnrichmentTable,
                                                     geneSets=geneSets,
                                                     rankings=motifRankings,
                                                     maxRank=geneErnMmaxRank,
                                                     plotCurve=FALSE, genesFormat="geneList",
                                                     method=geneErnMethod, nCores=nCores)
}
