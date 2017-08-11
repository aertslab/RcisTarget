.check_addMotifAnnotation <- function(motifs_AUC, motifAnnot_direct, motifAnnot_inferred)
{
  # No annot
  motifEnrichmentTable_noAnnot <- addMotifAnnotation(motifs_AUC)
  .check_motifEnrichmentTable(motifEnrichmentTable_noAnnot, annotDirect=NULL, annotInferred=NULL)

  # only direct
  motifEnrichmentTable <- addMotifAnnotation(motifs_AUC,
                                             motifAnnot_direct=motifAnnot_direct)
  .check_motifEnrichmentTable(motifEnrichmentTable, annotDirect=motifAnnot_direct, annotInferred=NULL)

  # only indirect
  motifEnrichmentTable <- addMotifAnnotation(motifs_AUC,
                                             motifAnnot_inferred=motifAnnot_inferred)
  .check_motifEnrichmentTable(motifEnrichmentTable, annotDirect=NULL, annotInferred=motifAnnot_inferred)

  # both annots, plus other options
  motifEnrichment_wIndirect <- suppressWarnings(addMotifAnnotation(motifs_AUC, nesThreshold=2,
                                                  motifAnnot_direct=motifAnnot_direct,
                                                  motifAnnot_inferred=motifAnnot_inferred,
                                                  highlightTFs = "HIF1A", digits=3))
  .check_motifEnrichmentTable(motifEnrichment_wIndirect, annotDirect=motifAnnot_direct, annotInferred=motifAnnot_inferred, highlightTFs="HIF1A")
  expect_true(min(motifEnrichment_wIndirect$NES)>=2)
}


.check_motifEnrichmentTable <- function(met, annotDirect, annotInferred, highlightTFs=NULL)
{
  nCols <- 4
  if(!is.null(annotDirect)) {
    nCols <- nCols+1
    expect_true("TF_direct" %in% colnames(met))
  }

  if(!is.null(annotInferred)) {
    nCols <- nCols+1
    expect_true("TF_inferred" %in% colnames(met))
  }
  if(!is.null(highlightTFs)) {
    nCols <- nCols+2
    expect_true("highlightedTFs" %in% colnames(met))
    expect_true(all(met$highlightedTFs=="HIF1A"))
  }


  expect_equal(class(met)[1], "data.table")
  expect_equal(ncol(met), nCols)
  expect_equal(colnames(met)[1:4], c("geneSet", "motif", "NES", "AUC"))
}






