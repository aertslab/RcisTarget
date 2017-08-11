

test_Workflow <- function()
{
  library(RcisTarget)

  # motifRankings <- fakeDatabase()
  # geneLists <- list(fakeGenes=sample(getRanking(motifRankings)$rn, 20))
  # library(RcisTarget.hg19.motifDatabases.20k)
  # data(hg19_direct_motifAnnotation)
  # data(hg19_inferred_motifAnnotation)

  ##### Example setup ##############################
  txtFile <- paste(file.path(system.file('examples', package='RcisTarget')),
                   "hypoxiaGeneSet.txt", sep="/")
  geneLists <- list(hypoxia=read.table(txtFile)[,1])

  # Load databases
  library(RcisTarget.hg19.motifDatabases.20k)
  data(hg19_10kbpAroundTss_motifRanking)
  motifRankings <- hg19_10kbpAroundTss_motifRanking

  # Fake example with 5000 random motifs
  library(data.table)
  set.seed(123)
  motifRankings <- subset(motifRankings,
                   sample(2:ncol(motifRankings@rankings), 5000), select="col")
  ##################################################

  ##################################################
  # Standard workflow
  # (To use as input for the different check functions)
  # Step 1.

  .check_calcAUC(motifRankings) # Checks different intput/ouputs for calcAUC
  motifs_AUC <- calcAUC(geneLists, motifRankings) # For input of next steps

  # Step 2.
  data(hg19_direct_motifAnnotation)
  data(hg19_inferred_motifAnnotation)
  .check_addMotifAnnotation(motifs_AUC, motifAnnot_direct=hg19_direct_motifAnnotation, motifAnnot_inferred=hg19_inferred_motifAnnotation)
  motifEnrichmentTable <- addMotifAnnotation(motifs_AUC,
                                             motifAnnot_direct=hg19_direct_motifAnnotation)

  # Step 3
  .check_addSignificantGenes(met=motifEnrichmentTable, motifRankings=motifRankings, geneLists=geneLists)
  motifEnrichmentTable_wGenes <- addSignificantGenes(resultsTable=motifEnrichmentTable,
                                                     geneSets=geneLists,
                                                     rankings=motifRankings,
                                                     method="aprox")
  ##################################################


  ##################################################
  # cisTarget (automated workflow)
  motifEnrichmentTable_wGenes <- cisTarget(geneLists, motifRankings,
                                           motifAnnot_direct=hg19_direct_motifAnnotation,
                                           nesThreshold=3.5, geneErnMethod="aprox", nCores=4)

  expect_equal(colnames(motifEnrichmentTable_wGenes),
               c("geneSet", "motif", "NES", "AUC", "TF_direct", "nEnrGenes", "rankAtMax", "enrichedGenes"))
  ##################################################


  motifEnrichmentTable <- addLogo(motifEnrichmentTable)
}

test_that("RcisTarget workflow tests", test_Workflow)



