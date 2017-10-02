.check_addSignificantGenes <- function(met, motifRankings, geneLists)
{
  motifEnrichmentTable_wGenes <- addSignificantGenes(met[1:3,],
                                                     plotCurve=TRUE,
                                                     geneSets=geneLists,
                                                     rankings=motifRankings,
                                                     method="aprox",
                                                     genesFormat="geneList")

  expect_equal(class(motifEnrichmentTable_wGenes)[1], "data.table")

  motifEnrichmentTable_wGenes <- as.data.frame(motifEnrichmentTable_wGenes)
  expect_equal(colnames(motifEnrichmentTable_wGenes),
               c("geneSet", "motif", "NES", "AUC", "TF_direct", "nEnrGenes", "rankAtMax", "enrichedGenes"))

  # As incidence matrix
  motifEnr_wIncidMat <- addSignificantGenes(met[1:3,],
                                            geneSets=geneLists,
                                            rankings=motifRankings,
                                            method="aprox",
                                            genesFormat = "incidMatrix")
  expect_true(ncol(motifEnr_wIncidMat)>7)

  ############# getSignificantGenes #########
  selectedMotif <- met$motif[1:2]
  onlyGenes <- getSignificantGenes(geneSet=geneLists$hypoxia,
                                   signifRankingNames=selectedMotif,
                                   genesFormat="incidMatrix",
                                   plotCurve=FALSE,
                                   rankings=motifRankings,
                                   method="aprox")


  ### Diferent input formats (geneSet)
  # a) Character vector (i.e. only one gene-set)
  fewGenes <- sample(getRanking(motifRankings)$rn, 10)
  gsOut <- getSignificantGenes(geneSet=fewGenes, signifRankingNames=selectedMotif[1],
             genesFormat="incidMatrix", plotCurve=FALSE, rankings=motifRankings, method="aprox")

  testthat::expect_equal(names(gsOut), c("enrStats", "incidMatrix"))

  # b) List
  otherGenes <- sample(getRanking(motifRankings)$rn, 5)
  geneSets <- list(geneSet1=fewGenes,
                   geneSet2=otherGenes)
  expect_error(getSignificantGenes(geneSet=geneSets, signifRankingNames=selectedMotif[1],
              genesFormat="incidMatrix", plotCurve=FALSE, rankings=motifRankings, method="aprox"))
  gsOut <- getSignificantGenes(geneSet=geneSets[1], signifRankingNames=selectedMotif[1],
              genesFormat="incidMatrix", plotCurve=FALSE, rankings=motifRankings, method="aprox")
  testthat::expect_equal(names(gsOut), c("enrStats", "incidMatrix"))

  # c) GeneSet object (from GSEABase)
  geneSetOne <- GSEABase::GeneSet(fewGenes, setName="geneSetOne")
  gsOut <- getSignificantGenes(geneSet=geneSetOne, signifRankingNames=selectedMotif[1],
                               genesFormat="incidMatrix", plotCurve=FALSE, rankings=motifRankings, method="aprox")

  testthat::expect_equal(names(gsOut), c("enrStats", "incidMatrix"))


  # d) GeneSetCollection object (from GSEABase)
  geneSetTwo <- GSEABase::GeneSet(otherGenes, setName="geneSetTwo")
  geneSets <- GSEABase::GeneSetCollection(geneSetOne, geneSetTwo)

  expect_error(getSignificantGenes(geneSet=geneSets, signifRankingNames=selectedMotif[1],
               genesFormat="incidMatrix", plotCurve=FALSE, rankings=motifRankings, method="aprox"))
  gsOut <- getSignificantGenes(geneSet=geneSets[1], signifRankingNames=selectedMotif[1],
              genesFormat="incidMatrix", plotCurve=FALSE, rankings=motifRankings, method="aprox")
  testthat::expect_equal(names(gsOut), c("enrStats", "incidMatrix"))


  ### Max aucMaxRank
  testthat::expect_error(getSignificantGenes(geneLists,motifRankings,
    signifRankingNames=selectedMotifs, plotCurve=TRUE, genesFormat="none",
    maxRank=50000, method="iCisTarget"))

  testthat::expect_error(getSignificantGenes(geneLists,motifRankings,
    signifRankingNames=selectedMotifs, plotCurve=TRUE, genesFormat="none",
    maxRank=50000, method="aprox"))

  testthat::expect_error(addSignificantGenes(motifEnrichmentTable,
    geneSets=geneLists, rankings=motifRankings, nCores=1,
    maxRank=50000, method="aprox"))
}
