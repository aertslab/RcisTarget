
.check_calcAUC <- function(motifRankings)
{
  motifRankings

  correctNcol <- ncol(motifRankings)-1

  # a) Character vector (i.e. only one gene-set)
  fewGenes <- sample(getRanking(motifRankings)$rn, 10)
  motifsAUC <- suppressWarnings(calcAUC(fewGenes, motifRankings, aucMaxRank=5))

  testthat::expect_equal(nrow(motifsAUC), 1)
  testthat::expect_equal(ncol(motifsAUC), correctNcol)

  # b) List
  otherGenes <- sample(getRanking(motifRankings)$rn, 5)
  geneSets <- list(geneSet1=fewGenes,
                   geneSet2=otherGenes)
  motifsAUC <- suppressWarnings(calcAUC(geneSets, motifRankings, aucMaxRank=5))

  testthat::expect_equal(nrow(motifsAUC), 2)
  testthat::expect_equal(ncol(motifsAUC), correctNcol)

  # c) GeneSet object (from GSEABase)
  geneSetOne <- GSEABase::GeneSet(fewGenes, setName="geneSetOne")
  motifsAUC <- suppressWarnings(calcAUC(geneSetOne, motifRankings, aucMaxRank=5))

  testthat::expect_equal(nrow(motifsAUC), 1)
  testthat::expect_equal(ncol(motifsAUC), correctNcol)

  # d) GeneSetCollection object (from GSEABase)
  geneSetTwo <- GSEABase::GeneSet(otherGenes, setName="geneSetTwo")
  geneSets <- GSEABase::GeneSetCollection(geneSetOne, geneSetTwo)
  motifsAUC <- suppressWarnings(calcAUC(geneSets, motifRankings, aucMaxRank=5))

  testthat::expect_equal(nrow(motifsAUC), 2)
  testthat::expect_equal(ncol(motifsAUC), correctNcol)

  testthat::expect_equal(class(motifsAUC)[1], "aucScores")
  testthat::expect_equal(SummarizedExperiment::assayNames(motifsAUC)[1], "AUC")

  ### Multicore
  motifsAUC_multicore <- suppressWarnings(calcAUC(geneSets, motifRankings, aucMaxRank=5, nCores=2))
  testthat::expect_equal(getAUC(motifsAUC), getAUC(motifsAUC_multicore))
}

