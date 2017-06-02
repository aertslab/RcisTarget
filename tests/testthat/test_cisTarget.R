## cisTarget()

test_cisTarget <- function()
{

  txtFile <- paste(file.path(system.file('examples', package='RcisTarget')),
                   "hypoxiaGeneSet.txt", sep="/")
  geneLists <- list(hypoxia=read.table(txtFile)[,1])

  # Load databases
  library(RcisTarget.hg19.motifDatabases.20k)
  data(hg19_direct_motifAnnotation)
  data(hg19_10kbpAroundTss_motifRanking)
  motifRankings <- hg19_10kbpAroundTss_motifRanking

  # Fake example with 5000 random motifs
  library(data.table)
  set.seed(123)
  motifRankings <- subset(motifRankings,
                   sample(2:ncol(motifRankings@rankings), 5000), select="col")

  ##################################################
  # Run (R)cisTarget
  motifEnrichmentTable_wGenes <- cisTarget(geneLists, motifRankings,
                           motifAnnot_direct=hg19_direct_motifAnnotation,
                           nesThreshold=3.5, geneErnMethod="aprox", nCores=2)

  ##################################################
  expect_equal(ncol(motifEnrichmentTable_wGenes),8)
}
