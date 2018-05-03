
# test_Workflow() calls:
## .check_calcAUC()
## .check_addSignificantGenes()
## .check_addMotifAnnotation()
## .check_motifEnrichmentTable()
.check_calcAUC <- function(motifRankings)
{
  motifRankings
  
  correctNmotifs <- nrow(motifRankings)
  
  # a) Character vector (i.e. only one gene-set)
  fewGenes <- sample(colnames(getRanking(motifRankings)), 10)
  motifsAUC <- suppressWarnings(calcAUC(fewGenes, motifRankings, aucMaxRank=5))
  
  testthat::expect_equal(nrow(motifsAUC), 1)
  testthat::expect_equal(ncol(motifsAUC), correctNmotifs)
  
  # b) List
  otherGenes <- sample(colnames(getRanking(motifRankings)), 5)
  geneSets <- list(geneSet1=fewGenes,
                   geneSet2=otherGenes)
  motifsAUC <- suppressWarnings(calcAUC(geneSets, motifRankings, aucMaxRank=5))
  
  testthat::expect_equal(nrow(motifsAUC), 2)
  testthat::expect_equal(ncol(motifsAUC), correctNmotifs)
  
  # c) GeneSet object (from GSEABase)
  geneSetOne <- GSEABase::GeneSet(fewGenes, setName="geneSetOne")
  motifsAUC <- suppressWarnings(calcAUC(geneSetOne, motifRankings, aucMaxRank=5))
  
  testthat::expect_equal(nrow(motifsAUC), 1)
  testthat::expect_equal(ncol(motifsAUC), correctNmotifs)
  
  # d) GeneSetCollection object (from GSEABase)
  geneSetTwo <- GSEABase::GeneSet(otherGenes, setName="geneSetTwo")
  geneSets <- GSEABase::GeneSetCollection(geneSetOne, geneSetTwo)
  motifsAUC <- suppressWarnings(calcAUC(geneSets, motifRankings, aucMaxRank=5))
  
  testthat::expect_equal(nrow(motifsAUC), 2)
  testthat::expect_equal(ncol(motifsAUC), correctNmotifs)
  
  testthat::expect_equal(class(motifsAUC)[1], "aucScores")
  testthat::expect_equal(SummarizedExperiment::assayNames(motifsAUC)[1], "AUC")
  
  ### Max aucMaxRank
  testthat::expect_error(calcAUC(geneSets, motifRankings, nCores=1, aucMaxRank=100000))
  
  ### Multicore
  motifsAUC_multicore <- suppressWarnings(calcAUC(geneSets, motifRankings, aucMaxRank=5, nCores=2))
  testthat::expect_equal(getAUC(motifsAUC), getAUC(motifsAUC_multicore))
}

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
               c("geneSet", "motif", "NES", "AUC", "TF_highConf","TF_lowConf", 
                 "nEnrGenes", "rankAtMax", "enrichedGenes"))
  
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
                                   maxRank=4500,
                                   rankings=motifRankings,
                                   method="aprox")
  
  
  ### Diferent input formats (geneSet)
  # a) Character vector (i.e. only one gene-set)
  fewGenes <- sample(colnames(getRanking(motifRankings)), 10)
  gsOut <- getSignificantGenes(geneSet=fewGenes, signifRankingNames=selectedMotif[1], maxRank=4500,
                               genesFormat="incidMatrix", plotCurve=FALSE, rankings=motifRankings, method="aprox")
  
  testthat::expect_equal(names(gsOut), c("enrStats", "incidMatrix"))
  
  # b) List
  otherGenes <- sample(colnames(getRanking(motifRankings)), 5)
  geneSets <- list(geneSet1=fewGenes,
                   geneSet2=otherGenes)
  expect_error(getSignificantGenes(geneSet=geneSets, signifRankingNames=selectedMotif[1], maxRank=4500,
                                   genesFormat="incidMatrix", plotCurve=FALSE, rankings=motifRankings, method="aprox"))
  gsOut <- getSignificantGenes(geneSet=geneSets[1], signifRankingNames=selectedMotif[1], maxRank=4500,
                               genesFormat="incidMatrix", plotCurve=FALSE, rankings=motifRankings, method="aprox")
  testthat::expect_equal(names(gsOut), c("enrStats", "incidMatrix"))
  
  # c) GeneSet object (from GSEABase)
  geneSetOne <- GSEABase::GeneSet(fewGenes, setName="geneSetOne")
  gsOut <- getSignificantGenes(geneSet=geneSetOne, signifRankingNames=selectedMotif[1], maxRank=4500,
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

.check_addMotifAnnotation <- function(motifs_AUC, motifAnnot)
{
  # No annot
  motifEnrichmentTable_noAnnot <- addMotifAnnotation(motifs_AUC)
  .check_motifEnrichmentTable(met=motifEnrichmentTable_noAnnot, motifAnnot=NULL)
  
  # provide annotation but confCat=NULL
  motifEnrichmentTable <- addMotifAnnotation(motifs_AUC,
                                             motifAnnot=motifAnnot,
                                             motifAnnot_highConfCat=NULL,
                                             motifAnnot_lowConfCat=NULL)
  .check_motifEnrichmentTable(met=motifEnrichmentTable,
                              motifAnnot=motifAnnot,
                              motifAnnot_highConfCat=NULL,
                              motifAnnot_lowConfCat=NULL)
  
  motifEnrichmentTable <- addMotifAnnotation(motifs_AUC,
                                             motifAnnot=motifAnnot,
                                             motifAnnot_highConfCat=NULL)
  .check_motifEnrichmentTable(met=motifEnrichmentTable,
                              motifAnnot=motifAnnot,
                              motifAnnot_highConfCat=NULL)
  
  # including annotation 
  motifEnrichmentTable <- addMotifAnnotation(motifs_AUC,
                                             motifAnnot=motifAnnot)
  .check_motifEnrichmentTable(met=motifEnrichmentTable, motifAnnot=motifAnnot)
  
  # some other options
  motifEnrichment_wIndirect <- suppressWarnings(addMotifAnnotation(motifs_AUC, 
                                                                   nesThreshold=2, motifAnnot=motifAnnot,
                                                                   highlightTFs = "HIF1A", digits=3))
  .check_motifEnrichmentTable(met=motifEnrichment_wIndirect, 
                              motifAnnot=motifAnnot,
                              highlightTFs="HIF1A")
  expect_true(min(motifEnrichment_wIndirect$NES)>=2)
}


.check_motifEnrichmentTable <- function(met, motifAnnot=NULL, 
                                        motifAnnot_highConfCat=c("directAnnotation", "inferredBy_Orthology"),
                                        motifAnnot_lowConfCat=c("inferredBy_MotifSimilarity",
                                                                "inferredBy_MotifSimilarity_n_Orthology"),
                                        highlightTFs=NULL)
{
  nCols <- 4
  
  if(!is.null(motifAnnot)) {
    
    if(!is.null(motifAnnot_highConfCat))
    {
      nCols <- nCols+1
      expect_true("TF_highConf" %in% colnames(met)) 
      
      # check content (not very exhaustive)
      expect_true(any(vapply(paste0(" (", motifAnnot_highConfCat,")."), 
                             function(x) any(grepl(x, met$TF_highConf, fixed=TRUE)), TRUE)))
      
      expect_false(any(vapply(paste0(" (", motifAnnot_lowConfCat,")."), 
                              function(x) any(grepl(x, met$TF_highConf, fixed=TRUE)), TRUE)))
    }
    
    if(!is.null(motifAnnot_lowConfCat))
    {
      nCols <- nCols+1
      expect_true("TF_lowConf" %in% colnames(met)) 
      
      # check content (not very exhaustive)
      expect_true(any(vapply(paste0(" (", motifAnnot_lowConfCat,")."), 
                             function(x) any(grepl(x, met$TF_lowConf, fixed=TRUE)), TRUE)))
      
      expect_false(any(vapply(paste0(" (", motifAnnot_highConfCat,")."), 
                              function(x) any(grepl(x, met$TF_lowConf, fixed=TRUE)), TRUE)))
    }
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

test_Workflow <- function()
{
  library(RcisTarget)
  
  ##### Example setup ##############################
  txtFile <- paste(file.path(system.file('examples', package='RcisTarget')),
                   "hypoxiaGeneSet.txt", sep="/")
  geneLists <- list(hypoxia=read.table(txtFile, stringsAsFactors = FALSE)[,1])

  # Load databases
  library(RcisTarget.hg19.motifDBs.cisbpOnly.500bp)
  data(hg19_500bpUpstream_motifRanking_cispbOnly)
  motifRankings <- hg19_500bpUpstream_motifRanking_cispbOnly
  
  ##################################################
  # Standard workflow
  # (To use as input for the different check functions)
  # Step 1.

  .check_calcAUC(motifRankings) # Checks different intput/ouputs for calcAUC
  motifs_AUC <- calcAUC(geneLists, motifRankings) # For input of next steps

  # Step 2.
  data(hg19_motifAnnotation_cisbpOnly)
  .check_addMotifAnnotation(motifs_AUC, 
                            motifAnnot=hg19_motifAnnotation_cisbpOnly)
  motifEnrichmentTable <- addMotifAnnotation(motifs_AUC,
                            motifAnnot=hg19_motifAnnotation_cisbpOnly)

  # Step 3
  .check_addSignificantGenes(met=motifEnrichmentTable, 
                             motifRankings=motifRankings, geneLists=geneLists)
  motifEnrichmentTable_wGenes <- addSignificantGenes(
                            resultsTable=motifEnrichmentTable,
                            geneSets=geneLists,
                            rankings=motifRankings,
                            method="aprox")
  ##################################################


  ##################################################
  # cisTarget (automated workflow)
  motifEnrichmentTable_wGenes <- cisTarget(geneLists, motifRankings,
                     motifAnnot=hg19_motifAnnotation_cisbpOnly,
                     nesThreshold=3.5, geneErnMethod="aprox", nCores=4)

  expect_equal(colnames(motifEnrichmentTable_wGenes),
               c("geneSet", "motif", "NES", "AUC", "TF_highConf", "TF_lowConf", 
                 "nEnrGenes", "rankAtMax", "enrichedGenes"))
  ##################################################


  motifEnrichmentTable_wLogo <- addLogo(motifEnrichmentTable)
}

test_that("RcisTarget workflow tests", test_Workflow())



