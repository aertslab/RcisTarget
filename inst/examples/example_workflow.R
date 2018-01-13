
# RcisTarget workflow for advanced users:
# Running the workflow steps individually

# This example is run using a (small) fake-database & gene-set:
load(system.file("examples", "fakeDb.RData", package="RcisTarget"))
motifRankings <- fakeDb$ranking
motifAnnot <- fakeDb$annotation

geneLists <- list(geneSet=sample(rownames(motifRankings), 100))

# Example of code for a real run:
##################################################
\dontrun{
  #### Setup gene sets
  # RcisTarget requires the gene sets to be stored in a list with this format:
  genes <- c("gene1", "gene2", "gene3")
  geneLists <- list(geneSet1=genes)
  # (The name 'geneSet1' can be changed)

  # In this example we will use an Hypoxia gene set included in the package:
  txtFile <- paste(file.path(system.file('examples', package='RcisTarget')),
                   "hypoxiaGeneSet.txt", sep="/")
  geneLists <- list(hypoxia=read.table(txtFile)[,1])

  #### Load databases
  # Select the package/database according to the organism and distance around TSS
  library(RcisTarget.hg19.motifDBs.20k)
  data(hg19_10kbpAroundTss_motifRanking)
  motifRankings <- hg19_10kbpAroundTss_motifRanking
  data(hg19_motifAnnotation)
  motifAnnot <- hg19_motifAnnotation
}
##################################################


#### Run RcisTarget

# Step 1. Calculate AUC
motifs_AUC <- calcAUC(geneLists, motifRankings)

# Step 2. Select significant motifs, add TF annotation & format as table
motifEnrichmentTable <- addMotifAnnotation(motifs_AUC,
                         motifAnnot=motifAnnot)

# Step 3 (optional). Identify genes that have the motif significantly enriched
# (i.e. genes from the gene set in the top of the ranking)
motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable,
                                                   geneSets=geneLists,
                                                   rankings=motifRankings,
                                                   method="aprox")
