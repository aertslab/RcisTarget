
# Example for running RcisTarget using cisTarget() function (workflow wrapper)

# A quick example can be run using a (small) fake-database & gene-set:
load(system.file("examples", "fakeDb.RData", package="RcisTarget"))
motifRankings <- fakeDb$ranking
motifAnnot_direct <- fakeDb$annotation
geneLists <- sample(rownames(motifRankings), 100)

# This is an example of code for a real run:
##################################################
\dontrun{
  #### Gene sets
  # As example, the package includes an Hypoxia gene set:
  txtFile <- paste(file.path(system.file('examples', package='RcisTarget')),
                   "hypoxiaGeneSet.txt", sep="/")
  geneLists <- list(hypoxia=read.table(txtFile)[,1])

  #### Databases
  # Select the package according to the organism and distance around TSS
  library(RcisTarget.hg19.motifDBs.20k)
  data(hg19_10kbpAroundTss_motifRanking)
  motifRankings <- hg19_10kbpAroundTss_motifRanking
  data(hg19_direct_motifAnnotation)
  motifAnnotation <- hg19_direct_motifAnnotation

  # Run (R)cisTarget
  motifEnrichmentTable_wGenes <- cisTarget(geneLists, motifRankings,
   motifAnnot_direct=hg19_direct_motifAnnotation,
   nesThreshold=3.5, geneErnMethod="aprox", nCores=2)
}

# Load results from analysis
load(paste(file.path(system.file('examples', package='RcisTarget')),
           "motifEnrichmentTable_wGenes.RData", sep="/"))
##################################################


# Exploring the output:
# Note: If using the fake-database, the results are not meaningful

# Number of enriched motifs (Over the given NES threshold)
nrow(motifEnrichmentTable_wGenes)

# Available info (columns)
colnames(motifEnrichmentTable_wGenes)

# The object returned is a data.table (for faster computation),
# which has a diferent syntax from the standard data.frame or matrix
# Feel free to convert it to a data.frame (as.data.frame())
class(motifEnrichmentTable_wGenes)
motifEnrichmentTable_wGenes[,1:5]

# Enriched genes
enrGenes <- as.character(motifEnrichmentTable_wGenes[1,"enrichedGenes"])
strsplit(enrGenes, ";")


# Interactive exploration
motifEnrichmentTable_wGenes <- addLogo(motifEnrichmentTable_wGenes)
DT::datatable(motifEnrichmentTable_wGenes[,1:9], escape = FALSE, filter="top",
              options=list(pageLength=5))
# Note: If using the fake database, the results of this analysis are meaningless

