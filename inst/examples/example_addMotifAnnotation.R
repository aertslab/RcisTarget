
##################################################
# Setup & previous steps in the workflow:

#### Gene sets
# As example, the package includes an Hypoxia gene set:
txtFile <- paste(file.path(system.file('examples', package='RcisTarget')),
                 "hypoxiaGeneSet.txt", sep="/")
geneLists <- list(hypoxia=read.table(txtFile)[,1])

#### Databases
# Select the package/database according to the organism and distance around TSS
# Recommended: Full database
# library(RcisTarget.hg19.motifDBs.20k)
# data(hg19_10kbpAroundTss_motifRanking)
# motifRankings <- hg19_10kbpAroundTss_motifRanking
# data(hg19_direct_motifAnnotation)
# motifAnnotation_direct <- hg19_direct_motifAnnotation
# data(hg19_inferred_motifAnnotation)
# motifAnnotation_inferred <- hg19_inferred_motifAnnotation

# For the example: Subset of database
library(RcisTarget.hg19.motifDBs.cisbpOnly.500bp)
data(hg19_500bpUpstream_motifRanking_cispbOnly)
  motifRankings <- hg19_500bpUpstream_motifRanking_cispbOnly
data(hg19_direct_motifAnnotation_cisbpOnly)
  motifAnnotation_direct <- hg19_direct_motifAnnotation_cisbpOnly
data(hg19_inferred_motifAnnotation_cisbpOnly)
  motifAnnotation_inferred <- hg19_inferred_motifAnnotation_cisbpOnly

### Run RcisTarget
# Step 1. Calculate AUC
motifs_AUC <- calcAUC(geneLists, motifRankings)

##################################################

### (This step: Step 2)
# Select significant motifs, add TF annotation & format as table
motifEnrichmentTable <- addMotifAnnotation(motifs_AUC,
                                           motifAnnot_direct=motifAnnotation_direct)

# Alternative: Adding indirect annotation and modifying some options
motifEnrichment_wIndirect <- addMotifAnnotation(motifs_AUC, nesThreshold=2,
                                                motifAnnot_direct=motifAnnotation_direct,
                                                motifAnnot_inferred=motifAnnotation_inferred,
                                                highlightTFs = "HIF1A", digits=3)

### Exploring the output:
# Note: Using the fake-database, these results are not meaningful.

# Number of enriched motifs (Over the given NES threshold)
nrow(motifEnrichmentTable)

# Interactive exploration
motifEnrichmentTable <- addLogo(motifEnrichmentTable)
DT::datatable(motifEnrichmentTable, filter="top", escape=FALSE,
              options=list(pageLength=50))
# Note: If using the fake database, the results of this analysis are meaningless

# The object returned is a data.table (for faster computation),
# which has a diferent syntax from the standard data.frame or matrix
# Feel free to convert it to a data.frame (as.data.frame())
motifEnrichmentTable[,1:6]


##################################################
# Next step (step 3, optional):
\dontrun{
  motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable,
                                                     geneSets=geneLists,
                                                     rankings=motifRankings,
                                                     method="aprox")
}
