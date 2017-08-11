
##################################################
# Setup & previous steps in the workflow:

#### Gene sets
# As example, the package includes an Hypoxia gene set:
txtFile <- paste(file.path(system.file('examples', package='RcisTarget')),
                 "hypoxiaGeneSet.txt", sep="/")
geneLists <- list(hypoxia=read.table(txtFile)[,1])

#### Databases
# Select the package/database according to the organism and distance around TSS
library(RcisTarget.hg19.motifDatabases.20k)
data(hg19_10kbpAroundTss_motifRanking)
motifRankings <- hg19_10kbpAroundTss_motifRanking

### Run RcisTarget
# Step 1. Calculate AUC
motifs_AUC <- calcAUC(geneLists, motifRankings)

##################################################

### (This step: Step 2)
# Select significant motifs, add TF annotation & format as table
data(hg19_direct_motifAnnotation)

motifEnrichmentTable <- addMotifAnnotation(motifs_AUC,
   motifAnnot_direct=hg19_direct_motifAnnotation)

# Alternative: Adding indirect annotation and modifying some options
data(hg19_inferred_motifAnnotation)
motifEnrichment_wIndirect <- addMotifAnnotation(motifs_AUC, nesThreshold=2,
                        motifAnnot_direct=hg19_direct_motifAnnotation,
                        motifAnnot_inferred=hg19_inferred_motifAnnotation,
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
