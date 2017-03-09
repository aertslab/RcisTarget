
##################################################
# Setup & previous steps in the workflow:
# Gene sets
txtFile <- paste(file.path(system.file('examples', package='RcisTarget')),
                 "hypoxiaGeneSet.txt", sep="/")
geneLists <- list(hypoxia=read.table(txtFile)[,1])

# Motif databases
# (Select the package/database according to the organism and distance around TSS)
library(RcisTarget.hg19.motifDatabases.20k)
data(hg19_10kbpAroundTss_motifRanking)
motifRankings <- hg19_10kbpAroundTss_motifRanking

# Fake-database with 5000 random motifs (to run the example faster)
# DO NOT use in real analyses!
set.seed(123)
motifRankings <- subset(hg19_10kbpAroundTss_motifRanking,
    sample(2:ncol(hg19_10kbpAroundTss_motifRanking@rankings), 5000), select="col")
motifRankings

# RcisTarget
# Step 1. Calculate AUC
motifs_AUC <- calcAUC(geneLists, motifRankings)

##################################################
# (This step: Step 2)
# Select significant motifs, add TF annotation & format as table
data(hg19_direct_motifAnnotation)

motifEnrichmentTable <- addMotifAnnotation(motifs_AUC,
   motifAnnot_direct=hg19_direct_motifAnnotation)

# Adding indirect annotation and modifying some options
data(hg19_bySimilarity_motifAnnotation)
motifEnrichment_wIndirect <- addMotifAnnotation(motifs_AUC, nesThreshold=2,
                        motifAnnot_direct=hg19_direct_motifAnnotation,
                        motifAnnot_bySimilarity=hg19_bySimilarity_motifAnnotation,
                        highlightTFs = "HIF1A", digits=3)

# Exploring the output:
# Note: Using the fake-database, these results are not meaningful.

# Number of enriched motifs (Over the given NES threshold)
nrow(motifEnrichmentTable)

# Interactive exploration
motifEnrichmentTable <- addLogo(motifEnrichmentTable)
DT::datatable(motifEnrichmentTable, filter="top", escape=FALSE, options=list(pageLength=50))
# Note: If using the fake database, the results of this analysis are meaningless

# The object returned is a data.table (for faster computation),
# which has a diferent syntax from the standard data.frame or matrix
# Feel free to convert it to a data.frame (as.data.frame())
motifEnrichmentTable[,1:6]


##################################################
# Next step (step 3, optional):
motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable,
                                   geneSets=geneLists,
                                   rankings=motifRankings,
                                   method="aprox")

