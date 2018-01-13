
##################################################
# Setup & previous steps in the workflow:

#### Gene sets
# As example, the package includes an Hypoxia gene set:
txtFile <- paste(file.path(system.file('examples', package='RcisTarget')),
                 "hypoxiaGeneSet.txt", sep="/")
geneLists <- list(hypoxia=read.table(txtFile, stringsAsFactors=FALSE)[,1])

#### Databases
# Select the package/database according to the organism and distance around TSS
# Recommended: Full database
# library(RcisTarget.hg19.motifDBs.20k)
# data(hg19_10kbpAroundTss_motifRanking)
# motifRankings <- hg19_10kbpAroundTss_motifRanking
# data(hg19_motifAnnotation)
# motifAnnotation <- hg19_motifAnnotation

# For the example: Subset of database
library(RcisTarget.hg19.motifDBs.cisbpOnly.500bp)
data(hg19_500bpUpstream_motifRanking_cispbOnly)
  motifRankings <- hg19_500bpUpstream_motifRanking_cispbOnly
data(hg19_motifAnnotation_cisbpOnly)
  motifAnnotation <- hg19_motifAnnotation_cisbpOnly

### Run RcisTarget
# Step 1. Calculate AUC
motifs_AUC <- calcAUC(geneLists, motifRankings)

##################################################

### (This step: Step 2)
# Select significant motifs, add TF annotation & format as table
motifEnrichmentTable <- addMotifAnnotation(motifs_AUC,
                                           motifAnnot=motifAnnotation)

# Alternative: Modifying some options
motifEnrichment_wIndirect <- addMotifAnnotation(motifs_AUC, nesThreshold=2, 
        motifAnnot=motifAnnotation,
        highlightTFs = "HIF1A", 
        motifAnnot_highConfCat=c("directAnnotation"), 
        motifAnnot_lowConfCat=c("inferredBy_MotifSimilarity",
                                "inferredBy_MotifSimilarity_n_Orthology",
                                "inferredBy_Orthology"),
        digits=3)

### Exploring the output:
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
