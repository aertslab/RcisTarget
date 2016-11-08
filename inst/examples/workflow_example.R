
# RcisTarget workflow (for advanced users: Running the workflow steps individually)

##################################################
# Format gene sets

# RcisTarget requires the gene sets to be stored in a list with this format:
genes <- c("gene1", "gene2", "gene3")
geneLists <- list(geneSet1=genes)
# (The name 'geneSet1' can be changed)

# In this example we will use an Hypoxia gene set included in the package:
txtFile <- paste(file.path(system.file('examples', package='RcisTarget')),
                 "hypoxiaGeneSet.txt", sep="/")
geneLists <- list(hypoxia=read.table(txtFile)[,1])


##################################################
# Load databases

# Select the package/database according to the organism and distance around TSS
# Load the motif rankings to calculate the AUC
library(RcisTarget.hg19.motifDatabases)
data(hg19_10kbpAroundTss_motifRanking)
motifRankings <- hg19_10kbpAroundTss_motifRanking

# This example is run using a fake-database with 5000 random motifs
# (for faster execution, only to explore the workflow & interface)
# DO NOT use in real analyses! To identify statistically significant motifs,
# RcisTarget needs the whole motif database.
set.seed(123)
motifRankings <- hg19_10kbpAroundTss_motifRanking[,c("rn",
        sample(colnames(hg19_10kbpAroundTss_motifRanking), 5000)), with=FALSE]

##################################################
# Run RcisTarget

# Step 1. Calculate AUC
motifs_AUC <- calcAUC(geneLists, motifRankings)

# Step 2. Select significant motifs, add TF annotation & format as table
data(hg19_direct_motifAnnotation)
motifEnrichmentTable <- addMotifAnnotation(motifs_AUC,
                              motifAnnot_direct=hg19_direct_motifAnnotation)

# Step 3 (optional). Identify the genes that have the motif significantly enriched
# (i.e. genes from the gene set in the top of the ranking)
motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable,
                                                      geneSets=geneLists,
                                                      rankings=motifRankings,
                                                      method="aprox")

# Note: Using the fake-database, these results are not meaningful.
