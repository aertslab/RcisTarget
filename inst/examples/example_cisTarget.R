
# Example for running RcisTarget using cisTarget() function (workflow wrapper)

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
# Load the motif databases
library(RcisTarget.hg19.motifDatabases)
data(hg19_10kbpAroundTss_motifRanking)
motifRankings <- hg19_10kbpAroundTss_motifRanking

data(hg19_direct_motifAnnotation)

# This example is run using a fake-database with 5000 random motifs
# (for faster execution, only to explore the workflow & interface)
# DO NOT use in real analyses! To identify statistically significant motifs,
# RcisTarget needs the whole motif database.
set.seed(123)
motifRankings <- hg19_10kbpAroundTss_motifRanking[,c("rn",
        sample(colnames(hg19_10kbpAroundTss_motifRanking), 5000)), with=FALSE]


##################################################
# Run (R)cisTarget

motifEnrichmentTable_wGenes <- cisTarget(geneLists, motifRankings,
                       motifAnnot_direct=hg19_direct_motifAnnotation,
                       highlightTFs="HIF1A", nesThreshold=5,
                       geneErnMethod="aprox", nCores=2)

##################################################
# Exploring the output:
# Note: Using the fake-database, these results are not meaningful.

# Number of enriched motifs (Over the given NES threshold)
nrow(motifEnrichmentTable_wGenes)

# Available info (columns)
colnames(motifEnrichmentTable_wGenes)

# The object returned is a data.table (for faster computation),
# which has a diferent syntax from the standard data.frame or matrix
class(motifEnrichmentTable_wGenes)
motifEnrichmentTable_wGenes[,1:9, with=FALSE]

# Feel free to convert it to a data.frame:
motifEnrichmentTable_wGenes <- as.data.frame(motifEnrichmentTable_wGenes)
motifEnrichmentTable_wGenes[,1:9]

# Enriched genes
enrGenes <- motifEnrichmentTable_wGenes[1,"enrichedGenes"]
strsplit(enrGenes, ";")

motifEnrichmentTable_wGenes <- addLogo(motifEnrichmentTable_wGenes)

# Interactive exploration
library(DT)
datatable(motifEnrichmentTable_wGenes[,1:9], escape = FALSE, filter="top", options=list(pageLength=5))
