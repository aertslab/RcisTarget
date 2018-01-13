
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
# Step 2. Select significant motifs, add TF annotation & format as table
motifEnrichmentTable <- addMotifAnnotation(motifs_AUC,
           motifAnnot=motifAnnotation)

##################################################

##################################################
# (This step: Step 3)
# Identify the genes that have the motif significantly enriched
# (i.e. genes from the gene set in the top of the ranking)
par(mfrow=c(1,2))
motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable,
                                       genesFormat="geneList",
                                       plotCurve=TRUE,
                                       geneSets=geneLists,
                                       rankings=motifRankings,
                                       method="aprox")

#### Exploring the output:
# Note: Using the fake-database, these results are not meaningful.

# The object returned is a data.table
# Feel free to convert it to a data.frame:
motifEnrichmentTable_wGenes <- as.data.frame(motifEnrichmentTable_wGenes)

# Enriched genes
enrGenes <- motifEnrichmentTable_wGenes[1,"enrichedGenes"]
enrGenes
strsplit(enrGenes, ";")

# As incidence matrix
motifEnr_wIncidMat <- addSignificantGenes(motifEnrichmentTable,
                geneSets=geneLists, rankings=motifRankings,
                method="aprox",
                genesFormat = "incidMatrix")

motifEnr_wIncidMat <- as.data.frame(motifEnr_wIncidMat)
which(colnames(motifEnr_wIncidMat) == "rankAtMax")

incidMat <- motifEnr_wIncidMat[,8:ncol(motifEnr_wIncidMat)]
rownames(incidMat) <- motifEnr_wIncidMat[,"motif"]
incidMat <- incidMat[, colSums(incidMat)>0, drop=FALSE]

# Plot as network
par(mfrow=c(1,1))
library(igraph)
plot(graph.incidence(incidMat))

###############################################################
# Alternative method: getSignificantGenes()
selectedMotif <- rownames(incidMat)
onlyGenes <- getSignificantGenes(geneSet=geneLists$hypoxia,
                            signifRankingNames=selectedMotif,
                            genesFormat="incidMatrix",
                            plotCurve=TRUE,
                            rankings=motifRankings,
                            method="aprox")


