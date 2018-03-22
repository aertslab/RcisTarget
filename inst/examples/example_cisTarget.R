

# Example for running RcisTarget using cisTarget() function (workflow wrapper)

\dontrun{

##################################################
### Load your gene sets
# As example, the package includes an Hypoxia gene set:
txtFile <- paste(file.path(system.file('examples', package='RcisTarget')),
                 "hypoxiaGeneSet.txt", sep="/")
geneLists <- list(hypoxia=read.table(txtFile, stringsAsFactors=FALSE)[,1])

### Load databases
# Motif rankings: Select according to organism and distance around TSS
# (See the vignette for URLs to download)
motifRankings <- importRankings("hg19-500bp-upstream-7species.mc9nr.feather")

# Motif - TF annotation:
data(motifAnnotations_hgnc) # human TFs (for motif collection 9)
motifAnnotation <- motifAnnotations_hgnc
##################################################

# Run (R)cisTarget
motifEnrichmentTable_wGenes <- cisTarget(geneLists, motifRankings,
  motifAnnot_direct=hg19_direct_motifAnnotation,
  nesThreshold=3.5, geneErnMethod="aprox", nCores=2)

}

# Load results from analysis
load(paste(file.path(system.file('examples', package='RcisTarget')),
           "motifEnrichmentTable_wGenes.RData", sep="/"))


### Exploring the output:
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

