##################################################
# Setup & previous steps in the workflow:

#### Gene sets
# As example, the package includes an Hypoxia gene set:
txtFile <- paste(file.path(system.file('examples', package='RcisTarget')),
                 "hypoxiaGeneSet.txt", sep="/")
geneLists <- list(hypoxia=read.table(txtFile, stringsAsFactors=FALSE)[,1])

#### Databases
## Motif rankings: Select according to organism and distance around TSS
## (See the vignette for URLs to download)
# motifRankings <- importRankings("~/databases/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")

## For this example we will use a SUBSET of the ranking/motif databases:
library(RcisTarget.hg19.motifDBs.cisbpOnly.500bp)
data(hg19_500bpUpstream_motifRanking_cispbOnly)
motifRankings <- hg19_500bpUpstream_motifRanking_cispbOnly

### Run RcisTarget
# Step 1. Calculate AUC
motifs_AUC <- calcAUC(geneLists, motifRankings)

##################################################

#Exploring the output:
motifs_AUC

class(motifs_AUC)

# Extracting the AUC matrix:
getAUC(motifs_AUC)[,1:5]

# Subsetting and regular manipulation methods are also available:
motifs_AUC[1,]
motifs_AUC[,3:4]

dim(motifs_AUC)
nrow(motifs_AUC)
ncol(motifs_AUC)
colnames(motifs_AUC)
rownames(motifs_AUC)
