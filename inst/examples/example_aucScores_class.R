############# Fake run of RcisTarget ########
set.seed(123)
motifRankings <- fakeDatabase()
motifRankings

geneLists <- list(geneSet=sample(rownames(motifRankings), 100))

motifs_AUC <- calcAUC(geneLists, motifRankings)
##############################################

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
