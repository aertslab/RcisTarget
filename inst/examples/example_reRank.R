
library(RcisTarget.hg19.motifDBs.cisbpOnly.500bp)
data(hg19_500bpUpstream_motifRanking_cispbOnly)
motifRankings <- hg19_500bpUpstream_motifRanking_cispbOnly


genes <- colnames(getRanking(motifRankings))[10:20]
reRank(motifRankings, columns=genes)
