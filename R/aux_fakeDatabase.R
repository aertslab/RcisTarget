#' @title Creates a small fake database to use in examples
#' @aliases fakeDatabase
#' @param nGenes Number of genes to include in the ranking
#' @param nMotifs Number of motifs to include in the ranking
#' @description
#' Creates a small fake database to use in examples
#' @export
fakeDatabase <- function(nGenes=1000, nMotifs=500)
{
  fakeRankings <- data.table(rn=paste0("gene", sprintf( "%04d",1:nGenes)),
             sapply(paste0("motif_", 1:nMotifs), function(x) sample(1:nGenes)))
  setkey(fakeRankings, "rn")

  new("rankingWrapper", rankings=fakeRankings, colType="motif", rowType="gene", org="fake", genome="-",description="Fake database, only for examples")
}

