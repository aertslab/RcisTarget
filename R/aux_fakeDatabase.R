#' @title Creates a small fake database to use in examples
#' @aliases fakeDatabase
#' @param nGenes Number of genes to include in the ranking
#' @param nMotifs Number of motifs to include in the ranking
#' @param incAnnotation Whether to also generate a fake motif annotation
#' @description
#' Creates a small fake database to use in examples
#' @return
#' Returns a small fake database, with random content but the right structure,
#' to run the examples.
#' @examples
#' fakeDatabase()
#' @export
fakeDatabase <- function(nGenes=1000, nMotifs=500, incAnnotation=FALSE)
{
  fakeRankings <- data.table(rn=paste0("gene", sprintf( "%04d",1:nGenes)),
             sapply(paste0("motif_", 1:nMotifs), function(x) sample(1:nGenes)))
  setkey(fakeRankings, "rn")

  fakeRanking <- new("rankingWrapper",
      rankings=fakeRankings,
      colType="motif",
      rowType="gene",
      org="fake",
      genome="-",
      description="Fake database, only for examples")

  ret <- fakeRanking
  if(incAnnotation)
  {
    fakeAnnotation <- setNames(lapply(seq_len(ncol(fakeRanking)-1),
                                      function(x){
      cbind(TF=sort(sample(getRanking(fakeRanking)$rn[1:30], sample(1:3,1))))
    }), colnames(getRanking(fakeRanking))[-1])
    ret <- list(ranking=fakeRanking, annotation=fakeAnnotation)
  }

  return(ret)
}

