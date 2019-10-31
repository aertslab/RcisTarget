# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)

#' @title Re-rank RcisTarget ranking
#' @description Re-ranks the genes/regions in the database for each motif. This allows to do motif enrichment over a background.
#' @param rankingsDb Results from RcisTarget (data.table)
#' @param columns Whether to add the HTML tag <img> around the URL or not
#' (boolean).
#' @param indexCol Name of the column (in the ranking) which contains the motif/feature ID.
#' @return Returns a new ranking database with the new ranking values.
#' @seealso See the "background" vignette for more examples:
#' \code{vignette("RcisTarget-withBackground")}
#' @example inst/examples/example_reRank.R
#' @export
reRank <- function(rankingsDb, columns=NULL, indexCol="features")
{
  rankingMat <- getRanking(rankingsDb)
  if(!is.null(columns)) rankingMat <- rankingMat[,unique(c(indexCol, columns))]
  
  # Re-rank the genes...
  featureNames <- unname(unlist(rankingMat[,1]))
  rankingMat <- t(apply(rankingMat[,-1], 1, rank, ties.method="random"))
  
  mode(rankingMat) <- "integer"
  rankingMat <- tibble::as.tibble(rankingMat)
  rankingMat <- tibble::add_column(rankingMat, features=featureNames, .before = 1)
  
  # Return
  reRanked <- new("rankingRcisTarget",
      rankings=rankingMat,
      colType=rankingsDb@colType,
      rowType=rankingsDb@rowType,
      org=rankingsDb@org,
      genome=rankingsDb@genome,
      nColsInDB=ncol(rankingMat)-1,
      maxRank=Inf,
      description=rankingsDb@description)
  return(reRanked)
}