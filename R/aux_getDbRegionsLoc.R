#' @title Gets the region location for the given database IDs.
#' @description
#' Gets the region location based on the region ID.
#' Only needed for drosophila\bold{drosophila} (fruit fly) regions.
#' 
#' For \bold{human/mouse} the region locations are stored in a separate object 
#' (i.e. \code{}). This function is not needed.
#' @param featherFilePath
#' @param spltChr Character(s) used to split the prefix from the region location.
#' The default is used for current Drosophila versions.
#' @return The region locations in a GRanges object, with the original region ID as name.
#'
#' See the package vignette for examples and more details:
#' \code{vignette("RcisTarget")}
#'
#' @examples
#' \dontrun{
#' featherFilePath <- "~/databases/dm6-regions-11species.mc9nr.feather"
#' dbRegionsLoc <- getDbRegionsLoc(featherFilePath)
#' }
#' @rdname getDbRegionsLoc
#' @importFrom GenomeInfoDb sortSeqlevels
#' @importFrom GenomicRanges GRanges
#' @export
getDbRegionsLoc <- function(featherFilePath, spltChr="__")
{
  dbCols <- getColumnNames(featherFilePath)[-1]
  dbCols <- setNames(dbCols, dbCols)
  
  dbCols <- sapply(strsplit(dbCols, spltChr), function(x) x[[2]])
  dbRegions <- GenomicRanges::GRanges(unname(dbCols), name=names(dbCols))
  dbRegions <- GenomeInfoDb::sortSeqlevels(dbRegions)
  return(dbRegions)
}
