#' @title convertToTargetRegions
#' @description Convert a set of input regions to the overlapping regions in the target set.
#' @param queryRegions List of regions to convert (normally the query region set).
#' @param targetRegions List of regions in the target set (normally the database).
#' @param minOverlap Minimum overlap to consider (in either direction, default: 0.40).
#' @param overlapType Parameter for \code{findOverlaps} (default: "any")
#' @param returnCorrespondence Returns a table containing the matches, 
#' or only a list of overlapping regions (default: FALSE).
#' @param verbose Print the number of regions selected?
#' @return IDs of the regions in the "target regions" overlapping with the query regions.
#' @seealso \code{\link{getDbRegionsLoc}}.
#'
#' See the package vignette for examples and more details:
#' \code{vignette("RcisTarget")}
#'
#' @examples
#' \dontrun{
#'  ## To apply on a list of regionSets:
#'  regionSets_db <- lapply(regionSets, function(x) 
#'     convertToTargetRegions(queryRegions=x, targetRegions=dbRegionsLoc))
#'  }
#' @rdname convertToTargetRegions
#' @importFrom GenomeInfoDb keepSeqlevels
#' @export
convertToTargetRegions <- function(queryRegions, targetRegions, minOverlap=0.4, overlapType="any", returnCorrespondence=FALSE, verbose=TRUE)
{
  ### Check types
  if(!is.numeric(minOverlap) || (minOverlap<0 && minOverlap>=1)) stop("minOverlap should be a number between 0 and 1 (percentage of overlap between the regions).")
  if(!isClass(targetRegions, "GRanges")) targetRegions <- GRanges(targetRegions)
  if(!isClass(queryRegions, "GRanges")) queryRegions <- GRanges(queryRegions)
  ###
  
  seqlvls <- intersect(seqlevels(queryRegions), seqlevels(targetRegions))
  queryRegions <- keepSeqlevels(queryRegions, seqlvls, pruning.mode = "coarse")
  overlapHits <- findOverlaps(queryRegions, targetRegions,
                              minoverlap=1,  
                              type=overlapType, select="all", ignore.strand=TRUE)
  
  if(minOverlap>0)
  {
    # In i-cisTarget, the default is 40% minimum overlap. Both ways: It takes the maximum percentage (of the peak or the ict region)
    # To reproduce those results:
    overlaps <- pintersect(queryRegions[queryHits(overlapHits)], targetRegions[subjectHits(overlapHits)])
    percentOverlapDb <- width(overlaps) / width(targetRegions[subjectHits(overlapHits)])
    percentOverlapQuery <- width(overlaps) / width(queryRegions[queryHits(overlapHits)])
    maxOverlap <- apply(cbind(percentOverlapDb, percentOverlapQuery), 1, max)
    overlapHits <- overlapHits[maxOverlap > 0.4]
  }
  
  # Get regions names
  dbHits <- targetRegions[subjectHits(overlapHits)]
  if("name" %in% colnames(elementMetadata(dbHits)))
  {
    dbHits <- as.character(dbHits$name)
  }else{
    dbHits <- as.character(dbHits)
  }
  ret <- unique(dbHits)
  if(verbose) message(paste("Number of regions selected: ", length(ret)))
  
  if(returnCorrespondence) 
  {
    queryHits <- queryRegions[queryHits(overlapHits)]
    ret <- cbind(query=unname(as.character(queryHits)), db=unname(dbHits))
  }
  
  return(ret)
}
