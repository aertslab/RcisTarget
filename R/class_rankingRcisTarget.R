#' @title Class to store the motif databases for RcisTarget.
#' @aliases rankingRcisTarget-class, rankingRcisTarget, getRanking,
#' ncol, nrow, show, subset
#' @param x [several methods] rankingRcisTarget object to apply the method
#' @param object [method: show] rankingRcisTarget object to show
#'
#' @description
#' This class is only meant as container for the motif rankings
#' (for internal use). Modify content at your own risk.
#'
#' Slots:
#' \itemize{
#' \item rankings: data.table with the rankings
#' \item colType: 'motif' or whatever other feature is stored (e.g. ChipSeq)
#' \item rowType: 'gene'or 'region'
#' \item org: human/mouse/fly
#' \item genome: hg19, mm9, ...
#' \item description: summary or any other information
#' }
#'
## Methods: See examples.
## @import BiocGenerics
#'
#' @importClassesFrom data.table data.table
#' @importFrom data.table setkey is.data.table rbindlist .SD
#   importMethodsFrom
#'
## @example inst/examples/example_class_rcos.R
#' @return
#' \itemize{
#' \item show: Prints a summary of the object
#' \item getRanking: Returns the data.frame containing the rankings
#' \item ncol, nrow, subset: Returns the number of columns, rows
#' or a subset of the ranking
#' }
#' @examples
#' library(RcisTarget.hg19.motifDBs.cisbpOnly.500bp)
#' data("hg19_500bpUpstream_motifRanking_cispbOnly")
#' hg19_500bpUpstream_motifRanking_cispbOnly
#' class(hg19_500bpUpstream_motifRanking_cispbOnly)

#' @rdname rankingRcisTarget-class
#' @export rankingRcisTarget
#' @exportClass rankingRcisTarget
rankingRcisTarget <- setClass(
  # Set the name for the class
  Class="rankingRcisTarget",

  contains="SummarizedExperiment",
  representation=representation(
    colType = "character", # motif or whatever feature is stored (e.g. ChipSeq)
    rowType = "character", # gene/region
    org = "character", # human/mouse
    genome = "character", # hg19, mm9 ...
    maxRank = "numeric", # Higher ranks are converted to Inf
    description = "character" # Other info, shown with "show"
  )
)




#' @rdname rankingRcisTarget-class
#' @aliases show,rankingRcisTarget-method
#' @export
setMethod("show",
  signature="rankingRcisTarget",
  definition = function(object) {
    message <- paste("Rankings for RcisTarget.", "\n",
    "  Organism: ", object@org, " (", object@genome,")","\n",
    "  Number of ", toupper(object@rowType),"S: ",nrow(getRanking(object)),"\n",
    "  Number of ", object@colType, "s: ", ncol(getRanking(object)),"\n\n",
    sep="")

    if(getMaxRank(object) < Inf)
    {
      message <- paste(message,
                      "** This database includes rankings up to ",
                      getMaxRank(object), "\n", sep="")
    }

    if(length(object@description)>0)
    {
      message <- paste(message, "\n",
                       object@description,
                       "\n", sep="")
    }


    cat(message)
  }
)

##### Access the rankings:

#' @rdname rankingRcisTarget-class
#' @aliases getRanking,rankingRcisTarget-method
#' @export
setGeneric(name="getRanking",
           def=function(x) standardGeneric("getRanking"))
setMethod("getRanking",
  signature="rankingRcisTarget",
  definition = function(x) {
    if("rankings" %in% assayNames(x)) {
      x@assays[["rankings"]] # assays() takes over 1 second...
    }else{
      stop("This object does not contain rankings.")
    }
  }
)

##### Access the maxRank:

#' @rdname rankingRcisTarget-class
#' @aliases getMaxRank,rankingRcisTarget-method
#' @export
setGeneric(name="getMaxRank",
           def=function(x) standardGeneric("getMaxRank"))
setMethod("getMaxRank",
          signature="rankingRcisTarget",
          definition = function(x) {
            # Previous versions didn't have @maxRank -> Contain all ranks
            ret <- tryCatch(
                  x@maxRank,
                  error=function(e) {Inf}
              )

            if(length(ret)==0)
              ret <- Inf

            return(ret)
          }
)
