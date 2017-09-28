#' @title Class to store the AUC scores for RcisTarget.
#' @aliases getAUC
#' @description
#'
#' Contains the AUC scores for each gene- or region-set.
#' They can be accessed through \code{getAUC()} and the regular methods
#' (i.e. nrow, rownames... ) available for SummarizedExperiment objects.
#'
#' @param object Results from \code{calcAUC}.
#' @return
#' \itemize{
#' \item show: Prints a summary of the object
#' \item getAUC: Returns the matrix containing the AUC scores
#' }
#' @method show aucScores
#' @method getAUC aucScores
#' @example inst/examples/example_aucScores_class.R
#' @import data.table
#' @import SummarizedExperiment
#' @importFrom utils head
# @importFrom R.utils capitalize
#' @rdname aucScores-class
#' @export aucScores
aucScores <- setClass("aucScores",
  contains="SummarizedExperiment",
  representation=representation(
    org = "character", # human/mouse
    genome = "character", # hg19, mm9 ...
    description = "character" # Other info, shown with "show"
  )
)

#' @rdname aucScores-class
#' @aliases show,aucScores-method
setMethod("show",
  signature="aucScores",
  definition = function(object) {
    message <- paste(R.utils::capitalize(assayNames(object)), " for ",
                 nrow(object)," ", names(dimnames(assay(object)))[1],
                 "s and ",
                 ncol(object)," ", names(dimnames(assay(object)))[2],
                 "s.\n", sep="")

    if(object@org!="")
      message <- paste(message, "  Organism: ", object@org,
                       " (", object@genome,")","\n  ",
                        object@description,"\n", sep="")

    message <- paste(message,
                       "\nAUC matrix preview:\n", sep="")
    cat(message)
    subsetToPrint <- head(assay(object)[,seq_len(min(5, ncol(object))),
                                        drop=FALSE])
    names(dimnames(subsetToPrint)) <- NULL
    show(subsetToPrint)

  }
)
##### Access the matrix:
#' @export
#' @rdname aucScores-class
#' @aliases getAUC,aucScores-method
setGeneric(name="getAUC",
           def=function(object) standardGeneric("getAUC"))
setMethod("getAUC",
  signature="aucScores",
  definition = function(object) {
    if("AUC" %in% assayNames(object)) {
      assays(object)[["AUC"]]
    }else{
      stop("This object does not contain an AUC matrix.")
    }
  }
)
