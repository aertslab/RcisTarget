
#' @title Class to store the motif databases for RcisTarget.
#' @description
#' This class is only meant for internal use. Modify at your own risk.
#'
## Methods: See examples.
#' @import BiocGenerics
## @example inst/examples/example_class_rankingWrapper.R
#' @export rankingWrapper
#' @exportClass rankingWrapper
rankingWrapper <- setClass(
  # Set the name for the class
  Class="rankingWrapper",

  # Define the slots
  slots = c(
    rankings = "data.table",
    colType = "character", # motif or whatever other feature is stored (e.g. ChipSeq)
    rowType = "character", # gene/region
    org = "character" # hg19, mm9 ...
  )
)

#' @export
setMethod("show",
          signature="rankingWrapper",
          definition = function(object) {
            message <- paste("Rankings for ", object@org, "containing", nrow(object@rankings)," ", object@rowType, "s (rows) and ", ncol(object@rankings)," ", object@colType, "s (columns).\n", sep="")
          }
)

#' @export
setMethod('[', signature(x="rankingWrapper"),
          definition=function(x, i, j, drop){
            x@rankings[i,j, drop=drop]
          })

#' @export
setMethod("subset",
          signature="rankingWrapper",
          definition = function(x, elements, select=x@rowType) {

            if(length(select) > 1) stop()

            if(grepl("cell", tolower(select))) {
              x@rankings <- x@rankings[,elements, drop=FALSE]
            }else{
              x@rankings <- x@rankings[elements,, drop=FALSE]
            }
            x
          }
)

#' @export
setMethod("dim",
          signature="rankingWrapper",
          definition = function(x) {
            dim(x@rankings)
          }
)

#' @export
setMethod("rownames",
          signature="rankingWrapper",
          definition = function(x) {
            rownames(x@rankings)
          }
)

#' @export
setMethod("colnames",
          signature="rankingWrapper",
          definition = function(x) {
            colnames(x@rankings)
          }
)

#' @export
setMethod("nrow",
          signature="rankingWrapper",
          definition = function(x) {
            nrow(x@rankings)
          }
)

#' @export
setMethod("ncol",
          signature="rankingWrapper",
          definition = function(x) {
            ncol(x@rankings)
          }
)


##### Access the rankings:
#' @export
setGeneric(name="getRanking", def=function(object) standardGeneric("getRanking"))
setMethod("getRanking",
          signature="rankingWrapper",
          definition = function(object) {
              object@rankings
          }
)

