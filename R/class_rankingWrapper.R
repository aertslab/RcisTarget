#' @title Class to store the motif databases for RcisTarget.
#' @aliases getRanking
#' @param rankings data.table with the rankings
#' @param colType 'motif' or whatever other feature is stored (e.g. ChipSeq)
#' @param rowType 'gene'or 'region'
#' @param org human/mouse/fly
#' @param genome hg19, mm9, ...
#' @param description summary or any other information
#' @description
#' This class is only meant for internal use. Modify at your own risk.
#'
## Methods: See examples.
## @import BiocGenerics
## @example inst/examples/example_class_rankingWrapper.R
#' @rdname rankingWrapper-class
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
    org = "character", # human/mouse
    genome = "character", # hg19, mm9 ...
    description = "character" # Other info, shown with "show"
  )
)

#' @rdname rankingWrapper
#' @aliases show,rankingWrapper-method
#' @export
setMethod("show",
          signature="rankingWrapper",
          definition = function(object) {
            message <- paste("Rankings for RcisTarget.", "\n",
                            "  Organism: ", object@org, " (", object@genome,")","\n",
                            "  Number of ", toupper(object@rowType), "S: ", nrow(object@rankings),"\n",
                            "  Number of ", object@colType, "s: ", ncol(object@rankings)-1,"\n\n",
                            object@description, sep="")

            cat(message)
          }
)

##### Access the rankings:

#' @rdname rankingWrapper
#' @aliases getRanking,rankingWrapper-method
#' @export
setGeneric(name="getRanking", def=function(object) standardGeneric("getRanking"))
setMethod("getRanking",
          signature="rankingWrapper",
          definition = function(object) {
            object@rankings
          }
)

##### Subset the object:
# #' @export
# setMethod('[', signature(x="rankingWrapper"),
#           definition=function(x, i, j, drop=FALSE){
#             x@ranking <- x@ranking[i,j, drop=drop]
#             x
#           })


#' @rdname rankingWrapper
#' @aliases subset,rankingWrapper-method
#' @export
setMethod("subset",
          signature="rankingWrapper",
          definition = function(x, elements, select=x@rowType) {

            if(length(select) > 1) stop()

            if(grepl("col", tolower(select))) {

              if(is.numeric(elements))
              {
                x@rankings <- x@rankings[, unique(c(1, elements)), with=FALSE]
              }else{
                x@rankings <- x@rankings[, unique(c("rn", elements)), with=FALSE]
              }

            }else{
              if(grepl("row", tolower(select))) {
                if(is.numeric(elements))
                {
                  x@rankings <- x@rankings[elements,]
                }else{
                  x@rankings <- x@rankings[rn %in% elements]
                }

            }}
            x
          }
)


# Regular "matrix" methods (rownames, []...) not behave as a data.table.
# Get from the @ranking slot manually...

#' @rdname rankingWrapper
#' @aliases nrow,rankingWrapper-method
#' @export
setMethod("nrow",
          signature="rankingWrapper",
          definition = function(x) {
            nrow(x@rankings)
          }
)

#' @rdname rankingWrapper
#' @aliases ncol,rankingWrapper-method
#' @export
setMethod("ncol",
          signature="rankingWrapper",
          definition = function(x) {
            ncol(x@rankings)
          }
)
