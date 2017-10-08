#' @title Class to store the motif databases for RcisTarget.
#' @aliases rankingWrapper-class, rankingWrapper, getRanking,
#' ncol, nrow, show, subset
#' @param x [several methods] rankingWrapper object to apply the method
#' @param object [method: show] rankingWrapper object to show
#' @param elements [method: subset]
#' @param select [method: subset]
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
## @example inst/examples/example_class_rankingWrapper.R
#' @return
#' \itemize{
#' \item show: Prints a summary of the object
#' \item getRanking: Returns the data.frame containing the rankings
#' \item ncol, nrow, subset: Returns the number of columns, rows
#' or a subset of the ranking
#' }
#' @examples
#' library(RcisTarget.hg19.motifDatabases.cisbpOnly.500bp)
#' data("hg19_500bpUpstream_motifRanking_cispbOnly")
#' hg19_500bpUpstream_motifRanking_cispbOnly
#' class(hg19_500bpUpstream_motifRanking_cispbOnly)
#' @rdname rankingWrapper-class
#' @export rankingWrapper
#' @exportClass rankingWrapper
rankingWrapper <- setClass(
  # Set the name for the class
  Class="rankingWrapper",

  # Define the slots
  slots = c(
    rankings = "data.table",
    colType = "character", # motif or whatever feature is stored (e.g. ChipSeq)
    rowType = "character", # gene/region
    org = "character", # human/mouse
    genome = "character", # hg19, mm9 ...
    maxRank = "numeric", # Higher ranks are converted to Inf
    description = "character" # Other info, shown with "show"
  )
)

#' @rdname rankingWrapper-class
#' @aliases show,rankingWrapper-method
#' @export
setMethod("show",
  signature="rankingWrapper",
  definition = function(object) {
    message <- paste("Rankings for RcisTarget.", "\n",
    "  Organism: ", object@org, " (", object@genome,")","\n",
    "  Number of ", toupper(object@rowType), "S: ", nrow(object@rankings),"\n",
    "  Number of ", object@colType, "s: ", ncol(object@rankings)-1,"\n\n",
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

#' @rdname rankingWrapper-class
#' @aliases getRanking,rankingWrapper-method
#' @export
setGeneric(name="getRanking",
           def=function(x) standardGeneric("getRanking"))
setMethod("getRanking",
  signature="rankingWrapper",
  definition = function(x) {
    x@rankings
  }
)

##### Access the maxRank:

#' @rdname rankingWrapper-class
#' @aliases getMaxRank,rankingWrapper-method
#' @export
setGeneric(name="getMaxRank",
           def=function(x) standardGeneric("getMaxRank"))
setMethod("getMaxRank",
          signature="rankingWrapper",
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

##### Subset the object:
# #' @export
# setMethod('[', signature(x="rankingWrapper"),
#           definition=function(x, i, j, drop=FALSE){
#             x@ranking <- x@ranking[i,j, drop=drop]
#             x
#           })


#' @rdname rankingWrapper-class
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
          x@rankings <- x@rankings[x@rankings$rn %in% elements]
        }

    }}
    x
  }
)


# Regular "matrix" methods (rownames, []...) not behave as a data.table.
# Get from the @ranking slot manually...

#' @rdname rankingWrapper-class
#' @aliases nrow,rankingWrapper-method
#' @export
setMethod("nrow",
  signature="rankingWrapper",
  definition = function(x) {
    nrow(x@rankings)
  }
)

#' @rdname rankingWrapper-class
#' @aliases ncol,rankingWrapper-method
#' @export
setMethod("ncol",
  signature="rankingWrapper",
  definition = function(x) {
    ncol(x@rankings)
  }
)
