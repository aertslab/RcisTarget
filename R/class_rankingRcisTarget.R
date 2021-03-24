#' @title Class to store the motif databases for RcisTarget.
#' @aliases rankingRcisTarget-class
#' @param x [several methods] rankingRcisTarget object to apply the method
#' @param object [method: show] rankingRcisTarget object to show
#' 
#' @description
#' This class contains the rankings used by RcisTarget. 
#' 
#' They are typically loaded from a .feather file with \code{importRankings()}.
#' 
#' If the associated .descr file is available, 
#' it will also load the description of the database.
#' 
#' Class slots:
#' 
#' \itemize{
#' \item rankings: data.frame (tbl_df) containing the rankings
#' \item colType: 'gene'or 'region'
#' \item rowType: 'motif' or the type of feature is stored (e.g. ChipSeq)
#' \item org: human/mouse/fly
#' \item genome: hg19, mm9, ...
#' \item nColsInDB: Number of columns (e.g. genes/regions) available 
#' in the database (.feather file). 
#' Note that not all might be loaded in the current object.
#' \item description: global description, summary, or any other information
#' \item maxRank: Maximum ranking included in the database, 
#' higher values are converted to Inf.
#' }
#' 
#' Note that the main slot is \code{@rankings}, 
#' which is the one used by RcisTarget 
#' (it can accessed with \code{getRanking()}).
#' The 'description' slots are mostly for user convenience. 
#' 
## @example inst/examples/example_class_rcos.R
#' @return
#' \itemize{
#' \item show: Prints a summary of the object
#' \item getRanking: Returns the rankings
#' \item ncol, nrow: Returns the number of columns or rows of the ranking
#' }

#' @examples
#' ## Loading from a .feather file:
#' # dbFile <- "hg19-500bp-upstream-7species.mc9nr.feather"
#' # motifRankings <- importRankings(dbFile)
#' # motifRankings
#' 
#' ## Loading a built object:
#' library(RcisTarget.hg19.motifDBs.cisbpOnly.500bp)
#' data("hg19_500bpUpstream_motifRanking_cispbOnly")
#' hg19_500bpUpstream_motifRanking_cispbOnly
#' class(hg19_500bpUpstream_motifRanking_cispbOnly)
#' @importFrom methods slotNames
#' @importClassesFrom tibble tbl_df
#' @rdname rankingRcisTarget-class
#' @export rankingRcisTarget
#' @exportClass rankingRcisTarget
rankingRcisTarget <- setClass(
  # Set the name for the class
  Class="rankingRcisTarget",
  
  slots = c(
    rankings = "tbl_df", # "data.frame",
    colType = "character", # gene/region
    rowType = "character", # motif or whatever feature is stored (e.g. ChipSeq)
    org = "character", # human/mouse
    genome = "character", # hg19, mm9 ...
    nColsInDB="numeric", # Number of columns (e.g. genes) available in DB
    maxRank = "numeric", # Higher ranks are converted to Inf
    description = "character", # Other info, shown with "show"
    rcistarget_version = "character" # Class/RcisTarget version
  )
)

#' @rdname rankingRcisTarget-class
#' @aliases show,rankingRcisTarget-method
#' @export
setMethod("show",
  signature="rankingRcisTarget",
  definition = function(object) {
    message <- paste("Rankings for RcisTarget.\n")
    
    if((object@org != "") || (object@genome != ""))
    {
      genome <- ""
      if((object@genome != ""))
        genome <- paste0(" (", object@genome,")")
      message <- paste0(message,
      "  Organism: ", object@org, genome,"\n")
    }
    
    nGenes <- ncol(getRanking(object))-1
    if(("nColsInDB" %in% slotNames(object)) && !is.na(object@nColsInDB)) nGenes <- paste0(nGenes," (",object@nColsInDB," available in the full DB)")
    message <- paste0(message,
     "  Number of ",object@colType, "s: ", nGenes, "\n",
     "  Number of ",toupper(object@rowType),"S: ",nrow(getRanking(object)),"\n"
     )
    
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
#' @importFrom AUCell getRanking
#' @name getRanking
#' @rdname rankingRcisTarget-class
#' @aliases getRanking,rankingRcisTarget-method
#' @exportMethod getRanking
setMethod("getRanking",
          signature="rankingRcisTarget",
          definition = function(object) {
            object@rankings
          }
)

#' @importFrom BiocGenerics nrow
#' @rdname rankingRcisTarget-class
#' @aliases nrow,rankingRcisTarget-method
#' @export
setMethod("nrow",
          signature="rankingRcisTarget",
          definition = function(x) {
            nrow(getRanking(x))
          }
)

#' @importFrom BiocGenerics ncol
#' @rdname rankingRcisTarget-class
#' @aliases ncol,rankingRcisTarget-method
#' @export
setMethod("ncol",
          signature="rankingRcisTarget",
          definition = function(x) {
            ncol(getRanking(x))
          }
)


#' @importFrom BiocGenerics colnames
#' @rdname rankingRcisTarget-class
#' @aliases colnames,rankingRcisTarget-method
#' @export
setMethod("colnames",
          signature="rankingRcisTarget",
          definition = function(x) {
            colnames(getRanking(x)[-1])
          }
)

#' @importFrom BiocGenerics rownames
#' @rdname rankingRcisTarget-class
#' @aliases rownames,rankingRcisTarget-method
#' @export
setMethod("rownames",
          signature="rankingRcisTarget",
          definition = function(x) {
            rownames(getRanking(x))
          }
)

##### Access the maxRank:
setGeneric(name="getMaxRank",
           def=function(x) standardGeneric("getMaxRank"))
#' @name getMaxRank
#' @rdname rankingRcisTarget-class
#' @aliases getMaxRank,rankingRcisTarget-method
#' @export
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

##### Access the nColsAvailable:
setGeneric(name="getNumColsInDB",
           def=function(x) standardGeneric("getNumColsInDB"))
#' @name getNumColsInDB
#' @rdname rankingRcisTarget-class
#' @aliases getNumColsInDB,rankingRcisTarget-method
#' @export
setMethod("getNumColsInDB",
          signature="rankingRcisTarget",
          definition = function(x) {
            # Previous versions didn't have @maxRank -> Contain all ranks
            ret <- tryCatch(
              x@nColsInDB,
              error=function(e) {NA}
            )
            
            if(length(ret)==0)
              ret <- NA
            
            return(ret)
          }
)
