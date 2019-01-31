#' @title Import the motif databases for RcisTarget.
#' @param dbFile .feather file containing the rankings
#' @param columns Columns to import from the .feather file 
#' (e.g. only selected genes or regions)
#' @param dbDescr The description fields are typically imported from the  
#' ".descr" file with the same filename as the \code{dbFile}. 
#' Otherwise they can be manually set e.g.: 
#' 
#' \code{dbDescr=list(colType="gene", rowType="motif", 
#' org="Human", genome="hg19", maxRank=Inf, description="")}
#' 
#' @description
#' The rankings are typically loaded from a .feather file 
#' with \code{importRankings()}.
#' 
#' If the associated .descr file is available, 
#' it will also load the description of the database.
#' @return
#' rankingRcisTarget object with the following slots:
#' #' \itemize{
#' \item rankings: data.frame containing the rankings
#' \item colType: 'gene'or 'region'
#' \item nColsInDB: Number of columns (e.g. genes/regions) available 
#' in the database (.feather file). 
#' Note that not all might be loaded in the current object.
#' \item rowType: 'motif' or the type of feature is stored (e.g. ChipSeq)
#' \item org: human/mouse/fly
#' \item genome: hg19, mm9, ...
#' \item description: global description, summary, or any other information
#' \item maxRank: Maximum ranking included in the database, 
#' higher values are converted to Inf.
#' }

#' @examples
#' ## Loading from a .feather file (the .descr file is read automatically):
#' #motifRankings<-importRankings("hg19-500bp-upstream-7species.mc9nr.feather")
#' 
#' ## The annotations for Motif collection 9 (sufix 'mc9nr') 
#' # are already included in RcisTarget, and can be loaded with:
#' data(motifAnnotations_hgnc)
#' 
#' ## For other versions, import the appropiate annotation. e.g.:
#' # annotDb <- importAnnotations("motifs-v9-nr.hgnc-m0.001-o0.0.tbl")
#' # optional: motifsInRanking <- getRanking(motifRankings)$features

##### Load/import the ranking from a feather file:
#' @rdname importRankings
#' @import feather
#' @import utils
#' @export
importRankings <- function(dbFile, columns=NULL, dbDescr=NULL)
{
  dbFile <- path.expand(dbFile)
  if(!file.exists(dbFile)) stop("File does not exist: ", dbFile)
  
  if(!is.null(columns)) columns <- unique(c("features", columns))
  rnks <- feather::read_feather(dbFile, columns=columns) # tibble
  #rnks <- data.frame... #to avoid replacing dash in names: check.names=FALSE
  nColsInDB <- feather::feather_metadata(dbFile)[["dim"]][2]
  
  dbFile_descr <- gsub(".feather",".descr", dbFile, fixed=TRUE)
  if(!is.null(dbDescr))
  {
    dbDescr <- as.matrix(dbDescr)
    if(file.exists(dbFile_descr)) 
      warning("Ignoring the DB file description (.descr)")
  } else {
    if(file.exists(dbFile_descr))
    {
      dbDescr <- utils::read.table(file=dbFile_descr,
                            sep = "\t", row.names=1, stringsAsFactors=FALSE)
      message("Imported description file:\n", 
              paste("\t", unname(sapply(rownames(dbDescr), 
            function(x) paste(x, dbDescr[x,1], sep=": "))), collapse="\n"))
    }else{
      # If not provided: keep empty
      dbDescr <- as.matrix(list(colType="column", 
                                rowType="row", 
                                org="", 
                                genome="", 
                                nColsAvailable=nColsInDB,
                                maxRank = Inf, 
                                description=""))
    }
  }
  
  dbDescr["nColsAvailable",] <- nColsInDB
  dbDescr["description",] <- paste0(dbDescr["description",], 
                                    " [Source file: ", basename(dbFile),"]")
  
  rownames(dbDescr) <- tolower(rownames(dbDescr))
  new("rankingRcisTarget",
      rankings=rnks,
      colType=as.character(dbDescr["coltype",]),
      rowType=as.character(dbDescr["rowtype",]),
      org=as.character(dbDescr["org",]),
      genome=as.character(dbDescr["genome",]),
      nColsInDB=as.numeric(dbDescr["ncolsavailable",]),
      maxRank = as.numeric(dbDescr["maxrank",]), 
      description=as.character(dbDescr["description",]))
}

