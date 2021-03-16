#' @title Import the motif databases for RcisTarget.
#' @param dbFile .feather or .parquet file containing the rankings
#' @param columns Columns to load from the .feather or .parquet file
#' (e.g. to read only specific genes or regions)
#' @param dbDescr Description fields (not used internally) e.g.:
#' \code{dbDescr=list(colType="gene", rowType="motif",
#' org="Human", genome="hg19", maxRank=Inf, description="")}
#' @param indexCol Column name containing the feature IDs (e.g. motif names or chip-seq tracks). 
#' @param warnMissingColumns If 'columns' is provided, warn if any ID is not available in the rankings? 
#' @description
#' The rankings are typically loaded from a .feather or .parquet file
#' with \code{importRankings()}.
#' @return
#' rankingRcisTarget object with the following slots:
#' #' \itemize{
#' \item rankings: data.frame containing the rankings
#' \item colType: 'gene'or 'region'
#' \item nColsInDB: Number of columns (e.g. genes/regions) available
#' in the database (.feather or .parquet file).
#' Note that not all might be loaded in the current object.
#' \item rowType: 'motif' or the type of feature is stored (e.g. ChipSeq)
#' \item org: human/mouse/fly
#' \item genome: hg19, mm9, ...
#' \item description: global description, summary, or any other information
#' \item maxRank: Maximum ranking included in the database,
#' higher values are converted to Inf.
#' }

#' @examples
#' ## Loading from a .feather or .parquet file:
#' #motifRankings<-importRankings("hg19-500bp-upstream-7species.mc9nr.feather")
#' #motifRankings<-importRankings("hg19-500bp-upstream-7species.mc9nr.parquet")
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
#' @import arrow
#' @importFrom utils read.table
#' @export
importRankings <- function(dbFile, columns=NULL, dbDescr=NULL, indexCol="features", warnMissingColumns=TRUE)
{
  dbFile <- path.expand(dbFile)
  if(!file.exists(dbFile)) stop("File does not exist: ", dbFile)

  if(!is.null(columns)){
    missingColumns <- columns[which(!columns %in% getColumnNames(dbFile))]
    if(length(missingColumns)>0 & warnMissingColumns)
    {
      warning("The following columns are missing from the database: ", paste(missingColumns, collapse=", "))
      columns <- columns[which(columns %in% getColumnNames(dbFile))]
    }
    columns <- unique(c(indexCol, columns))
  }
  extension <- strsplit(dbFile, "\\.") [[1]][length(strsplit(dbFile, "\\.") [[1]])]
  if (extension == 'feather'){
    rnks <- arrow::read_feather(dbFile, col_select=!!columns, mmap = TRUE); dim(rnks)
    nColsInDB <- length(names(arrow::FeatherReader$create(arrow::ReadableFile$create(dbFile))))-1
  }else if (extension == "parquet"){
    rnks <- read_parquet(dbFile, columns = !!columns)
    pq <- arrow::ParquetFileReader(dbFile)
    nColsInDB <- pq$GetSchema()$num_fields()-1
  }else{
    stop("Database format must be feather or parquet.")
  }

  dbFile_descr <- gsub(paste0(".", extension),".descr", dbFile, fixed=TRUE)
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

#' @rdname importRankings
#' @export
getRowNames <- function(dbFile)
{
  dbPath <- dbFile
  extension <- strsplit(dbPath, "\\.") [[1]][length(strsplit(dbPath, "\\.") [[1]])]
  if (extension == 'feather'){
    ret <- unname(unlist(arrow::read_feather(path.expand(dbPath), col_select=1, mmap = TRUE)))
  }
  else if (extension == "parquet"){
     stop("Not implemented") # TODO: add arrow
  }
  return(ret)
}

#' @rdname importRankings
#' @export
getColumnNames <- function(dbFile) # TODO: Check if they are really genes/regions
{
  dbPath <- dbFile
  extension <- strsplit(dbPath, "\\.") [[1]][length(strsplit(dbPath, "\\.") [[1]])]
  if (extension == 'feather'){
    ret <- names(arrow::FeatherReader$create(arrow::ReadableFile$create(dbPath)))[-1]
  }
  else if (extension == "parquet"){
    # Not tested!! (if it works, the if can be removed)
    ret <- names(arrow::FeatherReader$create(arrow::ReadableFile$create(dbPath)))[-1]
    # stop("Not implemented") # TODO: add arrow
  }
  return(ret)
}
# 
# getGeneNames <- getRegionNames <- getColumnNames 
