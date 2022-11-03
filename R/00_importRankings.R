#' @title Import the motif databases for RcisTarget.
#' @param dbFile .feather or .parquet file containing the rankings
#' @param indexCol Column name containing the feature IDs (e.g. motif names or chip-seq tracks). If NULL, it will try to use 'motifs', 'tracks', or 'features'.
#' @param colType Colum type. i.e.:'gene'or 'region'
#' @param columns Columns to load from the .feather or .parquet file
#' (e.g. to read only specific genes or regions)
#' @param warnMissingColumns If 'columns' is provided, warn if any ID is not available in the rankings? 
#' @param dbDescr Description fields (not used in the analysis, just for convenience of the user) 
#' e.g.:
#' \code{dbDescr=list(colType="gene", rowType="motif",
#' org="Human", genome="hg19", maxRank=Inf, description="")}
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
#' # dbFilePath = "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
#' # motifRankings<-importRankings(dbFilePath)
#'
#' ## The annotations for Motif collection 9 (sufix 'mc9nr')
#' # are included in RcisTarget, and can be loaded with:
#' data(motifAnnotations_hgnc)
#'
#' ## For other versions, import the appropiate annotation. e.g.:
#' # annotDb <- importAnnotations("motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl")
#' # optional: motifsInRanking <- getRanking(motifRankings)$motifs

##### Load/import the ranking from a feather file:
#' @rdname importRankings
#' @importFrom arrow read_feather read_parquet FeatherReader ParquetFileReader ReadableFile
#' @import dplyr
#' @importFrom utils read.table
#' @export
importRankings <- function(dbFile, indexCol=NULL, colType="gene", columns=NULL, warnMissingColumns=TRUE, dbDescr=NULL)
{
  dbFile <- path.expand(dbFile)
  if(!file.exists(dbFile)) stop("File does not exist: ", dbFile)
  
  allColumns <- getColumnNames(dbFile)
  indexCol <- .getIndexCol(allColumns, indexCol=indexCol, verbose=warnMissingColumns)
  
  ## If columns are subset: add the index column
  if(!is.null(columns)){
    missingColumns <- columns[which(!columns %in% allColumns)]
    if(length(missingColumns)>0 & warnMissingColumns)
    {
      warning("The following columns are missing from the database: ", paste(missingColumns, collapse=", "))
      columns <- columns[which(columns %in% allColumns)]
    }
    columns <- unique(c(indexCol, columns))
  }
  
  ## Read file (format chosen based on extension)
  extension <- strsplit(dbFile, "\\.") [[1]][length(strsplit(dbFile, "\\.") [[1]])]
  if (extension == 'feather'){
    rnks <- arrow::read_feather(dbFile, col_select=!!columns, mmap = TRUE); dim(rnks)
    nColsInDB <- length(allColumns)-1  #names(arrow::FeatherReader$create(arrow::ReadableFile$create(dbFile)))
  }else if (extension == "parquet"){
    rnks <- arrow::read_parquet(dbFile, columns = !!columns)
    nColsInDB <- length(allColumns)-1 # pq <- arrow::ParquetFileReader(dbFile); pq$GetSchema()$num_fields()-1
  }else{
    stop("Database format must be feather or parquet.")
  }
  
  # The index column should be the first one
  if(!is.null(indexCol)) 
  {
    rnks <- rnks %>% relocate(!!indexCol)
  }
  
  # Add description if available
  if(!is.null(dbDescr))
  {
    dbDescr <- as.matrix(dbDescr)
  } else {
    # If not provided: keep empty
    dbDescr <- as.matrix(list(colType=colType,
                              rowType=indexCol,
                              org="",
                              genome="",
                              nColsAvailable=nColsInDB,
                              maxRank = Inf,
                              description=""))
  }
  
  # Create object
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
getRowNames <- function(dbFile, indexCol=NULL)
{
  allColumns <- getColumnNames(path.expand(dbFile))
  indexCol <- .getIndexCol(allColumns, indexCol=indexCol, verbose=FALSE)
  
  dbPath <- dbFile
  extension <- strsplit(dbPath, "\\.") [[1]][length(strsplit(dbPath, "\\.") [[1]])]
  if (extension == 'feather'){
    ret <- unname(unlist(arrow::read_feather(path.expand(dbPath), col_select=indexCol, mmap = TRUE)))
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
  dbPath <- path.expand(dbFile)
  extension <- strsplit(dbPath, "\\.") [[1]][length(strsplit(dbPath, "\\.") [[1]])]
  if (extension == 'feather'){
    ret <- arrow::ReadableFile$create(dbPath)
    ret <- arrow::FeatherReader$create(ret)
    ret <- names(ret) #[-1]
  }
  else if (extension == "parquet"){
    # ret <- names(arrow::ParquetFileReader$create(arrow::ReadableFile$create(dbPath)))[-1]  #maybe?
    stop("Not implemented") # TODO: add arrow
  }
  return(ret)
}

.getIndexCol <-  function(allColumns, indexCol=NULL, verbose=TRUE)
{
  if(is.null(indexCol)) {
    indexCol <- intersect(allColumns, c('motifs', 'tracks', 'features'))#  [1]
    if(verbose) message("Using the column '", indexCol, "' as feature index for the ranking database.")
  }else{
    if(!indexCol %in% allColumns) stop(paste0("The index column '", indexCol,"' is not available in the file."))
  }
  
  return(indexCol)
}
