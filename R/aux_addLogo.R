# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)

#' @title Add motif logo to RcisTarget results table
#' @description Adds a column containing the logo URL to RcisTarget results
#' table. 
#' @param motifEnrDT Results from RcisTarget (data.table)
#' @param addHTML Whether to add the HTML tag <img> around the URL or not
#' (boolean).
#' @param dbVersion For current databases (mc9nr) use "v9"
#' @param motifCol Name of the column which contains the logo ID.
#' @return Returns the results table with a new column: 'logo'.
#' This column contains either a URL with the logo image, or the HTML code to
#' show the logo [e.g. with datatable()].
#' @seealso \code{vignette("showLogo")} to directly show the table as HTML. 
#' See the package vignette for more examples:
#' \code{vignette("RcisTarget")}
#' @example inst/examples/example_addLogo.R
#' @export
addLogo <- function(motifEnrDT, addHTML=TRUE, dbVersion="v9", motifCol="motif")
{
  isNA <- which(motifEnrDT[[motifCol]]=="")
  logos <- paste("http://motifcollections.aertslab.org/",
                dbVersion,"/logos/",
                motifEnrDT[[motifCol]],".png", sep="")
  
  if(addHTML)
  {
    logos <- paste('<img src="', logos,
                '" height="52" alt="',
                motifEnrDT[[motifCol]], '"></img>', sep="")
  }
  logos[isNA] <- ""
  
  data.table::data.table(logo=logos, motifEnrDT)
}

#' @title Show RcisTarget as HTML
#' @description Shows the results motif enrichment table as an HTML 
#' including the motif logos. Note that Transfac-Pro logos cannot be shown.
#' @param motifEnrDT Results from RcisTarget (data.table)
#' @param motifCol Name of the column which contains the logo ID.
#' @param dbVersion For current databases (mc9nr) use "v9"
#' @param nSignif Number of digits to show in numeric columns.
#' @param colsToShow Columns to show in the HTML 
#' (by default, the list of the enriched genes is hidden)
#' @param options List of options to pass to \code{DT::databable}
#' @param ... Other arguments to pass to \code{DT::databable}
#' @return Returns the DT::datatable() object which can be shown as HTML.
#' @seealso See the package vignette for more examples:
#' \code{vignette("RcisTarget")}
#' @example inst/examples/example_showLogo.R
#' @export
showLogo <- function(motifEnrDT, 
                       motifCol=c("motif", "bestMotif", "MotifID"), dbVersion="v9",
                       nSignif=3,
                       colsToShow=c(motifEnrichment=c("motifDb", "logo", "NES", "geneSet", "TF_highConf"),
                                    regulonTargets=c("TF", "gene", "nMotifs", "bestMotif", "logo", "NES", "highConfAnnot", "Genie3Weight")),
                       options=list(pageLength=50), 
                       ...)
{
  met <- motifEnrDT
  # Add motif logo
  if(!is.null(motifCol)){
    motifCol <- motifCol[which(motifCol %in% colnames(met))]
    if(length(motifCol) == 1){
      met <- RcisTarget::addLogo(met, motifCol=motifCol, dbVersion=dbVersion, addHTML=TRUE)
      if(!is.null(colsToShow)) colsToShow <- c("logo", colsToShow)
    }else{
      stop("Please indicate the column containing the motif id (argument 'motifCol') or set it to NULL.")
    }
    # TODO 
    met <- met[grep("transfac_pro__", met[[motifCol]], invert = T),]
  }
  
  # For numeric columns, show only the number of significant digits...
  for(i in which(sapply(met, is.numeric))){
    met[[i]] <- signif(met[[i]], nSignif)
  }
  
  # Keep only requested columns
  if(!is.null(colsToShow)) {
    colsToShow <- unique(unname(unlist(colsToShow)))
    colsToShow <- colsToShow[which(colsToShow %in% colnames(met))]
    met <- met[, colsToShow, with=F]
  }
  
  # Show...
  DT::datatable(met, 
                escape=FALSE, filter="top", options=options, ...)
}
