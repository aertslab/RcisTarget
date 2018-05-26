# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)

#' @title Add motif logo to RcisTarget results table
#' @description Adds a column containing the logo URL to RcisTarget results
#' table. Note that Transfac-Pro logos cannot be shown.
#' @param motifEnrDT Results from RcisTarget (data.table)
#' @param addHTML Whether to add the HTML tag <img> around the URL or not
#' (boolean).
#' @param dbVersion For current databases (mc9nr) use "v9"
#' @param motifCol Name of the column which contains the logo ID.
#' @return Returns the results table with a new column: 'logo'.
#' This column contains either a URL with the logo image, or the HTML code to
#' show the logo [e.g. with datatable()].
#' @seealso See the package vignette for more examples:
#' \code{vignette("RcisTarget")}
#' @example inst/examples/example_addLogo.R
#' @export
addLogo <- function(motifEnrDT, addHTML=TRUE, dbVersion="v9", motifCol="motif")
{
  logos <- paste("http://motifcollections.aertslab.org/",
                dbVersion,"/logos/",
                motifEnrDT[[motifCol]],".png", sep="")
  if(addHTML)
    logos <- paste('<img src="', logos,
                '") height="52" alt="',
                motifEnrDT[[motifCol]], '"></img>', sep="")

  data.table::data.table(logo=logos, motifEnrDT)
}
