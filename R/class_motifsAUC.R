
#' @title Class to store AUC
#' @description Class with only one slot: @AUC, a matrix which contains the AUC for the gene-sets & motifs.
#' @example
#' fakeAUC <- motifsAUC(AUC=matrix(1/1:30, nrow=3)) # gene-sets by motifs
#' fakeAUC
#' dim(fakeAUC@AUC)
#' class(fakeAUC@AUC)
#' fakeAUC@AUC[,1:5]
#' @export

motifsAUC <- setClass(
  # Set the name for the class
  Class="motifsAUC",

  # Define the slots
  slots = c(
    AUC = "matrix"
  )
)

setMethod("show",
          signature="motifsAUC",
          definition = function(object) {
            message(paste("AUC for",
                  nrow(object@AUC), "gene-sets and", ncol(object@AUC), "motifs."))
            if(!is.null(rownames(object@AUC))) message(paste("First gene-sets:", paste(head(rownames(object@AUC), 3), collapse=", ")))
            if(!is.null(colnames(object@AUC))) message(paste("First motifs:", paste(head(colnames(object@AUC), 3), collapse=", ")))
          }
)
