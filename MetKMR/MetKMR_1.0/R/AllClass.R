#' Main MetRKAT class
#'
#' @slot data data.frame.
#' @slot wsize numeric.
#' @slot wmethod character.
#' @slot gap numeric.
#' @slot intervals data.frame.
#' @slot results data.frame.
#' @slot annotation data.frame.
#' @slot distmethod character.
#' @slot max.na numeric.
#'
#' @return
#' @export
#'
#' @examples
setClass("MetRKAT",
         slots = c(
             data = "ANY",
             wsize = "numeric",
             wmethod = "character",
             gap = "numeric",
             intervals = "ANY",
             results = "ANY",
             annotation = "ANY",
             distmethod = "character",
             max.na = "numeric",
             dbsrc = "ANY"
         ),
         prototype = list(
            gap = 0,
            wmethod = "default",
            distmethod = "euclidean",
            max.na = 0.5,
            dbsrc = NULL,
            annotation = NULL,
            intervals = NULL
         ),

         validity = function(object) {
             if (! nrow(object@annotation)) {
                 # Check that provided sites on data are inside the annotation.
                 annot_loc <- IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Locations

                 if (! all(rownames(object@data) %in% rownames(annot_loc)))
                     message("Not all sites of the data provided are present on the annotation data frame. They will be discarded.")
             }

             if (object@wsize < 0)
                 stop("Invalid window size: must be greater than zero.")

             if (object@gap > object@wsize)
                 stop("The gap value cannot be greater than window size.")

             if (! all(object@distmethod %in% c("euclidean", "manhattan")))
                 stop("The distance matrix method must be 'euclidean' or 'manhattan'.")
         })
