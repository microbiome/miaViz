#' Plot tree with feature loadings after performing a reducedDim (e.g. PCA, PCoA, NMDS)
#'
#' Based on the \code{\link[diffTop:plotASVcircular]{plotASVcircular}} method using phyloseq 
#' and has been converted to use TreeSummarizedExperiment objects.
#'
#' @param object a
#'   \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'   object.
#' 
#' @param nfeatures An integer to specify the number of feature loadings to plot (i.e. number of topics).
#'   (default: \code{nfeatures = 2})
#' 
#' @param taxonomyLevel 
#' # ADD ALL PARAMS
#' 
#' @details
#' # TO DO 
#' @return 
#' A \code{ggplot2} object. A circular plot annotated with TreeSummarizedExperiment object.
#'
#' @name plotFeatureLoadings
#'
#' @examples
#' # TO DO
NULL

#' @rdname plotFeatureLoadings
setGeneric("plotFeatureLoadings", signature = c("object"),
           function(object, ...)
             standardGeneric("plotFeatureLoadings"))


#' @rdname plotTree
#' @export
setMethod("plotFeatureLoadings", signature = c(object = "TreeSummarizedExperiment"),
          function(object,
                   ...){
            
}