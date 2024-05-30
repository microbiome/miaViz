#' Wrapper for scater::plotReducedDim()
#' 
#' @param x a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object.
#'   
#' @name plotNMDS
NULL

#' @rdname plotNMDS
#' @export
plotNMDS <- function(x, ..., ncomponents = 2){
    plotReducedDim(x, ncomponents = ncomponents, dimred = "NMDS", ...)
}
