#' Wrapper for scater::plotReducedDim()
#'
#' @param x a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#'
#' @param ncomponents \code{Numeric scalar}. indicating the number of dimensions
#' to plot, starting from the first dimension. Alternatively, a numeric vector
#' specifying the dimensions to be plotted. (Default: \code{2})
#'
#' @param ... additional arguments passed to scater::plotReducedDim().
#'
#' @name plotNMDS
NULL

#' @rdname plotNMDS
#' @export
plotNMDS <- function(x, ..., ncomponents = 2){
    plotReducedDim(x, ncomponents = ncomponents, dimred = "NMDS", ...)
}
