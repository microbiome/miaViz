#' @name
#' plotOrdination
#'
#' @title
#' Create ordination plot
#'
#' @description
#' Ordinaton plotter
#'
#' @details
#' Creates ordination plot
#'
#' @return
#' A \code{ggplot2} object.
#'
#' @param x a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#'
#' @param ... Additional parameters for plotting.
#' \itemize{
#'   \item \code{colour.by}: \code{NULL} or \code{character scalar}. Specifies a
#'   variable from \code{colData(x)} or \code{rowData(x)} which is used to
#'   colour observations. (Default: \code{NULL})
#' }
#'
#' @examples
#' data("Tito2024QMP")
#' tse <- Tito2024QMP
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[scater:plotReducedDim]{scater::plotReducedDim}}
#' }
#'
NULL

#' @rdname plotOrdination
#' @export
setMethod("plotOrdination", signature = c(x = "SummarizedExperiment"),
    function(x, ...){
        temp <- .check_ordination_input(x, ...)
        df <- .get_ordination_data(x, ...)
        p <- .ordination_plotter(df, ...)
        return(p)
    }
)

################################ HELP FUNCTIONS ################################

.check_ordination_input <- function(x, ...){

}

.get_ordination_data <- function(x, ...){

}

.ordination_plotter <- function(df, ...){

}
