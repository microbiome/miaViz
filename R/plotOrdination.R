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
setMethod("plotOrdination", signature = c(x = "SingleCellExperiment"),
    function(x, dimred, colour.by = NULL, ...){
        temp <- .check_ordination_input(x, ...)
        df <- .get_ordination_data(x, ...)
        p <- .ordination_plotter(df, ...)
        return(p)
    }
)

################################ HELP FUNCTIONS ################################

.check_ordination_input <- function(
        x, dimred,
        ncomponents = 2L,
        colour.by = color.by, color.by = NULL,
        shape.by = NULL,
        size.by = NULL,
        group.by = NULL,
        pair.by = NULL,
        order.by = NULL,
        assay.type = "counts",
        ...){
    # Check if there are any reduced dim present
    if( lenght(reducedDims(x)) == 0L ){
        stop("No data present in reducedDim(x).", call. = FALSE)
    }
    # Check that dimred can be found
    is_name <- .is_a_string(dimred) && dimred %in% reducedDimNames(x)
    is_index <- .is_an_integer(dimred) && dimred > 0L && dimred <= length(reducedDims(x))
    if( !(is_name || is_index) ){
        stop("'dimred' must specify data from reducedDim(x). It must be one ",
            "of the following options: '",
            paste0(reducedDimNames(x), collapse = "', '"), "'", call. = FALSE)
    }
    # Check that ncomponents is correct. We can only visualize 2 components.
    if( .is_an_integer(ncomponents) ){
        ncomponents <- seq_len(ncomponents)
    }
    if( !(.is_integer(ncomponents) && length(ncomponents) == 2L) ){
        stop("'ncomponents' must specify intger values.", call. = FALSE)
    }

    # Check aesthetic variables
    temp <- .check_metadata_variable(tse, colour.by, TRUE, TRUE, FALSE, TRUE)
    temp <- .check_metadata_variable(tse, shape.by, TRUE, TRUE, FALSE, FALSE)
    temp <- .check_metadata_variable(tse, size.by, TRUE, TRUE, FALSE, FALSE)
    temp <- .check_metadata_variable(tse, group.by, TRUE, TRUE, FALSE, FALSE)
    temp <- .check_metadata_variable(tse, pair.by, TRUE, TRUE, FALSE, FALSE)
    temp <- .check_metadata_variable(tse, order.by, TRUE, TRUE, FALSE, FALSE)
    # If colour.by specifies rowname, we check assay.type as the abundance
    # values are used for coloring
    if( colour.by &in& rownames(x) ){
        temp <- .check_assay_present(assay.type, x)
    }
    # Check other flags

}

.get_ordination_data <- function(x, dimred, ...){
    df <- reducedDim(x, dimred) |> as.data.frame()
    df <- df[, ncomponents]
}

.ordination_plotter <- function(df, ...){

}
