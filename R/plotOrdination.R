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
        args <- .check_ordination_input(x, dimred, colour.by = colour.by, ...)
        df <- do.call(.get_ordination_data, args)
        return(df)
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
    if( length(reducedDims(x)) == 0L ){
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
    if( !(.is_integer(ncomponents) && length(ncomponents) == 2L && all(ncomponents > 0L & ncomponents <= ncol(reducedDim(x, dimred)))) ){
        stop("'ncomponents' must specify columns from reducedDim(x, dimred) with integer values.", call. = FALSE)
    }

    # Check aesthetic variables
    temp <- .check_metadata_variable(tse, colour.by, FALSE, TRUE, FALSE, TRUE)
    temp <- .check_metadata_variable(tse, shape.by, FALSE, TRUE, FALSE, FALSE)
    temp <- .check_metadata_variable(tse, size.by, FALSE, TRUE, FALSE, FALSE)
    temp <- .check_metadata_variable(tse, group.by, FALSE, TRUE, FALSE, FALSE)
    temp <- .check_metadata_variable(tse, pair.by, FALSE, TRUE, FALSE, FALSE)
    temp <- .check_metadata_variable(tse, order.by, FALSE, TRUE, FALSE, FALSE)
    # If colour.by specifies rowname, we check assay.type as the abundance
    # values are used for coloring
    if( !is.null(colour.by) && colour.by %in% rownames(x) ){
        temp <- .check_assay_present(assay.type, x)
    }
    # Check other flags

    # Put all the arguments into list that be fed to the data retrieval function
    args <- list(
        x = x, dimred = dimred,
        ncomponents = ncomponents,
        colour.by = colour.by,
        shape.by = shape.by,
        size.by = size.by,
        group.by = group.by,
        pair.by = pair.by,
        order.by = order.by,
        assay.type = "counts"
    )
    args <- c(args, list(...))
    return(args)
}

.get_ordination_data <- function(x, dimred, ncomponents, colour.by, shape.by, size.by, group.by, pair.by, order.by, assay.type, ...){
    df <- reducedDim(x, dimred)
    if( is.null(colnames(df)) ){
        colnames(df) <- paste0(dimred, seq_len(ncol(df)))
    }
    df <- df |> as.data.frame()
    df <- df[, ncomponents]
    x_var <- colnames(df)[[1L]]
    y_var <- colnames(df)[[2L]]

    cols <- c(shape.by, size.by, group.by, pair.by, order.by)
    if( !is.null(colour.by) && colour.by %in% colnames(colData(x)) ){
        cols <- c(cols, colour.by)
    } else if( !is.null(colour.by) ){
        df[[colour.by]] <- assay(x, assay.type)[colour.by, ]
    }
    cd <- colData(x)[, cols, drop = FALSE]
    df <- cbind(df, cd)

    attributes(df) <- c(
        attributes(df),
        x = x_var,
        y = y_var,
        colour.by = colour.by,
        shape.by = shape.by,
        size.by = size.by,
        group.by = group.by,
        pair.by = pair.by,
        order.by = order.by,
        assay.type = assay.type
    )
    return(df)
}

.ordination_plotter <- function(df, ...){

}
