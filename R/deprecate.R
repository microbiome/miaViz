#' These functions will be deprecated. Please use other functions instead.
#'
#' @param x -
#'
#' @param ... -
#'
#' @name deprecate
NULL

#' @rdname deprecate
#' @aliases plotRowPrevalence
#' @export
setGeneric("plotTaxaPrevalence", signature = c("x"),
    function(x, ...)
    standardGeneric("plotTaxaPrevalence"))

#' @rdname deprecate
#' @aliases plotRowPrevalence
#' @export
setMethod("plotTaxaPrevalence", signature = c(x = "ANY"), function(x, ...){
    .Deprecated(
        old ="plotTaxaPrevalence", new = "plotRowPrevalence",
        msg = paste0("The 'plotTaxaPrevalence' function is ",
            "deprecated. Use 'plotRowPrevalence' instead."))
    plotRowPrevalence(x, ...)
    }
)

#' @rdname deprecate
#' @aliases plotRowPrevalence
#' @export
setGeneric("plotFeaturePrevalence", signature = c("x"),
    function(x, ...)
    standardGeneric("plotFeaturePrevalence"))

#' @rdname deprecate
#' @aliases plotRowPrevalence
#' @export
setMethod("plotFeaturePrevalence", signature = c(x = "ANY"), function(x, ...){
    .Deprecated(
        old ="plotFeaturePrevalence", new = "plotRowPrevalence",
        msg = paste0("The 'plotFeaturePrevalence' function is ",
            "deprecated. Use 'plotRowPrevalence' instead."))
    plotRowPrevalence(x, ...)
    }
)
