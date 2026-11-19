#' These functions will be deprecated. Please use other functions instead.
#'
#' @param x -
#'
#' @param centering -
#'
#' @param ... -
#'
#' @name deprecate
NULL

#' @rdname deprecate
#' @export
setGeneric("plotTaxaPrevalence", signature = c("x"),
    function(x, ...)
    standardGeneric("plotTaxaPrevalence"))

#' @rdname deprecate
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
#' @export
setGeneric("plotFeaturePrevalence", signature = c("x"),
    function(x, ...)
    standardGeneric("plotFeaturePrevalence"))

#' @rdname deprecate
#' @export
setMethod("plotFeaturePrevalence", signature = c(x = "ANY"), function(x, ...){
    .Deprecated(
        old ="plotFeaturePrevalence", new = "plotRowPrevalence",
        msg = paste0("The 'plotFeaturePrevalence' function is ",
            "deprecated. Use 'plotRowPrevalence' instead."))
    plotRowPrevalence(x, ...)
    }
)

#' @rdname deprecate
#' @export
setGeneric("getNeatOrder", signature = c("x"),
    function(x, centering = "mean", ...)
    standardGeneric("getNeatOrder"))

#' @rdname deprecate
#' @export
setMethod("getNeatOrder", signature = c(x = "matrix"),
    function(x, centering = "mean", ...){
        .Deprecated(
            old = "getNeatOrder", new = "getHeatmapOrder",
            msg = paste0("The 'getNeatOrder' function is deprecated. Use ",
                "'getHeatmapOrder' with method = 'PCA_angle' instead."))
        # Check args
        .check_args(x, centering)
        # Get the theta values and order them
        theta_values <- .radial_theta(x, centering)
        ordering <- order(theta_values)
        return(ordering)
    }
)

# Checks the method arguments.
.check_args <- function(x, centering) {
    # Check data is a matrix
    if (!is.matrix(x)) {
        stop("Input data must be a matrix.", call. = FALSE)
    }
    # Check there is sufficient data
    if (nrow(x) == 0 || ncol(x) == 0) {
        stop(
            "No data to plot. Matrix must have at least one row and one ",
            "column.", call. = FALSE)
    }
    # Check there is sufficient data
    if (ncol(x) != 2) {
        stop("Matrix must have only 2 columns.", call. = FALSE)
    }
    # Check centering argument
    if ( !(is.null(centering) || (.is_a_string(centering) &&
            centering %in% c("mean", "median", NULL))) ){
        stop(
            "'centering' must be a single character value or NULL.",
            call. = FALSE)
    }
    return(NULL)
}

#' @importFrom stats median
# Computes the radial theta values for each row in the data matrix.
.radial_theta <- function(data, centering) {
    # Apply the centering if centering is specified
    if (!is.null(centering)) {
        # Choose the correct centering function based on the method
        center_fun <- switch(centering, "median" = median, "mean" = mean)
        center_vals <- apply(data, 2, center_fun)
        data <- scale(data, center = center_vals, scale = FALSE)
    }
    # Compute the radial theta values using the centered data
    theta <- atan2(data[, 2], data[, 1])
    # Set the names of theta values to the row names of the centered data and
    # return the theta values
    names(theta) <- rownames(data)
    return(theta)
}
