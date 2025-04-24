#' Sorting by radial theta angle
#'
#' @description \code{getNeatOrder} sorts already ordinated data by the radial
#' theta angle. This method is useful for organizing data points based on their
#' angular position in a 2D space, typically after an ordination technique such
#' as PCA or NMDS has been applied.
#'
#' The function takes in a matrix of ordinated data, optionally
#' centers the data using specified methods (\code{mean}, \code{median}, or
#' \code{NULL}), and then calculates the angle (theta) for each point relative
#' to the centroid. The data points are then sorted based on these theta values
#' in ascending order.
#'
#' One significant application of this sorting method is in plotting heatmaps.
#' By using radial theta sorting, the relationships between data points can be
#' preserved according to the ordination method's spatial configuration, rather
#' than relying on hierarchical clustering, which may distort these
#' relationships. This approach allows for a more faithful representation of the
#' data's intrinsic structure as captured by the ordination process.
#'
#' @param x A matrix containing the ordinated data to be sorted. Columns should
#' represent the principal components (PCs) and rows should represent the
#' entities being analyzed (e.g. features or samples). There should be 2 columns
#' only representing 2 PCs.
#'
#' @param centering \code{Character scalar}. Specifies the method to
#' center the data. Options are \code{"mean"}, \code{"median"}, or \code{NULL}
#' if your data is already centered. (Default: \code{"mean"})
#'
#' @param ... Additional arguments passed to other methods.
#'
#' @return A \code{character} vector of row indices in the sorted order.
#'
#' @details
#' It's important to note that the
#' [\pkg{sechm}](https://bioconductor.org/packages/3.18/bioc/vignettes/sechm/inst/doc/sechm.html#row-ordering)
#' package does actually have the functionality for plotting a heatmap using
#' this radial theta angle ordering, though only by using an  MDS ordination.
#'
#' That being said, the \code{getNeatOrder} function is more modular and
#' separate to the plotting, and can be applied to any kind of ordinated data
#' which can be valuable depending on the use case.
#'
#' [Rajaram & Oono (2010) NeatMap - non-clustering heat map alternatives in R](https://doi.org/10.1186/1471-2105-11-45)
#' outlines this in more detail.
#'
#' @name getNeatOrder
#'
#' @examples
#' # Load the required libraries and dataset
#' library(mia)
#' library(scater)
#' library(ComplexHeatmap) |> suppressPackageStartupMessages()
#' library(circlize) |> suppressPackageStartupMessages()
#' data(peerj13075)
#'
#' # Group data by taxonomic order
#' tse <- agglomerateByRank(peerj13075, rank = "order", onRankOnly = TRUE)
#'
#' # Transform the samples into relative abundances using CLR
#' tse <- transformAssay(
#'     tse, assay.type = "counts", method="clr", MARGIN = "cols",
#'     name="clr", pseudocount = TRUE)
#'
#' # Transform the features (taxa) into zero mean, unit variance
#' # (standardize transformation)
#' tse <- transformAssay(
#'     tse, assay.type="clr", method="standardize", MARGIN = "rows")
#'
#' # Perform PCA using calculatePCA
#' res <- calculatePCA(tse, assay.type = "standardize", ncomponents = 10)
#'
#' # Sort by radial theta and sort the original assay data
#' sorted_order <- getNeatOrder(res[, c(1,2)], centering = "mean")
#' tse <- tse[, sorted_order]
#'
#' # Define the color function and cap the colors at [-5, 5]
#' col_fun <- colorRamp2(c(-5, 0, 5), c("blue", "white", "red"))
#'
#' # Create the heatmap
#' heatmap <- Heatmap(assay(tse, "standardize"),
#'               name = "NeatMap",
#'               col = col_fun,
#'               cluster_rows = FALSE,  # Do not cluster rows
#'               cluster_columns = FALSE,  # Do not cluster columns
#'               show_row_dend = FALSE,
#'               show_column_dend = FALSE,
#'               row_names_gp = gpar(fontsize = 4),
#'               column_names_gp = gpar(fontsize = 6),
#'               heatmap_width = unit(20, "cm"),
#'               heatmap_height = unit(15, "cm")
#' )
#'
NULL

# Implementation for taking in a raw matrix.
#' @rdname getNeatOrder
#' @export
setMethod("getNeatOrder", signature = c("matrix"),
    function(x, centering = "mean", ...){
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
