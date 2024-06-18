#' Sorting by radial theta angle
#' 
#' @description \code{getNeatOrder} sorts already ordinated data by the radial theta angle.
#' This method is useful for organizing data points based on their angular 
#' position in a 2D space, typically after an ordination technique such as PCA or NMDS 
#' has been applied. 
#' 
#' The function takes in a matrix of ordinated data, optionally 
#' centers the data using specified methods (\code{mean}, \code{median}, or \code{none}), and then calculates 
#' the angle (theta) for each point relative to the centroid. The data points are then 
#' sorted based on these theta values in either ascending or descending order. 
#' 
#' One significant application of this sorting method is in plotting heatmaps. 
#' By using radial theta sorting, the relationships between data points can be preserved 
#' according to the ordination method's spatial configuration, rather than relying on 
#' hierarchical clustering, which may distort these relationships. This approach 
#' allows for a more faithful representation of the data's intrinsic structure as captured 
#' by the ordination process.
#' 
#' @param x A matrix containing the ordinated data to be sorted. Columns should represent the principal components (PCs) and rows should represent the entities being analyzed (e.g. features or samples).
#' @param dimensions A \code{character} or \code{integer} vector of length 2 specifying the columns of the matrix to use for the X and Y coordinates.
#' @param centering.method A single \code{character} value specifying the method to center the data. Options are \code{mean}, \code{median}, or \code{none} if your data is already centered. (default: method = \code{mean})
#' @param decreasing A \code{boolean} that when \code{TRUE} sorts the rows in a descending order by radial theta angle. (default: descending = \code{FALSE})
#' @param ... Additional arguments passed to other methods.
#' @return A \code{character} vector of row names in the sorted order.
#' 
#' @details 
#' It's important to note that the sechm package does actually have the functionality for plotting a heatmap using this radial theta angle ordering, though only by using an MDS ordination. 
#' 
#' This functionality can be found at:
#' 
#'  \url{https://bioconductor.org/packages/3.18/bioc/vignettes/sechm/inst/doc/sechm.html#row-ordering}. 
#'  
#' That being said, the \code{getNeatOrder} function is more modular and separate to the plotting, and 
#' can be applied to any kind of ordinated data which can be valuable depending on the use case.
#' 
#' @references
#' The below paper outlines the NeatMap method in more detail:
#' 
#' Rajaram, S., Oono, Y. NeatMap - non-clustering heat map alternatives in R. BMC Bioinformatics 11, 45 (2010). https://doi.org/10.1186/1471-2105-11-45.
#' 
#' It can be found at:
#' 
#'  \url{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-45}.
#' 
#' @name getNeatOrder
#' 
#' @examples
#' ## Load the required libraries and dataset
#' library(mia)
#' library(scater)
#' library(sechm)
#' data(peerj13075)
#' 
#' ## Group data by taxonomic order
#' tse <- agglomerateByRank(peerj13075, rank = "order", onRankOnly = TRUE)
#' 
#' ## Transform the samples into relative abundances
#' tse <- transformAssay(tse, assay.type = "counts", method="relabundance", MARGIN = "samples", name="relabundance")
#' 
#' ## Transform the features (taxa) into zero mean, unit variance (z transformation)
#' tse <- transformAssay(tse, assay.type="relabundance", method="z", MARGIN = "features", name="z")
#' 
#' ## Perform PCA using calculatePCA
#' pca_result <- calculatePCA(tse, ncomponents = 10, assay.type = "z")
#' 
#' ## Add PCA results to the TreeSE object
#' reducedDim(tse, "PCA") <- pca_result
#' 
#' ## Sort by radial theta and sort the original assay data
#' sorted_order <- getNeatOrder(reducedDim(tse, "PCA"), dimensions = c(1, 2), centering.method = "mean")
#' tse <- tse[, sorted_order]
#' 
#' ## Create the heatmap with sechm whilst retaining this radial theta ordering
#' features <- rownames(assay(tse, "z"))
#' sechm_plot <- sechm(tse, assayName = "z", features=features, do.scale=FALSE, cluster_rows=FALSE, 
#'                     sortRowsOn = NULL)
NULL

#' @rdname getNeatOrder
setGeneric("getNeatOrder", signature = c("x"),
           function(x, ...)
               standardGeneric("getNeatOrder"))


# Implementation for taking in a raw matrix.
#' @rdname getNeatOrder
#' @export
setMethod("getNeatOrder", signature = c("matrix"),
    function(x,
        dimensions = c(1, 2),
        centering.method = c("mean", "median", "none"),
        decreasing = FALSE,
        ...){
              
            # Check args
            .check_args(x, subset, dimensions, centering.method, decreasing)
            
            # Take the correct dimensions
            x <- x[, dimensions, drop = FALSE]

            # Get the theta values and order them
            theta_values <- .radial_theta(x, centering.method)
            ordering <- .get_sorted_rownames(theta_values, decreasing)
            
            return(ordering)
        }
    )


# Checks the method arguments.
.check_args <- function(x, subset, dimensions, centering.method, decreasing) {
    # Check data is a matrix
    if (!is.matrix(x)) {
        stop("Input data must be a matrix.", call. = FALSE)
    }
    
    # Check there is sufficient data
    if (nrow(x) == 0 || ncol(x) == 0) {
        stop("No data to plot. Matrix must have at least one row and one column.", call. = FALSE)
    }
    
    # Check dimensions are valid
    if (is.numeric(dimensions)) {
        if (any(dimensions > ncol(x) | dimensions < 1)) {
            stop("dimensions refer to columns that do not exist in the data.", call. = FALSE)
        }
    } else if (is.character(dimensions)) {
        if (any(!dimensions %in% colnames(x))) {
            stop("dimensions refer to column names that do not exist in the data.", call. = FALSE)
        }
    } else {
        stop("dimensions must be a vector of column indices or names.", call. = FALSE)
    }
    
    # Check dimension vector is of length 2
    if (length(dimensions) != 2) {
        stop("Exactly two dimensions must be specified.", call. = FALSE)
    }
    
    # Check centering.method
    centering.method <- match.arg(centering.method, c("mean", "median", "none"))
    
    # Check decreasing
    if (!is.logical(decreasing) || length(decreasing) != 1) {
        stop("decreasing must be a single boolean value.", call. = FALSE)
    }
    
    # Check for unique row names
    if (any(duplicated(rownames(x)))) {
        stop("Row names of the matrix must be unique.", call. = FALSE)
    }
    
    # Check for unique column names
    if (any(duplicated(colnames(x)))) {
        stop("Column names of the matrix must be unique.", call. = FALSE)
    }
    
    return(NULL)
}


# Computes the radial theta values for each row in the data matrix.
.radial_theta <- function(data, centering.method) {
    # Choose the correct centering function based on the centering.method
    center_fun <- switch(centering.method, "median" = median, "mean" = mean)
    
    # Apply the centering if there's a centering.method present
    if (!is.null(center_fun)) {
        center_vals <- apply(data, 2, center_fun)
        centered_data <- scale(data, center = center_vals, scale = FALSE)
    } else if (centering.method == "none") {
        centered_data <- data
    }
    
    # Compute the radial theta values using the centered data
    theta <- atan2(centered_data[, 2], centered_data[, 1])
    
    # Set the names of theta values to the row names of the centered data and return the theta values
    names(theta) <- rownames(centered_data)
    return(theta)
}


# Sorts the theta values and returns the ordered row names.
.get_sorted_rownames <- function(theta_values, decreasing) {
    sorted_indices <- order(theta_values, decreasing = decreasing)
    rownames <- names(theta_values)[sorted_indices]
    return(rownames)
}
