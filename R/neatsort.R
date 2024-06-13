#' Sorting by radial theta angle
#' 
#' @description \code{neatsort} sorts already ordinated data by the radial theta angle.
#' This method is useful for organizing data points based on their angular 
#' position in a 2D space, typically after an ordination technique such as PCA or NMDS 
#' has been applied. 
#' 
#' The function takes in a matrix of ordinated data, optionally 
#' centers the data using specified methods (mean, median, or none), and then calculates 
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
#' @param x A matrix containing the ordinated data to be sorted. Columns should represent the principal components (PCs) and rows should represent the entities being analyzed (e.g., features or samples).
#' @param subset A vector specifying a subset of rows to be used and retained. If NULL, all rows are used.
#' @param dimensions A vector of length 2 specifying the columns of the matrix to use for the X and Y coordinates.
#' @param centering_method A character string specifying the method to center the data. Options are "mean", "median", or "none" if your data is already centred.
#' @param sorting_order A character string specifying the order of sorting. Options are "ascending" or "descending".
#' @param ... Additional arguments passed to other methods.
#' @return A vector of row names in the sorted order.
#' 
#' @name neatsort
#' 
#' @examples
#' # Load required libraries
#' library(mia)
#' 
#' # Load the dataset
#' data(peerj13075)
#' 
#' # Agglomerate by Order and transform the data
#' tse_order <- mergeFeaturesByRank(peerj13075, rank = "order", onRankOnly = TRUE)
#' tse_order <- transformAssay(tse_order, assay.type = "counts", method="relabundance", MARGIN = "samples", name="relabundance")
#' tse_order <- transformAssay(tse_order, assay.type="relabundance", method="z", MARGIN = "features", name="z")
#' z_transformed_data <- assay(tse_order, "z")
#' 
#' # Get the top taxa and perform PCA
#' top_taxa <- getTopFeatures(tse_order, top = 10, assay.type="z")
#' top_feature_data <- z_transformed_data[top_taxa, ]
#' pca_results <- prcomp(top_feature_data, scale = TRUE)
#' scores_pca <- pca_results$x[, 1:2]
#' 
#' # Sort by radial theta and subset the transformed data
#' sorted_order <- neatsort(scores_pca, dimensions = c(1, 2), centering_method = "mean", sorting_order = "ascending")
#' ordered_transformed_data <- z_transformed_data[sorted_order, ]
NULL

#' @rdname neatsort
setGeneric("neatsort", signature = c("x"),
           function(x, ...)
               standardGeneric("neatsort"))


#' .check_neatsort_args
#' @description Checks if the input arguments for the neatsort function are valid.
#' @param x A matrix containing the ordinated data to be sorted.
#' @param subset A vector specifying a subset of rows to be used and retained.
#' @param dimensions A vector of length 2 specifying the columns of the matrix to use for the X and Y coordinates.
#' @param centering_method A character string specifying the method to center the data.
#' @param sorting_order A character string specifying the order of sorting.
.check_neatsort_args <- function(x, subset, dimensions, centering_method, sorting_order) {
    # Check data is a matrix
    if (!is.matrix(x)) {
        stop("Input data must be a matrix.")
    }
    
    # Check there is sufficient data
    if (nrow(x) == 0 || ncol(x) == 0) {
        stop("No data to plot. Matrix must have at least one row and one column.")
    }
    
    # Check subset validity
    if (!is.null(subset)) {
        if (is.numeric(subset) && any(subset > nrow(x))) {
            stop("Subset refers to rows that do not exist in the data.")
        } else if (is.character(subset) && any(!subset %in% rownames(x))) {
            stop("Subset refers to row names that do not exist in the data.")
        } else if (!is.numeric(subset) && !is.character(subset)) {
            stop("Subset must be a vector of row indices or names.")
        }
    }
    
    # Check dimensions are valid
    if (any(dimensions > ncol(x))) {
        stop("dimensions refer to columns that do not exist in the data.")
    }
    
    # Check dimension vector is of length 2
    if (length(dimensions) != 2) {
        stop("Exactly two dimensions must be specified.")
    }
    
    # Check centering_method
    centering_method <- match.arg(centering_method, c("mean", "median", "none"))
    
    # Check sorting_order
    sorting_order <- match.arg(sorting_order, c("ascending", "descending"))
    
    # Check for unique row names
    if (any(duplicated(rownames(x)))) {
        stop("Row names of the matrix must be unique.")
    }
    
    # Check for unique column names
    if (any(duplicated(colnames(x)))) {
        stop("Column names of the matrix must be unique.")
    }
}


#' @description Sorts a matrix by radial theta angle.
#' @param x A matrix containing the ordinated data to be sorted.
#' @param subset A vector specifying a subset of rows to be used and retained.
#' @param dimensions A vector of length 2 specifying the columns of the matrix to use for the X and Y coordinates.
#' @param centering_method A character string specifying the method to center the data.
#' @param sorting_order A character string specifying the order of sorting.
#' @return A vector of row names in the sorted order.
#' @rdname neatsort
#' @export
setMethod("neatsort", signature = c("matrix"),
          function(x,
                   subset = NULL,
                   dimensions = c(1, 2),
                   centering_method = c("mean", "median", "none"),
                   sorting_order = c("ascending", "descending"),
                   ...){
              
                # Check args
                .check_neatsort_args(x, subset, dimensions, centering_method, sorting_order)
              
                # Create subset if required
                if( !is.null(subset) ){
                    x <- .take_subset(x, subset)
                }

                # Take the correct dimensions
                x <- .take_dimensions(x, dimensions)

                # Get the theta values and order them
                theta_values <- .radial_theta(x, centering_method)
                ordering <- .get_sorted_rownames(theta_values, sorting_order)

                return(ordering)
          }
    )


#' .take_subset
#' @description Takes a subset of rows from the data matrix.
#' @param data The data matrix
#' @param subset The subset of rows to be retained
#' @return The subset of the data matrix
.take_subset <- function(data, subset) {
    data <- data[subset, , drop = FALSE]
    return(data)
}


#' .take_dimensions
#' @description Takes the specified columns (dimensions) from the data matrix.
#' @param data The data matrix
#' @param dimensions The columns to retain
#' @return The data matrix with only the specified dimensions
.take_dimensions <- function(data, dimensions) {
    data <- data[, dimensions, drop = FALSE]
    return(data)
}


#' .radial_theta
#' @description Computes the radial theta values for each row in the data matrix.
#' @param data The data matrix
#' @param centering_method The method used for centering the data
#' @return A named vector of theta values for each row
.radial_theta <- function(data, centering_method) {
    if (centering_method == "mean") {
        centered_data <- scale(data, center = TRUE, scale = FALSE)
    } else if (centering_method == "median") {
        centered_data <- scale(data, center = apply(data, 2, median), scale = FALSE)
    } else if (centering_method == "none") {
        centered_data <- data
    } else {
        stop("Unsupported centering method. Choose either 'mean', 'median', 'mode', or 'none'.")
    }
    
    theta <- atan2(centered_data[, 2], centered_data[, 1])
    names(theta) <- rownames(centered_data)
    
    return(theta)
}

#' .get_sorted_rownames
#' @description Sorts the theta values and returns the ordered row names.
#' @param theta_values A named vector of theta values
#' @param sorting_order The order of sorting (ascending or descending)
#' @return A vector of row names in the sorted order
.get_sorted_rownames <- function(theta_values, sorting_order) {
    sorted_indices <- order(theta_values, decreasing = (sorting_order == "descending"))
    rownames <- names(theta_values)[sorted_indices]
    return(rownames)
}

