#' Plot Scree Plot or Eigenvalues
#'
#' \code{plotScree} creates a scree plot or eigenvalues plot starting from a
#' TreeSummarizedExperiment object or a vector of eigenvalues. This visualization
#' shows how the eigenvalues decrease across components.
#'
#' @param x a
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-constructor]{TreeSummarizedExperiment}}
#' or a vector of eigenvalues.
#' 
#' @param dimred \code{Character scalar} or \code{integer scalar}. Determines
#' the reduced dimension to plot. This is used when x is a TreeSummarizedExperiment
#' to extract the eigenvalues from \code{reducedDim(x, dimred)}.
#' 
#' @param cumulative \code{Logical scalar}. Whether to show cumulative explained 
#' variance. (Default: \code{FALSE}).
#' 
#' @param names \code{Character vector}. Optional names for the components 
#' that will be displayed on the x-axis. If not provided, the components 
#' are labeled sequentially as 1, 2, 3, etc.
#' 
#' @param ... additional parameters for plotting
#' \describe{
#'   \item{\code{show.barplot}}{Logical scalar. Whether to show a barplot. 
#'   (Default: \code{TRUE}).}
#'   \item{\code{show.points}}{Logical scalar. Whether to show points. 
#'   (Default: \code{TRUE}).}
#'   \item{\code{show.line}}{Logical scalar. Whether to show a line connecting 
#'   points. (Default: \code{TRUE}).}
#'   \item{\code{show.labels}}{Logical scalar. Whether to show labels for each 
#'   point. (Default: \code{FALSE}).}
#' }
#'
#' @details
#' \code{plotScree} creates a scree plot or eigenvalues plot, which is useful
#' for visualizing the relative importance of components in dimensionality
#' reduction techniques like PCA, RDA, or CCA. When the input is a
#' TreeSummarizedExperiment, the function extracts eigenvalues from the specified
#' reduced dimension slot. When the input is a vector, it directly uses these
#' values as eigenvalues.
#' 
#' The plot can include a combination of barplot, points, connecting lines,
#' and labels, which can be controlled using the \code{show.*} parameters.
#' 
#' An option to show cumulative explained variance is also available by setting
#' \code{cumulative = TRUE}.
#' 
#' @return 
#' A \code{ggplot2} object 
#'
#' @name plotScree
#'
#' @examples
#' # Load necessary libraries
#' library(ggplot2)
#' 
#' # Load dataset
#' library(miaViz)
#' data("enterotype", package = "mia")
#' tse <- enterotype
#'  
#' # Run RDA and store results into TreeSE
#' tse <- addRDA(
#'     tse,
#'     formula = assay ~ ClinicalStatus + Gender + Age,
#'     FUN = getDissimilarity,
#'     distance = "bray",
#'     na.action = na.exclude
#'     )
#' 
#' # Plot scree plot
#' plotScree(tse, "RDA")
#' 
NULL

#' @rdname plotScree
#' @export
setGeneric("plotScree", signature = c("x"),
    function(x, ...) standardGeneric("plotScree"))

#' @rdname plotScree
#' @export
setMethod("plotScree", signature = c(x = "SingleCellExperiment"),
    function(x, dimred, cumulative = FALSE, ...) {
        # Check if dimred exists
        if (!dimred %in% reducedDimNames(x)) {
            stop("'dimred' must specify a valid reducedDim.", call. = FALSE)
        }
        # Get reducedDim
        reduced_dim <- reducedDim(tse, dimred)
      
        # Extract eigenvalues
        # Check if data is available
        ind <- names(attributes(reduced_dim)) %in% c("eig", "varExplained")
        if( any(ind) ){
            # Add explained variance
            eig <- attributes(reduced_dim)[ind][[1]]
        } else{
            stop("No eigenvalues found in the specified reducedDim.", 
                 call. = FALSE)
        }
      
        # Call the vector method
        plotScree(as.numeric(eig), names(eig), cumulative = cumulative, ...)
    }
)

#' @rdname plotScree
#' @export
setMethod("plotScree", signature = c(x = "vector"),
    function(x, names = NULL, cumulative = FALSE, ...) {
        # Ensure 'x' is numeric
        if (!is.numeric(x)) {
            stop("'x' must be a numeric vector.", call. = FALSE)
        }
        plot_data <- .prepare_data(x, names, cumulative)
        # plot vector
        .scree_plotter(plot_data, cumulative = cumulative, ...)
    }
)

.prepare_data <- function(x, names = NULL, cumulative = FALSE) {
    df <- data.frame(
        Component = if (!is.null(names)) names else seq_along(x),
        Eigenvalue = x
    )
    
    # Calculate cumulative proportion if needed
    if (cumulative) {
        df$CumulativeProportion <- cumsum(df$Eigenvalue) / sum(df$Eigenvalue)
    }
    
    return(df)
}

.scree_plotter <- function(df, show.barplot = TRUE, 
                           show.points = TRUE, 
                           show.line = TRUE, show.labels = FALSE, 
                           cumulative = FALSE, ...) {
    # Create base plot
    p <- ggplot(df, aes(x = Component, y = if (cumulative) 
        CumulativeProportion else Eigenvalue))
    
    # Add layers based on user preferences
    if (show.barplot) {
        p <- p + geom_col(fill = "lightblue", color = "black")
    }
    if (show.points) {
        p <- p + geom_point(size = 3)
    }
    if (show.line) {
        p <- p + geom_line()
    }
    if (show.labels) {
        p <- p + geom_text(aes(label = round(if (cumulative) 
            CumulativeProportion else Eigenvalue, 2)), 
            vjust = -0.5)
    }
    
    # Customize appearance
    p <- p + theme_minimal() +
        labs(x = "Component", 
             y = if (cumulative) "Cumulative Proportion of Variance" 
             else "Eigenvalue",
             title = if (cumulative) "Cumulative Explained Variance" 
             else "Scree Plot")
    
    return(p)
}
