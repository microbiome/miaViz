#' @name
#' plotScree
#'
#' @title
#' Create a scree plot
#'
#' @description
#' \code{plotScree} generates a scree plot to visualize the eigenvalues.
#' The eigenvalues can be provided either as a part of a
#' \code{TreeSummarizedExperiment} object or as a separate \code{vector}.
#' This plot illustrates the decline in eigenvalues across components,
#' helping to assess the importance of each component.
#'
#' @details
#' \code{plotScree} generates a scree plot to visualize the relative importance
#' of components in dimensionality reduction techniques such as Principal
#' Component Analysis (PCA) or Principal Coordinate Analysis (PCoA). If the
#' input is a \code{TreeSummarizedExperiment} object, the function extracts
#' eigenvalues from the specified reduced dimension slot, which requires that
#' dimensionality reduction has been performed beforehand using a dedicated
#' function. Alternatively, if the input is a \code{vector} or an
#' \code{eigenvals} object, these values are directly used as eigenvalues for
#' the plot.
#'
#' The plot can include a combination of barplot, points, connecting lines,
#' and labels, which can be controlled using the \code{show.*} parameters.
#'
#' An option to show cumulative explained variance is also available by setting
#' \code{add.cumulative = TRUE}.
#'
#' @return
#' A \code{ggplot2} object
#'
#' @param x a
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-constructor]{TreeSummarizedExperiment}}
#' \code{\link[vegan:eigenvals]{eigenvals}} or a vector.
#'
#' @param dimred \code{Character scalar} or \code{integer scalar}. Determines
#' the reduced dimension to plot. This is used when \code{x} is a
#' \code{TreeSummarizedExperiment} to extract the eigenvalues from
#' \code{reducedDim(x, dimred)}.
#'
#' @param ... additional parameters for plotting
#' \itemize{
#'   \item \code{show.barplot}: \code{Logical scalar}. Whether to show a
#'   barplot. (Default: \code{TRUE})
#'
#'   \item \code{show.points}: \code{Logical scalar}. Whether to show a
#'   points. (Default: \code{TRUE})
#'
#'   \item \code{show.line}: \code{Logical scalar}. Whether to show a
#'   line. (Default: \code{TRUE})
#'
#'   \item \code{show.labels}: \code{Logical scalar}. Whether to show a
#'   labels for each point. (Default: \code{FALSE})
#'
#'   \item \code{add.proportion}: \code{Logical scalar}. Whether to show
#'   proportion of explained variance, i.e., raw eigenvalues.
#'   (Default: \code{TRUE})
#'
#'   \item \code{add.cumulative}: \code{Logical scalar}. Whether to show
#'   cumulative explained variance calculated from eigenvalues.
#'   (Default: \code{FALSE})
#'
#'   \item \code{n}: \code{Integer scalar}. Number of eigenvalues to plot.
#'   If \code{NULL}, all eigenvalues are plotted. (Default: \code{NULL})
#'
#'   \item \code{show.names}: \code{Logical scalar}. Whether to show names of
#'   components in x-axis. If \code{FALSE}, the index of component is shown
#'   instead of names. (Default: \code{FALSE})
#'
#'   \item \code{eig.name}: \code{Character scalar}. The name of the attribute
#'   in \code{reducedDim(x, dimred)} that contains the eigenvalues.
#'   (Default: \code{c("eig", "varExplained")})
#' }
#'
#' @examples
#'
#' library(miaViz)
#' library(scater)
#'
#' data("enterotype", package = "mia")
#' tse <- enterotype
#'
#' # Run  PCA and store results into TreeSE
#' tse <- transformAssay(tse, method = "clr", pseudocount = TRUE)
#' tse <- runPCA(tse, assay.type = "clr")
#'
#' # Plot scree plot
#' plotScree(tse, "PCA", add.cumulative = TRUE)
#'
NULL

#' @rdname plotScree
#' @export
setMethod("plotScree", signature = c(x = "SingleCellExperiment"),
    function(x, dimred, ...){
        eig <- .get_eigenvalues(x, dimred, ...)
        p <- plotScree(eig, ...)
        return(p)
    }
)

#' @rdname plotScree
#' @export
setMethod("plotScree", signature = c(x = "ANY"),
    function(x, ...){
        # Check that the the values are in correct format
        is_correct <- length(x) > 0L &&
            ((is.vector(x) && is.numeric(x)) || is(x, "eigenvals"))
        if( !is_correct ){
            stop("'x' must be a numeric vector or class 'eigenvals'.",
                call. = FALSE)
        }
        # Prepare data for plotting
        plot_data <- .prepare_data(x, ...)
        # Create a plot
        p <- .scree_plotter(plot_data, ...)
        return(p)
    }
)

################################ HELP FUNCTIONS ################################

# This function retrieves the eigenvalues from reducedDim. The ordination must
# be calculated with dedicaded function in mia or scater so that the eigenvalues
# are stored in correct place.
#' @importFrom SingleCellExperiment reducedDim
.get_eigenvalues <- function(
        x, dimred, eig.name = c("eig", "varExplained"), ...){
    # Get reducedDim
    if( !((.is_a_string(dimred) && dimred %in% reducedDimNames(x)) ||
            (.is_an_integer(dimred) && dimred > 0 &&
            dimred <= length(reducedDims(x)))) ){
        stop("'dimred' must specify a valid reducedDim.", call. = FALSE)
    }
    reduced_dim <- reducedDim(x, dimred)
    # Get eigenvalues
    eig.name <- eig.name[ eig.name %in% names(attributes(reduced_dim)) ]
    if( length(eig.name) != 1L ){
        stop("'eig.name' must specify a name of attributes from ",
            "reducedDim(x, dimred).", call. = FALSE)
    }
    eig <- attributes(reduced_dim)[[ eig.name ]]
    return(eig)
}

# This function creates a data.frame from the input vector. The output is ready
# for plotter.
.prepare_data <- function(
        x, add.proportion = TRUE, add.cumulative = FALSE, n = NULL,
        show.names = FALSE, ...){
    # Input check
    if( !.is_a_bool(add.proportion) ){
        stop("'add.proportion' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(add.cumulative) ){
        stop("'add.cumulative' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !(is.null(n) || .is_an_integer(n) ) ){
        stop("'n' must be NULL or integer.", call. = FALSE)
    }
    if( !.is_a_bool(show.names) ){
        stop("'show.names' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    # Create a data.frame with eigenvalues
    df <- data.frame(y = x)
    df[["x"]] <- factor(rownames(df), levels = unique(rownames(df)))
    df[["type"]] <- "proportion"
    # Calculate cumulative proportion
    df_cum <- df
    df_cum[["y"]] <- cumsum(df_cum[["y"]]) /
        sum(df_cum[["y"]], na.rm = TRUE)
    df_cum[["type"]] <- "cumulative"
    df <- rbind(df, df_cum)

    # Based on user preference, keep proportion or/and cumulative values
    if( !add.proportion ){
        df <- df[df[["type"]] != "proportion", ]
    }
    if( !add.cumulative ){
        df <- df[df[["type"]] != "cumulative", ]
    }
    # If user has specified, take only n first eigenvalues
    if( !is.null(n) ){
        n <- levels(df[["x"]])[ seq_len(n) ]
        df <- df[ df[["x"]] %in% n,  ]
    }
    # Replace names with integers to keep the x-axis of plot tidier
    if( !show.names ){
        df[["x"]] <- as.integer(df[["x"]])
    }
    return(df)
}

# This function creates a scree plot. The input is data.frame that includes
# 2 columns: one for eigenvalues and other for principal component name.
#' @importFrom scales pretty_breaks
.scree_plotter <- function(
        df, show.points = TRUE, show.line = TRUE, show.barplot = FALSE,
        show.labels = FALSE, ...){
    # To disable "no visible binding for global variable" message in cmdcheck
    x <- y <- type <- NULL
    # Input check
    if( !.is_a_bool(show.points) ){
        stop("'show.points' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(show.line) ){
        stop("'show.line' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(show.barplot) ){
        stop("'show.barplot' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(show.labels) ){
        stop("'show.labels' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    # If user wants to add proportion and cumulative values to same plot, we
    # scale cumulative values into same scale as proportion.
    if( length(unique(df[["type"]])) > 1L && !(show.labels || show.barplot) ){
        ind <- df[["type"]] == "cumulative"
        df[ind, "y"] <- df[ind, "y"] * max(df[!ind, "y"]) # Scale
    }

    # Create base plot
    p <- ggplot(df, aes(
        x = x,
        y = y,
        group = type,
        colour = if(length(unique(.data[["type"]])) > 1L &&
            !(show.labels || show.barplot) ) type
        ))
    # Add layers based on user preference
    if( show.points ){
        p <- p + geom_point()
    }
    if( show.line ){
        p <- p + geom_line()
    }
    if( show.barplot ){
        p <- p + geom_col(width = 0.5)
    }
    if( show.labels ){
        p <- p + geom_label(aes(label = round(y, 2)))
    }

    # If user wants to add barplots or labels with both cumulative and
    # propotion values, the plot is splitted into two facets. Otherwise the
    # the plot would be too messy to read.
    if( length(unique(df[["type"]])) > 1L && (show.labels || show.barplot) ){
        p <- p + facet_grid(rows = vars(type), scales = "free_y")
    }
    # If user wants to plot both cumulative and proportion, but with points
    # and lines only, we can add tem to same plot. Now the plot has two y-axis;
    # proportion at left and cumulative at right
    if( length(unique(df[["type"]])) > 1L && !(show.labels || show.barplot) ){
        p <- p + scale_y_continuous(
            name = "Proportion",
            sec.axis = sec_axis(
                ~ . / max(df[["y"]]), name = "Cumulative proportion"))
    }
    # Adjust labels in a case where either proportion or cumulative was plotted
    if( length(unique(df[["type"]])) == 1L ){
        p <- p + labs(x = "PC", y = "Eigenvalue")
    }
    # Adjust the x-axis to display a subset of evenly spaced values for
    # improved readability
    if( is.numeric(df[["x"]]) ){
        p <- p +
            scale_x_continuous(breaks = pretty_breaks())
    }
    # Adjust theme and remove legend
    p <- p + theme_classic() +
        theme(legend.position = "none")
    return(p)
}
