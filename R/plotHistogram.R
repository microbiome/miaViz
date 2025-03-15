#' @name
#' plotHistogram
#'
#' @title
#' Create histogram or barplot of \code{assay}, \code{rowData} or \code{colData}
#'
#' @description
#' This methods visualizes abundances or variables from \code{rowData} or
#' \code{colData}.
#'
#' @details
#' Histogram and bar plot are a basic visualization techniques in quality
#' control. It helps to visualize the distribution of data. \code{plotAbundance}
#' allows researcher to visualise the abundance from \code{assay}, or variables
#' from \code{rowData} or \code{colData}. For visualizing categorical values,
#' one can utilize \code{plotBarplot}.
#'
#' \code{\link[=plotAbundanceDensity]{plotAbundanceDensity}} function is related
#' to \code{plotHistogram}. However, the former visualizes the most prevalent
#' features, while the latter can be used more freely to explore the
#' distributions.
#'
#' @return
#' A \code{ggplot2} object.
#'
#' @inheritParams plotAbundance
#'
#' @param assay.type \code{NULL} or \code{character scalar}. Specifies the
#' abundace table to plot. (Default: \code{NULL})
#'
#' @param features \code{NULL} or \code{character vector}. If \code{assay.type}
#' is specified, this specifies rows to visualize in different facets. If
#' \code{NULL}, whole data is visualized as a whole. (Default: \code{NULL})
#'
#' @param row.var \code{NULL} or \code{character vector}. Specifies a variable
#' from \code{rowData(x)} to visualize. (Default: \code{NULL})
#'
#' @param col.var \code{NULL} or \code{character vector} Specifies a variable
#' from \code{colData(x)} to visualize. (Default: \code{NULL})
#'
#' @param ... Additional parameters for plotting.
#' \itemize{
#'   \item \code{layout}: \code{Character scalar}. Specifies the layout of plot.
#'   Must be either \code{"histogram"} or \code{"density"}.
#'   (Default: \code{"histogram"})
#' }
#'
#' @examples
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#'
#' # Visualize the counts data. There are lots of zeroes.
#' plotHistogram(tse, assay.type = "counts")
#'
#' # Apply transformation
#' tse <- transformAssay(tse, method = "clr", pseudocount = TRUE)
#' # And plot specified rows
#' plotHistogram(tse, assay.type = "clr", features = rownames(tse)[1:10])
#'
#' # Calculate shannon diversity and visualize its distribution with density
#' # plot
#' tse <- addAlpha(tse, index = "shannon")
#' plotHistogram(tse, col.var = "shannon", layout = "density")
#'
#' # For categorical values, one can utilize a bar plot
#' plotBarplot(tse, col.var = "SampleType")
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[=plotAbundanceDensity]{plotAbundanceDensity}}
#'   \item \code{\link[scater:plotExpression]{scater::plotExpression}}
#'   \item \code{\link[scater:plotRowData]{scater::plotRowData}}
#'   \item \code{\link[scater:plotColData]{scater::plotColData}}
#' }
#'
NULL

#' @rdname plotHistogram
#' @export
setMethod("plotHistogram", signature = c(x = "SummarizedExperiment"),
    function(x, assay.type = NULL, features = NULL, row.var = NULL,
            col.var = NULL, ...){
        # Check input
        args <- list(
            x = x, assay.type = assay.type, features = features,
            row.var = row.var, col.var = col.var)
        temp <- do.call(.check_input_for_histogram, args)
        # Get the data from the object
        args <- c(args, list(...))
        args[["mode"]] <- "histogram"
        df <- do.call(.get_histogram_data, args)
        # Create a histogram
        p <- .plot_histogram(df, ...)
        return(p)
    }
)

#' @rdname plotHistogram
#' @export
setMethod("plotBarplot", signature = c(x = "SummarizedExperiment"),
    function(x, assay.type = NULL, features = NULL, row.var = NULL,
            col.var = NULL, ...){
        # Check input
        args <- list(
            x = x, assay.type = assay.type, features = features,
            row.var = row.var, col.var = col.var)
        temp <- do.call(.check_input_for_histogram, args)
        # Get the data from the object
        args <- c(args, list(...))
        args[["mode"]] <- "barplot"
        df <- do.call(.get_histogram_data, args)
        # Create a barplot
        p <- .plot_barplot(df, ...)
        return(p)
    }
)
################################ HELP FUNCTIONS ################################

# This function harmonizes the input check for histogram and bar plot
.check_input_for_histogram <- function(
        x, assay.type, features, row.var, col.var){
    # Either assay.type. row.var or col.var must be specified
    if( sum(c(is.null(assay.type), is.null(row.var), is.null(col.var))) !=
            2L ){
        stop("Please specify either 'assay.type', 'row.var', or 'col.var'.",
            call. = FALSE)
    }
    # features cannot be specified if row.var or col.var is specified
    if( is.null(assay.type) && !is.null(features) ){
        stop("'features' can be specified only when 'assay.type is ",
            "specified.", call. = FALSE)
    }
    # Check that the values are correct
    if( !is.null(features) && is.null(rownames(x)) ){
        stop("'x' must have rownames.", call. = FALSE)
    }
    if( !(is.null(assay.type) || .is_a_string(assay.type)) ){
        stop("'assay.type' must be NULL or single character value.",
            call. = FALSE)
    }
    if( !(is.null(row.var) || is.character(row.var)) ){
        stop("'row.var' must be NULL or single character value.",
            call. = FALSE)
    }
    if( !(is.null(col.var) || is.character(col.var)) ){
        stop("'col.var' must be NULL or single character value.",
            call. = FALSE)
    }
    if( !(is.null(features) || is.character(features)) ){
        stop("'features' must be NULL or single character value.",
            call. = FALSE)
    }
    # If parameters are specified, check that they can be found
    if( !is.null(assay.type) ){
        .check_assay_present(assay.type, x)
    }
    if( !is.null(row.var) && !all(row.var %in% colnames(rowData(x))) ){
        stop("'row.var' must be from the following options: '",
            paste0(colnames(rowData(x)), collapse = "', '"), "'",
            call. = FALSE)
    }
    if( !is.null(col.var) && !all(col.var %in% colnames(colData(x))) ){
        stop("'col.var' must be from the following options: '",
            paste0(colnames(colData(x)), collapse = "', '"), "'",
            call. = FALSE)
    }
    if( !is.null(features) && !all(features %in% rownames(x)) ){
        stop("'feature' must specify features from rownames(x).",
            call. = FALSE)
    }
    return(NULL)
}

# This function retrieves the data from TreeSE and returns a data.frame, ready
# for inputting it to plotting function.
#' @importFrom tidyr pivot_longer
.get_histogram_data <- function(
        x, assay.type, features, row.var, col.var, mode = "histogram", ...){
    #
    if( !(.is_a_string(mode) && mode %in% c("histogram", "barplot")) ){
        stop("'mode' must be 'histogram' or 'barplot'.", call. = FALSE)
    }
    # If assay.type is specified, get melted data
    if( !is.null(assay.type) ){
        df <- meltSE(x, assay.type = assay.type, col.var = "id")
        colnames(df)[ colnames(df) == assay.type ] <- "value"
    }
    # If features wre specified, subset data
    if( !is.null(features) ){
        df <- df[ df[["FeatureID"]] %in% features, , drop = FALSE]
        colnames(df)[ colnames(df) == "FeatureID" ] <- "facet_by"
    }
    # If row.var was specified, get the data from rowData
    if( !is.null(row.var) ){
        df <- rowData(x)[, row.var, drop = FALSE]
    }
    # If col.var was specified, get the data from colData
    if( !is.null(col.var) ){
        df <- colData(x)[, col.var, drop = FALSE]
    }
    # If either row.var or col.var was specified, convert data into long format
    if( !is.null(row.var) || !is.null(col.var) ){
        cols <- colnames(df)
        df[["id"]] <- rownames(df)
        df <- df |> as.data.frame() |>
            pivot_longer(
                cols = all_of(cols),
                names_to = "facet_by",
                values_to = "value"
                )
    }
    # If there is single facetting value, disable facetting
    if( length(unique(df[["facet_by"]])) == 1L ){
        df[["facet_by"]] <- NULL
    }
    # Check that values are numeric for histogram and categorical for barplot
    if( mode == "histogram" && !is.numeric(df[["value"]]) ){
        stop("Values must be numeric.", call. = FALSE)
    }
    if( mode == "barplot" &&
            !(is.character(df[["value"]]) || is.factor(df[["value"]])) ){
        stop("Values must be categorical.", call. = FALSE)
    }
    return(df)
}

# This function gets data.frame and creates a plot.
.plot_histogram <- function(
        df, layout = "histogram", color = colour, colour = "black",
        fill = "white", ...){
    # To disable "no visible binding for global variable" message in cmdcheck
    value <- facet_by <- NULL
    # Check layout
    supported_layouts <- c("histogram", "density")
    if( !(.is_a_string(layout) && all(layout %in% supported_layouts) ) ){
        stop("'layout' must be from the following options: '",
            paste0(supported_layouts, collapse = "', '"), "'", call. = FALSE)
    }
    # Initialize a plot
    p <- ggplot(df, aes(x = value))
    # Either create histogram or density
    if( layout == "density" ){
        p <- p + geom_density(...)
    } else{
        p <- p + geom_histogram(color = color, fill = fill, ...)
    }
    # If there are multiple features and user wants to plot them separately,
    # apply facetting
    if( "facet_by" %in% colnames(df) ){
        p <- p + facet_wrap(vars(facet_by))
    }
    # Adjust theme
    p <- p + theme_classic()
    return(p)
}

# This function gets data.frame and creates a plot.
.plot_barplot <- function(
        df, color = colour, colour = "black", fill = "white", ...){
    # To disable "no visible binding for global variable" message in cmdcheck
    value <- facet_by <- NULL
    # Initialize a plot
    p <- ggplot(df, aes(x = value))
    # Either create barplot
    p <- p + geom_bar(color = color, fill = fill, ...)
    # If there are multiple features and user wants to plot them separately,
    # apply facetting
    if( "facet_by" %in% colnames(df) ){
        p <- p + facet_wrap(vars(facet_by))
    }
    # Adjust theme
    p <- p + theme_classic()
    return(p)
}