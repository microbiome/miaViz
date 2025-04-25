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
#'   \item \code{facet.by}: \code{Character vector}. Specifies variables from
#'   \code{colData(x)} or \code{rowData(x)} used for facetting.
#'   (Default: \code{NULL})
#'   \item \code{fill.by}: \code{Character scalar}. Specifies variable from
#'   \code{colData(x)} or \code{rowData(x)} used for coloring.
#'   (Default: \code{NULL})
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
#' plotHistogram(
#'     tse,
#'     assay.type = "clr",
#'     features = rownames(tse)[1:5],
#'     facet.by = "rownames"
#' )
#'
#' # Calculate shannon diversity and visualize its distribution with density
#' # plot. Different sample types are separated with color.
#' tse <- addAlpha(tse, index = "shannon")
#' plotHistogram(
#'     tse,
#'     col.var = "shannon",
#'     layout = "density",
#'     fill.by = "SampleType"
#' )
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
        # Add rownames to rowData so that they are available for plotting.
        rowData(x)[["rownames"]] <- rownames(x)
        # Check input
        args <- list(
            x = x, assay.type = assay.type, features = features,
            row.var = row.var, col.var = col.var)
        args <- c(args, list(...))
        temp <- do.call(.check_input_for_histogram, args)
        # Get the data from the object
        args[["mode"]] <- "histogram"
        df <- do.call(.get_histogram_data, args)
        # Create a histogram
        args <- c(list(df = df), list(...))
        args <- args[ !names(args) %in% c("fill.by", "facet.by") ]
        p <- do.call(.plot_histogram, args)
        return(p)
    }
)

#' @rdname plotHistogram
#' @export
setMethod("plotBarplot", signature = c(x = "SummarizedExperiment"),
    function(x, assay.type = NULL, features = NULL, row.var = NULL,
            col.var = NULL, ...){
        # Add rownames to rowData so that they are available for plotting.
        rowData(x)[["rownames"]] <- rownames(x)
        # Check input
        args <- list(
            x = x, assay.type = assay.type, features = features,
            row.var = row.var, col.var = col.var)
        args <- c(args, list(...))
        temp <- do.call(.check_input_for_histogram, args)
        # Get the data from the object
        args[["mode"]] <- "barplot"
        df <- do.call(.get_histogram_data, args)
        # Create a barplot
        args <- c(list(df = df), list(...))
        args <- args[ !names(args) %in% c("fill.by", "facet.by") ]
        p <- do.call(.plot_barplot, args)
        return(p)
    }
)

################################ HELP FUNCTIONS ################################

# This function harmonizes the input check for histogram and bar plot
.check_input_for_histogram <- function(
        x, assay.type, features, row.var, col.var,
        facet.by = NULL, fill.by = NULL, ...){
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
    if( !(is.null(row.var) || .is_a_string(row.var)) ){
        stop("'row.var' must be NULL or single character value.",
            call. = FALSE)
    }
    if( !(is.null(col.var) || .is_a_string(col.var)) ){
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
    # Check that facetting and cooring variables can be found correctly from
    # row or column metadata
    temp <- .check_metadata_variable(
        x, fill.by,
        row = length(c(row.var, assay.type))>0,
        col = length(c(col.var, assay.type))>0,
        multiple = TRUE
    )
    temp <- .check_metadata_variable(
        x, facet.by,
        row = length(c(row.var, assay.type))>0,
        col = length(c(col.var, assay.type))>0,
        multiple = TRUE
    )
    # User cannot specify more than 2 variables for facetting
    if( length(facet.by) > 2L ){
        stop("'facet.by' cannot specify more than 2 variables.", call. = FALSE)
    }
    # It does not make sense to visualize same variable as used for facetting or
    # coloring
    if( !is.null(col.var) && col.var %in% c(facet.by, fill.by) ){
        stop("'col.var' must not equal to 'fill.by' or 'facet.by'.",
            call. = FALSE)
    }
    if( !is.null(row.var) && row.var %in% c(facet.by, fill.by) ){
        stop("'row.var' must not equal to 'fill.by' or 'facet.by'.",
            call. = FALSE)
    }
    return(NULL)
}

# This function retrieves the data from TreeSE and returns a data.frame, ready
# for inputting it to plotting function.
#' @importFrom tidyr pivot_longer
.get_histogram_data <- function(
        x, assay.type, features, row.var, col.var,
        fill.by = NULL, facet.by = NULL, mode = "histogram", ...){
    #
    if( !(.is_a_string(mode) && mode %in% c("histogram", "barplot")) ){
        stop("'mode' must be 'histogram' or 'barplot'.", call. = FALSE)
    }
    # If assay.type is specified, get melted data
    all_vars <- c(fill.by, facet.by)
    if( !is.null(assay.type) ){
        # Specify whether to retrieve data from rowData or colData
        row_vars <- vapply(all_vars, function(var){
            var %in% colnames(rowData(x))
        }, logical(1L))
        col_vars <- all_vars[ !row_vars ]
        row_vars <- all_vars[ row_vars ]
        #
        df <- meltSE(
            x, assay.type = assay.type,
            col.var = "id",
            add.col = col_vars,
            add.row = row_vars
        )
        colnames(df)[ colnames(df) == assay.type ] <- "value"
    }
    # If features wre specified, subset data
    if( !is.null(features) ){
        df <- df[ df[["FeatureID"]] %in% features, , drop = FALSE]
    }
    # If row.var was specified, get the data from rowData
    if( !is.null(row.var) ){
        df <- rowData(x)[, c(row.var, fill.by, facet.by), drop = FALSE]
        colnames(df)[ colnames(df) == row.var ] <- "value"
    }
    # If col.var was specified, get the data from colData
    if( !is.null(col.var) ){
        df <- colData(x)[, c(col.var, fill.by, facet.by), drop = FALSE]
        colnames(df)[ colnames(df) == col.var ] <- "value"
    }
    # Check that values are numeric for histogram and categorical for barplot
    if( mode == "histogram" && !is.numeric(df[["value"]]) ){
        stop("Values must be numeric.", call. = FALSE)
    }
    if( mode == "barplot" &&
            !(is.character(df[["value"]]) || is.factor(df[["value"]])) ){
        stop("Values must be categorical.", call. = FALSE)
    }
    # Facetting and coloring values must be categorical
    are_correct <- vapply(df[, facet.by, drop = FALSE], function(col){
        is.character(col) || is.factor(col)
    }, logical(1L))
    if( !all(are_correct) ){
        stop("'facet.by' must specify categorical values.", call. = FALSE)
    }
    are_correct <- vapply(df[, fill.by, drop = FALSE], function(col){
        is.character(col) || is.factor(col)
    }, logical(1L))
    if( !all(are_correct) ){
        stop("'fill.by' must specify categorical values.", call. = FALSE)
    }
    # Add x-axis title, facetting, and coloring info to attributes so that we
    # can use it in plotting function.
    attributes(df)[["x"]] <- c(assay.type, col.var, row.var)
    attributes(df)[["facet.by"]] <- facet.by
    attributes(df)[["fill.by"]] <- fill.by
    return(df)
}

# This function gets data.frame and creates a plot.
.plot_histogram <- function(
        df, layout = "histogram", color = colour, colour = "black",
        fill = "grey35", alpha = 0.4, scales = "fixed",
        position = ifelse(
            !is.null(attributes(df)[["fill.by"]]), "dodge2", "identity"),
        ...){
    # To disable "no visible binding for global variable" message in cmdcheck
    value <- facet_by <- NULL
    # Check layout
    supported_layouts <- c("histogram", "density")
    if( !(.is_a_string(layout) && all(layout %in% supported_layouts) ) ){
        stop("'layout' must be from the following options: '",
            paste0(supported_layouts, collapse = "', '"), "'", call. = FALSE)
    }
    # Initialize a plot
    p <- ggplot(df, aes(
        x = value,
        fill = if(!is.null(attributes(df)[["fill.by"]]))
            .data[[attributes(df)[["fill.by"]]]] else fill
        ))
    # Either create histogram or density
    if( layout == "density" ){
        p <- p + geom_density(color = color, alpha = alpha, ...)
    } else{
        p <- p + geom_histogram(
            color = color, alpha = alpha, position = position, ...)
    }
    # Apply facetting
    if( length(attributes(df)[["facet.by"]]) > 0L ){
        p <- p + facet_grid(attributes(df)[["facet.by"]], scales = scales)
    }
    # Adjust theme
    p <- p + theme_classic()
    # Adjust titles
    p <- p + labs(x = attributes(df)[["x"]])
    if( !is.null(attributes(df)[["fill.by"]]) ){
        p <- p + labs(fill = attributes(df)[["fill.by"]])
    } else{
        p <- p + guides(fill = "none")
    }
    return(p)
}

# This function gets data.frame and creates a plot.
.plot_barplot <- function(
        df, color = colour, colour = "black", fill = "grey35", alpha = 0.4,
        scales = "fixed",
        position = ifelse(
            !is.null(attributes(df)[["fill.by"]]), "dodge2", "identity"),
        ...){
    # To disable "no visible binding for global variable" message in cmdcheck
    value <- facet_by <- NULL
    # Initialize a plot
    p <- ggplot(df, aes(
        x = value,
        fill = if(!is.null(attributes(df)[["fill.by"]]))
            .data[[attributes(df)[["fill.by"]]]] else fill
        ))
    # Either create barplot
    p <- p + geom_bar(color = color, alpha = alpha, position = position, ...)
    # Apply facetting
    if( length(attributes(df)[["facet.by"]]) > 0L ){
        p <- p + facet_grid(attributes(df)[["facet.by"]], scales = scales)
    }
    # Adjust theme
    p <- p + theme_classic()
    # Adjust titles
    p <- p + labs(x = attributes(df)[["x"]])
    if( !is.null(attributes(df)[["fill.by"]]) ){
        p <- p + labs(fill = attributes(df)[["fill.by"]])
    } else{
        p <- p + guides(fill = "none")
    }
    return(p)
}
