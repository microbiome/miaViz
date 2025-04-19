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
        p <- .plot_histogram(df, ...)
        return(p)
    }
)

#' @rdname plotHistogram
#' @export
setMethod("plotBarplot", signature = c(x = "SummarizedExperiment"),
    function(x, assay.type = NULL, features = NULL, row.var = NULL,
            col.var = NULL, ...){
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
        p <- .plot_barplot(df, ...)
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
    temp <- .check_metadata_variable(x, fill.by, row = TRUE, col = TRUE, enable.multi = TRUE)
    temp <- .check_metadata_variable(x, facet.by, row = TRUE, col = TRUE, enable.multi = TRUE)
    if( length(facet.by) > 2L ){
        stop("'facet.by' cannot specify more than 2 variables.", call. = FALSE)
    }
    return(NULL)
}

# This function retrieves the data from TreeSE and returns a data.frame, ready
# for inputting it to plotting function.
#' @importFrom tidyr pivot_longer
.get_histogram_data <- function(
        x, assay.type, features, row.var, col.var, mode = "histogram",
        fill.by = NULL, facet.by = NULL, ...){
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
        # colnames(df)[ colnames(df) == "FeatureID" ] <- "facet_by"
    }
    # If row.var was specified, get the data from rowData
    if( !is.null(row.var) && is.null(assay.type) ){
        df <- rowData(x)[, row.var, drop = FALSE]
    }
    # If col.var was specified, get the data from colData
    if( !is.null(col.var) && is.null(assay.type) ){
        df <- colData(x)[, col.var, drop = FALSE]
    }
    # If either row.var or col.var was specified, convert data into long format
    if( (!is.null(row.var) || !is.null(col.var)) && is.null(assay.type) ){
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
    attributes(df)[["facet.by"]] <- c("facet_by", facet.by)[ c("facet_by", facet.by) %in% colnames(df) ]
    attributes(df)[["fill.by"]] <- fill.by
    attributes(df)[["x"]] <- c(assay.type, col.var, row.var)
    return(df)
}

# This function gets data.frame and creates a plot.
.plot_histogram <- function(
        df, layout = "histogram", color = colour, colour = "black", alpha = 0.4,
        position = ifelse(!is.null(attributes(df)[["fill.by"]]), "dodge2", "identity"), ...){
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
            .data[[attributes(df)[["fill.by"]]]]
        ))
    # Either create histogram or density
    if( layout == "density" ){
        p <- p + geom_density(alpha = alpha, ...)
    } else{
        p <- p + geom_histogram(color = color, alpha = alpha, position = position, ...)
    }
    # If there are multiple features and user wants to plot them separately,
    # apply facetting
    # if( "facet_by" %in% colnames(df) ){
    #     p <- p + facet_wrap(vars(facet_by))
    # }
    if( length(attributes(df)[["facet.by"]]) > 0L ){
        p <- p + facet_grid( attributes(df)[["facet.by"]] )
    }
    # Adjust theme
    p <- p + theme_classic()
    # Otherwise, we add correct title from original variable name
    p <- p + labs(
        x = attributes(df)[["x"]]
    )
    # Add correct titles for aesthetics. The titles are the original variable
    # names.
    if( !is.null(attributes(df)[["fill.by"]]) ){
        p <- p + labs(fill = attributes(df)[["fill.by"]])
    }
    return(p)
}

# This function gets data.frame and creates a plot.
.plot_barplot <- function(
        df, color = colour, colour = "black", ...){
    # To disable "no visible binding for global variable" message in cmdcheck
    value <- facet_by <- NULL
    # Initialize a plot
    p <- ggplot(df, aes(
        x = value,
        fill = if(!is.null(attributes(df)[["fill.by"]]))
            .data[[attributes(df)[["fill.by"]]]]
        ))
    # Either create barplot
    p <- p + geom_bar(color = color, ...)
    # If there are multiple features and user wants to plot them separately,
    # apply facetting
    if( "facet_by" %in% colnames(df) ){
        p <- p + facet_wrap(vars(facet_by))
    }
    # Adjust theme
    p <- p + theme_classic()
    # Otherwise, we add correct title from original variable name
    p <- p + labs(
        x = attributes(df)[["x"]]
    )
    # Add correct titles for aesthetics. The titles are the original variable
    # names.
    if( !is.null(attributes(df)[["fill.by"]]) ){
        p <- p + labs(fill = attributes(df)[["fill.by"]])
    }
    return(p)
}

# This function checks whether variable can be found from colData or rowData.
.check_metadata_variable <- function(
        tse, var, row = FALSE, col = FALSE,
        var.name = .get_name_in_parent(var), enable.multi = FALSE){
    # If the variable is not NULL
    if( !is.null(var) ){
        # It must be a string and found from colData/rowData
        is_string <- ifelse(enable.multi, is.character(var), .is_a_string(var))
        check_values <- c()
        check_values <- c(check_values, if(col) colnames(colData(tse)))
        check_values <- c(check_values, if(row) colnames(rowData(tse)))
        var_found <- all( var %in% check_values )
        if( !(is_string && var_found) ){
            stop("'", var.name, "' must be", ifelse(enable.multi, "", "a single "), "character value from the ",
                 "following options: '",
                 paste0(check_values, collapse = "', '"), "'", call. = FALSE)
        }
    }
    return(NULL)
}
