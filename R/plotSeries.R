#' Plot Series
#'
#' This function plots time series data.
#'
#' @param x a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#'
#' @param assay.type \code{Character scalar}. Specifies the
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}} to be
#' plotted. (Default: \code{"counts"})
#'
#' @param time.col \code{Character scalar}. selecting the column from
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{colData}} that
#' will specify values of x-axis.
#'
#' @param features \code{Character scalar}. Selects the taxa from
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{rownames}}.
#' This parameter specifies taxa whose abundances will be plotted.
#'
#' @param group \code{Character scalar}. Specifies a sample grouping. Must be
#' value from
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{colData}}. If
#' \code{NULL}, grouping is not applied. (Default: \code{NULL})
#'
#' @param ... additional parameters for plotting.
#' \itemize{
#'   \item \code{rank} \code{Character scalar}. A taxonomic rank, that is used
#'   to agglomerate the data. (Default: \code{NULL})
#'
#'   \item \code{colour.by} \code{Character scalar}. A column name from
#'   \code{rowData(x)} or \code{colData(x)}, that is used to divide observations
#'   to different colors. If \code{NULL}, this is not applied.
#'   (Default: \code{NULL})
#'
#'   \item \code{linetype.by} \code{Character scalar}. A column name from
#'   \code{rowData(x)} or \code{colData(x)}, that is used to divide observations
#'   to different line types. If \code{NULL}, this is not applied.
#'   (Default: \code{NULL})
#'
#'   \item \code{size.by}: \code{Character scalar}. A column name from
#'   \code{rowData(x)} or \code{colData(x)}, that is used to divide observations
#'   to different size types. If \code{NULL}, this is not applied.
#'   (Default: \code{NULL})
#'
#'   \item \code{facet.cols}: \code{Logical scalar}. Should the sample groups
#'   specified by \code{group} be spitted into facets? If \code{group} is
#'   specified and \code{facet.cols = FALSE}, features are facetted instead.
#'   (Default: \code{FALSE})
#'
#'   \item \code{ncol}: \code{Numeric scalar}. if facets are applied,
#'   \code{ncol} defines many columns should be for plotting the different
#'   facets. (Default: \code{1L})
#'
#'   \item \code{scales} \code{Character scalar}. Defines the behavior of the
#'   scales of each facet. The value is passed into
#'   \code{\link[ggplot2:facet_wrap]{facet_wrap}}. (Default: \code{"fixed"})
#' }
#' See \code{\link{mia-plot-args}} for more details i.e. call
#' \code{help("mia-plot-args")}
#'
#' @details
#' This function creates series plot, where x-axis includes e.g. time points,
#' and y-axis abundances of selected taxa. If there are multiple observations
#' for single system and time point, mean and standard deviation is plotted.
#'
#' @return
#' A \code{ggplot2} object
#'
#' @name plotSeries
#'
#' @examples
#' \dontrun{
#' library(mia)
#' # Load data from miaTime package
#' library("miaTime")
#' data(SilvermanAGutData)
#' tse <- SilvermanAGutData
#'
#' # Plots 2 most abundant taxa, which are colored by their family
#' plotSeries(
#'     tse,
#'     assay.type = "counts",
#'     time.col = "DAY_ORDER",
#'     features = getTop(tse, 2),
#'     colour.by = "Family"
#' )
#'
#' # Counts relative abundances
#' tse <- transformAssay(tse, method = "relabundance")
#'
#' # Selects taxa
#' taxa <- c("seq_1", "seq_2", "seq_3", "seq_4", "seq_5")
#'
#' # Plots relative abundances of phylums
#' plotSeries(
#'     tse[taxa,],
#'     time.col = "DAY_ORDER",
#'     colour.by = "Family",
#'     linetype.by = "Phylum",
#'     assay.type = "relabundance"
#' )
#'
#' # In addition to 'colour.by' and 'linetype.by', 'size.by' can also be used
#' # to group taxa.
#' plotSeries(
#'     tse,
#'     time.col = "DAY_ORDER",
#'     features = getTop(tse, 5),
#'     colour.by = "Family",
#'     size.by = "Phylum",
#'     assay.type = "counts"
#' )
#'
#' # If the data includes multiple systems, e.g., patients or bioreactors,
#' # one can plot each system separately
#' plotSeries(
#'     tse,
#'     time.col = "DAY_ORDER",
#'     assay.type = "relabundance",
#'     features = getTop(tse, 5),
#'     group = "Vessel",
#'     linetype.by = "Pre_Post_Challenge",
#'     scales = "free"
#' )
#'
#' # One can visualize colData variables by specifying col.var
#' # First calculate alpha diversity index to visualize.
#' tse <- addAlpha(tse, index = "shannon")
#' # Then create a plot
#' plotSeries(
#'     tse,
#'     col.var = "shannon",
#'     time.col = "DAY_ORDER",
#'     group = "Vessel"
#' )
#'
#' }
NULL

#' @rdname plotSeries
#' @importFrom SummarizedExperiment colData rowData assay
#' @importFrom mia meltSE
#' @importFrom stats sd
#' @export
setMethod("plotSeries", signature = c(x = "SummarizedExperiment"),
    function(
        x,
        time.col,
        assay.type = NULL,
        col.var = NULL,
        features = NULL,
        group = NULL,
        ...){
        ###################### Input check #######################
        if( sum(c(is.null(assay.type), is.null(col.var))) != 1L ){
            stop("Either 'assay.type' or 'col.var' must be specified.",
                call. = FALSE)
        }
        if( !is.null(col.var) && !is.null(features) ){
            stop("'features' can be specified only when 'assay.type' is ",
                "specified.", call. = FALSE)
        }
        # Checks assay.type
        if( !is.null(assay.type) ){
            .check_assay_present(assay.type, x)
        }
        # Check col.var
        if( !(is.null(col.var) || (is.character(col.var) &&
                all(col.var %in% colnames(colData(x))))) ){
            stop("'col.var' must be NULL or specify a column(s) from ",
                "colData(x).", call. = FALSE)
        }
        # Checks X
        if( !(.is_a_string(time.col) && time.col %in% names(colData(x))) ){
            stop("'time.col' must be a name of column of colData(x)",
                call. = FALSE)
        }
        # Agglomerate data if specified
        x <- .merge_features(x, ...)
        # Checks features
        if( !(is.null(features) || (is.character(features) &&
                all(features %in% rownames(x)))) ){
            stop("'y' must be in rownames(x). \n If 'rank' was used, ",
                "check that 'y' matches agglomerated data.", call. = FALSE)
        }
        # Select taxa that user has specified
        if (!is.null(features) ){
            x <- x[features,]
        }
        # Checks group
        if( !(is.null(group) ||
                (.is_a_string(group) && group %in% names(colData(x)))) ){
            stop("'group' must be a name of column of colData(x)",
                call. = FALSE)
        }
        # Gets warning or error if too many taxa are selected. Too many taxa
        # cannot be plotted since otherwise the plot is too crowded.
        if( !is.null(assay.type) && length(rownames(x)) > 20 ){
            stop("Over 20 taxa selected. 20 or under allowed.", call. = FALSE)
        } else if ( !is.null(assay.type) && length(rownames(x)) > 10 ){
            warning("Over 10 taxa selected.", call. = FALSE)
        }
        ###################### Input check end ####################
        # Get the data
        args <- .get_series_data(x, assay.type, col.var, time.col, group, ...)
        # Create the plot
        p <- do.call(.series_plotter, args)
        return(p)
    }
)

################## HELP FUNCTIONS ##########################

# This function fetches data from SE object. It outputs data in a format that
# can directly be plotted with .series_plotter().
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr group_by mutate ungroup
#' @importFrom stats sd
#' @importFrom mia meltSE
#' @importFrom SummarizedExperiment rowData<-
.get_series_data <- function(
        x, assay.type, col.var, time.col, group,
        colour.by = color.by, color.by = colour_by, colour_by = color_by,
        color_by = NULL,
        size.by = size_by, size_by = NULL,
        linetype.by = linetype_by, linetype_by = NULL,
        facet.cols = FALSE,
        ...){
    # To disable "no visible binding for global variable" message in cmdcheck
    Y <- NULL

    # Check that the styling parameters can be found from rowData or colData
    col_names <- colData(x) |> colnames()
    row_names <- rowData(x) |> colnames()
    if( !(is.null(colour.by) || (.is_a_string(colour.by) &&
            colour.by %in% c(row_names, col_names)) ) ){
        stop("'colour.by' must be a string from rowData(x) or colData(x).",
            call. = FALSE)
    }
    if( !(is.null(size.by) || (.is_a_string(size.by) &&
            size.by %in% c(row_names, col_names)) ) ){
        stop("'size.by' must be a string from rowData(x) or colData(x).",
            call. = FALSE)
    }
    if( !(is.null(linetype.by) || (.is_a_string(linetype.by) &&
            linetype.by %in% c(row_names, col_names)) ) ){
        stop("'linetype.by' must be a string from rowData(x) or colData(x).",
            call. = FALSE)
    }
    # Check that facet.cols is boolean value
    if( !.is_a_bool(facet.cols) ){
        stop("'facet.cols' must be TRUE or FALSE.", call. = FALSE)
    }

    # Check from where these parameters can be found
    cols <- c(
        colour_by = colour.by, size_by = size.by, linetype_by = linetype.by)
    ind <- match(cols, col_names)
    col_vars <- setNames(col_names[ind], names(cols)) |> na.omit()
    ind <- match(cols, row_names)
    row_vars <- setNames(row_names[ind], names(cols)) |> na.omit()

    # If user wanted to plot variable from colData, user cannot specify
    # variables from rowData
    if( !is.null(col.var) && length(row_vars)> 0L ){
        stop("If 'col.var' is specified, variables from rowData(x) cannot be ",
            "visualized.", call. = FALSE)
    }

    # Melt SE object or take variables from colData
    if( !is.null(assay.type) ){
        plot_data <- meltSE(
            x, assay.type = assay.type,
            row.name = "feature",
            add.row = row_vars,
            add.col = c(time.col, group, col_vars)
        )
    } else{
        plot_data <- colData(x)[
            , c(col.var, time.col, group, col_vars), drop = FALSE] |>
            as.data.frame()
        # Put to long format so that it matches with data when assay.type is
        # specified
        plot_data <- plot_data |>
            pivot_longer(
                cols = col.var, names_to = "feature", values_to = "values")
        assay.type <- "values"
    }

    # If time point replicates are present, calculate sd and mean for each
    # timepoint
    cols <- setNames(c("feature", time.col, group), c("feature", "X", group))
    cols <- cols[ cols %in% colnames(plot_data) ]
    if( anyDuplicated(plot_data[, cols, drop = FALSE]) ){
        # Summarize the data to mean and sd
        plot_data <- plot_data |>
            group_by(across(all_of(cols))) |>
            mutate(
                sd = sd(.data[[assay.type]], na.rm = TRUE),
                !!assay.type := mean(.data[[assay.type]], na.rm = TRUE)
            ) |>
            ungroup()
    }

    # Select only certain columns and remove duplicates
    cols <- c(cols, Y = assay.type, sd = "sd", row_vars, col_vars)
    cols <- cols[ cols %in% colnames(plot_data) ]
    plot_data <- plot_data[, cols, drop = FALSE]
    plot_data <- plot_data[!duplicated(plot_data), , drop = FALSE]
    # Rename columns to harmonize input of plotter function
    colnames(plot_data) <- names(cols)

    # If user did not specify coloring, different features are colored by
    # default. If grouping was specified, we color features (if we facet by
    # groups) or groups (if we facet by features).
    if( is.null(colour.by) && (is.null(group) ||
            (!is.null(group) && facet.cols)) ){
        colour.by <- "feature"
        plot_data[["colour_by"]] <-  plot_data[[colour.by]]
    } else if( is.null(colour.by) && !is.null(group) && !facet.cols ){
        colour.by <- group
        plot_data[["colour_by"]] <-  as.factor( plot_data[[colour.by]] )
    }
    # If grouping was specified, user can specify if facetting is applied to
    # features or sample group
    if( !is.null(group) && facet.cols ){
        plot_data[["facet_by"]] <-  as.factor( plot_data[[group]] )
    } else if( !is.null(group) && !facet.cols ){
        plot_data[["facet_by"]] <-  plot_data[["feature"]]
    }

    # Create an argument list to feed to plotter
    args <- list(
        plot_data = plot_data,
        colour_by = colour.by,
        linetype_by = linetype.by,
        size_by = size.by,
        xlab = as.character(time.col),
        ylab = as.character(assay.type)
    )
    args <- c(args, list(...))
    return(args)
}

# This function gets time series data as an input and creates a plot from it.
.series_plotter <- function(
        plot_data,
        xlab,
        ylab,
        colour_by,
        linetype_by,
        size_by,
        add_legend = add.legend,
        add.legend = TRUE,
        line_alpha = line.alpha,
        line.alpha = 1,
        line_type = line.type,
        line.type = NULL,
        line_width = line.width,
        line.width = 1,
        line_width_range = line.width.range,
        line.width.range =  c(0.5,3),
        ribbon_alpha = ribbon.alpha,
        ribbon.alpha = 0.3,
        ncol = 1L,
        scales = "fixed",
        ...){
    #
    if( !.is_an_integer(ncol) ){
        stop("'ncol' must be an integer.", call. = FALSE)
    }
    if( !.is_a_string(scales) ){
        stop("'scales' must be an integer.", call. = FALSE)
    }
    #
    # Initialize a plot a plot
    plot_out <- ggplot(plot_data, aes(x = .data[["X"]], y = .data[["Y"]])) +
        labs(x = xlab, y = ylab)

    # if sd column is present add a ribbon
    if(!is.null(plot_data[["sd"]])){
        ribbon_args <- .get_ribbon_args(
            colour_by = colour_by, alpha = ribbon_alpha)
        plot_out <- plot_out +
            do.call(geom_ribbon, ribbon_args$args)
    }

    # Fetches arguments for geom_line to plot mean
    line_args <- .get_line_args(
        colour_by = colour_by,
        linetype_by = linetype_by,
        size_by = size_by,
        alpha = line_alpha,
        linetype = line_type,
        linewidth = line_width)
    # Adds arguments to the plot
    plot_out <- plot_out +
        do.call(geom_line, line_args$args)
    # apply line_width_range
    if( !is.null(size_by) ){
        if(is.numeric(plot_data$size_by)){
            SIZEFUN <- scale_size_continuous
        } else {
            SIZEFUN <- scale_size_discrete
        }
        plot_out <- plot_out +
            SIZEFUN(range = line_width_range)
    }

    # Resolves the colours
    plot_out <- .resolve_plot_colours(
        plot_out, plot_data$colour_by, colour_by, fill = FALSE)
    if( !is.null(plot_data[["sd"]]) ){
        plot_out <- .resolve_plot_colours(
            plot_out, plot_data$colour_by, colour_by, fill = TRUE)
    }

    # If facetting is specified, create separate panels
    if( "facet_by" %in% colnames(plot_data) ){
        plot_out <- plot_out +
            facet_wrap(~facet_by, ncol = ncol, scales = scales)
    }

    # add additional guides
    plot_out <- .add_extra_line_guide(plot_out, linetype_by, size_by)
    # To choose if legend is kept, and its position
    plot_out <- .add_legend(plot_out, add_legend)
    # Set a theme
    plot_out <- plot_out +
        theme_classic()

    return(plot_out)
}

# This function ensures that we add legends for styling parameters.
.add_extra_line_guide <- function(plot_out, linetype_by, size_by) {
    guide_args <- list()
    if (!is.null(linetype_by)) {
        guide_args$linetype <- guide_legend(title = linetype_by)
    }
    if (!is.null(size_by)) {
        guide_args$linewidth <- guide_legend(title = size_by)
    }
    if (length(guide_args)) {
        plot_out <- plot_out + do.call(guides, guide_args)
    }
    return(plot_out)
}
