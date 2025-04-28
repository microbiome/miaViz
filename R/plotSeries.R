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
#' plotted.
#'
#' @param col.var \code{Character scalar}. Selecting the column from
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{colData}} that
#' will be plotted. This can be used instead of \code{assay.type} for
#' visualizing temporal changes in sample metadata variable.
#'
#' @param time.col \code{Character scalar}. Selecting the column from
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{colData}} that
#' will specify values of x-axis.
#'
#' @param features \code{Character scalar}. Selects the taxa from
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{rownames}}.
#' This parameter specifies taxa whose abundances will be plotted.
#'
#' @param facet.by \code{Character scalar}. Specifies a sample grouping. Must be
#' value from
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{rowData}} or
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
#'     facet.by = "Vessel",
#'     colour.by = "rownames", colour.lab = "Feature",
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
#'     facet.by = "Vessel",
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
        facet.by = NULL,
        ...){
        # Add rownames to rowData so that they are available for plotting
        rowData(x)[["rownames"]] <- rownames(x)
        # Check input and optionally filter data
        x <- .check_series_data(
            x, time.col, assay.type, col.var, features, facet.by, ...)
        # Get the data
        args <- .get_series_data(
            x, assay.type, col.var, time.col, facet.by, ...)
        # Create the plot
        p <- do.call(.series_plotter, args)
        return(p)
    }
)

################## HELP FUNCTIONS ##########################

# This function validates the input
.check_series_data <- function(
        x, time.col, assay.type = NULL, col.var = NULL, features = NULL,
        facet.by = NULL, colour.by = color.by, color.by = colour_by,
        colour_by = color_by, color_by = NULL,
        size.by = size_by, size_by = NULL,
        linetype.by = linetype_by, linetype_by = NULL,
        ...){
    # Agglomerate data if specified
    x <- .merge_features(x, ...)
    # Check assay.type and col.var
    if( sum(c(is.null(assay.type), is.null(col.var))) != 1L ){
        stop("Either 'assay.type' or 'col.var' must be specified.",
            call. = FALSE)
    }
    if( !is.null(assay.type) ){
        .check_assay_present(assay.type, x)
    }
    temp <-  .check_metadata_variable(x, col.var, FALSE, TRUE)
    # Checks time.col
    if( !(.is_a_string(time.col) && time.col %in% names(colData(x))) ){
        stop("'time.col' must be a name of column of colData(x)",
            call. = FALSE)
    }
    # Checks features
    if( !(is.null(features) || (is.character(features) &&
            all(features %in% rownames(x)))) ){
        stop("'y' must be in rownames(x). \n If 'rank' was used, ",
            "check that 'y' matches agglomerated data.", call. = FALSE)
    }
    if( !is.null(col.var) && !is.null(features) ){
        stop("'features' can be specified only when 'assay.type' is ",
            "specified.", call. = FALSE)
    }
    # Select taxa that user has specified
    if( !is.null(features) ){
        x <- x[features,]
    }
    # Gets warning or error if too many taxa are selected. Too many taxa
    # cannot be plotted since otherwise the plot is too crowded.
    if( !is.null(assay.type) && length(rownames(x)) > 20 ){
        stop("Over 20 taxa selected. 20 or under allowed.", call. = FALSE)
    } else if ( !is.null(assay.type) && length(rownames(x)) > 10 ){
        warning("Over 10 taxa selected.", call. = FALSE)
    }
    # Check that the styling parameters can be found from rowData or colData
    temp <- .check_metadata_variable(x, facet.by, !is.null(assay.type), TRUE)
    temp <- .check_metadata_variable(x, colour.by, !is.null(assay.type), TRUE)
    temp <- .check_metadata_variable(x, size.by, !is.null(assay.type), TRUE)
    temp <- .check_metadata_variable(x, linetype.by, !is.null(assay.type), TRUE)
    return(x)
}

# This function fetches data from SE object. It outputs data in a format that
# can directly be plotted with .series_plotter().
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr group_by mutate ungroup
#' @importFrom stats sd
#' @importFrom mia meltSE
#' @importFrom SummarizedExperiment rowData<-
.get_series_data <- function(
        x, assay.type, col.var, time.col, facet.by,
        colour.by = color.by, color.by = colour_by, colour_by = color_by,
        color_by = NULL,
        size.by = size_by, size_by = NULL,
        linetype.by = linetype_by, linetype_by = NULL,
        ...){
    # To disable "no visible binding for global variable" message in cmdcheck
    Y <- NULL

    # Melt SE object or take variables from colData
    cols <- c(col.var, time.col, facet.by, colour.by, size.by, linetype.by)
    if( !is.null(assay.type) ){
        # Check from where these parameters can be found
        col_names <- colData(x) |> colnames()
        row_names <- rowData(x) |> colnames()
        col_vars <- cols[ cols %in% colnames(colData(x)) ]
        row_vars <- cols[ cols %in% colnames(rowData(x)) ]
        # Get data from SE
        df <- meltSE(
            x, assay.type = assay.type,
            row.name = "feature",
            add.row = row_vars,
            add.col = col_vars
        )
    } else{
        df <- colData(x)[, cols, drop = FALSE]
        assay.type <- col.var
    }
    df <- df |> as.data.frame()

    # Check that data is numeric. Both x and y axis must be numeric.
    if( !is.numeric(df[[assay.type]]) ){
        stop("'", ifelse(!is.null(col.var), "col.var", "assay.type"),
            "' must specify numeric data.", call. = FALSE)
    }
    if( !is.numeric(df[[time.col]]) ){
        stop("'time.col' must specify numeric data.", call. = FALSE)
    }

    # If time point replicates are present, calculate sd and mean for each
    # timepoint
    cols <- c("feature", time.col, facet.by)
    cols <- cols[ cols %in% colnames(df) ]
    sd_col <- NULL
    if( anyDuplicated(df[, cols, drop = FALSE]) ){
        sd_col <- "sd"
        df <- df |>
            group_by(across(all_of(cols))) |>
            mutate(
                !!sd_col := sd(.data[[assay.type]], na.rm = TRUE),
                !!assay.type := mean(.data[[assay.type]], na.rm = TRUE),
            ) |>
            distinct()
    }

    # Ribbon can only have one color per time point and group, otherwise it lead
    # to an error. Of course, it does not make sense to color lines with
    # multiple values either, but it does not lead to an error.
    colour.ribbon <- colour.by
    if( !is.null(colour.by) ){
        sum_of_colors <- df |>
            group_by(across(all_of(cols))) |>
            summarise(n = n_distinct(.data[[colour.by]]))
        if( any(sum_of_colors[["n"]] > 1L) ){
            colour.ribbon <- NULL
            message("Multiple values per facet and time point detected. ",
                "Coloring might not make sense.")
        }
    }

    # Add all the necessary information to attributes
    attributes(df) <- c(
        attributes(df),
        y = assay.type,
        sd = sd_col,
        x = time.col,
        facet.by = facet.by,
        colour.by = colour.by,
        colour.ribbon = colour.ribbon,
        size.by = size.by,
        linetype.by = linetype.by
    )
    # Return argument list with all the passed arguments
    args <- c(list(df = df), list(...))
    return(args)
}

# This function gets time series data as an input and creates a plot from it.
.series_plotter <- function(
        df,
        add.legend = add_legend, add_legend = TRUE,
        line.alpha = line_alpha, line_alpha = 1,
        line.type = line_type, line_type = NULL,
        line.width = line_width, line_width = 1,
        line.width.range = line_width_range, line_width_range = c(0.5,3),
        ribbon.alpha = ribbon_alpha, ribbon_alpha = 0.3,
        ncol = 1L,
        scales = "fixed",
        xlab = NULL, ylab = NULL, colour.lab = color.lab, color.lab = NULL,
        ...){
    #
    if( !.is_an_integer(ncol) ){
        stop("'ncol' must be an integer.", call. = FALSE)
    }
    if( !.is_a_string(scales) ){
        stop("'scales' must be a string.", call. = FALSE)
    }
    if( !(is.null(xlab) || .is_a_string(xlab)) ){
        stop("'xlab' must be a string.", call. = FALSE)
    }
    if( !(is.null(ylab) || .is_a_string(ylab)) ){
        stop("'ylab' must be a string.", call. = FALSE)
    }
    if( !(is.null(colour.lab) || .is_a_string(colour.lab)) ){
        stop("'colour.lab' must be a string.", call. = FALSE)
    }
    #
    # Initialize a plot a plot
    plot_out <- ggplot(df, aes(
        x = .data[[attributes(df)[["x"]]]],
        y = .data[[attributes(df)[["y"]]]])
        )

    # Fetches arguments for geom_line to plot mean
    line_args <- .get_line_args(
        colour_by = if(!is.null(attributes(df)[["colour.by"]]))
            df[[attributes(df)[["colour.by"]]]],
        linetype_by = if(!is.null(attributes(df)[["linetype.by"]]))
            df[[attributes(df)[["linetype.by"]]]],
        size_by = if(!is.null(attributes(df)[["size.by"]]))
            df[[attributes(df)[["size.by"]]]],
        alpha = line.alpha,
        linetype = line.type,
        linewidth = line.width)
    # Adds arguments to the plot
    plot_out <- plot_out +
        do.call(geom_line, line_args$args)

    # If data included multiple observations per timepoint, standard deviation
    # is calculated. Create a ribbon that shows the SD.
    if( !is.null(attributes(df)[["sd"]]) ){
        ribbon_args <- .get_ribbon_args(
            colour_by = attributes(df)[["colour.ribbon"]],
            mean_col = attributes(df)[["y"]], sd_col = attributes(df)[["sd"]],
            alpha = ribbon.alpha)
        plot_out <- plot_out +
            do.call(geom_ribbon, ribbon_args$args)
    }

    # Apply line_width_range, i.e., how much line gets bigger and narrower
    if( !is.null(attributes(df)[["size.by"]]) ){
        if( is.numeric(df[[attributes(df)[["size.by"]]]]) ){
            SIZEFUN <- scale_size_continuous
        } else {
            SIZEFUN <- scale_size_discrete
        }
        plot_out <- plot_out +
            SIZEFUN(range = line_width_range)
    }

    # Resolves the colours
    if( !is.null(attributes(df)[["colour.by"]]) ){
        plot_out <- .resolve_plot_colours(
            plot_out, colour_by = df[[attributes(df)[["colour.by"]]]],
            colour_by_name = attributes(df)[["colour.by"]], fill = FALSE)
        if( !is.null(attributes(df)[["sd"]]) ){
            plot_out <- .resolve_plot_colours(
                plot_out, colour_by = df[[attributes(df)[["colour.by"]]]],
                colour_by_name = attributes(df)[["colour.by"]], fill = TRUE)
        }
    }

    # If facetting is specified, create separate panels
    if( !is.null(attributes(df)[["facet.by"]]) ){
        plot_out <- plot_out + facet_wrap(
            ~ .data[[attributes(df)[["facet.by"]]]],
            ncol = ncol, scales = scales)
    }

    # Add additional guides
    plot_out <- .add_extra_line_guide(
        plot_out, linetype_by = attributes(df)[["linetype.by"]],
        size_by = attributes(df)[["size.by"]])
    # To choose if legend is kept, and its position
    plot_out <- .add_legend(plot_out, add_legend)
    # Set a theme
    plot_out <- plot_out +
        theme_classic()

    # Adjust labels if they are specified by user
    if( !is.null(xlab) ){
        plot_out <- plot_out + labs(x = xlab)
    }
    if( !is.null(ylab) ){
        plot_out <- plot_out + labs(y = ylab)
    }
    if( !is.null(colour.lab) ){
        plot_out <- plot_out + guides(
            fill = guide_legend(title = colour.lab),
            colour = guide_legend(title = colour.lab)
            )
    }
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
