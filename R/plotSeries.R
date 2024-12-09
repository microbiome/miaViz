#' Plot Series
#'
#' This function plots series data.
#'
#' @param object a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#'
#' @param assay.type \code{Character scalar}. selecting the
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}} to be
#' plotted. (Default: \code{"counts"})
#'
#' @param assay_name Deprecated. Use \code{assay.type} instead.
#'
#' @param x \code{Character scalar}. selecting the column from
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{ColData}} that
#' will specify values of x-axis.
#'
#' @param y \code{Character scalar}. Selects the taxa from
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{rownames}}.
#'   This parameter specifies taxa whose abundances will be plotted.
#'
#' @param rank \code{Character scalar}. A taxonomic rank, that is used
#'   to agglomerate the data. Must be a value of \code{taxonomicRanks()}
#'   function. (Default: \code{NULL})
#'
#' @param colour.by \code{Character scalar}. A taxonomic rank, that is used to
#' color plot. Must be a value of \code{taxonomicRanks()} function.
#' (Default: \code{NULL})
#'
#' @param colour_by Deprecated. Use \code{colour.by} instead.
#'
#' @param linetype.by \code{Character scalar}. A taxonomic rank, that
#' is used to divide taxa to different line types. Must be a value of
#' \code{taxonomicRanks()} function. (Default: \code{NULL})
#'
#' @param linetype_by Deprecated. Use \code{linetype.by} instead.
#'
#' @param size.by \code{Character scalar}. A taxonomic rank, that is
#' used to divide taxa to different line size types. Must be a value of
#' \code{taxonomicRanks()} function. (Default: \code{NULL})
#'
#' @param size_by Deprecated. Use \code{size.by} instead.
#'
#' @param ... additional parameters for plotting. See
#' \code{\link{mia-plot-args}} for more details i.e. call
#' \code{help("mia-plot-args")}
#'
#' @details
#' This function creates series plot, where x-axis includes e.g. time points,
#' and y-axis abundances of selected taxa.
#'
#' @return
#' A \code{ggplot2} object
#'
#' @name plotSeries
#'
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#'
#' @examples
#' \dontrun{
#' library(mia)
#' # Load data from miaTime package
#' library("miaTime")
#' data(SilvermanAGutData)
#' object <- SilvermanAGutData
#'
#' # Plots 2 most abundant taxa, which are colored by their family
#' plotSeries(object,
#'            x = "DAY_ORDER",
#'            y = getTop(object, 2),
#'            colour.by = "Family")
#'
#' # Counts relative abundances
#' object <- transformAssay(object, method = "relabundance")
#'
#' # Selects taxa
#' taxa <- c("seq_1", "seq_2", "seq_3", "seq_4", "seq_5")
#'
#' # Plots relative abundances of phylums
#' plotSeries(object[taxa,],
#'            x = "DAY_ORDER",
#'            colour.by = "Family",
#'            linetype.by = "Phylum",
#'            assay.type = "relabundance")
#'
#' # In addition to 'colour.by' and 'linetype.by', 'size.by' can also be used
#' # to group taxa.
#' plotSeries(object,
#'            x = "DAY_ORDER",
#'            y = getTop(object, 5),
#'            colour.by = "Family",
#'            size.by = "Phylum",
#'            assay.type = "counts")
#' }
NULL

#' @rdname plotSeries
#' @importFrom SummarizedExperiment colData rowData assay
#' @importFrom mia meltSE
#' @importFrom stats sd
#' @export
setMethod("plotSeries", signature = c(object = "SummarizedExperiment"),
    function(
        object,
        x,
        y = NULL,
        rank = NULL,
        colour.by = colour_by,
        colour_by = NULL,
        size.by = size_by,
        size_by = NULL,
        linetype.by = linetype_by,
        linetype_by = NULL,
        assay.type = assay_name, assay_name = "counts",
        ...){
        ###################### Input check #######################
        # Checks assay.type
        .check_assay_present(assay.type, object)
        # Checks X
        if( !.is_a_string(x) ||
            !(x %in% names(colData(object))) ){
            stop("'x' must be a name of column of colData(object)",
                call. = FALSE)
        }
        # If rank is not null, data will be agglomerated by rank
        if( !is.null(rank) ){
            # Check rank
            .check_taxonomic_rank(rank, object)

            # Agglomerates the object
            object <- agglomerateByRank(object, rank = rank)
        }
        # Checks Y
        # If Y is not null, user has specified it
        if (!is.null(y)){
            if(!is.character(y) || !all( y %in% rownames(object))){
                stop("'y' must be in rownames(x). \n If 'rank' was used, ",
                    "check that 'y' matches agglomerated data.", call. = FALSE)
            }
            # Select taxa that user has specified
            object <- object[y,]
        }
        # Gets warning or error if too many taxa are selected. Too many taxa
        # cannot be plotted since otherwise the plot is too crowded.
        if( length(rownames(object)) > 20 ){
            stop("Over 20 taxa selected. 20 or under allowed.", call. = FALSE)
        } else if ( length(rownames(object)) > 10 ){
            warning("Over 10 taxa selected.", call. = FALSE)
        }
        ###################### Input check end ####################
        # Get the data
        plot_data <- .get_series_data(
            object, assay.type, x, colour.by, size.by, linetype.by)
        # Adjust labels
        xlab <- paste0(x)
        ylab <- paste0(assay.type)
        # Create the plot
        p <- .series_plotter(
            plot_data,
            xlab = xlab,
            ylab = ylab,
            colour_by = colour.by,
            linetype_by = linetype.by,
            size_by = size.by,
            ...)
        return(p)
    }
)

################## HELP FUNCTIONS ##########################

# This function fetches data from SE object. It outputs data in a format that
# can directly be plotted with .series_plotter().
#' @importFrom dplyr group_by summarize ungroup
#' @importFrom stats sd
#' @importFrom mia meltSE
#' @importFrom SummarizedExperiment rowData<-
.get_series_data <- function(
        object, assay.type, x, colour.by, size.by, linetype.by){
    # To disable "no visible binding for global variable" message in cmdcheck
    Y <- NULL
    # Get variables that can be found from rowData
    row_vars <- c(
        colour_by = colour.by, size_by = size.by, linetype_by = linetype.by)
    # Rename rowData columns. If we do not do this, this might cause
    # problems in melting step if "x" is named equally to colour.by, size.by or
    # linetype.by. Duplicate colnames get suffix, and variable names are not
    # then detected correctly.
    colnames(rowData(object))[ match(row_vars, colnames(rowData(object))) ] <-
        names(row_vars)
    row_vars <- names(row_vars)

    # Melt SE object. If value is not found from rowData/colData, user get
    # informative error message.
    plot_data <- meltSE(
        object, assay.type = assay.type,
        row.name = "feature",
        add.row = row_vars,
        add.col = x)
    # Rename abundance value and timepoint columns
    colnames(plot_data)[ colnames(plot_data) == assay.type ] <- "Y"
    colnames(plot_data)[ colnames(plot_data) == x ] <- "X"
    # If time point replicates are present calculate sd and mean for each
    # timepoint
    if( anyDuplicated(plot_data[["X"]]) ){
        # Step 1: Summarize the data
        summary_data <- plot_data %>%
            group_by(!!sym("X"), !!sym("feature")) %>%
            summarize(
                sd = sd(.data[["Y"]], na.rm = TRUE),
                Y = mean(.data[["Y"]], na.rm = TRUE)
            ) %>%
            ungroup()
        # Step 2: Join the summarized data back to the original data
        plot_data <- plot_data %>%
            select(-Y) %>%
            dplyr::left_join(
                summary_data, by = c("X" = "X", "feature" = "feature"))
    }
    return(plot_data)
}

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
        ribbon.alpha = 0.3){
    # fall back for feature grouping
    if(is.null(colour_by)){
        colour_by <- "feature"
        plot_data$colour_by <-  plot_data$feature
    }
    # Creates a "draft" of a plot
    plot_out <- ggplot(plot_data, aes(x = .data[["X"]], y = .data[["Y"]])) +
        labs(x = xlab, y = ylab)
    # if sd column is present add a ribbon
    if(!is.null(plot_data$sd)){
        ribbon_args <- .get_ribbon_args(
            colour_by = colour_by, alpha = ribbon_alpha)
        plot_out <- plot_out +
            do.call(geom_ribbon, ribbon_args$args)
    }
    # Fetches arguments for geom_line
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
    if (!is.null(size_by)) {
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
    if(!is.null(plot_data$sd)){
        plot_out <- .resolve_plot_colours(
            plot_out, plot_data$colour_by, colour_by, fill = TRUE)
    }
    # add additional guides
    plot_out <- .add_extra_line_guide(plot_out, linetype_by, size_by)
    # To choose if legend is kept, and its position
    plot_out <- .add_legend(plot_out, add_legend)
    plot_out <- plot_out +
        theme_classic()
    plot_out
}

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
