#' Plot Rarefaction Curves
#'
#' Plot rarefaction curves of observed richness 
#' \code{\link[vegan:rarecurve]{vegan::rarecurve}}.
#'
#' @param x a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#'
#' @param assay.type \code{Character scalar}. Specifies which assay to use.
#'   Must contain count data with singletons.
#'  (Default: \code{"counts"})
#'
#' @param colour.by \code{Character scalar}. Specifies a column from
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{colData}} to
#'   colour the rarefaction curves by. (Default: \code{NULL})
#'
#' @param linetype.by \code{Character scalar}. Specifies a column from
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{colData}} to
#'   set the line type of the rarefaction curves by. (Default: \code{NULL})
#'
#' @param size.by \code{Character scalar}. Specifies a column from
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{colData}} to
#'   set the line width of the rarefaction curves by. (Default: \code{NULL})
#'
#' @param ... additional arguments
#' \itemize{
#'   \item passed to \code{\link[vegan:rarecurve]{vegan::rarecurve}}, e.g.
#'   \code{step} and \code{sample}.
#'
#'   \item additional plotting arguments. See \code{\link{mia-plot-args}} for
#'   more details i.e. call \code{help("mia-plot-args")}
#' }
#'
#' @details
#' \code{plotRarefaction} generates rarefaction curves for each sample showing 
#' the number of observed species as a function of read count. Uses 
#' \code{\link[vegan:rarecurve]{vegan::rarecurve}}.
#'
#' The \code{other.fields} argument allows extra
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{colData}}
#' columns to be attached to the underlying plot data, so they are accessible
#' when customizing the returned \code{ggplot2} object.
#'
#' @return
#' A \code{ggplot2} object.
#'
#' @seealso
#' \code{\link[vegan:rarecurve]{vegan::rarecurve}}
#'
#' @name plotRarefaction
#'
#' @examples
#' data(esophagus, package = "mia")
#' tse <- esophagus
#' 
#' # Create colData column for example
#' tse$Sample_ID <- colnames(tse)
#' 
#' # Basic rarefaction curve
#' plotRarefaction(tse, assay.type = "counts")
#' 
#' plotRarefaction(tse, assay.type = "counts", linetype.by = "Sample_ID",
#' colour.by = "Sample_ID", size.by = "Sample_ID", line.width.range = c(0, 5))
#'
NULL

################################################################################
# plotRarefaction

#' @rdname plotRarefaction
#' @importFrom SummarizedExperiment assay colData
#' @export
setMethod("plotRarefaction", signature = c(x = "SummarizedExperiment"),
    function(
            x,
            assay.type = "counts",
            colour.by = NULL,
            linetype.by = NULL,
            size.by = NULL,
            other.fields = NULL,
            ...){
        # Input checks
        .check_assay_present(assay.type, x)
        .check_metadata_variable(x, colour.by, col = TRUE)
        .check_metadata_variable(x, linetype.by, col = TRUE)
        .check_metadata_variable(x, size.by, col = TRUE)
        .check_metadata_variable(x, other.fields, col = TRUE, multiple = TRUE)
        # Get data for plotting
        plot_data <- .get_rarefaction_plot_data(
            x, assay.type, colour.by, linetype.by, size.by, other.fields, ...)
        # Create the plot
        p <- .rarefaction_plotter(plot_data, colour.by, 
                                  linetype.by, size.by, ...)
        return(p)
    }
)

################################################################################
# Helper functions

# Computes rarefaction curves via vegan and merges relevant colData columns
.get_rarefaction_plot_data <- function(
        x, assay.type, colour.by, linetype.by, size.by, other.fields, ...) {
    .require_package("vegan")
    # Get rarefaction data
    curve_data <- vegan::rarecurve(t(assay(x, assay.type)), tidy = TRUE, ...)
    # Merge with colData for styling variables and any other.fields
    fields <- unique(c(colour.by, linetype.by, size.by, other.fields))
    if (length(fields) > 0L) {
        col_data <- as.data.frame(colData(x))[, fields, drop = FALSE]
        col_data[["Site"]] <- rownames(col_data)
        curve_data <- merge(
            curve_data, col_data, by = "Site", all.x = TRUE)
    }
    return(curve_data)
}

# Builds the ggplot2 rarefaction curve plot from the prepared data.frame.
.rarefaction_plotter <- function(
        plot_data,
        colour.by = NULL,
        linetype.by = NULL,
        size.by = NULL,
        add.legend = add_legend,
        add_legend = TRUE,
        line.alpha = line_alpha,
        line_alpha = 1,
        line.type = line_type,
        line_type = NULL,
        line.width = line_width,
        line_width = 0.5,
        line.width.range = line_width_range,
        line_width_range = c(0.5, 3), ...) {
    # Initialise plot
    plot_out <- ggplot(
        plot_data,
        aes(x = .data[["Sample"]], y = .data[["Species"]],
            group = .data[["Site"]])
    )
    # Build line geom arguments
    line_args <- .get_line_args(
        colour_by = if (!is.null(colour.by)) plot_data[[colour.by]],
        linetype_by = if (!is.null(linetype.by)) plot_data[[linetype.by]],
        size_by = if (!is.null(size.by)) plot_data[[size.by]],
        alpha = line.alpha,
        linetype = line.type,
        linewidth = line.width
    )
    plot_out <- plot_out + do.call(geom_line, line_args$args)
    # Resolve colour scale and legend title
    if (!is.null(colour.by)) {
        plot_out <- .resolve_plot_colours(
            plot_out,
            plot_data[[colour.by]],
            colour.by,
            fill = FALSE)
    }
    # Scale line width range if size.by is used
    if (!is.null(size.by)) {
        if (is.numeric(plot_data[[size.by]])) {
            plot_out <- plot_out + scale_size_continuous(
                range = line.width.range)
        } else {
            plot_out <- plot_out + scale_size_discrete(
                range = line.width.range)
        }
    }
    
    plot_out <- .add_extra_line_guide(plot_out,
      linetype_by = linetype.by,
      size_by = size.by)
        
    plot_out <- .add_legend(plot_out, add_legend)
    
    plot_out <- plot_out + theme_classic()
    
    return(plot_out)
}