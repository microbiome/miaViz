#' Plot abundance density
#'
#' This function plots abundance of the most abundant taxa.
#'
#' @inheritParams plotAbundance
#'
#' @param layout \code{Character scalar}. Selects the layout of the plot.
#'   There are three different options: \code{jitter}, \code{density}, and
#'   \code{point} plot. (default: \code{layout = "jitter"})
#'
#' @param n \code{Integer scalar}. Specifies the number of the most abundant
#' taxa to show. (Default: \code{min(nrow(x), 25L)})
#'
#' @param colour.by \code{Character scalar}. Defines a column from
#'   \code{colData}, that is used to color plot. Must be a value of
#'   \code{colData()} function. (Default: \code{NULL})
#'
#' @param colour_by Deprecated. Use \code{colour.by} instead.
#'
#' @param shape.by \code{Character scalar}. Defines a column from
#'   \code{colData}, that is used to group observations to different point shape
#'   groups. Must be a value of \code{colData()} function. \code{shape.by} is
#'   disabled when \code{layout = "density"}. (Default: \code{NULL})
#'
#' @param shape_by Deprecated. Use \code{shape.by} instead.
#'
#' @param size.by \code{Character scalar}. Defines a column from
#'   \code{colData}, that is used to group observations to different point size
#'   groups. Must be a value of \code{colData()} function. \code{size.by} is
#'   disabled when \code{layout = "density"}. (Default: \code{NULL})
#'
#' @param size_by Deprecated. Use \code{size.by} instead.
#'
#' @param decreasing \code{Logical scalar}. Indicates whether the results should
#' be ordered in a descending order or not. If \code{NA} is given the order
#'   as found in \code{x} for the \code{n} most abundant taxa is used.
#'   (Default: \code{TRUE})
#'
#' @param order_descending Deprecated. Use \code{order.descending} instead.
#'
#' @param ... additional parameters for plotting.
#' \itemize{
#'   \item \code{xlab} \code{Character scalar}. Selects the x-axis label.
#'   (Default: \code{assay.type})
#'
#'   \item \code{ylab} \code{Character scalar}. Selects the y-axis label.
#'   \code{ylab} is disabled when \code{layout = "density"}.
#'   (Default: \code{"Taxa"})
#'
#'   \item \code{point.alpha} \code{Numeric scalar}. From range 0 to 1. Selects
#'   the transparency of
#'   colour in \code{jitter} and \code{point} plot. (Default: \code{0.6})
#'
#'   \item \code{point.shape} \code{Positive integer scalar}. Value selecting
#'   the shape of point in
#'   \code{jitter} and \code{point} plot. (Default: \code{21})
#'
#'   \item \code{point.size} \code{Positive integer scalar}. Selects the size of
#'   point in
#'   \code{jitter} and \code{point} plot. (Default: \code{2})
#'
#'   \item \code{add_legend} \code{Logical scalar}. Determines if legend is
#'   added. (Default: \code{TRUE})
#'
#'   \item \code{flipped}: \code{Logical scalar}. Determines if the orientation
#'   of plot is changed so that x-axis and y-axis are swapped.
#'   (Default: \code{FALSE})
#'
#'   \item \code{add_x_text} \code{Logical scalar}. Determines if text that
#'   represents values is included in x-axis. (Default: \code{TRUE})
#' }
#' See \code{\link{mia-plot-args}} for more details i.e. call
#' \code{help("mia-plot-args")}
#'
#' @details
#' This function plots abundance of the most abundant taxa. Abundance can be
#' plotted as a jitter plot, a density plot, or a point plot. By default, x-axis
#' represents abundance and y-axis taxa. In a jitter and point plot, each point
#' represents abundance of individual taxa in individual sample. Most common
#' abundances are shown as a higher density.
#'
#' A density plot can be seen as a smoothened bar plot. It visualized
#' distribution of abundances where peaks represent most common abundances.
#'
#' @return
#' A \code{ggplot2} object
#'
#' @name plotAbundanceDensity
#'
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#'
#' @seealso
#' \code{\link[scater:plotExpression]{scater::plotExpression}}
#'
#' @examples
#' data("peerj13075", package = "mia")
#' tse <- peerj13075
#'
#' # Plots the abundances of 25 most abundant taxa. Jitter plot is the default
#' # option.
#' plotAbundanceDensity(tse, assay.type = "counts")
#'
#' # Counts relative abundances
#' tse <- transformAssay(tse, method = "relabundance")
#'
#' # Plots the relative abundance of 10 most abundant taxa.
#' # "nationality" information is used to color the points. X-axis is
#' # log-scaled.
#' plotAbundanceDensity(
#'     tse, layout = "jitter", assay.type = "relabundance", n = 10,
#'     colour.by = "Geographical_location") +
#'     scale_x_log10()
#'
#' # Plots the relative abundance of 10 most abundant taxa as a density plot.
#' # X-axis is log-scaled
#' plotAbundanceDensity(
#'     tse, layout = "density", assay.type = "relabundance", n = 10 ) +
#'     scale_x_log10()
#'
#' # Plots the relative abundance of 10 most abundant taxa as a point plot.
#' # Point shape is changed from default (21) to 41.
#' plotAbundanceDensity(
#'     tse, layout = "point", assay.type = "relabundance", n = 10,
#'     point.shape = 41)
#'
#' # Plots the relative abundance of 10 most abundant taxa as a point plot.
#' # In addition to colour, groups can be visualized by size and shape in point
#' # plots, and adjusted for point size
#' plotAbundanceDensity(
#'     tse, layout = "point", assay.type = "relabundance", n = 10,
#'     shape.by = "Geographical_location", size.by = "Age", point.size=1)
#'
#' # Ordering via decreasing
#' plotAbundanceDensity(
#'     tse, assay.type = "relabundance", decreasing = FALSE)
#'
#' # for custom ordering set decreasing = NA and order the input object
#' # to your wishes
#' plotAbundanceDensity(
#'     tse, assay.type = "relabundance", decreasing = NA)
#'
#' # Box plots and violin plots are supported by scater::plotExpression.
#' # Plots the relative abundance of 5 most abundant taxa as a violin plot.
#' library(scater)
#' top <- getTop(tse, top = 5)
#' plotExpression(tse, features = top, assay.type = "relabundance") +
#'     ggplot2::coord_flip()
#'
#' # Plots the relative abundance of 5 most abundant taxa as a box plot.
#' plotExpression(tse, features = top, assay.type = "relabundance",
#'     show_violin = FALSE, show_box = TRUE) + ggplot2::coord_flip()
#'
NULL

#' @rdname plotAbundanceDensity
#' @export
setMethod("plotAbundanceDensity", signature = c(x = "SummarizedExperiment"),
    function(
        x,
        layout = c("jitter", "density", "point"),
        assay.type = assay_name, assay_name = "counts",
        n = min(nrow(x), 25L), colour.by = colour_by,
        colour_by = NULL,
        shape.by = shape_by,
        shape_by = NULL,
        size.by = size_by,
        size_by = NULL,
        decreasing = order_descending,
        order_descending = TRUE,
        ...){
        ############################# Input Check ##############################
        # Check layout
        layout <- match.arg(layout, c("jitter", "density", "point"))
        # Checks assay.type
        .check_assay_present(assay.type, x)
        # Checks n
        if( !(length(n)==1 && is.numeric(n) && n%%1==0 && n>0) ){
            stop("'n' must be a positive integer.", call. = FALSE)
        }
        # Checks colour.by
        if( !is.null(colour.by) && !.is_a_string(colour.by)){
            stop("'colour.by' must be a single character value.", call. = FALSE)
        }
        # Checks shape.by
        if( !is.null(shape.by) && !.is_a_string(shape.by)){
            stop("'shape.by' must be a single character value.", call. = FALSE)
        }
        # Checks shape.by
        if( !is.null(shape.by) && !.is_a_string(shape.by)){
            stop("'shape.by' must be a single character value.", call. = FALSE)
        }
        # Checks decreasing
        if( !is.na(decreasing) && !.is_a_bool(decreasing)){
            stop("'decreasing' must be TRUE, FALSE or NA.", call. = FALSE)
        }
        ########################### Input Check end ############################
        # Gets data that will be plotted. Gets a list
        density_data_list <- .incorporate_density_data(
            object = x,
            assay.type = assay.type,
            n = n,
            colour_by = colour.by,
            shape_by = shape.by,
            size_by = size.by,
            order_descending = decreasing)
        # Extracts the density data and aesthetic from the list
        density_data <- density_data_list$density_data
        colour_by <- density_data_list$colour_by
        shape_by <- density_data_list$shape_by
        size_by <- density_data_list$size_by

        # Gets the plot from plotter
        plot_out <- .density_plotter(
            density_data = density_data,
            layout = layout,
            xlab = assay.type,
            colour_by = colour_by,
            shape_by = shape_by,
            size_by = size_by,
            ...)
        return(plot_out)
    }
)

################################ HELP FUNCTIONS ################################

.incorporate_density_data <- function(
        object, assay.type, n, colour_by, shape_by, size_by,
        order_descending = TRUE){
    # Gets the assay
    mat <- assay(object, assay.type, withDimnames = TRUE)
    # Gets the most abundant taxa
    top_taxa <- getTop(object, top = n, assay.type = assay.type)
    # Subsets abundance table  by taking taxa of highest abundance
    mat <- mat[top_taxa, , drop=FALSE]
    # enable conversion to data.frame for non-matrix assays, e.g. sparseMatrices
    if(!is.matrix(mat)){
        mat <- as.matrix(mat)
    }
    # melt the data
    density_data <- t(mat) %>%
        as.data.frame() %>%
        rownames_to_column("Sample")
    # Gets coloring information if 'colour_by' is not NULL
    if (!is.null(colour_by)) {
        # Gets information from colData
        colour_out <- retrieveCellInfo(object, colour_by)
        # Mirrors back the variable name, if a partial match was used
        colour_by <- colour_out$name
        # Adds information to the data frame that includes all density data
        density_data$colour_by <- colour_out$value
    }
    # Gets shape information if 'shape_by' is not NULL
    if (!is.null(shape_by)) {
        shape_out <- retrieveCellInfo(object, shape_by)
        shape_by <- shape_out$name
        density_data$shape_by <- shape_out$value
    }
    # Gets size information if 'size_by' is not NULL
    if (!is.null(size_by)) {
        size_out <- retrieveCellInfo(object, size_by)
        size_by <- size_out$name
        density_data$size_by <- size_out$value
    }
    cols <- c("Sample","colour_by","shape_by","size_by")
    cols <- cols[cols %in% colnames(density_data)]
    density_data <- density_data %>%
        pivot_longer(cols = !cols, names_to = "Y", values_to = "X")
    # Converts taxa to factor. Order of levels is the opposite than in
    # 'top_taxa' so that taxa with highest abundance is on top
    if(is.na(order_descending)){
        lvls <- rownames(object[rownames(object) %in% top_taxa,])
    } else if(order_descending) {
        lvls <- rev(top_taxa)
    } else {
        lvls <- top_taxa
    }
    density_data$Y <- factor( density_data$Y, lvls )
    res <- list(
        density_data = density_data,
        colour_by = colour_by,
        shape_by = shape_by,
        size_by = size_by)
    return(res)
}

.density_plotter <- function(
        density_data,
        layout,
        add_legend = TRUE,
        xlab,
        ylab = NULL,
        colour_by = NULL,
        shape_by = NULL,
        size_by = NULL,
        point_shape = point.shape,
        point.shape = 21,
        point_size = point.size,
        point.size = 2,
        point_alpha = point.alpha,
        point.alpha = 0.6,
        point_colour = point.colour,
        point.colour = "grey70",
        flipped = FALSE,
        scales_free = scales.free,
        scales.free = TRUE,
        angle_x_text = angle.x.text,
        angle.x.text = TRUE){
    # start plotting
    plot_out <- ggplot(density_data, aes(x=.data[["X"]])) +
        xlab(xlab) +
        ylab(ylab)
    # Layout can be "density", "jitter", or "point"
    if (layout == "density"){
        plot_out$data$Y <- factor(
            plot_out$data$Y, levels = rev(levels(plot_out$data$Y)) )
        point_args <- .get_density_args(colour_by, alpha = point_alpha)
        # density specific options for flipping
        grid_args <- list(
            switch = ifelse(flipped, "x", "y"),
            scales = ifelse(scales_free, "free", "fixed"))
        if(flipped){
            grid_args$cols <- vars(!!sym("Y"))
        } else {
            grid_args$rows <- vars(!!sym("Y"))
        }
        #
        plot_out <- plot_out +
            do.call(geom_density, point_args$args) +
            do.call(facet_grid, grid_args)
        shape_by <- NULL
        size_by <- NULL
        angle_x_text <- FALSE
    } else if (layout %in% c("point","jitter")) {
        point_args <- .get_point_args(
            colour_by,
            shape_by = shape_by,
            size_by = size_by,
            alpha = point_alpha,
            shape = point_shape,
            size = point_size,
            colour = point_colour)
        point_args$args$mapping$y <- sym("Y")
        if (layout == "point"){
            plot_out <- plot_out +
                do.call(geom_point, point_args$args)
        } else if (layout == "jitter") {
            point_args$args$height <- 0.25
            plot_out <- plot_out +
                do.call(geom_jitter, point_args$args)
        } else {
            stop(".")
        }
    } else{
        stop("Unsupported layout option: '",layout,"'.", call. = FALSE)
    }
    # If colour_by is specified, colours are resolved
    if (!is.null(colour_by)) {
        plot_out <- .resolve_plot_colours(
            plot_out,
            density_data$colour_by,
            colour_by,
            fill = point_args$fill,
            na.translate = FALSE)
        if(layout == "density"){
            plot_out <- .resolve_plot_colours(
                plot_out,
                density_data$colour_by,
                colour_by,
                fill = !point_args$fill,
                na.translate = FALSE)
        }
    }
    # set the theme
    plot_out <- plot_out +
        theme_classic()
    # fine tuning for density layout
    if( layout == "density" ){
        plot_out <- plot_out +
            theme(strip.background = element_blank())
        if(flipped){
            plot_out <- plot_out +
                # Removes label grid, horizontal labels
                theme(strip.text.x.bottom = element_text(angle = 90, hjust = 1),
                    axis.ticks.x = element_blank(),
                    axis.text.x = element_blank(), # Removes x-axis
                    axis.title.x = element_blank(),
                    axis.line.x = element_blank()) # Removes x-axis
        } else {
            plot_out <- plot_out +
                # Removes label grid, horizontal labels
                theme(strip.text.y.left = element_text(angle = 0, hjust = 1),
                    axis.ticks.y = element_blank(),
                    axis.text.y = element_blank(), # Removes y-axis
                    axis.title.y = element_blank(),
                    axis.line.y = element_blank()) # Removes y-axis
        }
    }
    # add additional guides
    plot_out <- .add_extra_guide(plot_out, shape_by, size_by)
    # add legend
    plot_out <- .add_legend(plot_out, add_legend)
    # flip
    plot_out <- .flip_plot(
        plot_out, flipped = flipped, add_x_text = TRUE,
        angle_x_text = angle_x_text)
    return(plot_out)
}
