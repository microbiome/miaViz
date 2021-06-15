#' Plot abundance density
#'
#' This function plots abundance of the most abundant taxa. 
#'
#' @param object a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#' 
#' @param layout a single character value for selecting the layout of the plot. 
#' There are three different options: \code{jitter}, \code{density}, and \code{point} plot.
#' (default: \code{layout = "jitter"})
#'
#' @param abund_values a single character value for selecting the
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}} to be
#' plotted. (default: \code{abund_values = "counts"})
#'
#' @param n a positive integer specifying the number of the most abundant taxa to show.
#' (default: \code{n = 50})
#'  
#' @param colour_by a single character value defining a column from \code{colData}, that is used to
#' color plot. Must be a value of \code{colData()} function. \code{colour_by} is disabled when
#' \code{layout = "density"}.
#'   
#' @param ... additional parameters for plotting. 
#' \itemize{
#'   \item{xlab}{ a single character value for selecting the x-axis label.}
#'   \item{ylab}{ a single character value for selecting the y-axis label.}
#'   \item{point_alpha}{ a numeric value from range [0,1] selecting the transparency of colour in
#'   \code{jitter} and \code{point} plot. (default \code{point_alpha = 0.6}) }
#'   \item{point_shape}{ a positive integer value selecting the shape of point in
#'   \code{jitter} and \code{point} plot. (default \code{point_shape = 21}) }
#'   \item{point_size}{ a positive numeric value selecting the size of point in
#'   \code{jitter} and \code{point} plot. (default \code{point_size = 2}) }
#'   \item{add_legend}{ a boolean value selecting if legend is added. 
#'   (default \code{add_legend = TRUE}) }
#'   \item{flipped}{ a boolean value selecting if the orientation of plot is changed 
#'   so that x-axis and y-axis are swapped. \code{flipped} is disabled when 
#'   \code{layout = "density"}. (default \code{flipped = FALSE}) }
#'   \item{add_x_text}{ a boolean value selecting if text that represents values is included in x-axis. 
#'   (default \code{add_x_text = TRUE}) }
#'   \item{log_scale_x_axis}{ a boolean value selecting if x-axis is log-scaled. 
#'   (default \code{log_scale_x_axis = FALSE}) }
#' }
#' See \code{\link{mia-plot-args}} for more details
#'
#' @details
#' This function plots abundance of the most abundant taxa. Abundance is plotted as
#' a point plot, where x-axis represents abundance and y-axis taxa. Each point represents
#' abundance of individual taxa in individual sample. Most common abundances are shown
#' as a higher density.
#'
#' @return 
#' A \code{ggplot2} object 
#'
#' @name plotAbundanceDensity
#'
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#'
#' @examples
#' tse <- microbiomeDataSets::atlas1006()
#' # Plots the abundances of 100 most abundant taxa
#' plotAbundanceDensity(tse, abund_values = "counts", n = 100)
#' # Counts relative abundances
#' tse <- transformSamples(tse, method = "relabundance")
#' # Plots the relative abundance of 50 msot abundant taxa. 
#' # "nationality" information is used to color the points. 
#' plotAbundanceDensity(tse, abund_values = "relabundance", colour_by = "nationality")
#' 
NULL

#' @rdname plotAbundanceDensity
#' @export
setGeneric("plotAbundanceDensity", signature = c("object"),
           function(object,
                    layout = c("jitter", "density", "point"),
                    abund_values = "counts",
                    n = 50,
                    colour_by = NULL, 
                    ...)
             standardGeneric("plotAbundanceDensity"))

#' @rdname plotAbundanceDensity
#' @export
setMethod("plotAbundanceDensity", signature = c(object = "SummarizedExperiment"),
    function(object,
             layout = c("jitter", "density", "point"),
             abund_values = "counts",
             n = 50, 
             colour_by = NULL, 
             ...){
        ############################# Input Check ##############################
        # Check layout
        layout <- match.arg(layout, several.ok = FALSE)
        # Checks abund_values
        .check_assay_present(abund_values, object)
        # Checks n
        if( !(length(n)==1 && is.numeric(n) && n%%1==0 && n>0) ){
            stop("'n' must be a positive integer.", call. = FALSE)
        }
        # Checks colour_by
        if( !is.null(colour_by) &&
            (!.is_a_string(colour_by) ||
            !(colour_by %in% names(colData(object)))) ){
            stop("'colour_by' must be a name of column of colData(object)",
                 call. = FALSE)
        }
        ########################### Input Check end ############################
        # Gets data that will be plotted. Gets a list
        density_data_list <- .incorporate_density_data(object = object,
                                                       abund_values = abund_values,
                                                       n = n,
                                                       colour_by = colour_by)
        # Extracts the density data and aesthetic from the list
        density_data <- density_data_list$density_data
        colour_by <- density_data_list$colour_by
        
        # Gets the plot from plotter
        plot_out <- .density_plotter(density_data = density_data, 
                                     layout = layout,
                                     xlab = abund_values,
                                     ylab = "Taxa",
                                     colour_by = colour_by,
                                     ...)
        return(plot_out)
    }
)

################################ HELP FUNCTIONS ################################

.incorporate_density_data <- function(object, abund_values, n, colour_by){
    # Gets the assay
    mat <- assay(object, abund_values)
    # Gets the most abundant taxa
    top_taxa <- getTopTaxa(object, n)
    # Subsets abundance table  by taking taxa of highest abundance
    mat <- mat[top_taxa, , drop=FALSE]
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
    cols <- c("Sample","colour_by")[c("Sample","colour_by") %in% colnames(density_data)]
    density_data <- density_data %>%
        pivot_longer(cols = !cols, names_to = "Y", values_to = "X")
    # Converts taxa to factor. Order of levels is the opposite than in 'top_taxa'
    # so that taxa with highest abundance is on top
    density_data$Y <- factor( density_data$Y, rev(top_taxa) )
    return(list(density_data = density_data, 
                colour_by = colour_by))
}

.density_plotter <- function(density_data, 
                             layout,
                             xlab = NULL,
                             ylab = NULL,
                             colour_by = NULL,
                             point_alpha = 0.6,
                             point_shape = 21,
                             point_size = 2,
                             add_legend = TRUE,
                             flipped = FALSE,
                             add_x_text = TRUE,
                             log_scale_x_axis = FALSE){
    # start plotting
    # Density plot needs different kind of structure
    if (layout == "density"){
        # # Reorders the levels to reverse order so that the taxa with highest abundance
        # is on the top
        density_data$Y <- factor(density_data$Y , levels = rev(levels(density_data$Y)) )
        plot_out <- ggplot(density_data, aes_string(x="X", colour = "Y")) +
            xlab(xlab) +
            ylab(ylab)
    }
    else {
        plot_out <- ggplot(density_data, aes_string(x="X", y="Y")) +
            xlab(xlab) +
            ylab(ylab)
    }
    
    # Layout can be "density", "jitter", or "point"
    if (layout == "point"){
        density_out <- .get_point_args(colour_by,
                                       shape_by = NULL,
                                       size_by = NULL,
                                       alpha = point_alpha,
                                       shape = point_shape,
                                       size = point_size,
                                       position = "identity")
        plot_out <- plot_out +
            do.call(geom_point, density_out$args)
    }
    else if (layout == "density"){
        density_out <- .get_density_args(alpha = 0.65)
        plot_out <- plot_out +
            do.call(geom_density, density_out$args) + 
            facet_grid(Y~., switch = "y", scales="free") + 
            theme_classic() +
            theme(strip.background = element_blank(), strip.text.y.left = element_text(angle = 0), # Removes label grid, horizontal labels
                  axis.ticks.y = element_blank(), axis.text.y = element_blank(), # Removes y-axis
                  axis.title.y = element_blank(), axis.line.y = element_blank()) # Removes y-axis
        # colour_by and flipped are disabled in density_plot
        colour_by <- NULL
        flipped <- FALSE
    }
    else {
        density_out <- .get_point_args(colour_by,
                                       shape_by = NULL,
                                       size_by = NULL,
                                       alpha = point_alpha,
                                       shape = point_shape,
                                       size = point_size,
                                       position = "jitter")
        plot_out <- plot_out +
            do.call(geom_point, density_out$args)
    }
    if (log_scale_x_axis){
        plot_out <- plot_out + scale_x_log10() 
    }
    # If colour_by is specified, colours are added
    if (!is.null(colour_by)) {
        # resolve the colours
        plot_out <- .resolve_plot_colours(plot_out,
                                          density_data$colour_by,
                                          colour_by,
                                          fill = FALSE,
                                          na.translate = FALSE)
    }
    if ( !layout == "density" ){
        plot_out <- plot_out +
            theme_classic()
    }
    # add legend
    plot_out <- .add_legend(plot_out, add_legend)
    # flip
    plot_out <- .flip_plot(plot_out, flipped, add_x_text)
    return(plot_out)
}
