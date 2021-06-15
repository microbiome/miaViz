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
#' (default: \code{n = 25})
#'  
#' @param colour_by a single character value defining a column from \code{colData}, that is used to
#' color plot. Must be a value of \code{colData()} function. (default: \code{colour_by = NULL})
#' 
#' @param shape_by a single character value defining a column from \code{colData}, that is used to
#' group observations to different point shape groups. Must be a value of \code{colData()} function.
#' \code{shape_by} is disabled when \code{layout = "density"}. (default: \code{shape_by = NULL})
#' 
#' @param size_by a single character value defining a column from \code{colData}, that is used to
#' group observations to different point size groups. Must be a value of \code{colData()} function.
#' \code{size_by} is disabled when \code{layout = "density"}. (default: \code{size_by = NULL})
#'   
#' @param ... additional parameters for plotting. 
#' \itemize{
#'   \item{xlab}{ a single character value for selecting the x-axis label. 
#'   (default: \code{xlab = abund_values}) }
#'   
#'   \item{ylab}{ a single character value for selecting the y-axis label. 
#'   \code{ylab} is disabled when \code{layout = "density"}. 
#'   (default: \code{ylab = "Taxa")} }
#'   
#'   \item{point_alpha}{ a numeric value from range 0 to 1 selecting the transparency of colour in
#'   \code{jitter} and \code{point} plot. (default: \code{point_alpha = 0.6}) }
#'   
#'   \item{point_shape}{ a positive integer value selecting the shape of point in
#'   \code{jitter} and \code{point} plot. (default: \code{point_shape = 21}) }
#'   
#'   \item{point_size}{ a positive numeric value selecting the size of point in
#'   \code{jitter} and \code{point} plot. (default: \code{point_size = 2}) }
#'   
#'   \item{add_legend}{ a boolean value selecting if legend is added. 
#'   (default: \code{add_legend = TRUE}) }
#'   
#'   \item{flipped}{ a boolean value selecting if the orientation of plot is changed 
#'   so that x-axis and y-axis are swapped. (default \code{flipped = FALSE}) }
#'   
#'   \item{add_x_text}{ a boolean value selecting if text that represents values is included in x-axis. 
#'   (default: \code{add_x_text = TRUE}) }
#' }
#' See \code{\link{mia-plot-args}} for more details
#'
#' @details
#' This function plots abundance of the most abundant taxa. Abundance can be plotted as
#' a jitter plot, a density plot, or a point plot. By default, x-axis represents abundance 
#' and y-axis taxa. In a jitter and point plot, each point represents abundance of individual taxa 
#' in individual sample. Most common abundances are shown as a higher density. 
#' 
#' A density plot can be seen as a smoothened bar plot. It visualized distribution of 
#' abundances where peaks represent most common abundances.
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
#' 
#' # Plots the abundances of 25 most abundant taxa. Jitter plot is the default option.
#' plotAbundanceDensity(tse, abund_values = "counts")
#' 
#' # Counts relative abundances
#' tse <- transformSamples(tse, method = "relabundance")
#' 
#' # Plots the relative abundance of 10 most abundant taxa. 
#' # "nationality" information is used to color the points. X-axis is log-scaled.
#' plotAbundanceDensity(tse, layout = "jitter", abund_values = "relabundance", 
#'                      n = 10, colour_by = "nationality") +
#'     scale_x_log10() 
#'                      
#' # Plots the relative abundance of 10 most abundant taxa as a density plot.
#' # X-axis is log-scaled
#' plotAbundanceDensity(tse, layout = "density", abund_values = "relabundance",
#'                      n = 10 ) +
#'     scale_x_log10()
#'                      
#' # Plots the relative abundance of 10 most abundant taxa as a point plot.
#' # Point shape is changed from default (21) to 41.
#' plotAbundanceDensity(tse, layout = "point", abund_values = "relabundance", n = 10,
#'                      point_shape = 41)
#'                      
#' # Plots the relative abundance of 10 most abundant taxa as a point plot.
#' # In addition to colour, groups can be visualized by size and sahep in point plots.
#' plotAbundanceDensity(tse, layout = "point", abund_values = "relabundance", n = 10,
#'                      shape_by = "sex", size_by = "time")
#' 
NULL

#' @rdname plotAbundanceDensity
#' @export
setGeneric("plotAbundanceDensity", signature = c("object"),
           function(object, ...)
             standardGeneric("plotAbundanceDensity"))

#' @rdname plotAbundanceDensity
#' @export
setMethod("plotAbundanceDensity", signature = c(object = "SummarizedExperiment"),
    function(object,
             layout = c("jitter", "density", "point"),
             abund_values = "counts",
             n = 25, 
             colour_by = NULL, 
             shape_by = NULL, 
             size_by = NULL, 
             ...){
        ############################# Input Check ##############################
        # Check layout
        layout <- match.arg(layout, c("jitter", "density", "point"))
        # Checks abund_values
        .check_assay_present(abund_values, object)
        # Checks n
        if( !(length(n)==1 && is.numeric(n) && n%%1==0 && n>0) ){
            stop("'n' must be a positive integer.", call. = FALSE)
        }
        # Checks colour_by
        if( !is.null(colour_by) && !.is_a_string(colour_by)){
            stop("'colour_by' must be a single character value.",
                 call. = FALSE)
        }
        # Checks shape_by
        if( !is.null(shape_by) && !.is_a_string(shape_by)){
            stop("'shape_by' must be a single character value.",
                 call. = FALSE)
        }
        # Checks size_by
        if( !is.null(size_by) && !.is_a_string(size_by)){
            stop("'size_by' must be a single character value.",
                 call. = FALSE)
        }
        ########################### Input Check end ############################
        # Gets data that will be plotted. Gets a list
        density_data_list <- .incorporate_density_data(object = object,
                                                       abund_values = abund_values,
                                                       n = n,
                                                       colour_by = colour_by,
                                                       shape_by = shape_by,
                                                       size_by = size_by)
        # Extracts the density data and aesthetic from the list
        density_data <- density_data_list$density_data
        colour_by <- density_data_list$colour_by
        shape_by <- density_data_list$shape_by
        size_by <- density_data_list$size_by
        
        # Gets the plot from plotter
        plot_out <- .density_plotter(density_data = density_data, 
                                     layout = layout,
                                     xlab = abund_values,
                                     colour_by = colour_by,
                                     shape_by = shape_by,
                                     size_by = size_by,
                                     ...)
        return(plot_out)
    }
)

################################ HELP FUNCTIONS ################################

.incorporate_density_data <- function(object, abund_values, n,
                                      colour_by,
                                      shape_by,
                                      size_by){
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
    # Converts taxa to factor. Order of levels is the opposite than in 'top_taxa'
    # so that taxa with highest abundance is on top
    density_data$Y <- factor( density_data$Y, rev(top_taxa) )
    return(list(density_data = density_data, 
                colour_by = colour_by,
                shape_by = shape_by,
                size_by = size_by))
}

.density_plotter <- function(density_data, 
                             layout,
                             add_legend = TRUE,
                             xlab,
                             ylab = NULL,
                             colour_by = NULL,
                             shape_by = NULL,
                             size_by = NULL,
                             point_shape = 21,
                             point_size = 2,
                             alpha = 0.6,
                             flipped = FALSE,
                             scales_free = TRUE,
                             angle_x_text = TRUE){
    # start plotting
    plot_out <- ggplot(density_data, aes_string(x="X")) +
        xlab(xlab) +
        ylab(ylab)
    # Layout can be "density", "jitter", or "point"
    if (layout == "density"){
        plot_out$data$Y <- factor(plot_out$data$Y,
                                  levels = rev(levels(plot_out$data$Y)) )
        point_args <- .get_density_args(colour_by,
                                        alpha = alpha)
        # density specific options for flipping
        grid_args <- list(switch = ifelse(flipped, "x", "y"),
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
        point_args <- .get_point_args(colour_by,
                                      shape_by = shape_by,
                                      size_by = size_by,
                                      alpha = alpha,
                                      shape = point_shape,
                                      size = point_size)
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
    # If colour_by is specified, colours are resolve
    if (!is.null(colour_by)) {
        plot_out <- .resolve_plot_colours(plot_out,
                                          density_data$colour_by,
                                          colour_by,
                                          fill = point_args$fill,
                                          na.translate = FALSE)
        if(layout == "density"){
            plot_out <- .resolve_plot_colours(plot_out,
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
                theme(strip.text.x.bottom = element_text(angle = 90, hjust = 1), # Removes label grid, horizontal labels
                      axis.ticks.x = element_blank(),
                      axis.text.x = element_blank(), # Removes x-axis
                      axis.title.x = element_blank(),
                      axis.line.x = element_blank()) # Removes x-axis
        } else {
            plot_out <- plot_out +
                theme(strip.text.y.left = element_text(angle = 0, hjust = 1), # Removes label grid, horizontal labels
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
    plot_out <- .flip_plot(plot_out,
                           flipped = flipped,
                           add_x_text = TRUE,
                           angle_x_text = angle_x_text)
    return(plot_out)
}
