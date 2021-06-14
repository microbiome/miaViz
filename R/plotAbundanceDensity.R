#' Plot abundance density
#'
#' This function plots abundance of the most abundant taxa. 
#'
#' @param object a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'    object.
#'
#' @param abund_values a single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}} to be
#'   plotted. (default: \code{abund_values = "counts"})
#'
#' @param n a positive integer specifying the number of the most abundant taxa to show. 
#'  
#' @param colour_by a single character value defining a column from \code{colData}, that is used to
#'   color plot. Must be a value of \code{colData()} function.
#'   
#' @param ... additional parameters for plotting. See 
#'   \code{\link{mia-plot-args}} for more details
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
                    abund_values = "counts",
                    n = 50,
                    colour_by = NULL, 
                    ...)
             standardGeneric("plotAbundanceDensity"))

#' @rdname plotAbundanceDensity
#' @export
setMethod("plotAbundanceDensity", signature = c(object = "SummarizedExperiment"),
    function(object,
             abund_values = "counts",
             n = 50, 
             colour_by = NULL, 
             ...){
        ############################# Input Check ##############################
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
                             xlab = NULL,
                             ylab = NULL,
                             colour_by = NULL,
                             layout = "point",
                             point_alpha = 0.6,
                             point_shape = 124,
                             point_size = 2,
                             text_size = 8,
                             add_legend = TRUE,
                             flipped = FALSE,
                             add_x_text = FALSE){
    # start plotting
    plot_out <- ggplot(density_data, aes_string(x="X", y="Y")) +
        xlab(xlab) +
        ylab(ylab)
    
    # Layout "density" or "point", "density" will be added
    if (layout == "point"){
        density_out <- .get_point_args(colour_by,
                                       shape_by = NULL,
                                       size_by = NULL,
                                       alpha = point_alpha,
                                       shape = point_shape,
                                       size = point_size)
        plot_out <- plot_out +
            do.call(geom_point, density_out$args)
        # resolve the colours
        plot_out <- .resolve_plot_colours(plot_out,
                                          density_data$colour_by,
                                          colour_by,
                                          fill = FALSE,
                                          na.translate = FALSE)
    }
    plot_out <- plot_out +
        theme_classic()
    # add legend
    plot_out <- .add_legend(plot_out, add_legend)
    # flip
    plot_out <- .flip_plot(plot_out, flipped, add_x_text)
    return(plot_out)
}
