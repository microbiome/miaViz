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
#' @param n A positive integer specifying the number of the most abundant taxa to show. 
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
        if( !(n%%1==0 && n>0) ){
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
        aesthetic <- density_data_list$aesthetic
        
        # Gets the plot from plotter
        plot_out <- .density_plotter(density_data = density_data, 
                                     aesthetic = aesthetic,
                                     abund_values = abund_values, 
                                     colour_by = colour_by)
        return(plot_out)
    }
)

################################ HELP FUNCTIONS ################################

.incorporate_density_data <- function(object, abund_values, n, colour_by){
    # Gets the assay
    abundance_table <- assay(object, abund_values, withDimnames=FALSE)
    # Gets the most abundant taxa
    top_taxa <- getTopTaxa(object, n)
    # Subsets abundance table  by taking taxa of highest abundance
    abundance_table_subset <- abundance_table[top_taxa, , drop=FALSE]
    
    # Converts top taxa names to factor
    top_taxa <- factor(top_taxa, rev(top_taxa))
    
    # Creates a data frame that includes sample names, taxa names as factors,
    # and abundance values
    density_data <- data.frame(
        # Sample data = sample 1 is repeated n(taxa) times, sample 2 repeated n(taxaa) times...
        sample = rep(seq_len(ncol(abundance_table_subset)), each = nrow(abundance_table_subset)), 
        # taxa data = taxa 1, taxa 2,.. taxa n, taxa 1, taxa 2,... taxa n, as many times as there are samples
        taxa = rep(top_taxa, ncol(abundance_table_subset)),
        # Values grom abundance table
        value = as.numeric(abundance_table_subset)
    )
    
    # Gets coloring information if 'colour_by' is not NULL, otherwise specifies just aesthetic
    if (!is.null(colour_by)) {
        # Gets information from colData
        colour_out <- retrieveCellInfo(object, colour_by)
        # Mirrors back the variable name, if a partial match was used
        colour_by <- colour_out$name
        # Adds information to the data frame that includes all density data
        density_data$colour_by <- colour_out$val[density_data$sample]
        # Specifies aesthetic
        aesthetic <- aes_string(y = "taxa", x = "value", colour = "colour_by")
    } else {
        # Specifies aesthetic, when colour_by is not used
        aesthetic <- aes_string(y = "taxa", x = "value")
    }
    
    # Creates a list that includes the data and aesthetic
    density_data_list <- list(density_data = density_data, 
                              aesthetic = aesthetic)
    return(density_data_list)
}

.density_plotter <- function(density_data, 
                             aesthetic,
                             abund_values,
                             colour_by,
                             point_alpha = 0.6,
                             point_shape = 124,
                             text_size = 8,
                             text_colour = "gray35"){
    # Creates the plot and annotations. 
    plot_out <- ggplot(density_data, aesthetic) + geom_point(alpha = point_alpha, shape = point_shape)
    plot_out <- plot_out + xlab(abund_values) + ylab("Feature")  + theme_bw(text_size) +
        theme(axis.text.x = element_text(colour = text_colour),
              axis.text.y = element_text(colour = text_colour),
              axis.title.x = element_text(colour = text_colour),
              axis.title.y = element_text(colour = text_colour),
              title = element_text(colour = text_colour))
    
    if (!is.null(density_data$colour_by)) {
        plot_out <- .resolve_plot_colours(plot_out, density_data$colour_by, colour_by)
    }
    return(plot_out)
}