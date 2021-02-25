#' Plot Series
#'
#' This function plots series data.
#'
#' @param object
#' A \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#'
#' @param abund_values
#' A single character value for selecting the
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}}
#' to be plotted.
#'
#' @param X
#' A single character value for selecting the column from 
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{ColData}} 
#' that will specify values of x-axis. 
#'  
#' @param Y
#' A single character value for selecting the taxa from 
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{rownames}}. 
#' This parameter specifies taxa whose abundances will be plotted.
#'  
#' @param rank a single character value defining a taxonomic rank, that is used to
#' agglomerate the data. Must be a value of \code{taxonomicRanks()} function.
#'  
#' @param colour_by
#' a single character value defining a taxonomic rank, that is used to color plot. 
#' Must be a value of \code{taxonomicRanks()} function.
#' 
#' @param linetype_by
#' a single character value defining a taxonomic rank, that is used to divide taxa to
#' different line types. Must be a value of \code{taxonomicRanks()} function.
#'
#' @details
#' This function creates series plot, where x-axis includes e.g. time points, and
#' y-axis abundances of selected taxa.
#'
#' @return
#' A \code{ggplot2} object 
#'
#' @name plotSeries
#' @export
#'
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#'
#' @examples
#' x <- microbiomeDataSets::SilvermanAGutData()
#' # Plots 2 most abudant taxa, which are colore by their family
#' plotSeries(x, abund_values = "counts", X = "DAY_ORDER", Y = mia::getTopTaxa(x, 2), colour_by = "Family")
#' 
#' # Counts relative abundances
#' x <- mia::transformCounts(x, method = "relabundance")
#' 
#' # Plots relative abundances of phylums
#' plotSeries(x, abund_values = "relabundance", X = "DAY_ORDER", rank = "Phylum", colour_by = "Phylum")
#' 
#'
NULL

#' @rdname plotSeries
#' @export
setGeneric("plotSeries", signature = c("object"),
           function(object,
                    abund_values,
                    X,
                    Y = NULL,
                    rank = NULL,
                    colour_by = NULL,
                    size_by = NULL,
                    linetype_by = NULL)
               standardGeneric("plotSeries"))


#' @rdname plotSeries
#' @export
setMethod("plotSeries", signature = c(object = "SummarizedExperiment"),
          function(object,
                   abund_values,
                   X,
                   Y = NULL,
                   rank = NULL,
                   colour_by = NULL,
                   size_by = NULL,
                   linetype_by = NULL){
              
              ###################### Input check #######################
              # Check abund_values
              .check_abund_values(abund_values, object)
              
              # Check X
              if( !(X %in% names(colData(object))) ){
                  stop("'X' must be a name of column of colData(object)", call. = FALSE)
              }
              
              # If rank is not null, data will be agglomerated by rank
              if( !is.null(rank) ){
                  # Check rank
                  .check_taxonomic_rank(rank, object)
                  
                  # Agglomerates the object
                  object <- agglomerateByRank(object, rank = rank)
                  
              }
              
              # Check Y
              # If Y is not null, user has specified it
              if (!is.null(Y)){
                  if(!all( is.element(Y, rownames(object)) ) ){
                      stop("'Y' must be in rownames(x). If 'rank' was used,
                           check that 'Y' matches agglomerated data.", call. = FALSE)
                  }
                  # Select taxa that user has specified
                  object <- object[Y]
              }
              
              # Check colour_by
              if( !is.null(colour_by) ){
                  .check_taxonomic_rank(colour_by, object)
              }
              
              # Check linetype_by
              if( !is.null(linetype_by) ){
                  .check_taxonomic_rank(linetype_by, object)
              }
              
              # Get assay data
              assay <- .get_assay_data(object, abund_values)
              
              # Fetch series and features data as a list. 
              series_and_features_data <- .incorporate_series_vis(object, X, colour_by, linetype_by, size_by)
              # Divides it to series and feature data
              series_data <- series_and_features_data$series_data$value
              feature_data <- series_and_features_data$feature_data
              
              # Melts the data
              melted_data <- .melt_series_data(assay, series_data, feature_data)
              
              # Creates variables for series_plotter
              plot_data <- melted_data
              colour_by_title <- colour_by
              colour_by <- melted_data$colour_by
              linetype_by_title <- linetype_by
              linetype_by <- melted_data$linetype_by
              xlab <- paste0(X)
              ylab <- paste0(abund_values)
              
              # Plots the data
              plot <- .series_plotter(plot_data, 
                                      xlab = xlab,
                                      ylab = ylab,
                                      colour_by_title = colour_by_title,
                                      colour_by = colour_by,
                                      linetype_by_title = linetype_by_title,
                                      linetype_by = linetype_by)
              
              return(plot)
              
          }
)

################## HELP FUNCTIONS ##########################

.get_assay_data <- function(object, abund_values){
    
    # Gets warning or error if too many taxa are selected. 
    if( length(rownames(object)) > 20 ){
        stop("Over 20 taxa selected. 20 or under allowed.", call. = FALSE)
    } else if ( length(rownames(object)) > 10 ){
        warning("Over 10 taxa selected.", call. = FALSE)
    }
    
    # Retrieves the assay
    assay <- assay(object, abund_values)
    
    return(assay)
}

.incorporate_series_vis <- function(object, X, colour_by, linetype_by, size_by){
    
    # Stores the variables
    variables <- c(X = X, 
                   colour_by = colour_by,
                   linetype_by = linetype_by,
                   size_by = size_by)
    
    # If variables there are variables
    if(!is.null(variables)){
        
        # Loops through variables
        for(i in seq_along(variables)){
            
            # If the variable is in colData, it is the series data
            if( variables[i] %in% names(colData(object)) ){
                
                # Retrieve series x-axis points from colData
                series_data <- retrieveCellInfo(object, variables[i], search = "colData")
                # mirror back variable name, if a partial match was used
                series_data$name <- names(variables)[i]
            }
            else{
                # get data from rowData
                feature_info <- retrieveFeatureInfo(object, variables[i], search = "rowData")
                # mirror back variable name, if a partial match was used
                feature_info$name <- names(variables)[i]

                # If feature_data dataframe does not exist, create one
                if(!exists("feature_data")){
                    
                    # Store values
                    feature_data <- data.frame(feature_info$value)
                    # Name the column by parameter name, like "colour_by"
                    names(feature_data)[names(feature_data) == "feature_info.value"] <- feature_info$name
                }
                else{
                    # Store values to data frame that already exist
                    feature_data <- cbind(feature_data, feature_info$value)
                    # Name the column by parameter name, like "colour_by"
                    names(feature_data)[names(feature_data) == "feature_info$value"] <- feature_info$name
                }
            }
        }
    }

    # If feature_data exists, add feature_data in addition to series_data
    if(exists("feature_data")){
        returned_list <- list(series_data = series_data, feature_data = feature_data)
    }# If it does not exist, just add the series_data
    else{
        returned_list <- list(series_data = series_data)
    }
    return(returned_list)
}

.melt_series_data <- function(assay, series_data, feature_data){
    
    # Melt assay table 
    melted_data <- as.data.frame(assay) %>% pivot_longer(colnames(assay), names_to = "sample", values_to = "Y")
    
    # Add series data to the data frame
    melted_data <- cbind(melted_data, X = series_data)
    
    # If feature_data is not null
    if( !is.null(feature_data) ){
        # Loops through feature_data
        for( i in 1:ncol(feature_data) ){
            # Assigns feature data to data points
            melted_data <- cbind(melted_data, feature = rep(feature_data[[i]], nrow(melted_data)/length(feature_data[[i]])))
            # Renames the column
            names(melted_data)[names(melted_data) == "feature"] <- names(feature_data)[i]
        }
    }
    
    return(melted_data)
}

.series_plotter <- function(plot_data,
                            xlab = NULL,
                            ylab = NULL,
                            colour_by = NULL,
                            colour_by_title = NULL,
                            linetype_by = NULL,
                            linetype_by_title = NULL,
                            add_legend = TRUE,
                            point_alpha = 1,
                            point_size = 2,
                            line_alpha = 1,
                            line_size = 1,
                            ...){
    
    # Creates a "draft" of a plot
    plot_out <- ggplot(plot_data, aes_string(x = "X", y = "Y")) +
        labs(x = xlab, y = ylab)
    
    # Fetch arguments of line
    line_args <- .get_line_args(colour_by = colour_by,
                                linetype_by = linetype_by,
                                size_by = NULL,
                                alpha = line_alpha,
                                size = line_size)
    
    line_args$args$mapping$group <- sym("colour_by")
    
    # Adds arguments to the plot
    plot_out <- plot_out +
        do.call(geom_line, line_args$args)
    
    # resolve the colours
    plot_out <- .resolve_plot_colours(plot_out,
                                      plot_data$colour_by,
                                      colour_by,
                                      fill = TRUE)
    plot_out <- .resolve_plot_colours(plot_out,
                                      plot_data$colour_by,
                                      colour_by,
                                      fill = FALSE)
    
    # Changes theme
    plot_out <- plot_out +
        theme_classic()
    # Adds legend
    plot_out <- .add_legend(plot_out, add_legend)
    
    return(plot_out)
    
}
