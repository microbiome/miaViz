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
#' a single character value defining a feature, that is used to color plot. 
#' Must be a value of \code{taxonomicRanks()} function or 
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{ColData}}.
#' 
#' @param linetype_by
#' a single character value defining a taxonomic rank, that is used to divide taxa to
#' different line types. Must be a value of \code{taxonomicRanks()} function.
#' 
#' @param size_by
#' a single character value defining a taxonomic rank, that is used to divide taxa to
#' different line size types. Must be a value of \code{taxonomicRanks()} function.
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
#' plotSeries(x, abund_values = "counts", X = "DAY_ORDER", 
#'            Y = mia::getTopTaxa(x, 2), colour_by = "Family")
#' 
#' # Counts relative abundances
#' x <- mia::transformCounts(x, method = "relabundance")
#' 
#' # Selects taxa
#' taxa <- c("seq_1", "seq_2", "seq_3", "seq_4", "seq_5")
#' 
#' # Plots relative abundances of phylums
#' plotSeries(x[taxa], abund_values = "relabundance", X = "DAY_ORDER", 
#'            colour_by = "Family", linetype_by = "Phylum")
#' 
#' # In addition to 'colour_by' and 'linetype_by', 'size_by' can also be used to group taxa.
#' plotSeries(x, abund_values = "counts", X = "DAY_ORDER", Y = mia::getTopTaxa(x, 5), 
#'            colour_by = "Family", size_by = "Phylum") +
#'            scale_size_discrete(range=c(0.5, 2))
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
              # Checks abund_values
              .check_abund_values(abund_values, object)
              
              # Checks X
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
              
              # Checks Y
              # If Y is not null, user has specified it
              if (!is.null(Y)){
                  if(!all( is.element(Y, rownames(object)) ) ){
                      stop("'Y' must be in rownames(x). 
                      If 'rank' was used, check that 'Y' matches agglomerated data.", 
                           call. = FALSE)
                  }
                  # Select taxa that user has specified
                  object <- object[Y]
              }
              
              # Checks colour_by
              if( !is.null(colour_by) ){
                  if( !(colour_by %in% taxonomyRanks(object) || colour_by %in% names(colData(object))) ){
                      stop("'colour_by' must be value of taxonomyRanks(x) or colData(x).", 
                           call. = FALSE)
                  }
              }
              
              # Checks linetype_by
              if( !is.null(linetype_by) ){
                  if( !(linetype_by %in% taxonomyRanks(object)) ){
                      stop("'linetype_by' must be value of taxonomyRanks(x).", 
                           call. = FALSE)
                  }
              }
              
              # Checks size_by
              if( !is.null(size_by) ){
                  if( !(size_by %in% taxonomyRanks(object)) ){
                      stop("'size_by' must be value of taxonomyRanks(x).", 
                           call. = FALSE)
                  }
              }
              
              # Gets assay data
              assay <- .get_assay_data(object, abund_values)
              
              # Fetches sample and features data as a list. 
              sample_and_features_data <- .incorporate_series_vis(object, X, colour_by, linetype_by, size_by)
              # Divides it to sample and feature data
              sample_data <- sample_and_features_data$sample_data
              feature_data <- sample_and_features_data$feature_data
              # Gets rownames
              rownames <- rownames(object)
              
              # Melts the data
              melted_data <- .melt_series_data(assay, sample_data, feature_data, rownames)
              
              # Creates variables for series_plotter
              plot_data <- melted_data
              colour_by <- melted_data$colour_by
              linetype_by <- melted_data$linetype_by
              size_by <- melted_data$size_by
              xlab <- paste0(X)
              ylab <- paste0(abund_values)
              
              # Plots the data
              plot <- .series_plotter(plot_data, 
                                      xlab = xlab,
                                      ylab = ylab,
                                      colour_by = colour_by,
                                      linetype_by = linetype_by,
                                      size_by = size_by)
              
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
    
    # If there are variables
    if(!is.null(variables)){
        
        # Loops through variables
        for(i in seq_along(variables)){
            
            # If the variable is in colData
            if( variables[i] %in% names(colData(object)) ){
                
                # Retrieves series x-axis points from colData
                sample_info <- retrieveCellInfo(object, variables[i], search = "colData")
                # Mirrors back variable name, if a partial match was used
                sample_info$name <- names(variables)[i]
                
                # If sample_data data frame does not exist, create one
                if(!exists("sample_data")){
                    # Stores values
                    sample_data <- data.frame(sample_info$value)
                    # Names the column by parameter name, like "colour_by"
                    names(sample_data)[names(sample_data) == "sample_info.value"] <- sample_info$name
                } # If sample_data data frame already exists
                else{
                    # Stores values to data frame that already exist
                    sample_data <- cbind(sample_data, sample_info$value)
                    # Names the column by parameter name, like "colour_by"
                    names(sample_data)[names(sample_data) == "sample_info$value"] <- sample_info$name
                }
            }
            else{
                # Gets data from rowData
                feature_info <- retrieveFeatureInfo(object, variables[i], search = "rowData")
                # Mirrors back variable name, if a partial match was used
                feature_info$name <- names(variables)[i]
                
                # If feature_data data frame does not exist, create one
                if(!exists("feature_data")){
                    # Store values
                    feature_data <- data.frame(feature_info$value)
                    # Name the column by parameter name, like "colour_by"
                    names(feature_data)[names(feature_data) == "feature_info.value"] <- feature_info$name
                } # If feature_data data frame already exists
                else{
                    # Store values to data frame that already exist
                    feature_data <- cbind(feature_data, feature_info$value)
                    # Name the column by parameter name, like "colour_by"
                    names(feature_data)[names(feature_data) == "feature_info$value"] <- feature_info$name
                }
            }
        }
    }
    
    # If feature_data exists, add feature_data + sample_data
    if(exists("feature_data")){
        returned_list <- list(sample_data = sample_data, feature_data = feature_data)
    }# If it does not exist, just add the sample_data
    else{
        returned_list <- list(sample_data = sample_data)
    }
    return(returned_list)
}

.melt_series_data <- function(assay, sample_data, feature_data, rownames){
    
    # Melts assay table
    melted_data <- as.data.frame(assay) %>% pivot_longer(colnames(assay), names_to = "sample", values_to = "Y")
    
    # Adds rowname information. Repeats as many times there are samples. Repeats 1st element x times, then 2nd element x times...
    melted_data <- cbind(melted_data, feature = rep(rownames, each = nrow(melted_data)/length(rownames)))
    
    # Loops through sample_data
    for( i in 1:ncol(sample_data) ){
        # Assigns sample data to data points. When there are more sample-taxa combinations than different sample data values,
        # sample data is repeated as many times there are taxa. Repeats whole list x times.
        melted_data <- cbind(melted_data, temp_name = rep(sample_data[[i]], nrow(melted_data)/length(sample_data[[i]])))
        # Renames the column
        names(melted_data)[names(melted_data) == "temp_name"] <- names(sample_data)[i]
    }
    
    # If feature_data is not null
    if( !is.null(feature_data) ){
        # Loops through feature_data
        for( i in 1:ncol(feature_data) ){
            # Assigns feature data to data points. When there are more sample-taxa combinations than feature values,
            # features are repeated as many times there are samples. Repeats 1st element x times, then 2nd element x times...
            melted_data <- cbind(melted_data, temp_name = rep(feature_data[[i]], each = nrow(melted_data)/length(feature_data[[i]])))
            # Renames the column
            names(melted_data)[names(melted_data) == "temp_name"] <- names(feature_data)[i]
        }
    }
    
    return(melted_data)
}

.series_plotter <- function(plot_data,
                            xlab = NULL,
                            ylab = NULL,
                            colour_by = NULL,
                            linetype_by = NULL,
                            size_by = NULL,
                            add_legend = TRUE,
                            point_alpha = 1,
                            point_size = 2,
                            line_alpha = 1,
                            line_type = NULL,
                            line_size = 1,
                            ...){
    
    # Creates a "draft" of a plot
    plot_out <- ggplot(plot_data, aes_string(x = "X", y = "Y")) +
        labs(x = xlab, y = ylab)
    
    # Fetches arguments of line
    line_args <- .get_line_args(colour_by = colour_by,
                                linetype_by = linetype_by,
                                size_by = size_by,
                                alpha = line_alpha,
                                linetype = line_type,
                                size = line_size)
    
    # Adds information, what column is used to group observations. Feature includes rownames.
    line_args$args$mapping$group <- sym("feature")
    
    # Adds arguments to the plot
    plot_out <- plot_out +
        do.call(geom_line, line_args$args)
    
    # Resolves the colours
    plot_out <- .resolve_plot_colours(plot_out,
                                      plot_data$colour_by,
                                      colour_by,
                                      fill = FALSE)
    
    # Changes theme
    plot_out <- plot_out +
        theme_classic()
    
    # Changes legend titles
    plot_out<- plot_out + guides(col=guide_legend("Colour"),
                              linetype=guide_legend("Linetype"),
                              size=guide_legend("Size"))
    
    # To choose if legend is kept, and its position
    plot_out <- .add_legend(plot_out, add_legend)
    
    return(plot_out)
    
}
