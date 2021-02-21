#' Plot Series
#'
#' Function that plots core taxa against time. #####################################
#'
#' @param x
#' A \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#'
#' @param abund_values
#' A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}}
#'   to be plotted.
#'
#' @param x
#'  A single character value for selecting the column from 
#'  \code{\link[SummarizedExperiment:SummarizedExperiment-class]{ColData}} 
#'  that will specify values of x-axis. 
#'  
#' @param y
#'  A single character value for selecting the taxa from 
#'  \code{\link[SummarizedExperiment:SummarizedExperiment-class]{rownames}}. 
#'  This parameter specifies taxa whose abundances will be plotted.
#'  
#' @param rank a single character defining a taxonomic rank, that is used to
#'  agglomerate the data. Must be a value of \code{taxonomicRanks()} function.
#'  
#' @param colour_by
#'  
#' @param size_by
#'  
#' @param linetype_by
#'  
#'
#' @details
#' Creates a plot where core taxa is presented against time. ###############################
#'
#'
#' @references
#' Add here.#####################################
#'
#' @return
#' Returns plot###############################
#'
#' @name plotSeries
#' @export
#'
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#'
#' @examples
#' data(esophagus)
#' x <- esophagus
#'
NULL

#' @rdname plotSeries
#' @export
setGeneric("plotSeries", signature = c("object"),
           function(object,
                    abund_values,
                    x,
                    y = NULL,
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
                   x,
                   y = NULL,
                   rank = NULL,
                   colour_by = NULL,
                   size_by = NULL,
                   linetype_by = NULL){
              
              ###################### Input check #######################
              # Check abund_values
              .check_abund_values(abund_values, object)
              
              # Check x_axis
              if( !(x %in% names(colData(object))) ){
                  stop("'x' must be a name of column of colData(object)", call. = FALSE)
              }
              
              # If rank is not null, data will be agglomerated by rank
              if( !is.null(rank) ){
                  # Check rank
                  .check_taxonomic_rank(rank, object)
                  
                  object <- agglomerateByRank(object, rank = rank)
                  
                  
                  ######################## should be removed?
                  # Rownames are now changed, so y-axis must be checked again
                  # Check y_axis
                  if (!is.null(y)){
                      if(!(y %in% rownames(object)) ){
                          stop("'y' does not match with rownames(x). 
                          Data is agglomerated by 'rank', and rownames(object) are
                          also changed in relation to 'rank'. Check that 'y'
                          match with agglomerated data.", call. = FALSE)
                      }
                  }##########################################
              }
              
              # Check y_axis
              # If y_axis is not null, user has specified it
              if (!is.null(y)){
                  if(!(y %in% rownames(object)) ){
                      stop("'y' must be in rownames(x).", call. = FALSE)
                  }
              }# If it is null, assign rownames
              else{
                  y <- rownames(object)
              }
              
              # Get assay data
              assay <- .get_assay_data(object, abund_values, y)
              
              # Fetch series and features data as a list. 
              series_and_features_data <- .incorporate_series_vis(object, x, colour_by, linetype_by, size_by)
              # Splits it to individual datas
              series_data <- series_and_features_data$series_data$value
              feature_data <- series_and_features_data$feature_data
              
              # Melts the data
              melted_data <- .melt_series_data(assay, series_data, feature_data)
              
              plot_data <- data.frame(x = melted_data$x, y = melted_data$y)
              colour_by <- melted_data$colour_by
              size_by <- melted_data$size_by
              size_by <- melted_data$size_by
              xlab <- paste0(x)
              ylab <- paste0("Abundance ", abund_values)
              
              plot <- .series_plotter(plot_data, 
                                          xlab = xlab,
                                          ylab = ylab,
                                          colour_by = colour_by,
                                          size_by = size_by)
              
              return(plot)
              
          }
)

################## HELP FUNCTIONS ##########################

.get_assay_data <- function(object, abund_values, y){
    
    # Gets warning or error if too many taxa are selected. 
    if(length(y) > 20 ){
        stop("Over 20 taxa selected. 20 or under allowed.", call. = FALSE)
    } else if (length(y) > 10 ){
        warning("Over 10 taxa selected.", call. = FALSE)
    }
    
    # Take only those taxa that are specified by 'y'
    sub_object <- object[y]
    
    assay <- assay(sub_object, abund_values)
    
    return(assay)
}

.incorporate_series_vis <- function(object, x, colour_by, linetype_by, size_by){
    
    # Stores the variables
    variables <- c(x = x, 
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
                    names(feature_data)[names(feature_data) == "feature_info.value"] <- feature_info$name
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
    melted_data <- as.data.frame(assay) %>% pivot_longer(colnames(assay), names_to = "sample", values_to = "y")
    
    # Add series data to the data frame
    melted_data <- cbind(melted_data, x = series_data)
    
    # If feature_data is not null
    if( !is.null(feature_data) ){
        # Loops through feature_data
        for( i in ncol(feature_data) ){
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
                            size_by = NULL,
                            
                            flipped = FALSE,
                            add_legend = TRUE,
                            point_alpha = 1,
                            point_size = 2,
                            line_alpha = 1,
                            line_type = NULL,
                            line_size = 1,
                            ...){
    
    plot_out <- ggplot(plot_data) + geom_line(aes(x = x, y = y, color = colour_by)) +
        labs(x = xlab, y = ylab) 
    
    
}
