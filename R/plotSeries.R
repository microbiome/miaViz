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
#' @param x_axis
#'  A single character value for selecting the column from 
#'  \code{\link[SummarizedExperiment:SummarizedExperiment-class]{ColData}} 
#'  that will specify values of x-axis. 
#'  
#' @param x_axis
#'  A single character value for selecting the taxa from 
#'  \code{\link[SummarizedExperiment:SummarizedExperiment-class]{rownames}}. 
#'  This parameter specifies taxa whose abundances will be plotted.
#'  
#'  @param rank a single character defining a taxonomic rank, that is used to
#'  agglomerate the data. Must be a value of \code{taxonomicRanks()} function.
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
setGeneric("plotSeries", signature = c("x"),
           function(x,
                    abund_values = NULL,
                    x_axis = NULL,
                    y_axis = NULL,
                    rank = NULL)
               standardGeneric("plotSeries"))


#' @rdname plotSeries
#' @export
setMethod("plotSeries", signature = c(x = "SummarizedExperiment"),
          function(x,
                   abund_values = NULL,
                   x_axis = NULL,
                   y_axis = NULL,
                   rank = NULL){
              
              ###################### Input check #######################
              # Check abund_values
              .check_abund_values(abund_values, x)
              
              # Check x_axis
              if( !(x_axis %in% names(colData(x))) ){
                  stop("'x_axis' must be a name of column of colData(x)", call. = FALSE)
              }
              
              # If rank is not null, data will be agglomerated by rank
              if( !is.null(rank) ){
                  # Check rank
                  .check_taxonomic_rank(rank, x)
                  
                  x <- agglomerateByRank(object, rank = rank)
                  
                  
                  ######################## should be removed?
                  # Rownames are now changed, so y-axis must be checked again
                  # Check y_axis
                  if (!is.null(y)){
                      if(!(y_axis %in% rownames(object)) ){
                          stop("'y' does not match with rownames(x). 
                          Data is agglomerated by 'rank', and rownames(object) are
                          also changed in relation to 'rank'. Check that 'y'
                          match with agglomerated data.", call. = FALSE)
                      }
                  }##########################################
              }
              
              # Check y_axis
              # If y_axis is not null, user has specified it
              if (!is.null(y_axis)){
                  if(!(y_axis %in% rownames(x)) ){
                      stop("'y_axis' must be in rownames(x).", call. = FALSE)
                  }
              }# If it is null, assign rownames
              else{
                  y <- rownames(object)
              }
              
              # Get assay data
              assay <- .get_assay_data(object, abund_values, y)
              
              data <- .incorporate_series_vis(assay, object, X, colour_by, linetype_by, size_by)
              
              
          }
)

################## HELP FUNCTIONS ##########################

.get_assay_data <- function(object, abund_values, y){
    
    # Gets warning or error if too many taxa are selected. 
    if(length(y) > 10 ){
        warning("Over 10 taxa selected.", call. = FALSE)
    }
    if(length(y) > 20 ){
        stop("Over 20 taxa selected. 20 or under allowed.", call. = FALSE)
    }
    
    # Take only those taxa that are specified by 'y'
    sub_object <- object[y]
    
    assay <- assay(sub_object, abund_values)
    
    return(assay)
}

.incorporate_series_vis <- function(assay, se, X, colour_by, linetype_by, size_by){
    
    # Stores the variables, se is the object
    
    variables <- c(X = X, 
                   colour_by = colour_by,
                   linetype_by = linetype_by,
                   size_by = size_by)
    
    # If variables there are variables
    if(!is.null(variables)){
        
        # Loops through variables
        for(i in seq_along(variables)){
            
            # If the variable is in colData
            if( variables[i] %in% names(colData(se)) ){
                # Retrieve series x-axis points from colData
                feature_info <- retrieveCellInfo(se, variables[1], search = "colData")
            }
            else{
                # get data from rowData
                feature_info <- retrieveFeatureInfo(se, variables[i],
                                                    search = "rowData")
            }
            # mirror back variable name, if a partial match was used
            var_name <- names(variables)[i]
            
        }
    }
    
    return(list(df = assay,
                X = X,
                colour_by = colour_by,
                linetype_by = linetype_by,
                size_by = size_by))
}



.melt_series_data <- function(){
    
    
}

