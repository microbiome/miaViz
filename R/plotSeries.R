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
              
              # Check y_axis
              # If y_axis is not null, user has specified it
              if (!is.null(y_axis)){
                  if(!(y_axis %in% rownames(x)) ){
                      stop("'y_axis' must be in rownames(x).", call. = FALSE)
                  }
              }
              
              
              # Check rank
              # IF rank is not null, it has been specified by user
              if( !is.null(rank) ){
                  if(!.is_non_empty_string(rank)){
                      stop("'rank' must be an non empty single character value,
                           and it must be one of taxonomyRanks(x).",
                           call. = FALSE)
                  }
                  .check_taxonomic_rank(rank, x)
              }
              ##############################################
              
              # If rank is not null, data will be agglomerated by rank
              if( !is.null(rank) ){
                  x <- agglomerateByRank(x, rank = rank)
                  
                  # Rownames are now changed, so y-axis must be checked again
                  # Check y_axis
                  if (!is.null(y_axis)){
                      if(!(y_axis %in% rownames(x)) ){
                          stop("'y-axis' does not match to rownames(x). 
                          Data is agglomerated by 'rank', and rownames(x) are
                          also changed in relation to 'rank'. Check that 'y-axis'
                          match to agglomerated data.", call. = FALSE)
                      }
                  }
              }
              
              # If y_axis is not null, take only those taxa into account that user
              # is specified
              if( !is.null(y_axis) ){
                  x <- x[y_axis]
              }
              
              # Melt data
              melted_data <- meltAssay(x, add_col_data = TRUE, add_row_data = TRUE)
              color_by = "Kingdom"
              line_by = "FeatureID"
              
              # # Take only those columns that are needed
              # melted_data <- melted_data[, c("FeatureID", abund_values, x_axis)]
              
              ggplot2::ggplot(melted_data) + 
                  ggplot2::geom_line(
                      ggplot2::aes_string(x = x_axis, y = abund_values, 
                                          color = color_by, linetype = line_by))
              
              
                  #     ) +
                  # theme_bw() + 
                  # theme(legend.position="top") + 
                  # xlab(x_axis) + 
                  # ylab(abund_values) + 
                  # scale_color_brewer("Core ASVs",palette = "Paired") +
                  # guides(col = guide_legend(ncol = 3, nrow = 3))
              
              
              
              
          }
)

