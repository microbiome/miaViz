#' Plot Core Taxa Time Series
#'
#' Function that plots cor taxa against time.
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
#' @param timepoints_column
#'  A single character value for selecting the column from 
#'  \code{\link[SummarizedExperiment:SummarizedExperiment-class]{ColData}} 
#'  that includes time points. 
#'  
#' @param core_taxa_size
#' An integer value specifying the number of species in core taxa. By default, it 
#' is 5.
#'
#' @details
#' Creates a plot where core taxa is presented against time. 
#'
#'
#' @references
#' Add here.#####################################
#'
#' @return
#' Returns plot###############################
#'
#' @name plotCoreTaxaTimeSeries
#' @export
#'
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#'
#' @examples
#' data(esophagus)
#' x <- esophagus
#'
NULL

#' @rdname plotCoreTaxaTimeSeries
#' @export
setGeneric("plotCoreTaxaTimeSeries", signature = c("x"),
           function(x,
                    abund_values = "counts",
                    timepoints_column = NULL,
                    core_taxa_size = 5)
               standardGeneric("plotCoreTaxaTimeSeries"))


#' @rdname plotCoreTaxaTimeSeries
#' @export
setMethod("plotCoreTaxaTimeSeries", signature = c(x = "SummarizedExperiment"),
          function(x,
                   abund_values = "counts",
                   timepoints_column = NULL,
                   core_taxa_size = 5){
              
              # Input check
              # Check abund_values
              .check_abund_values(abund_values, x)
              # Other checks
              ###############################################
              
              
              # Nothing ready yet #############################
              
              # Load moving pictures, it includes time points
              
              #devtools::load_all()
              
              mp <- moving_pictures
              # Create TSE from it
              TSE <- mia::makeTreeSummarizedExperimentFromphyloseq(mp)
              colData(TSE)
              ########################################################
              
              # Transforms abundance table to relative (there could be an option for relative or absolute values)
              TSE <- mia::transformCounts(TSE, method = "relabundance")
              
              # Get the names of most abundant taxa
              core_taxa <- mia::getTopTaxa(TSE, top = 5)
              core_taxa
              
              
              # Melt data to get a data frame PROBLEM: when creating molted df, it does not include colData, all values are NA. Problem in meltAssay?
              # I am not yet 100 % sure that TSE[core_taxa] works as intended. The purpose of this is to take only those taxa that are in the core taxa. 
              # It seems to work, but I have to check more carefully and with other data
              molten_se <- mia::meltAssay(TSE[core_taxa],
                                          add_row_data = TRUE,
                                          add_col_data = TRUE,
                                          abund_values = "relabundance")
              
              # Create a plot. days in x-axis and abundance in y-axis
              ggplot2::ggplot(na.omit(as.data.frame(molten_se))) + 
                  ggplot2::geom_line(ggplot2::aes(days_since_experiment_start, 
                                                  relabundance, color = FeatureID)) +
                  theme_bw() + 
                  theme(legend.position="top") + 
                  xlab("Days since experiment start") + 
                  ylab("Relative abundance") + 
                  scale_color_brewer("Core ASVs",palette = "Paired") +
                  guides(col = guide_legend(ncol = 3, nrow = 3))
              
              
              
              
          }
)

