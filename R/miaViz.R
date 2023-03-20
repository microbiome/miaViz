#' miaViz - Microbiome Analysis Plotting and Visualization
#'
#' The scope of this package is the plotting and visualization of microbiome 
#' data. The main class for interfacing is the \code{TreeSummarizedExperiment}
#' class.
#'
#' @name miaViz-package
#' @docType package
#' @seealso
#' \link[mia:mia-package]{mia} class
NULL

#' @import methods
#' @import TreeSummarizedExperiment
#' @import mia
#' @import ggplot2
#' @import ggraph
#' @importFrom rlang sym !! :=
#' @importFrom dplyr %>%
#' @importFrom BiocGenerics ncol nrow
NULL

#' @title miaViz example data
#'
#' @description
#' These example data objects were prepared to serve as examples. See the 
#' details for more information.
#' 
#' @details
#' For \code{*_graph} data:
#' 
#' 1. \dQuote{Jaccard} distances were calculated via 
#' \code{calculateDistance(genus, FUN = vegan::vegdist, method = "jaccard",
#' exprs_values = "relabundance")}, either using transposed assay data or not
#' to calculate distances for samples or features. NOTE: the function
#' mia::calculateDistance is now deprecated.
#' 
#' 2. \dQuote{Jaccard} dissimilarites were converted to similarities and values
#' above a threshold were used to construct a graph via 
#' \code{graph.adjacency(mode = "lower", weighted = TRUE)}.
#' 
#' 3. The \code{igraph} object was converted to \code{tbl_graph} via 
#' \code{as_tbl_graph} from the \code{tidygraph} package.
#' 
#' 
#' @name mia-datasets
#' @docType data
#' @keywords datasets
#' @usage data(col_graph)
"col_graph"
#' @name mia-datasets
#' @usage data(row_graph)
"row_graph"
#' @name mia-datasets
#' @usage data(row_graph_order)
"row_graph_order"
