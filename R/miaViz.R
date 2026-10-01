#' miaViz - Microbiome Analysis Plotting and Visualization
#'
#' The scope of this package is the plotting and visualization of microbiome
#' data. The main class for interfacing is the \code{TreeSummarizedExperiment}
#' class.
#'
#' @name miaViz-package
#' @seealso
#' \link[mia:mia-package]{mia} class
"_PACKAGE"
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
#' \code{getDissimilarity(genus, FUN = vegan::vegdist, method = "jaccard",
#' assay.type = "relabundance")}, either using transposed assay data or not
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

#' HMP_2019_ibdmdb: IBDMDB TreeSummarizedExperiment data
#'
#' An example gut microbiome dataset derived from the
#' \pkg{curatedMetagenomicData} \code{HMP_2019_ibdmdb} resource (iHMP/HMP2 IBD
#' cohort; Lloyd-Price et al., 2019).
#'
#' The object is provided as a \code{TreeSummarizedExperiment}.
#'
#' @format A \code{TreeSummarizedExperiment} with 585 microbial species (rows)
#' and 44 stool samples (columns) from 22 subjects with inflammatory bowel
#' disease (two time points per subject: baseline \code{T0} and follow-up
#' \code{T1}).
#'
#' Assays include relative abundance (\code{relabundance}) and a scaled version
#' of the same data (\code{scaled}). Taxonomic annotation is stored in
#' \code{rowData}, and a phylogenetic tree is available via \code{rowTree}.
#' Sample metadata (\code{colData}) includes subject identifiers, time point,
#' IBD subtype (CD/UC), and additional technical and demographic variables.
#'
#' @source
#' Derived from \code{curatedMetagenomicData} (Lloyd-Price et al., 2019).
#'
#' @references
#' Lloyd-Price J, Arze C, Ananthakrishnan AN, et al. (2019).
#' Multi-omics of the gut microbial ecosystem in inflammatory bowel diseases.
#' \emph{Nature} 569, 655--662. \doi{10.1038/s41586-019-1237-9}
#'
#' @name HMP_2019_ibdmdb
#' @docType data
#' @keywords datasets
#' @usage data(HMP_2019_ibdmdb)
#'
#' @examples
#' data("HMP_2019_ibdmdb", package = "miaViz")
#' HMP_2019_ibdmdb
#' SummarizedExperiment::assayNames(HMP_2019_ibdmdb)
#' TreeSummarizedExperiment::rowTree(HMP_2019_ibdmdb)
"HMP_2019_ibdmdb"
