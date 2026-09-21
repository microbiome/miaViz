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

#' @title Multimedia example datasets
#'
#' @description
#' Example microbiome data objects used for demonstrations and examples across
#' the miaViz package, especially those that involve multi-view / multi-omics
#' visualization.
#'
#' Two related objects are provided:
#' \itemize{
#'   \item \code{se_relative}: a \code{SummarizedExperiment} with relative
#'         abundance data at the taxonomic level.
#'   \item \code{tse_relative}: a \code{TreeSummarizedExperiment} with the same
#'         relative abundance data as \code{se_relative}, together with a
#'         phylogenetic tree.
#' }
#'
#' Both objects share the same sample metadata columns (see \code{colData}) so
#' they can be used together to demonstrate paired visualizations at different
#' levels of structure.
#'
#' @details
#' \strong{se_relative} — \code{SummarizedExperiment}. Assays:
#' \code{relative_abundance}, \code{relabundance}, \code{scaled}.
#'
#' \strong{tse_relative} — \code{TreeSummarizedExperiment}. Same assays as
#' \code{se_relative} plus an attached row tree.
#'
#' Shared \code{colData} columns include: \code{study_name},
#' \code{subject_id}, \code{body_site}, \code{study_condition},
#' \code{disease}, \code{age}, \code{age_category}, \code{gender},
#' \code{country}, \code{BMI}, \code{visit_number}, \code{disease_subtype},
#' \code{disease_binary}, \code{dysbiosis}, \code{Time_point},
#' \code{treatment}, and others.
#'
#' @source
#' Derived from the \code{HMP_2019_ibdmdb.relative_abundance} dataset in the
#' \pkg{curatedMetagenomicData} R/Bioconductor package (Lloyd-Price et al.,
#' 2019). The full preparation pipeline — download, dysbiosis-score
#' computation, subsetting to paired baseline / post-treatment IBD samples,
#' treatment-factor definition and assay standardization — is shipped with
#' miaViz at \code{inst/scripts/multimedia_relative_data_prep.R}.
#'
#' Lloyd-Price J. et al. (2019). Multi-omics of the gut microbial ecosystem
#' in inflammatory bowel diseases. \emph{Nature} 569, 655-662.
#'
#' @name multimedia-datasets
#' @docType data
#' @keywords datasets
#' @examples
#' data(se_relative)
#' data(tse_relative)
#' @usage data(se_relative)
"se_relative"

#' @name multimedia-datasets
#' @usage data(tse_relative)
"tse_relative"
