#' @title Multimedia example datasets
#'
#' @description
#' Example microbiome data objects bundled together as `multimedia_data_demo`
#' for use in demonstrations and examples across the miaViz package,
#' especially those that involve multi-view / multi-omics visualization.
#'
#' Loading `data(multimedia_data_demo)` makes three objects available:
#' \itemize{
#'   \item \code{se_relative}: a \code{SummarizedExperiment} with relative
#'         abundance data at the taxonomic level.
#'   \item \code{tse_relative}: a \code{TreeSummarizedExperiment} with the same
#'         relative abundance data as \code{se_relative}, together with a
#'         phylogenetic tree.
#'   \item \code{tse_pathway}: a \code{TreeSummarizedExperiment} containing
#'         pathway-level functional abundance data.
#' }
#'
#' All three objects share the same sample metadata columns (see
#' \code{colData}) so they can be used together to demonstrate paired
#' taxonomy / pathway visualizations.
#'
#' @details
#' \strong{se_relative} — \code{SummarizedExperiment} with 579 features x
#' 44 samples. Assays: \code{relative_abundance}, \code{relabundance},
#' \code{scaled}. Row names are taxa; column names are sample IDs.
#'
#' \strong{tse_relative} — \code{TreeSummarizedExperiment} with 579 features x
#' 44 samples. Same assays as \code{se_relative} plus an attached
#' \code{rowTree}.
#'
#' \strong{tse_pathway} — \code{TreeSummarizedExperiment} with 476 features x
#' 36 samples. Assays: \code{pathway_abundance}, \code{scaled}. Row names are
#' MetaCyc pathway identifiers.
#'
#' Shared \code{colData} columns include: \code{study_name},
#' \code{subject_id}, \code{body_site}, \code{study_condition},
#' \code{disease}, \code{age}, \code{age_category}, \code{gender},
#' \code{country}, \code{BMI}, \code{visit_number}, \code{disease_subtype},
#' \code{disease_binary}, \code{dysbiosis}, \code{Time_point},
#' \code{treatment}, and others.
#'
#' @name multimedia_data_demo
#' @aliases se_relative tse_relative tse_pathway
#' @docType data
#' @keywords datasets
#' @usage data(multimedia_data_demo)
#' @format
#' A single \code{.rda} file that loads three objects: \code{se_relative}
#' (\code{SummarizedExperiment}, 579 x 44), \code{tse_relative}
#' (\code{TreeSummarizedExperiment}, 579 x 44), and \code{tse_pathway}
#' (\code{TreeSummarizedExperiment}, 476 x 36).
#' @examples
#' data(multimedia_data_demo)
#' se_relative
#' tse_relative
#' tse_pathway
"se_relative"

#' @rdname multimedia_data_demo
"tse_relative"

#' @rdname multimedia_data_demo
"tse_pathway"
