#' Order rows or columns for heatmap visualization
#'
#' @description
#' \code{getHeatmapOrder} determines an order for the rows or columns of a
#' matrix or an assay so that similar entities are placed next to each other.
#' The order can be used to sort heatmaps, and it is used internally by
#' plotting functions such as \code{\link{plotHeatmap}} and
#' \code{\link{plotLoadings}}.
#'
#' @details
#' The ordering method is specified with \code{method}:
#'
#' \itemize{
#'   \item \strong{Single variable}: If \code{method} is a name of a column
#'   (when \code{MARGIN = "rows"}) or a row (when \code{MARGIN = "cols"}),
#'   entities are ordered based on the values of that variable. For instance,
#'   features of a loadings matrix can be ordered based on a single principal
#'   component.
#'
#'   \item \strong{Seriation}: Other \code{character} values specify a
#'   seriation method that is passed to
#'   \code{\link[seriation:seriate]{seriation::seriate()}}. See
#'   \code{seriation::list_seriation_methods("matrix")} for available
#'   methods. For example, \code{"PCA_angle"} orders entities by their angle
#'   in the space of the first two principal components (NeatMap ordering,
#'   see Rajaram & Oono 2010), and \code{"Heatmap"} uses hierarchical
#'   clustering with optimal leaf ordering.
#'
#'   \item \strong{Clustering}: If \code{method} is a \code{BlusterParam}
#'   object, entities are clustered with
#'   \code{\link[bluster:clusterRows]{bluster::clusterRows()}}, and entities
#'   belonging to the same cluster are placed next to each other.
#' }
#'
#' If \code{x} is a \code{dist} object, seriation methods for distances (see
#' \code{seriation::list_seriation_methods("dist")}), such as
#' \code{"MDS_angle"} (NeatMap-like angular ordering in the space of classical
#' multidimensional scaling) or \code{"OLO"} (hierarchical clustering with
#' optimal leaf ordering), are applied to the distances directly. Other methods
#' are applied to the full distance matrix, i.e., each entity is described by
#' its distances to all entities.
#'
#' @param x A \code{matrix}, a \code{dist} object or a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#'
#' @param assay.type \code{Character scalar}. Specifies the assay to use when
#' \code{x} is a \code{SummarizedExperiment}. (Default: \code{"counts"})
#'
#' @param method \code{Character scalar} or \code{BlusterParam} object.
#' Specifies the ordering method. See details. (Default: \code{"PCA_angle"},
#' or \code{"MDS_angle"} when \code{x} is a \code{dist} object)
#'
#' @param MARGIN \code{Character scalar}. Specifies whether to order
#' \code{"rows"} or \code{"cols"}. Ignored when \code{x} is a \code{dist}
#' object. (Default: \code{"rows"})
#'
#' @param decreasing \code{Logical scalar}. When ordering based on a single
#' variable, should the order be decreasing? (Default: \code{TRUE})
#'
#' @param control \code{NULL} or \code{list}. Control parameters passed to
#' \code{\link[seriation:seriate]{seriation::seriate()}}.
#' (Default: \code{NULL})
#'
#' @param ... Additional arguments. Currently not used.
#'
#' @return An \code{integer} vector of row or column indices in the sorted
#' order.
#'
#' @references
#' Rajaram S & Oono Y (2010) NeatMap - non-clustering heat map alternatives
#' in R. \emph{BMC Bioinformatics} 11, 45. \doi{10.1186/1471-2105-11-45}
#'
#' Hahsler M, Hornik K & Buchta C (2008) Getting Things in Order: An
#' Introduction to the R Package seriation. \emph{Journal of Statistical
#' Software} 25(3), 1-34. \doi{10.18637/jss.v025.i03}
#'
#' @seealso
#' \code{\link{plotHeatmap}}, \code{\link{plotLoadings}}
#'
#' @name getHeatmapOrder
#'
#' @examples
#' library(mia)
#' data(GlobalPatterns)
#' tse <- agglomerateByPrevalence(GlobalPatterns, rank = "Phylum")
#' tse <- transformAssay(tse, method = "clr", pseudocount = TRUE)
#'
#' # Order features with NeatMap-like angular ordering
#' row_order <- getHeatmapOrder(tse, assay.type = "clr")
#' # Order samples with hierarchical clustering
#' col_order <- getHeatmapOrder(
#'     tse, assay.type = "clr", method = "Heatmap", MARGIN = "cols")
#' tse <- tse[row_order, col_order]
#'
#' # Order based on clustering
#' row_order <- getHeatmapOrder(
#'     tse, assay.type = "clr", method = bluster::KmeansParam(3))
#'
#' # Order features of a matrix based on a single column
#' mat <- assay(tse, "clr")
#' getHeatmapOrder(mat, method = "CL3")
#'
#' # Order samples based on dissimilarities
#' diss <- getDissimilarity(tse, method = "euclidean", assay.type = "clr")
#' getHeatmapOrder(diss)
#'
NULL

#' @rdname getHeatmapOrder
#' @export
setMethod("getHeatmapOrder", signature = c(x = "SummarizedExperiment"),
    function(x, assay.type = "counts", ...){
        .check_assay_present(assay.type, x)
        mat <- assay(x, assay.type)
        getHeatmapOrder(mat, ...)
    }
)

#' @rdname getHeatmapOrder
#' @export
#' @importClassesFrom BiocGenerics dist
setMethod("getHeatmapOrder", signature = c(x = "dist"),
    function(x, method = "MDS_angle", control = NULL, ...){
        if( !(.is_a_string(method) || is(method, "BlusterParam")) ){
            stop("'method' must be a single character value or a ",
                "BlusterParam object.", call. = FALSE)
        }
        if( !(is.null(control) || is.list(control)) ){
            stop("'control' must be NULL or a list.", call. = FALSE)
        }
        # Seriation methods for distances are applied to the dist object. Other
        # methods are applied to the distance matrix.
        if( .is_a_string(method) && attr(x, "Size") > 1L ){
            .require_package("seriation")
            if( method %in% seriation::list_seriation_methods("dist") ){
                res <- seriation::seriate(
                    x, method = method, control = control) |>
                    seriation::get_order()
                return(unname(as.integer(res)))
            }
        }
        .get_heatmap_order(as.matrix(x), method = method, control = control)
    }
)

#' @rdname getHeatmapOrder
#' @export
setMethod("getHeatmapOrder", signature = c(x = "ANY"),
    function(
        x, method = "PCA_angle", MARGIN = "rows", decreasing = TRUE,
        control = NULL, ...){
        .get_heatmap_order(
            x, method = method, MARGIN = MARGIN, decreasing = decreasing,
            control = control)
    }
)

############################### HELPER FUNCTIONS ###############################

# This function determines order of rows or columns based on a single variable,
# seriation or clustering.
.get_heatmap_order <- function(
        x, method = "PCA_angle", MARGIN = "rows", decreasing = TRUE,
        control = NULL){
    # Check arguments
    if( !(.is_a_string(method) || is(method, "BlusterParam")) ){
        stop("'method' must be a single character value or a BlusterParam ",
            "object.", call. = FALSE)
    }
    MARGIN <- .check_MARGIN(MARGIN)
    if( !.is_a_bool(decreasing) ){
        stop("'decreasing' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !(is.null(control) || is.list(control)) ){
        stop("'control' must be NULL or a list.", call. = FALSE)
    }
    x <- as.matrix(x)
    if( !is.numeric(x) ){
        stop("Values must be numeric.", call. = FALSE)
    }
    # Order is always calculated for rows
    if( MARGIN == 2L ){
        x <- t(x)
    }
    # With fewer than two entities, there is nothing to order
    if( nrow(x) < 2L ){
        return(seq_len(nrow(x)))
    }
    if( .is_a_string(method) && method %in% colnames(x) ){
        # Order based on single variable
        res <- order(x[, method], decreasing = decreasing)
    } else if( .is_a_string(method) ){
        # Order based on seriation
        .require_package("seriation")
        res <- seriation::seriate(
            x, method = method, control = control, margin = 1L) |>
            seriation::get_order()
    } else{
        # Order based on clustering. Entities in same cluster are placed next
        # to each other.
        .require_package("bluster")
        res <- bluster::clusterRows(x, BLUSPARAM = method) |>
            order()
    }
    res <- unname(as.integer(res))
    return(res)
}
