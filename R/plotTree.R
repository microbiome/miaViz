#' plotting tree information enriched with information
#'
#' Based on the stored data in a \code{TreeSummarizedExperiment} a tree can
#' be plotted. From the \code{rowData}, the \code{assays} as well as the
#' \code{colData} information can be taken for enriching the tree plots with
#' additional information.
#'
#' @param x a
#'   \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'   object.
#'
#' @return a \code{\link{ggtree}} plot
#'
#' @name plotTree
#'
#' @examples
#'
NULL

#' @rdname plotTree
setGeneric("plotRowTree", signature = c("x"),
           function(x, ...)
               standardGeneric("plotRowTree"))
#' @rdname plotTree
setGeneric("plotColTree", signature = c("x"),
           function(x, ...)
               standardGeneric("plotColTree"))

#' @importFrom ape keep.tip
.trim_tree <- function(tree, x){
    rl <- rowLinks(x)
    tree <- ape::keep.tip(tree, rl$nodeNum)
}

#' @importFrom tibble as_tibble
.get_tree_data <- function(tree, x){
    tree_data <- as_tibble(tree)
}


.add_data_to_tree <- function(tree, x){

}

#' @importFrom ape as.phylo
#' @importFrom ggtree ggtree
.plot_tree <- function(tree_data, layout = "circular"){
    p <- ggtree(as.phylo(tree_data), layout = layout)
    p
}

#' @rdname plotTree
#'
#'
#'
#' @export
setMethod("plotRowTree", signature = c(x = "TreeSummarizedExperiment"),
    function(x, ..., additional_data = NULL){
        # input check
        if(is.null(rowTree(x))){
            stop("rowTree() is empty.", call. = FALSE)
        }
        #
        tree <- rowTree(x)
        tree <- .trim_tree(tree, x)
        tree_data <- .get_tree_data(tree, x)
        tree_data <- .add_data_to_tree(tree, x, additional_data)
        .plot_tree(tree)
    }
)
