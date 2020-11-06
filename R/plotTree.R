#' Plotting tree information enriched with information
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


#' @rdname plotTree
#' @export
setMethod("plotColTree", signature = c(x = "TreeSummarizedExperiment"),
    function(x, ..., add_traits = NULL, trait_data = NULL,
             by_exprs_values = "counts"){
        # input check
        if(is.null(colTree(x))){
          stop("colTree() is empty.", call. = FALSE)
        }

        #
        tree <- .get_trimed_tree(x, type = "column", relabel = FALSE,
                                 colnames(x))
        tree_data <- .get_tree_data(tree, x)
        #
        trait_data <- .norm_trait_data(trait_data, tree_data)
        trait_data <- .combine_trait_data(x, add_traits, trait_data,
                                          by_exprs_values = by_exprs_values,
                                          type = "column")
        tree_data <- .add_trait_data_to_tree(tree_data, x, trait_data)
        #
        .plot_tree(tree_data, add_traits)
    }
)

#' @rdname plotTree
#' @export
setMethod("plotRowTree", signature = c(x = "TreeSummarizedExperiment"),
    function(x, ..., add_traits = NULL, trait_data = NULL,
             by_exprs_values = "counts", relabel_tree = FALSE){
        browser()
        # input check
        if(is.null(rowTree(x))){
          stop("rowTree() is empty.", call. = FALSE)
        }
        if(!is.logical(relabel_tree) || length(relabel_tree) != 1L){
            stop("'relabel_tree' must be either TRUE or FALSE.", call. = FALSE)
        }
        #
        tree <- .get_trimed_tree(x, type = "row", relabel = relabel_tree,
                                 rownames(x))
        tree_data <- .get_tree_data(tree, x)
        #
        trait_data <- .norm_trait_data(trait_data, rownames(x))
        trait_data <- .combine_trait_data(x, add_traits, trait_data,
                                          by_exprs_values = by_exprs_values,
                                          type = "row")
        tree_data <- .add_trait_data_to_tree(tree_data, x, trait_data)
        #
        .plot_tree(tree_data, add_traits)
    }
)

#' @importFrom tibble rownames_to_column
.norm_trait_data <- function(trait_data, dimnames,
                             dimnames_name = .get_name_in_parent(dimnames)){
    if(is.null(trait_data)){
        return(NULL)
    }
    if(is(trait_data,"DataFrame") || !is.data.frame(trait_data)){
        trait_data <- as.data.frame(trait_data)
    }
    if(!is.data.frame(trait_data)){
        stop("'trait_data' must be a data.frame or coercible to one.",
             call. = FALSE)
    }
    if(!is.null(trait_data$label)){
        if(!all(trait_data$label %in% dimnames)){
            warning("Not all values of trait_data$label match entries in ",
                    dimnames_name,".",
                    call. = FALSE)
        } else if(!any(trait_data$label %in% dimnames)){
            stop("No overlapp of trait_data$label and ",
                 dimnames_name,".",
                 call. = FALSE)
        }
    } else {
        if(!all(rownames(trait_data) %in% dimnames)){
            warning("Not all rownames(trait_data) match entries in ",
                    dimnames_name,".",
                    call. = FALSE)
        } else if(is.null(rownames(trait_data)) ||
           !any(rownames(trait_data) %in% dimnames)){
            stop("rownames of 'trait_data' must be set and contain values of ",
                 dimnames_name,".",
                 call. = FALSE)
        }
        trait_data <- trait_data %>%
            rownames_to_column("label")
    }
    trait_data <- trait_data[trait_data$label %in% dimnames,]
    trait_data
}

.get_info <- function(by, x, FUN, exprs_values){
    feature_info <- try(type_FUN(x, by = add_traits,
                                 exprs_values = by_exprs_values),
                        silent = TRUE)
    if(is(feature_info,"try-error")){
        return(NULL)
    }
    tibble(!!sym(feature_info$name) := feature_info$value)
}

#' @importFrom scater retrieveFeatureInfo retrieveCellInfo
.combine_trait_data <- function(x, add_traits, trait_data,
                                by_exprs_values = "counts",
                                type = c("row","column")){
    if(!is.null(trait_data)){
        add_traits <- add_traits[!(add_traits %in% colnames(trait_data))]
    }
    if(length(add_traits) > 0L){
        type <- match.arg(type)
        type_FUN <- switch(type,
                           row = scater::retrieveFeatureInfo,
                           column = scater::retrieveCellInfo,
                           stop("."))
        feature_info <- lapply(add_traits, .get_info, x = x, FUN = type_FUN,
                               exprs_values = by_exprs_values)
        f <- !vapply(feature_info,is.null,logical(1))
        if(all(!f)){
            warning("No data for values of 'add_traits' found in 'x'.",
                    call. = FALSE)
        } else if(any(!f)) {
            warning("No data for values '",
                    paste(add_traits, collapse = "', '"),
                    "' found in 'x'.",
                    call. = FALSE)
        }
        feature_info <- feature_info[f]
        if(length(feature_info) != 0L){
            feature_info <- do.call(cbind,feature_info) %>%
                rownames_to_column("label")
            if(!is.null(trait_data)){
                trait_data <- feature_info %>%
                    left_join(trait_data, by = "label")
            } else {
                trait_data <- feature_info
            }
        }
    }
    trait_data
}


#' @importFrom ape keep.tip
.get_trimed_tree <- function(x, type = c("row","columns"),
                             relabel = FALSE, dimnames){
    type <- match.arg(type)
    tree_FUN <- switch(type, row = rowTree, column = colTree, stop("."))
    links_FUN <- switch(type, row = rowLinks, column = colLinks, stop("."))
    tree <- tree_FUN(x)
    links <- links_FUN(x)
    tree <- ape::keep.tip(tree, unique(links$nodeNum))
    m <- match(unique(links$nodeNum),links$nodeNum)
    if(relabel){
        if(type == "column"){
            stop(".") # Now taxonomic info on the cols
        }
        new_tip_labels <- getTaxonomyLabels(x[m,], with_type = TRUE)
    } else {
        new_tip_labels <- dimnames[m]
    }
    tree$tip.label <- new_tip_labels
    tree
}

#' @importFrom tibble as_tibble
.get_tree_data <- function(tree, x){
    tree_data <- as_tibble(tree)
    tree_data
}

#' @importFrom dplyr left_join
.add_trait_data_to_tree <- function(tree_data, x, trait_data){
    if(!is.null(trait_data)){
        tree_data <- tree_data %>%
            left_join(trait_data, by = "label")
    }
    tree_data
}

#' @importFrom ape as.phylo
#' @importFrom ggtree ggtree
.plot_tree <- function(tree_data, layout = "circular"){
    p <- ggtree(as.phylo(tree_data), layout = layout)
    p
}
