#' Adding information to tree data in \code{TreeSummarizedExperiment}
#'
#' To facilitate the dressing of the tree data stored in a
#' \code{TreeSummarizedExperiment} object, \code{rowTreeData} and
#' \code{colTreeData} can be used.
#'
#' @param x a
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object.
#'
#' @param other.fields,value a \code{data.frame} or coercible to one, with at
#' least one type of id information. See details.(Default: \code{list()})
#'
#' @param other_fields Deprecated. Use \code{other.fields} instead.
#'
#' @param tree.name \code{Character scalar}. Specifies a rowTree/colTree from
#'   \code{x}. (Default: \code{"phylo"})
#'
#' @param tree_name  Deprecated. Use \code{tree.name} instead.
#'
#' @param ... additional arguments, currently not used.
#'
#' @details
#' To match information to nodes, the id information in \code{other.fields} are
#' used. These can either be a column, named \sQuote{node} or \sQuote{label}
#' (\sQuote{node} taking precedent), or rownames. If all rownames can be coerced
#' to \code{integer}, they are considered as \sQuote{node} values, otherwise as
#' \sQuote{label} values. The id information must be unique and match available
#' values of \code{rowTreeData(c)}
#'
#' The result of the accessors, \code{rowTreeData} and \code{colTreeData},
#' contain at least a \sQuote{node} and \sQuote{label} column.
#'
#' @return a \code{data.frame} for the accessor and the modified
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object
#'
#' @name treeData
#'
#' @examples
#' data(GlobalPatterns)
#' td <- rowTreeData(GlobalPatterns)
#' td
#' td$test <- rnorm(nrow(td))
#' rowTreeData(GlobalPatterns) <- td
#' rowTreeData(GlobalPatterns)
#' combineTreeData(rowTree(GlobalPatterns), td)
NULL

#' @importFrom tidytree as_tibble
.get_tree_data <- function(tree){
    tree %>%
        as_tibble()
}


#' @rdname treeData
#' @importFrom dplyr last_col
#' @export
setMethod("colTreeData", signature = c(x = "TreeSummarizedExperiment"),
    function(x, tree.name = tree_name, tree_name = "phylo"){
        # Check tree.name
        if( !.is_a_string(tree.name) ){
            stop("'tree.name' must be a single character value specifying a ",
                "colTree.", call. = FALSE)
        }
        if(is.null(colTree(x, tree.name))){
            return(NULL)
        }
        .get_tree_data(colTree(x, tree.name)) %>%
            select(c("node","label":last_col()))
    }
)
#' @rdname treeData
#' @importFrom dplyr last_col
#' @importFrom dplyr all_of
#' @export
setMethod("rowTreeData", signature = c(x = "TreeSummarizedExperiment"),
    function(x, tree.name = tree_name, tree_name = "phylo"){
        # Check tree.name
        if( !.is_a_string(tree.name) ){
            stop("'tree.name' must be a single character value specifying a ",
                "rowTree.", call. = FALSE)
        }
        if(is.null(rowTree(x, tree.name))){
            return(NULL)
        }
        .get_tree_data(rowTree(x, tree.name)) %>%
            select(c("node","label":last_col()))
    }
)

DEFAULT_TREE_DATA_COLS <- c("parent","node","branch.length","label")
.clean_tree_data <- function(tree_data){
    tree_data %>%
        select(all_of(DEFAULT_TREE_DATA_COLS))
}

#' @rdname treeData
#' @export
setReplaceMethod("colTreeData", signature = c(x = "TreeSummarizedExperiment"),
    function(x, tree.name = tree_name, tree_name = "phylo", value){
        # Check tree.name
        if( !.is_a_string(tree.name) ){
            stop("'tree.name' must be a single character value specifying a ",
                "colTree.", call. = FALSE)
        }
        tree <- colTree(x, tree.name)
        # input check
        if(is.null(tree)){
            stop("'colTree(x, tree.name)' is NULL.", call. = FALSE)
        }
        # this is just temporary solution since phylo does not support data
        tree <- tidytree::as.phylo(combineTreeData(tree, value))
        x@colTree[[tree.name]] <- tree
        return(x)
    }
)

#' @rdname treeData
#' @importFrom tidytree as.phylo
#' @export
setReplaceMethod("rowTreeData", signature = c(x = "TreeSummarizedExperiment"),
    function(x, tree.name = tree_name, tree_name = "phylo", value){
        # Check tree.name
        if( !.is_a_string(tree.name) ){
            stop("'tree.name' must be a single character value specifying a ",
                "rowTree.", call. = FALSE)
        }
        tree <- rowTree(x)
        # input check
        if(is.null(tree)){
            stop("'rowTree(x)' is NULL.", call. = FALSE)
        }
        # this is just temporary solution since phylo does not support data
        tree <- tidytree::as.phylo(combineTreeData(tree, value))
        x@rowTree[[tree.name]] <- tree
        return(x)
    }
)

#' @importFrom tibble rownames_to_column
#' @importFrom dplyr select relocate
.norm_other_fields <- function(other_fields){
    if(is.null(other_fields) || length(other_fields) == 0L){
        return(NULL)
    }
    if(is(other_fields,"DataFrame") || !is.data.frame(other_fields)){
        other_fields <- as.data.frame(other_fields)
    }
    if(!is.data.frame(other_fields)){
        stop("'other_fields' must be a data.frame or coercible to one.",
            call. = FALSE)
    }
    if(nrow(other_fields) == 0L){
        return(NULL)
    }
    # check for node or label column and rownames
    rn <- rownames(other_fields)
    cn <- colnames(other_fields)
    f <- c("node","label") %in% cn
    if(!any(f)){
        # populate if necessary
        if(is.null(rn)){
            stop("Neither one of the following columns 'node'/'label' nor ",
                "rownames set for 'other_fields'.", call. = FALSE)
        }
        rn_i <- suppressWarnings(as.integer(rn))
        if(!anyNA(rn_i)){
            other_fields <- other_fields %>%
                rownames_to_column("node")
        } else {
            other_fields <- other_fields %>%
                rownames_to_column("label")
        }
    } else {
        other_fields <- other_fields %>%
            relocate(c("node","label")[f])
    }
    other_fields
}

.norm_id_col_of_other_fields <- function(other_fields, tree_data){
    if(is.null(other_fields)){
        return(other_fields)
    }
    if(missing(tree_data)){
        stop(".")
    }
    cn <- colnames(other_fields)
    # select one id column and reorder columns
    f <- which(cn %in% c("label","node"))
    if(length(f) == 2L){
        other_fields <- other_fields %>%
            select(!sym(cn[f[2L]]))
    }
    other_fields <- other_fields %>%
        relocate(sym(cn[f[1L]]))
    #
    by_col_name <- colnames(other_fields)[1L]
    if(anyDuplicated(other_fields[[by_col_name]])){
        stop("'",by_col_name,"' contains duplicate entries.", call. = FALSE)
    }
    if(!all(other_fields[[by_col_name]] %in% tree_data[[by_col_name]])){
        warning("Not all '",by_col_name,"' values found in tree data.",
                call. = FALSE)
    } else if(!any(other_fields[[by_col_name]] %in% tree_data[[by_col_name]])){
        stop("No overlap between '",by_col_name,"'and tree data.",
            call. = FALSE)
    }
    other_fields
}

.combine_tree_data_and_other_fields <- function(tree_data, other_fields){
    if(!is.null(other_fields)){
        by_col_name <- colnames(other_fields)[1L]
        tree_data <- tree_data %>%
            dplyr::left_join(other_fields, by = by_col_name)
    }
    tree_data
}

#' @importFrom tidytree as.treedata
.combine_tree_and_other_fields <- function(tree, other_fields = list()){
    tree_data <- .get_tree_data(tree)
    other_fields <- .norm_other_fields(other_fields)
    if(is.null(other_fields)){
        tree_data <- .clean_tree_data(tree_data)
    } else {
        other_fields <- .norm_id_col_of_other_fields(other_fields, tree_data)
        tree_data <- .combine_tree_data_and_other_fields(
            tree_data, other_fields)
    }
    tidytree::as.treedata(tree_data)
}

#' @rdname treeData
#' @export
setMethod("combineTreeData", signature = c(x = "phylo"),
    function(x, other.fields = other_fields, other_fields = list()){
        .combine_tree_and_other_fields(x, other.fields)
    }
)

#' @rdname treeData
#' @export
setMethod("combineTreeData", signature = c(x = "treedata"),
    function(x, other.fields = other_fields, other_fields = list()){
        .combine_tree_and_other_fields(x, other.fields)
    }
)
