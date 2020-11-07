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
setGeneric("plotRowTree", signature = c("object"),
           function(object, ...)
               standardGeneric("plotRowTree"))
#' @rdname plotTree
setGeneric("plotColTree", signature = c("object"),
           function(object, ...)
               standardGeneric("plotColTree"))


#' @rdname plotTree
#' @export
setMethod("plotColTree", signature = c(object = "TreeSummarizedExperiment"),
    function(object,
             y = NULL,
             layout = "circular",
             show_label = FALSE,
             colour_by = NULL,
             shape_by = NULL,
             size_by = NULL,
             by_exprs_values = "counts",
             other_fields = list(),
             ...){
        # input check
        if(is.null(colTree(object))){
          stop("colTree(object) is empty.", call. = FALSE)
        }
        .check_plotting_options(layout = layout,
                                show_label = show_label,
                                colour_by = colour_by,
                                shape_by = shape_by,
                                size_by = size_by)
        #
        tree <- .get_trimed_tree(object, type = "column", relabel = FALSE,
                                 colnames(object))
        tree_data <- .get_tree_data(tree)
        #
        variables <- .norm_variables_for_tree_plotting(c(y,colour_by,shape_by,
                                                         size_by))
        fields <- .norm_fields(other_fields, tree_data)
        fields <- .add_object_data_to_fields(object, variables, fields,
                                             by_exprs_values = by_exprs_values,
                                             type = "column")
        fields <- .norm_id_col_of_fields(fields, tree_data)
        tree_data <- .combine_tree_data_and_fields(tree_data, fields)
        #
        .plot_tree(tree_data, y,
                   layout = layout,
                   show_label = show_label,
                   colour_by = colour_by,
                   shape_by = shape_by,
                   size_by = size_by)
    }
)

#' @rdname plotTree
#' @export
setMethod("plotRowTree", signature = c(object = "TreeSummarizedExperiment"),
    function(object,
             y = NULL,
             relabel_tree = FALSE,
             layout = "circular",
             show_label = FALSE,
             colour_by = NULL,
             shape_by = NULL,
             size_by = NULL,
             by_exprs_values = "counts",
             other_fields = list(),
             ...){
        browser()
        # input check
        if(is.null(rowTree(object))){
          stop("rowTree(object) is empty.", call. = FALSE)
        }
        if(!is.logical(relabel_tree) || length(relabel_tree) != 1L){
            stop("'relabel_tree' must be either TRUE or FALSE.", call. = FALSE)
        }
        .check_plotting_options(layout = layout,
                                show_label = show_label,
                                colour_by = colour_by,
                                shape_by = shape_by,
                                size_by = size_by)
        #
        tree <- .get_trimed_tree(object, type = "row", relabel = relabel_tree,
                                 rownames(object))
        tree_data <- .get_tree_data(tree)
        #
        variables <- .norm_variables_for_tree_plotting(c(y,colour_by,shape_by,
                                                         size_by))
        fields <- .norm_fields(other_fields, rownames(object))
        fields <- .add_object_data_to_fields(object, variables, fields,
                                             by_exprs_values = by_exprs_values,
                                             type = "row")
        fields <- .norm_id_col_of_fields(fields, tree_data)
        tree_data <- .combine_tree_data_and_fields(tree_data, fields)
        #
        .plot_tree(tree_data, y,
                   layout = layout,
                   show_label = show_label,
                   colour_by = colour_by,
                   shape_by = shape_by,
                   size_by = size_by)
    }
)

.check_plotting_options <- function(...){

}

.norm_variables_for_tree_plotting <- function(x){
    if(is.null(x)){
        return(NULL)
    }
    names_x <- x
    f <- x %in% c("node","label")
    if(any(f)){
        x <- paste0(x[f], "_tmp")
    }
    names(x) <- names_x
    x
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

#' @importFrom tibble tibble
.get_feature_info <- function(by, object, FUN, exprs_values){
    feature_info <- try(FUN(object, by = by, exprs_values = exprs_values),
                        silent = TRUE)
    if(is(feature_info,"try-error")){
        return(NULL)
    }
    tibble(!!sym(feature_info$name) := feature_info$value)
}

#' @importFrom scater retrieveFeatureInfo retrieveCellInfo
#' @importFrom dplyr bind_cols mutate relocate
#' @importFrom tibble rownames_to_column
.add_object_data_to_fields <- function(object, variables, fields,
                                       by_exprs_values = "counts",
                                       type = c("row","column")){
    if(is.null(variables)){
        return(fields)
    }
    type <- match.arg(type)
    # remove any variables values, which are already available
    cn <- colnames(fields)
    cn_data <- cn[!(cn %in% c("node","label"))]
    if(!is.null(fields) && length(cn_data) > 0L){
        variables <- variables[!(names(variables) %in% cn_data)]
    }
    # for remaining variables try to get data
    if(length(variables) > 0L){
        type_FUN <- switch(type,
                           row = scater::retrieveFeatureInfo,
                           column = scater::retrieveCellInfo,
                           stop("."))
        feature_info <- vector(mode = "list", length = length(variables))
        for(i in seq_along(variables)){
            feature_info[[i]] <-
                .get_feature_info(names(variables)[i], object = object,
                                  FUN = type_FUN, exprs_values = by_exprs_values)
        }
        f <- !vapply(feature_info,is.null,logical(1))
        if(any(!f)) {
            stop("No data for values '", paste(variables[!f],collapse = "', '"),
                 "' found in 'object'.",
                 call. = FALSE)
        }
        feature_info <- feature_info[f]
        if(length(feature_info) != 0L){
            feature_info <- bind_cols(feature_info)
            colnames(feature_info) <- names(variables)[f]
            rn <- rownames(feature_info)
            rn_i <- suppressWarnings(as.integer(rn))
            if(!anyNA(rn_i)){
                feature_info <- feature_info %>%
                    rownames_to_column("node") %>%
                    mutate(node = as.integer(node)) %>%
                    relocate("node")
            } else {
                feature_info <- feature_info %>%
                    rownames_to_column("label") %>%
                    relocate("label")
            }
            if(!is.null(fields)){
                fields <- feature_info %>%
                    left_join(fields, by = colnames(feature_info)[1L])
            } else {
                fields <- feature_info
            }
        }
    }
    fields
}

#' @importFrom tidytree as.treedata
#' @importFrom ggtree ggtree geom_tiplab geom_tippoint
.plot_tree <- function(tree_data,
                       y,
                       layout = "circular",
                       show_label = TRUE,
                       colour_by = NULL,
                       shape_by = NULL,
                       size_by = NULL){
    # assemble arg list
    args <- list(y = y,
                 layout = layout,
                 show_label = show_label,
                 colour_by = colour_by,
                 shape_by = shape_by,
                 size_by = size_by)
    # start plotting
    p <- ggtree(tidytree::as.treedata(tree_data), layout = layout)
    if(show_label){
        p <- p + geom_tiplab(offset = .1)
    }
    point_aes <- .get_point_aes(args)
    colour <- .get_colour(tree_data, args)
    p <- p +
        geom_tippoint(point_aes, shape = 21) +
        geom_nodepoint(point_aes, shape = 21) +
        colour
    p
}

.get_point_aes <- function(args){
    colour <- args$colour_by
    shape <- args$shape_by
    size <- args$size_by
    aes <- aes_string(fill = colour,
                      shape = shape,
                      size = size)
    aes
}

.get_colour <- function(tree_data, args){
    if(is.null(args$colour_by)){
        return(NULL)
    }
    if(is.null(args$shape_by)){
        if(is.numeric(tree_data[[args$colour_by]])){
            scale_fill_viridis_c()
        } else {
            scale_fill_viridis_d()
        }
    } else {
        if(is.numeric(tree_data[[args$colour_by]])){
            scale_colour_viridis_c()
        } else {
            scale_colour_viridis_d()
        }
    }
}
