#' Plotting tree information enriched with information
#'
#' Based on the stored data in a \code{TreeSummarizedExperiment} a tree can
#' be plotted. From the \code{rowData}, the \code{assays} as well as the
#' \code{colData} information can be taken for enriching the tree plots with
#' additional information.
#'
#' @param object a
#'   \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'   object.
#'
#' @param relabel_tree logical scalar, Should the tip labels be relabeled using 
#'   the output of \code{getTaxonomyLabels(object, with_type = TRUE)}? This
#'   has consequences on how data can be merged from \code{other_fields}. 
#'   (default: \code{relabel_tree = FALSE})
#' @param show_label logical scalar or named logical vector. Should tip labels
#'   be plotted or if a logical vector is provided, which labels should be 
#'   shown? Only names corresponding to \code{TRUE}, will be plotted. 
#'   (default: \code{show_label = FALSE})
#' @param add_legend logical scalar. Should legens be plotted? 
#'   (default: \code{add_legend = TRUE})
#' @param layout layout for the plotted tree. See 
#'   \code{\link[ggtree:ggtree]{ggtree}} for details.
#' @param edge_colour_by Specification of a column metadata field or a feature 
#'   to colour tree edges by, see the by argument in 
#'   \code{\link[scater:retrieveCellInfo]{?retrieveCellInfo}} for possible 
#'   values.
#' @param edge_size_by Specification of a column metadata field or a feature 
#'   to size tree edges by, see the by argument in 
#'   \code{\link[scater:retrieveCellInfo]{?retrieveCellInfo}} for possible 
#'   values.
#' @param tip_colour_by Specification of a column metadata field or a feature to
#'   colour tree tips by, see the by argument in 
#'   \code{\link[scater:retrieveCellInfo]{?retrieveCellInfo}} for possible 
#'   values.
#' @param tip_shape_by Specification of a column metadata field or a feature to
#'   shape tree tips by, see the by argument in 
#'   \code{\link[scater:retrieveCellInfo]{?retrieveCellInfo}} for possible 
#'   values.
#' @param tip_size_by Specification of a column metadata field or a feature to
#'   size tree tips by, see the by argument in 
#'   \code{\link[scater:retrieveCellInfo]{?retrieveCellInfo}} for possible 
#'   values.
#' @param node_colour_by Specification of a column metadata field or a feature to
#'   colour tree nodes by. Must be a field from \code{other_fields}.
#' @param node_shape_by Specification of a column metadata field or a feature to
#'   shape tree nodes by. Must be a field from \code{other_fields}.
#' @param node_size_by Specification of a column metadata field or a feature to
#'   size tree nodes by. Must be a field from \code{other_fields}.
#' @param by_exprs_values A string or integer scalar specifying which assay to
#'   obtain expression values from, for use in point aesthetics - see the 
#'   \code{exprs_values} argument in 
#'   \code{\link[scater:retrieveCellInfo]{?retrieveCellInfo}}.
#' @param other_fields a \code{data.frame} or coercible to one, with at least 
#'   one type of id information. See details in
#'   \code{\link[=treeData]{treeData}}.
#' @param ... additional argument, currently not used.
#'
#' @return a \code{\link{ggtree}} plot
#' 
#' @seealso
#' \code{\link[mia:splitByRanks]{splitByRanks}}
#'
#' @name plotTree
#'
#' @examples
#' library(scater)
#' library(mia)
#' # preparation of some data
#' data(GlobalPatterns)
#' altExps(GlobalPatterns) <- splitByRanks(GlobalPatterns)
#' altExp(GlobalPatterns,"Genus") <- addPerFeatureQC(altExp(GlobalPatterns,"Genus"))
#' rowData(altExp(GlobalPatterns,"Genus"))$log_mean <- 
#'   log(rowData(altExp(GlobalPatterns,"Genus"))$mean)
#' rowData(altExp(GlobalPatterns,"Genus"))$detected <- 
#'   rowData(altExp(GlobalPatterns,"Genus"))$detected / 100
#' top_taxa <- getTopTaxa(altExp(GlobalPatterns,"Genus"),
#'                        method="mean",
#'                        top=100L,
#'                        abund_values="counts")
#' #
#' plotRowTree(altExp(GlobalPatterns,"Genus")[top_taxa,],
#'             tip_colour_by = "log_mean",
#'             tip_size_by = "detected")
#'
#' # plot with tip labels
#' plotRowTree(altExp(GlobalPatterns,"Genus")[top_taxa,],
#'             tip_colour_by = "log_mean",
#'             tip_size_by = "detected",
#'             show_label = TRUE)
#' # plot with selected labels
#' labels <- c("Genus:Providencia" = TRUE, "Genus:Morganella" = FALSE,
#'             "0.961.60" = TRUE)
#' plotRowTree(altExp(GlobalPatterns,"Genus")[top_taxa,],
#'             tip_colour_by = "log_mean",
#'             tip_size_by = "detected",
#'             show_label = labels,
#'             layout="rectangular")
#'             
#' # plot with labeled edges
#' plotRowTree(altExp(GlobalPatterns,"Genus")[top_taxa,],
#'             edge_colour_by = "Kingdom",
#'             edge_size_by = "detected",
#'             tip_colour_by = "log_mean")
#' 
#' # aggregating data over the taxonomic levels for plotting a taxonomic tree
#' # please note that the original tree of GlobalPatterns is dropped by
#' # unsplitByRanks
#' altExps(GlobalPatterns) <- splitByRanks(GlobalPatterns)
#' altExps(GlobalPatterns) <- lapply(altExps(GlobalPatterns), addPerFeatureQC)
#' altExps(GlobalPatterns) <-
#'    lapply(altExps(GlobalPatterns),
#'           function(y){
#'               rowData(y)$log_mean <- log(rowData(y)$mean)
#'               rowData(y)$detected <- rowData(y)$detected / 100
#'               y
#'           })
#' x <- unsplitByRanks(GlobalPatterns)
#' x <- addTaxonomyTree(x)
#' plotRowTree(x,
#'             edge_colour_by = "Kingdom",
#'             edge_size_by = "detected",
#'             tip_colour_by = "log_mean",
#'             node_colour_by = "log_mean")
NULL

#' @rdname plotTree
setGeneric("plotRowTree", signature = c("object"),
           function(object, ...)
               standardGeneric("plotRowTree"))
#' @rdname plotTree
setGeneric("plotColTree", signature = c("object"),
           function(object, ...)
               standardGeneric("plotColTree"))

.check_tree_plot_switches <- function(relabel_tree, show_label, add_legend){
    if(!.is_a_bool(relabel_tree)){
        stop("'relabel_tree' must be either TRUE or FALSE.", call. = FALSE)
    }
    if(!.is_a_bool(show_label)){
        if(!is.logical(show_label) || is.null(names(show_label))){
            stop("'show_label' must be either TRUE or FALSE or named logical ",
                 "vector. Names should match the label of the tree.",
                 call. = FALSE)
        }
    }
    if(!.is_a_bool(add_legend)){
        stop("'add_legend' must be either TRUE or FALSE.", call. = FALSE)
    }
}

#' @rdname plotTree
#' @export
setMethod("plotColTree", signature = c(object = "TreeSummarizedExperiment"),
    function(object,
             relabel_tree = FALSE,
             show_label = FALSE,
             add_legend = TRUE,
             layout = "circular",
             edge_colour_by = NULL,
             edge_size_by = NULL,
             tip_colour_by = NULL,
             tip_shape_by = NULL,
             tip_size_by = NULL,
             node_colour_by = NULL,
             node_shape_by = NULL,
             node_size_by = NULL,
             by_exprs_values = "counts",
             other_fields = list(),
             ...){
        # input check
        if(is.null(colTree(object))){
          stop("colTree(object) is empty.", call. = FALSE)
        }
        .check_tree_plot_switches(relabel_tree = relabel_tree,
                                  show_label = show_label,
                                  add_legend = add_legend)
        #
        tree <- .get_trimed_tree(object, type = "column")
        relabel_out <- relabel_object_tree(object, tree, type = "column",
                                           relabel = relabel_tree)
        object <- relabel_out$object
        tree <- relabel_out$tree
        tree_data <- .get_tree_data(combineTreeData(tree, other_fields))
        #
        vis_out <- .incorporate_tree_vis(tree_data,
                                         se = object,
                                         edge_colour_by = edge_colour_by,
                                         edge_size_by = edge_size_by,
                                         tip_colour_by = tip_colour_by,
                                         tip_shape_by = tip_shape_by,
                                         tip_size_by = tip_size_by,
                                         node_colour_by = node_colour_by,
                                         node_shape_by = node_shape_by,
                                         node_size_by = node_size_by,
                                         by_exprs_values = by_exprs_values)
        tree_data <- vis_out$df
        edge_colour_by <- vis_out$edge_colour_by
        edge_size_by <- vis_out$edge_size_by
        colour_by <- vis_out$colour_by
        shape_by <- vis_out$shape_by
        size_by <- vis_out$size_by
        show_tips <- any(!vapply(c(tip_colour_by, tip_shape_by, tip_size_by),
                                 is.null, logical(1)))
        show_nodes <- any(!vapply(c(node_colour_by, node_shape_by, node_size_by),
                                  is.null, logical(1)))
        #
        .tree_plotter(tree_data,
                      layout = layout,
                      add_legend = add_legend,
                      show_label = show_label,
                      show_tips = show_tips,
                      show_nodes = show_nodes,
                      edge_colour_by = edge_colour_by,
                      edge_size_by = edge_size_by,
                      colour_by = colour_by,
                      shape_by = shape_by,
                      size_by = size_by,
                      ...)
    }
)
#' @rdname plotTree
#' @export
setMethod("plotRowTree", signature = c(object = "TreeSummarizedExperiment"),
    function(object,
             relabel_tree = FALSE,
             show_label = FALSE,
             add_legend = TRUE,
             layout = "circular",
             edge_colour_by = NULL,
             edge_size_by = NULL,
             tip_colour_by = NULL,
             tip_shape_by = NULL,
             tip_size_by = NULL,
             node_colour_by = NULL,
             node_shape_by = NULL,
             node_size_by = NULL,
             by_exprs_values = "counts",
             other_fields = list(),
             ...){
        # input check
        if(is.null(rowTree(object))){
          stop("rowTree(object) is empty.", call. = FALSE)
        }
        .check_tree_plot_switches(relabel_tree = relabel_tree,
                                  show_label = show_label,
                                  add_legend = add_legend)
        #
        tree <- .get_trimed_tree(object, type = "row")
        relabel_out <- relabel_object_tree(object, tree, type = "row",
                                           relabel = relabel_tree)
        object <- relabel_out$object
        tree <- relabel_out$tree
        tree_data <- .get_tree_data(combineTreeData(tree, other_fields))
        #
        vis_out <- .incorporate_tree_vis(tree_data,
                                         se = object,
                                         edge_colour_by = edge_colour_by,
                                         edge_size_by = edge_size_by,
                                         tip_colour_by = tip_colour_by,
                                         tip_shape_by = tip_shape_by,
                                         tip_size_by = tip_size_by,
                                         node_colour_by = node_colour_by,
                                         node_shape_by = node_shape_by,
                                         node_size_by = node_size_by,
                                         by_exprs_values = by_exprs_values)
        tree_data <- vis_out$df
        edge_colour_by <- vis_out$edge_colour_by
        edge_size_by <- vis_out$edge_size_by
        colour_by <- vis_out$colour_by
        shape_by <- vis_out$shape_by
        size_by <- vis_out$size_by
        show_tips <- any(!vapply(c(tip_colour_by, tip_shape_by, tip_size_by),
                                 is.null, logical(1)))
        show_nodes <- any(!vapply(c(node_colour_by, node_shape_by, node_size_by),
                                  is.null, logical(1)))
        #
        .tree_plotter(tree_data,
                      layout = layout,
                      add_legend = add_legend,
                      show_label = show_label,
                      show_tips = show_tips,
                      show_nodes = show_nodes,
                      edge_colour_by = edge_colour_by,
                      edge_size_by = edge_size_by,
                      colour_by = colour_by,
                      shape_by = shape_by,
                      size_by = size_by,
                      ...)
    }
)

#' @importFrom ape keep.tip
.get_trimed_tree <- function(object, type = c("row","columns")){
    type <- match.arg(type)
    tree_FUN <- switch(type, row = rowTree, column = colTree, stop("."))
    links_FUN <- switch(type, row = rowLinks, column = colLinks, stop("."))
    tree <- tree_FUN(object)
    links <- links_FUN(object)
    #
    tree <- ape::keep.tip(tree, unique(links$nodeNum[links$isLeaf]))
    tree
}

relabel_object_tree <- function(object, tree, type = c("row","columns"),
                                relabel = FALSE){
    type <- match.arg(type)
    links_FUN <- switch(type, row = rowLinks, column = colLinks, stop("."))
    dimnames_FUN <- switch(type, row = rownames, column = colnames, stop("."))
    dimnames_rep_FUN <- switch(type, row = rownames, column = colnames, stop("."))
    links <- links_FUN(object)
    dimnames <- dimnames_FUN(object)
    # change labels
    leaf_nodes <- unique(links$nodeNum[links$isLeaf])
    # in case leafs are not the first in the list
    leaf_seq <- seq_along(leaf_nodes)
    leaf_nodes <- leaf_seq[leaf_seq[order(leaf_nodes)]]
    if(relabel ||
       any(tree$tip.label[leaf_nodes] != dimnames[links$isLeaf]) ||
       anyDuplicated(dimnames[links$isLeaf])){
        if(type == "column"){
            stop(".") # Now taxonomic info on the cols
        }
        new_tip_labels <- getTaxonomyLabels(object, with_type = TRUE,
                                            resolve_loops = TRUE)
    } else {
        new_tip_labels <- dimnames
    }
    tree$tip.label[leaf_nodes] <- new_tip_labels[links$isLeaf]
    if(type == "row"){
        rownames(object) <- new_tip_labels
    } else {
        colnames(object) <- new_tip_labels
    }
    list(object = object, tree = tree)
}


#' @importFrom tibble tibble
.get_feature_info <- function(by, se, FUN, exprs_values){
    feature_info <- FUN(se, by = by, exprs_values = exprs_values)
    feature_info <- tibble(!!sym(feature_info$name) := feature_info$value)
    feature_info
}

TIP_VARIABLES <- c("tip_colour_by", "tip_shape_by", "tip_size_by")
NODE_VARIABLES <- c("node_colour_by", "node_shape_by", "node_size_by")

#' @importFrom scater retrieveFeatureInfo retrieveCellInfo
#' @importFrom dplyr bind_cols mutate relocate
#' @importFrom tibble rownames_to_column
.incorporate_tree_vis <- function(tree_data,
                                  se,
                                  edge_colour_by,
                                  edge_size_by,
                                  tip_colour_by,
                                  tip_shape_by,
                                  tip_size_by,
                                  node_colour_by,
                                  node_shape_by,
                                  node_size_by,
                                  by_exprs_values = "counts",
                                  type = c("row","column")){
    
    variables <- c(edge_colour_by = edge_colour_by,
                   edge_size_by = edge_size_by,
                   tip_colour_by = tip_colour_by,
                   tip_shape_by = tip_shape_by,
                   tip_size_by = tip_size_by,
                   node_colour_by = node_colour_by,
                   node_shape_by = node_shape_by,
                   node_size_by = node_size_by)
    edge_colour_by <- NULL
    edge_size_by <- NULL
    colour_by <- NULL
    shape_by <- NULL
    size_by <- NULL
    if(!is.null(variables)){
        # remove any variables values, which are already available and
        # rename columns by their usage
        cn <- colnames(tree_data)
        cn_data <- cn[!(cn %in% c(DEFAULT_TREE_DATA_COLS))]
        if(length(cn_data) > 0L){
            f <- variables %in% cn_data
            if(any(f)){
                tree_data <- tree_data[,c(DEFAULT_TREE_DATA_COLS,variables[f])]
                # rename columns by their usage and merge by node type
                colnames(tree_data) <- c(DEFAULT_TREE_DATA_COLS,names(variables)[f])
                # mirror back variable name
                for(i in variables[f]){
                    var_name <- gsub("tip_|node_","",names(variables)[f][i])
                    assign(var_name,
                           paste(get(var_name),
                                 variables[f][i],
                                 collapse = " & ", sep = ""))
                }
                variables <- variables[!f]
            }
        }
        if(length(variables) > 0L){
            type <- match.arg(type)
            type_FUN <- switch(type,
                               row = scater::retrieveFeatureInfo,
                               column = scater::retrieveCellInfo)
            feature_info <- vector(mode = "list", length = length(variables))
            for(i in seq_along(variables)){
                # get data
                feature_info[[i]] <-
                    .get_feature_info(variables[i], se = se,
                                      FUN = type_FUN, exprs_values = by_exprs_values)
                # mirror back variable name, if a partial match was used
                var_name <- gsub("tip_|node_","",names(variables)[i])
                assign(var_name,
                       paste(get(var_name),
                             colnames(feature_info[[i]]),
                             collapse = " & ", sep = ""))
            }
            feature_info <- bind_cols(feature_info)
            # rename columns by their usage
            colnames(feature_info) <- names(variables)
            feature_info <- feature_info %>%
                mutate(label = rownames(se)) %>%
                relocate("label")
            tree_data <- .merge_tree_vis_data(tree_data, feature_info)
        }
        tree_data <- .merge_tip_node_tree_data(tree_data)
    }
    return(list(df = tree_data,
                edge_colour_by = edge_colour_by,
                edge_size_by = edge_size_by,
                colour_by = colour_by,
                shape_by = shape_by,
                size_by = size_by))
}

.merge_tip_node_tree_data <- function(tree_data){
    # setup variables for ordering and order tree_data
    is_leaf <- !(tree_data$node %in% unique(tree_data$parent))
    bak_o <- tree_data$node
    o <- order(is_leaf)
    tree_data <- tree_data[o,]
    is_leaf_o <- !(tree_data$node %in% unique(tree_data$parent))
    # default values
    edge_colour_by <- NULL
    colour_by <- NULL
    shape_by <- NULL
    size_by <- NULL
    #
    cn <- colnames(tree_data)
    if(all(c("tip_colour_by","node_colour_by") %in% cn)){
        colour_by <- c(tree_data$node_colour_by[!is_leaf_o],
                       tree_data$tip_colour_by[is_leaf_o])
    } else if("tip_colour_by" %in% cn) {
        colour_by <- tree_data$tip_colour_by
    } else if("node_colour_by" %in% cn) {
        colour_by <- tree_data$node_colour_by
    }
    if(all(c("tip_shape_by","node_shape_by") %in% cn)){
        shape_by <- c(tree_data$node_shape_by[!is_leaf_o],
                      tree_data$tip_shape_by[is_leaf_o])
    } else if("tip_shape_by" %in% cn) {
        shape_by <- tree_data$tip_shape_by
    } else if("node_shape_by" %in% cn) {
        shape_by <- tree_data$node_shape_by
    }
    if(all(c("tip_size_by","node_size_by") %in% cn)){
        size_by <- c(tree_data$node_size_by[!is_leaf_o],
                     tree_data$tip_size_by[is_leaf_o])
    } else if("tip_size_by" %in% cn) {
        size_by <- tree_data$tip_size_by
    } else if("node_size_by" %in% cn) {
        size_by <- tree_data$node_size_by
    }
    #
    tree_data <- tree_data[,cn[!grepl("tip_|node_",cn)]]
    tree_data$colour_by <- colour_by
    tree_data$shape_by <- shape_by
    tree_data$size_by <- size_by
    # return tree_data with original ordering
    tree_data[match(bak_o,tree_data$node),]
}

.merge_tree_vis_data <- function(tree_data, feature_info){
    if(anyDuplicated(tree_data$label) || anyDuplicated(feature_info$label)){
        stop(".")
    }
    tree_data <- tree_data %>%
        left_join(feature_info, by = "label")
    tree_data
}

#' @importFrom tidytree as.treedata
#' @importFrom ggplot2 scale_size_identity
#' @importFrom ggtree ggtree geom_tree geom_tippoint geom_nodepoint groupOTU
.tree_plotter <- function(object,
                          layout = "circular",
                          add_legend = TRUE,
                          show_label = TRUE,
                          show_tips = FALSE,
                          show_nodes = FALSE,
                          edge_colour_by = NULL,
                          edge_size_by = NULL,
                          colour_by = NULL,
                          shape_by = NULL,
                          size_by = NULL,
                          tree_edge_alpha = 1,
                          tree_edge_size = NULL,
                          point_alpha = 1,
                          point_size = 2){
    # assemble arg list
    point_out <- .get_point_args(colour_by,
                                 shape_by,
                                 size_by,
                                 alpha = point_alpha,
                                 size = point_size)
    edge_out <- .get_edge_args(edge_colour_by,
                               edge_size_by,
                               alpha = tree_edge_alpha,
                               size = tree_edge_size)
    object <- .na_replace_from_plot_data(object,
                                         edge_size_by,
                                         shape_by,
                                         size_by)
    # start plotting
    if (!is.null(edge_colour_by)) {
        object <- groupOTU(object, split(object$node, object$edge_colour_by))
    }
    plot_out <- ggtree(tidytree::as.treedata(object), layout = layout) +
        do.call(geom_tree, edge_out$args)
    if (!is.null(edge_size_by)) {
        plot_out <- plot_out + 
            scale_size_identity(guide="legend")
    }
    tip_point_FUN <- geom_tippoint
    node_point_FUN <- geom_nodepoint
    if(show_tips){
        plot_out <- plot_out +
            do.call(tip_point_FUN, point_out$args)
    }
    if(show_nodes){
        plot_out <- plot_out +
            do.call(node_point_FUN, point_out$args)
    }
    plot_out <- .add_tree_labels(plot_out, show_label)
    if(!is.null(colour_by)){
        plot_out <- .resolve_plot_colours(plot_out,
                                          object$colour_by,
                                          colour_by,
                                          fill = point_out$fill)
    }
    plot_out <- .add_extra_guide(plot_out, shape_by, size_by)
    if (!add_legend) {
        plot_out <- plot_out + theme(legend.position = "none")
    }
    plot_out
}

#' @importFrom ggtree geom_tiplab geom_nodelab
.add_tree_labels <- function(plot_out, show_label){
    if(length(show_label) == 1L && show_label){
        plot_out <- plot_out +
            geom_tiplab(offset = .1)
    } else if(length(show_label) > 1L && any(show_label)) {
        is_leaf <- !(plot_out$data$node %in% unique(plot_out$data$parent))
        label <- plot_out$data$label
        m <- match(label,names(show_label))
        m <- m[!is.na(m)]
        m <- m[show_label[m]]
        tmp <- rep("",length(label))
        names(tmp) <- label
        label <- tmp
        label[names(show_label)[m]] <- names(show_label)[m]
        plot_out$data$label <- label
        # add tip and node labels
        plot_out <- plot_out +
            geom_tiplab(offset = .1) +
            #geom_nodelab(offset = .1)
            geom_nodelab()
    }
    plot_out
}
