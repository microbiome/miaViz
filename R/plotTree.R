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
#'   the output of \code{getTaxonomyLabels(object, with_rank = TRUE)}? This
#'   has consequences on how data can be merged from \code{other_fields}. 
#'   (default: \code{relabel_tree = FALSE})
#'   
#' @param show_label logical scalar or character vector. Should tip labels
#'   be plotted or if a logical vector is provided, which labels should be 
#'   shown? Only values corresponding to actual labels will be plotted. 
#'   (default: \code{show_label = FALSE})
#'   
#' @param show_highlights logical scalar or character vector. Should highlights
#'   around nodes and its decendents be plotted? If set \code{TRUE} this will be
#'   limit do direct decendents of the rootnode. If a character vector is
#'   provided, which nodes should be shown? Only values corresponding to actual
#'   node labels will be plotted. (default: \code{show_highlights = FALSE})
#'   
#' @param colour_highlights Should the highlights be colour differently?
#'   If \code{show_highlights = TRUE}, \code{colour_highlights} will be set to
#'   \code{TRUE} as default. (default: \code{colour_highlights = FALSE})
#'   
#' @param add_legend logical scalar. Should legends be plotted? 
#'   (default: \code{add_legend = TRUE})
#'   
#' @param layout layout for the plotted tree. See 
#'   \code{\link[ggtree:ggtree]{ggtree}} for details.
#'   
#' @param edge_colour_by Specification of a column metadata field or a feature 
#'   to colour tree edges by, see the by argument in 
#'   \code{\link[scater:retrieveCellInfo]{?retrieveCellInfo}} for possible 
#'   values.
#'   
#' @param edge_size_by Specification of a column metadata field or a feature 
#'   to size tree edges by, see the by argument in 
#'   \code{\link[scater:retrieveCellInfo]{?retrieveCellInfo}} for possible 
#'   values.
#'   
#' @param tip_colour_by Specification of a column metadata field or a feature to
#'   colour tree tips by, see the by argument in 
#'   \code{\link[scater:retrieveCellInfo]{?retrieveCellInfo}} for possible 
#'   values.
#'   
#' @param tip_shape_by Specification of a column metadata field or a feature to
#'   shape tree tips by, see the by argument in 
#'   \code{\link[scater:retrieveCellInfo]{?retrieveCellInfo}} for possible 
#'   values.
#'   
#' @param tip_size_by Specification of a column metadata field or a feature to
#'   size tree tips by, see the by argument in 
#'   \code{\link[scater:retrieveCellInfo]{?retrieveCellInfo}} for possible 
#'   values.
#'   
#' @param node_colour_by Specification of a column metadata field or a feature to
#'   colour tree nodes by. Must be a field from \code{other_fields}.
#'   
#' @param node_shape_by Specification of a column metadata field or a feature to
#'   shape tree nodes by. Must be a field from \code{other_fields}.
#'   
#' @param node_size_by Specification of a column metadata field or a feature to
#'   size tree nodes by. Must be a field from \code{other_fields}.
#'   
#' @param by_exprs_values A string or integer scalar specifying which assay to
#'   obtain expression values from, for use in point aesthetics - see the 
#'   \code{exprs_values} argument in 
#'   \code{\link[scater:retrieveCellInfo]{?retrieveCellInfo}}.
#'   
#' @param other_fields a \code{data.frame} or coercible to one, with at least 
#'   one type of id information. See details in
#'   \code{\link[=treeData]{treeData}}.
#'   
#' @param ... additional arguments for plotting.
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
#'    rowData(altExp(GlobalPatterns,"Genus"))$detected / 100
#' top_genus <- getTopTaxa(altExp(GlobalPatterns,"Genus"),
#'                         method="mean",
#'                         top=100L,
#'                         abund_values="counts")
#' #
#' x <- altExp(GlobalPatterns,"Genus")
#' plotRowTree(x[rownames(x) %in% top_genus,],
#'             tip_colour_by = "log_mean",
#'             tip_size_by = "detected")
#' 
#' # plot with tip labels
#' plotRowTree(x[rownames(x) %in% top_genus,],
#'             tip_colour_by = "log_mean",
#'             tip_size_by = "detected",
#'             show_label = TRUE)
#' # plot with selected labels
#' labels <- c("Genus:Providencia" = TRUE, "Genus:Morganella" = FALSE,
#'             "0.961.60" = TRUE)
#' plotRowTree(x[rownames(x) %in% top_genus,],
#'             tip_colour_by = "log_mean",
#'             tip_size_by = "detected",
#'             show_label = labels,
#'             layout="rectangular")
#' 
#' # plot with labeled edges
#' plotRowTree(x[rownames(x) %in% top_genus,],
#'             edge_colour_by = "Phylum",
#'             tip_colour_by = "log_mean")
#' # if edges are sized, colours might disappear depending on plotting device
#' plotRowTree(x[rownames(x) %in% top_genus,],
#'             edge_colour_by = "Phylum",
#'             edge_size_by = "detected",
#'             tip_colour_by = "log_mean")
#' 
#' # aggregating data over the taxonomic levels for plotting a taxonomic tree
#' # please note that the original tree of GlobalPatterns is dropped by
#' # unsplitByRanks
#' altExps(GlobalPatterns) <- splitByRanks(GlobalPatterns)
#' top_phyla <- getTopTaxa(altExp(GlobalPatterns,"Phylum"),
#'                         method="mean",
#'                         top=10L,
#'                         abund_values="counts")
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
#' plotRowTree(x[rowData(x)$Phylum %in% top_phyla,],
#'             edge_colour_by = "Phylum",
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

.check_tree_plot_switches <- function(relabel_tree,
                                      show_label, 
                                      show_highlights,
                                      colour_highlights,
                                      add_legend){
    if(!.is_a_bool(relabel_tree)){
        stop("'relabel_tree' must be either TRUE or FALSE.", call. = FALSE)
    }
    if(!.is_a_bool(show_label)){
        if( (!is.logical(show_label) && !is.character(show_label)) ||
            is.null(show_label)){
            stop("'show_label' must be either TRUE or FALSE or character ",
                 "vector. Values should match the label of the tree.",
                 call. = FALSE)
        }
    }
    if(!.is_a_bool(show_highlights)){
        if( (!is.logical(show_highlights) && !is.character(show_highlights)) ||
            is.null(show_highlights)){
            stop("'show_highlights' must be either TRUE or FALSE or character ",
                 "vector. Values should match the label of the tree.",
                 call. = FALSE)
        }
    }
    if(!.is_a_bool(colour_highlights)){
        stop("'colour_highlights' must be either TRUE or FALSE.",
             call. = FALSE)
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
             show_highlights = FALSE,
             colour_highlights = FALSE,
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
                                  show_highlights = show_highlights,
                                  colour_highlights = colour_highlights,
                                  add_legend = add_legend)
        #
        tree_out <- .get_trimed_object_and_tree(object, type = "column",
                                                relabel = relabel_tree)
        object <- tree_out$object
        tree <- tree_out$tree
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
                                         by_exprs_values = by_exprs_values,
                                         type = "column")
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
                      show_highlights = show_highlights,
                      colour_highlights = colour_highlights,
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
             show_highlights = FALSE,
             colour_highlights = FALSE,
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
                                  show_highlights = show_highlights,
                                  colour_highlights = colour_highlights,
                                  add_legend = add_legend)
        #
        tree_out <- .get_trimed_object_and_tree(object, type = "row",
                                                relabel = relabel_tree)
        object <- tree_out$object
        tree <- tree_out$tree
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
                                         by_exprs_values = by_exprs_values,
                                         type = "row")
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
                      show_highlights = show_highlights,
                      colour_highlights = colour_highlights,
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

#' @importFrom ape keep.tip as.phylo
#' @importFrom tidytree as_tibble 
.get_trimed_object_and_tree <- function(object, type = c("row","columns"),
                                        relabel = FALSE){
    type <- match.arg(type)
    tree_FUN <- switch(type, row = rowTree, column = colTree, stop("."))
    links_FUN <- switch(type, row = rowLinks, column = colLinks, stop("."))
    dimnames_FUN <- switch(type, row = rownames, column = colnames, stop("."))
    tree <- tree_FUN(object)
    links <- links_FUN(object)
    #
    tips <- sort(setdiff(tree$edge[, 2], tree$edge[, 1]))
    drop_tip <- tips[!(tips %in% unique(links$nodeNum[links$isLeaf]))]
    oldTree <- tree
    newTree <- ape::drop.tip(oldTree, tip = drop_tip, collapse.singles = FALSE)
    track <- trackNode(oldTree)
    track <- ape::drop.tip(track, tip = drop_tip, collapse.singles = FALSE)
    #
    oldAlias <- links$nodeLab_alias
    newNode <- convertNode(tree = track, node = oldAlias)
    newAlias <- convertNode(tree = newTree, node = newNode)
    if(type == "row"){
        object <- changeTree(x = object, rowTree = newTree, rowNodeLab = newAlias)
    } else {
        object <- changeTree(x = object, colTree = newTree, colNodeLab = newAlias)
    }
    #
    tree <- tree_FUN(object)
    links <- links_FUN(object)
    dimnames <- dimnames_FUN(object)
    #
    tree_data <- as_tibble(newTree)
    m <- match(links$nodeNum,tree_data$node)
    node_labels <- tree_data$label[m]
    if(relabel || 
       !all(node_labels %in% dimnames)){
        new_node_labels <- getTaxonomyLabels(object, with_rank = TRUE,
                                             resolve_loops = TRUE)
        if(type == "row"){
            rownames(object) <- new_node_labels
        } else {
            colnames(object) <- new_node_labels
        }
        tree_data$label[m] <- new_node_labels
        tree <- as.phylo(tree_data)
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

.get_new_var_name_value <- function(var_name_value, add){
    if(!is.null(var_name_value) &&
       add != var_name_value){
        new_var_name_value <-
            paste0(var_name_value,
                   ifelse(is.null(var_name_value),"", " & "),
                   add)
    } else {
        new_var_name_value <- add
    }
    new_var_name_value
}

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
                           .get_new_var_name_value(get(var_name),
                                                   variables[f][i]))
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
                       .get_new_var_name_value(get(var_name),
                                               colnames(feature_info[[i]])))
                # rename columns by their usage
                colnames(feature_info[[i]]) <- names(variables[i])
            }
            feature_info <- bind_cols(feature_info)
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
#'   theme_tree
.tree_plotter <- function(object,
                          layout = "circular",
                          add_legend = TRUE,
                          show_label = FALSE,
                          show_highlights = FALSE,
                          colour_highlights = FALSE,
                          show_tips = FALSE,
                          show_nodes = FALSE,
                          edge_colour_by = NULL,
                          edge_size_by = NULL,
                          colour_by = NULL,
                          shape_by = NULL,
                          size_by = NULL,
                          line_alpha = 1,
                          line_width = NULL,
                          line_width_range = c(0.5,3),
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
                               alpha = line_alpha,
                               size = line_width)
    #
    if (!is.null(edge_colour_by) && anyNA(object$edge_colour_by)) {
        object <- groupOTU(object, split(object$node, object$edge_colour_by),
                           group_name = "group")
        f_zero <- object$group != 0
        f_zero <- f_zero[!is.na(f_zero)]
        object$edge_colour_by[f_zero] <- as.character(object$group[f_zero])
    }
    object <- .na_replace_from_plot_data(object,
                                         edge_size_by,
                                         shape_by,
                                         size_by)
    # start plotting
    plot_out <- ggtree(tidytree::as.treedata(object), layout = layout)
    # add highlights
    plot_out <- .add_tree_highlights(plot_out, show_highlights,
                                     colour_highlights)
    # add tree and adjust edges
    plot_out <- plot_out +
        do.call(geom_tree, edge_out$args) + 
        theme_tree()
    if (!is.null(edge_size_by)) {
        if(is.numeric(object$edge_size_by)){
            SIZEFUN <- scale_size_continuous
        } else {
            SIZEFUN <- scale_size_discrete
        }
        plot_out <- .add_extra_guide_tree(plot_out, edge_size_by) +
            SIZEFUN(range = line_width_range)
    }
    # add tip and node points
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
    # add tip and node labels
    plot_out <- .add_tree_labels(plot_out, show_label)
    # adjust edge colours
    if(!is.null(edge_colour_by)){
        plot_out <- .resolve_plot_colours(plot_out,
                                          object$edge_colour_by,
                                          edge_colour_by,
                                          na.translate = FALSE)
    }
    # adjust point colours
    if(!is.null(colour_by)){
        plot_out <- .resolve_plot_colours(plot_out,
                                          object$colour_by,
                                          colour_by,
                                          fill = point_out$fill,
                                          na.translate = FALSE)
    }
    # add additional guides
    plot_out <- .add_extra_guide(plot_out, shape_by, size_by)
    # add theme
    plot_out <- .theme_plotTree(plot_out)
    # optionally hide legends
    if (!add_legend) {
        plot_out <- plot_out +
            theme(legend.position = "none")
    }
    plot_out
}

#' @importFrom ggtree geom_highlight
#' @importFrom ggnewscale new_scale_fill new_scale_colour
#' @importFrom tidytree rootnode
.add_tree_highlights <- function(plot_out, show_highlights = FALSE,
                                 colour_highlights = FALSE){
    if(length(show_highlights) == 1L && is.logical(show_highlights) &&
       show_highlights){
        rootnode <- tidytree::rootnode(plot_out$data)
        f <- plot_out$data$parent == rootnode$node &
            plot_out$data$node != rootnode$node
        show_highlights <- plot_out$data$label[f]
        colour_highlights <- TRUE
    }
    if(length(show_highlights) >= 1L) {
        m <- plot_out$data$label %in% show_highlights
        highlights <- plot_out$data$node[m]
        if(length(highlights) > 0L){
            # add tip and node labels
            hl_args <- .get_hightlight_args(highlights,
                                            colour_highlights)
            plot_out <- plot_out +
                do.call(geom_highlight, hl_args$args)
            plot_out <- .resolve_plot_colours(plot_out,
                                              plot_out$data$label[m],
                                              "",
                                              fill = TRUE)
            plot_out <- plot_out + 
                new_scale_fill() +
                new_scale_colour()
        }
    }
    plot_out
}

#' @importFrom ggtree geom_tiplab geom_nodelab
.add_tree_labels <- function(plot_out, show_label){
    if(length(show_label) == 1L && is.logical(show_label) &&
       show_label){
        plot_out <- plot_out +
            geom_tiplab(offset = .02)
    } else if(length(show_label) > 1L) {
        m <- plot_out$data$label %in% show_label
        nodes <- plot_out$data$node[m]
        f_tip <- plot_out$data$node %in% nodes & plot_out$data$isTip
        f_node <- plot_out$data$node %in% nodes & !plot_out$data$isTip
        if(any(f_tip)){
            # add tip labels
            plot_out <- plot_out +
                geom_tiplab(mapping = aes_string(subset = f_tip),
                            offset = .02)
        }
        if(any(f_node)){
            # add node labels
            plot_out <- plot_out +
                geom_nodelab(mapping = aes_string(subset = f_node))
        }
    }
    plot_out
}

.theme_plotTree <- function(plot){
    plot + 
        theme(legend.background = element_rect(fill = "transparent",colour = NA),
              legend.box.background = element_rect(fill = "transparent",colour = NA),
              panel.background = element_rect(fill = "transparent",colour = NA),
              plot.background = element_rect(fill = "transparent",colour = NA))
}
