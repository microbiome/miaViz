#' Plotting tree information enriched with information
#'
#' Based on the stored data in a \code{TreeSummarizedExperiment} a tree can
#' be plotted. From the \code{rowData}, the \code{assays} as well as the
#' \code{colData} information can be taken for enriching the tree plots with
#' additional information.
#'
#' @param x a
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}.
#'
#' @param tree.name \code{Character scalar}. Specifies a rowTree/colTree from
#' \code{x}. (Default: \code{tree.name = "phylo"})
#'
#' @param tree_name Deprecated. Use \code{tree.name} instead.
#'
#' @param relabel.tree \code{Logical scalar}. Should the tip labels be relabeled
#' using the output of \code{getTaxonomyLabels(x, with_rank = TRUE)}?
#' (Default: \code{FALSE})
#'
#' @param relabel_tree Deprecated. Use \code{relavel.tree} instead.
#'
#' @param order.tree \code{Logical scalar}. Should the tree be ordered based on
#'   alphabetic order of taxonomic levels? (Default: \code{FALSE})
#'
#' @param order_tree Deprecated. Use \code{order.tree} instead.
#'
#' @param levels.rm \code{Logical scalar}. Should taxonomic level information
#'   be removed from labels? (Default: \code{FALSE})
#'
#' @param remove_levels Deprecated. Use \code{levels.rm} instead.
#'
#' @param show.label,show.highlights,show.highlight.label,abbr.label
#'   \code{logical} (scalar), \code{integer} or \code{character} vector. If a
#'   \code{logical} scalar is given, should tip labels be plotted or if a
#'   logical vector is provided, which labels should be shown? If an
#'   \code{integer} or \code{character} vector is provided, it will be converted
#'   to a logical vector. The \code{integer} values must be in the range of 1
#'   and number of nodes, whereas the values of a \code{character} vector must
#'   match values of the \code{label} column in the node data. In case of a
#'   \code{character} vector only values corresponding to actual labels will be
#'   plotted and if no labels are provided no labels will be shown. (Default:
#'   \code{FALSE})
#'
#' @param show_label,show_highlights,show_highlight_label,abbr_label Deprecated.
#' Use \code{show.label, show.highlights, show.highlight.label, abbr_label}
#' instead.
#'
#' @param add.legend \code{Logical scalar}. Should legends be plotted?
#'   (Default: \code{TRUE})
#'
#' @param add_legend Deprecated. Use \code{add.legend} instead.
#'
#' @param layout layout for the plotted tree. See
#'   \code{\link[ggtree:ggtree]{ggtree}} for details.
#'
#' @param edge.colour.by \code{Character scalar}. Specification of a column
#' metadata field or a feature to colour tree edges by, see the by argument in
#' \code{\link[scater:retrieveCellInfo]{?retrieveCellInfo}} for possible
#' values.
#'
#' @param edge_colour_by Deprecated. Use \code{edge.colour.by} instead.
#'
#' @param edge.size.by \code{Character scalar}. Specification of a column
#' metadata field or a feature to size tree edges by, see the by argument in
#' \code{\link[scater:retrieveCellInfo]{?retrieveCellInfo}} for possible
#' values. (Default: \code{NULL})
#'
#' @param edge_size_by Deprecated. Use \code{edge.size.by} instead.
#'
#' @param tip.colour.by \code{Character scalar}. Specification of a column
#' metadata field or a feature to colour tree tips by, see the by argument in
#' \code{\link[scater:retrieveCellInfo]{?retrieveCellInfo}} for possible
#' values. (Default: \code{NULL})
#'
#' @param tip_colour_by Deprecated. Use \code{tip.colour.by} instead.
#'
#' @param tip.shape.by \code{Character scalar}. Specification of a column
#' metadata field or a feature to shape tree tips by, see the by argument in
#' \code{\link[scater:retrieveCellInfo]{?retrieveCellInfo}} for possible
#' values. (Default: \code{NULL})
#'
#' @param tip_shape_by Deprecated. Use \code{tip.shape.by} isntead.
#'
#' @param tip.size.by \code{Character scalar}. Specification of a column
#' metadata field or a feature to size tree tips by, see the by argument in
#' \code{\link[scater:retrieveCellInfo]{?retrieveCellInfo}} for possible
#' values. (Default: \code{NULL})
#'
#' @param tip_size_by Deprecated. Use \code{tip.size.by} instead.
#'
#' @param node.colour.by \code{Character scalar}. Specification of a column
#' metadata field or a feature to colour tree nodes by. Must be a field from
#' \code{other.fields}. (Default: \code{NULL})
#'
#' @param node_colour_by Deprecated. Use \code{node.colour.by} instead.
#'
#' @param node.shape.by \code{Character scalar}. Specification of a column
#' metadata field or a feature to shape tree nodes by. Must be a field from
#' \code{other.fields}. (Default: \code{NULL})
#'
#' @param node_shape_by Deprecated. Use \code{node.shape.by} instead.
#'
#' @param node.size.by \code{Character scalar}. Specification of a column
#' metadata field or a feature to size tree nodes by. Must be a field from
#' \code{other.fields}. (Default: \code{NULL})
#'
#' @param node_size_by Deprecated. Use \code{node.size.by} instead.
#'
#' @param colour.highlights.by \code{Logical scalar}. Should the highlights be
#' colour differently? If \code{show.highlights = TRUE},
#' \code{colour_highlights} will be set to \code{TRUE} as default.
#' (Default: \code{FALSE})
#'
#' @param colour_highlights_by Deprecated. Use \code{colour.highlights.by}
#' instead.
#'
#' @param assay.type \code{Character scalar}. or \code{integer scalar}.
#' Specifies which assay to obtain expression values from, for use in point
#' aesthetics - see the \code{exprs_values} argument in
#' \code{\link[scater:retrieveCellInfo]{?retrieveCellInfo}}.
#' (Default: \code{"counts"})
#'
#' @param by_exprs_values Deprecated. Use \code{assay.type} instead.
#'
#' @param other.fields \code{Character vector}. Additional fields to include in
#' the node information without plotting them. (Default: \code{list()})
#'
#' @param other_fields Deprecated. Use \code{other.fields} instead.
#'
#' @param ... additional arguments for plotting. See
#' \code{\link{mia-plot-args}} for more details i.e. call
#' \code{help("mia-plot-args")}
#'
#' @details
#' If \code{show.label} or \code{show.highlight.label} have the same length
#' as the number of nodes, the vector will be used to relabel the nodes.
#'
#' @return a \code{\link{ggtree}} plot
#'
#' @seealso
#' \code{\link[mia:agglomerate-methods]{agglomerateByRanks}}
#'
#' @name plotTree
#'
#' @examples
#' library(scater)
#' library(mia)
#' # preparation of some data
#' data(GlobalPatterns)
#' GlobalPatterns <- agglomerateByRanks(GlobalPatterns)
#' altExp(GlobalPatterns,"Genus") <- addPerFeatureQC(
#'   altExp(GlobalPatterns,"Genus"))
#' rowData(altExp(GlobalPatterns,"Genus"))$log_mean <-
#'   log(rowData(altExp(GlobalPatterns,"Genus"))$mean)
#' rowData(altExp(GlobalPatterns,"Genus"))$detected <-
#'    rowData(altExp(GlobalPatterns,"Genus"))$detected / 100
#' top_genus <- getTop(altExp(GlobalPatterns,"Genus"),
#'                         method="mean",
#'                         top=100L,
#'                         assay.type="counts")
#' #
#' x <- altExp(GlobalPatterns,"Genus")
#' plotRowTree(x[rownames(x) %in% top_genus,],
#'             tip.colour.by = "log_mean",
#'             tip.size.by = "detected")
#'
#' # plot with tip labels
#' plotRowTree(x[rownames(x) %in% top_genus,],
#'             tip.colour.by = "log_mean",
#'             tip.size.by = "detected",
#'             show.label = TRUE)
#' # plot with selected labels
#' labels <- c("Genus:Providencia", "Genus:Morganella", "0.961.60")
#' plotRowTree(x[rownames(x) %in% top_genus,],
#'             tip.colour.by = "log_mean",
#'             tip.size.by = "detected",
#'             show.label = labels,
#'             layout="rectangular")
#'
#' # plot with labeled edges
#' plotRowTree(x[rownames(x) %in% top_genus,],
#'             edge.colour.by = "Phylum",
#'             tip.colour.by = "log_mean")
#' # if edges are sized, colours might disappear depending on plotting device
#' plotRowTree(x[rownames(x) %in% top_genus,],
#'             edge.colour.by = "Phylum",
#'             edge.size.by = "detected",
#'             tip.colour.by = "log_mean")
#'
#' # aggregating data over the taxonomic levels for plotting a taxonomic tree
#' # please note that the original tree of GlobalPatterns is dropped by
#' # unsplitByRanks
#' altExps(GlobalPatterns) <- splitByRanks(GlobalPatterns)
#' top_phyla <- getTop(altExp(GlobalPatterns,"Phylum"),
#'                         method="mean",
#'                         top=10L,
#'                         assay.type="counts")
#' altExps(GlobalPatterns) <- lapply(altExps(GlobalPatterns), addPerFeatureQC)
#' altExps(GlobalPatterns) <-
#'    lapply(altExps(GlobalPatterns),
#'           function(y){
#'               rowData(y)$log_mean <- log(rowData(y)$mean)
#'               rowData(y)$detected <- rowData(y)$detected / 100
#'               y
#'           })
#' x <- unsplitByRanks(GlobalPatterns)
#' x <- addHierarchyTree(x)
#'
#' highlights <- c("Phylum:Firmicutes","Phylum:Bacteroidetes",
#'                 "Family:Pseudomonadaceae","Order:Bifidobacteriales")
#' plotRowTree(x[rowData(x)$Phylum %in% top_phyla,],
#'             tip.colour.by = "log_mean",
#'             node.colour.by = "log_mean",
#'             show.highlights = highlights,
#'             show.highlight.label = highlights,
#'             colour.highlights.by = "Phylum")
#'
#' plotRowTree(x[rowData(x)$Phylum %in% top_phyla,],
#'             edge.colour.by = "Phylum",
#'             edge.size.by = "detected",
#'             tip.colour.by = "log_mean",
#'             node.colour.by = "log_mean")
NULL

#' @rdname plotTree
#' @export
setMethod("plotColTree", signature = c(x = "TreeSummarizedExperiment"),
    function(x,
        tree.name = tree_name,
        tree_name = "phylo",
        relabel.tree = relabel_tree,
        relabel_tree = FALSE,
        order.tree = order_tree,
        order_tree = FALSE,
        levels.rm = remove_levels,
        remove_levels = FALSE,
        show.label = show_label,
        show_label = FALSE,
        show.highlights = show_highlights,
        show_highlights = FALSE,
        show.highlight.label = show_highlight_label,
        show_highlight_label = FALSE,
        abbr.label = abbr_label,
        abbr_label = FALSE,
        add.legend = add_legend,
        add_legend = TRUE,
        layout = "circular",
        edge.colour.by = edge.colour.by,
        edge_colour_by = NULL,
        edge.size.by = edge_size_by,
        edge_size_by = NULL,
        tip.colour.by = tip_colour_by,
        tip_colour_by = NULL,
        tip.shape.by = tip_shape_by,
        tip_shape_by = NULL,
        tip.size.by = tip_size_by,
        tip_size_by = NULL,
        node.colour.by = node_colour_by,
        node_colour_by = NULL,
        node.shape.by = node_shape_by,
        node_shape_by = NULL,
        node.size.by = node_size_by,
        node_size_by = NULL,
        colour.highlights.by = colour_highlights_by,
        colour_highlights_by = NULL,
        assay.type = by_exprs_values,
        by_exprs_values = "counts",
        other.fields = other_fields,
        other_fields = list(),
        ...){
        .plot_row_column_tree(x,
            tree_name = tree.name,
            relabel_tree = relabel.tree,
            order_tree = order.tree,
            remove_levels = levels.rm,
            show_label = show.label,
            show_highlights = show.highlights,
            show_highlight_label = show.highlight.label,
            abbr_label = abbr.label,
            add_legend = add.legend,
            layout = layout,
            edge_colour_by = edge.colour.by,
            edge_size_by = edge.size.by,
            tip_colour_by = tip.colour.by,
            tip_shape_by = tip.shape.by,
            tip_size_by = tip.size.by,
            node_colour_by = node.colour.by,
            node_shape_by = node.shape.by,
            node_size_by = node.size.by,
            colour_highlights_by = colour.highlights.by,
            by_exprs_values = assay.type,
            other_fields = other.fields,
            type = "column",
            ...)
    }
)
#' @rdname plotTree
#' @export
setMethod("plotRowTree", signature = c(x = "TreeSummarizedExperiment"),
    function(x,
        tree.name = tree_name,
        tree_name = "phylo",
        relabel.tree = relabel_tree,
        relabel_tree = FALSE,
        order.tree = order_tree,
        order_tree = FALSE,
        levels.rm = remove_levels,
        remove_levels = FALSE,
        show.label = show_label,
        show_label = FALSE,
        show.highlights = show_highlights,
        show_highlights = FALSE,
        show.highlight.label = show_highlight_label,
        show_highlight_label = FALSE,
        abbr.label = abbr_label,
        abbr_label = FALSE,
        add.legend = add_legend,
        add_legend = TRUE,
        layout = "circular",
        edge.colour.by = edge_colour_by,
        edge_colour_by = NULL,
        edge.size.by = edge_size_by,
        edge_size_by = NULL,
        tip.colour.by = tip_colour_by,
        tip_colour_by = NULL,
        tip.shape.by = tip_shape_by,
        tip_shape_by = NULL,
        tip.size.by = tip_size_by,
        tip_size_by = NULL,
        node.colour.by = node_colour_by,
        node_colour_by = NULL,
        node.shape.by = node_shape_by,
        node_shape_by = NULL,
        node.size.by = node_size_by,
        node_size_by = NULL,
        colour.highlights.by = colour_highlights_by,
        colour_highlights_by = NULL,
        assay.type = by_exprs_values,
        by_exprs_values = "counts",
        other.fields = other_fields,
        other_fields = list(),
        ...){
        #
        .plot_row_column_tree(x,
            tree_name = tree.name,
            relabel_tree = relabel.tree,
            order_tree = order.tree,
            remove_levels = levels.rm,
            show_label = show.label,
            show_highlights = show.highlights,
            show_highlight_label = show.highlight.label,
            abbr_label = abbr.label,
            add_legend = add.legend,
            layout = layout,
            edge_colour_by = edge.colour.by,
            edge_size_by = edge.size.by,
            tip_colour_by = tip.colour.by,
            tip_shape_by = tip.shape.by,
            tip_size_by = tip.size.by,
            node_colour_by = node.colour.by,
            node_shape_by = node.shape.by,
            node_size_by = node.size.by,
            colour_highlights_by = colour.highlights.by,
            by_exprs_values = assay.type,
            other_fields = other.fields,
            type = "row",
            ...)
    }
)

.check_tree_plot_switches <- function(layout,
        relabel_tree,
        remove_levels,
        order_tree,
        show_label,
        show_highlights,
        show_highlight_label,
        abbr_label,
        add_legend){
    if(!.is_a_string(layout)){
        stop("'layout' must be a single character value.", call. = FALSE)
    }
    if(!.is_a_bool(relabel_tree)){
        stop("'relabel.tree' must be either TRUE or FALSE.", call. = FALSE)
    }
    if(!.is_a_bool(remove_levels)){
        stop("'level.rm' must be either TRUE or FALSE.", call. = FALSE)
    }
    if(!.is_a_bool(order_tree)){
        stop("'order.tree' must be either TRUE or FALSE.", call. = FALSE)
    }
    if(!.is_a_bool(show_label)){
        if( (!is.logical(show_label) && !is.character(show_label) &&
                !is.numeric(show_label)) || is.null(show_label)){
            stop("'show.label' must be either TRUE or FALSE or logical, ",
                "integer or character ",
                "vector. Character alues should match the label of the tree.",
                call. = FALSE)
        }
    }
    if(!.is_a_bool(show_highlights)){
        if( (!is.logical(show_highlights) && !is.character(show_highlights) &&
                !is.numeric(show_highlights)) || is.null(show_highlights)){
            stop("'show.label' must be either TRUE or FALSE or logical, ",
                "integer or character ",
                "vector. Character alues should match the label of the tree.",
                call. = FALSE)
        }
    }
    if(!.is_a_bool(show_highlight_label)){
        if( (!is.logical(show_highlight_label) &&
                !is.character(show_highlight_label) &&
                !is.numeric(show_highlight_label)) ||
                is.null(show_highlight_label)){
            stop("'show.highlight.label' must be either TRUE or FALSE or ",
                "logical, integer or character ",
                "vector. Character alues should match the label of the tree.",
                call. = FALSE)
        }
    }
    if(!.is_a_bool(abbr_label)){
        if( (!is.logical(abbr_label) && !is.character(abbr_label) &&
                !is.numeric(abbr_label)) || is.null(abbr_label)){
            stop("'abbr.label' must be either TRUE or FALSE or logical, ",
                "integer or character ",
                "vector. Character alues should match the label of the tree.",
                call. = FALSE)
        }
    }
    if(!.is_a_bool(add_legend)){
        stop("'add.legend' must be either TRUE or FALSE.", call. = FALSE)
    }
}

.plot_row_column_tree <- function(object,
        tree_name = "phylo",
        relabel_tree = FALSE,
        order_tree = FALSE,
        remove_levels = FALSE,
        show_label = FALSE,
        show_highlights = FALSE,
        show_highlight_label = FALSE,
        abbr_label = FALSE,
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
        colour_highlights_by = NULL,
        by_exprs_values = "counts",
        other_fields = list(),
        type = c("row","column"),
        ...){
    type <- match.arg(type)
    # input check
    # Check tree_name
    if( !.is_a_string(tree_name) ){
        stop("'tree.name' must be a single character value specifying a ",
            "colTree.", call. = FALSE)
    }
    FUN <- switch(
        type,
        row = "rowTree",
        column = "colTree")
    if(is.null(do.call(FUN,list(x = object, whichTree = tree_name)))){
        stop(FUN,"(object, tree.name) is empty.", call. = FALSE)
    }
    .check_tree_plot_switches(layout = layout,
        relabel_tree = relabel_tree,
        remove_levels = remove_levels,
        order_tree = order_tree,
        show_label = show_label,
        show_highlights = show_highlights,
        show_highlight_label = show_highlight_label,
        abbr_label = abbr_label,
        add_legend = add_legend)
    #
    tree_out <- .get_object_and_trimmed_tree(
        object,
        tree_name = tree_name,
        type = type,
        relabel = relabel_tree,
        order = order_tree)
    object <- tree_out$object
    tree <- tree_out$tree
    tree_data <- .get_tree_data(tree)
    label_out <- .add_tree_node_labels(tree_data, show_label, remove_levels)
    tree_data <- label_out$df
    show_label <- label_out$show_label
    label_out <- .add_tree_highlights(tree_data, show_highlights)
    tree_data <- label_out$df
    show_highlights <- label_out$show_highlights
    label_out <- .add_tree_highlight_labels(
        tree_data, show_highlight_label, remove_levels)
    tree_data <- label_out$df
    show_highlight_label <- label_out$show_highlight_label
    #
    vis_out <- .incorporate_tree_vis(
        tree_data,
        se = object,
        edge_colour_by = edge_colour_by,
        edge_size_by = edge_size_by,
        tip_colour_by = tip_colour_by,
        tip_shape_by = tip_shape_by,
        tip_size_by = tip_size_by,
        node_colour_by = node_colour_by,
        node_shape_by = node_shape_by,
        node_size_by = node_size_by,
        colour_highlights_by = colour_highlights_by,
        by_exprs_values = by_exprs_values,
        other_fields = other_fields,
        type = type)
    tree_data <- vis_out$df
    edge_colour_by <- vis_out$edge_colour_by
    edge_size_by <- vis_out$edge_size_by
    colour_by <- vis_out$colour_by
    shape_by <- vis_out$shape_by
    size_by <- vis_out$size_by
    colour_highlights_by <- vis_out$colour_highlights_by
    show_tips <- any(!vapply(
        c(tip_colour_by, tip_shape_by, tip_size_by),
        is.null, logical(1)))
    show_nodes <- any(!vapply(
        c(node_colour_by, node_shape_by, node_size_by),
        is.null, logical(1)))
    #
    object <- .create_treedata_for_plotting(
        tree_data, tree, edge_colour_by, edge_size_by, shape_by, size_by)
    .tree_plotter(
        object,
        layout = layout,
        add_legend = add_legend,
        show_label = show_label,
        show_highlights = show_highlights,
        show_highlight_label = show_highlight_label,
        abbr_label = abbr_label,
        show_tips = show_tips,
        show_nodes = show_nodes,
        edge_colour_by = edge_colour_by,
        edge_size_by = edge_size_by,
        colour_by = colour_by,
        shape_by = shape_by,
        size_by = size_by,
        colour_highlights_by = colour_highlights_by,
        order_tree = order_tree,
        ...)
}

#' @importFrom ape keep.tip as.phylo drop.tip
#' @importFrom tidytree as_tibble
.get_object_and_trimmed_tree <- function(
        object,
        tree_name = "phylo",
        type = c("row","column"),
        relabel = FALSE,
        order = FALSE){
    # Check type
    type <- match.arg(type)
    # Get correct functions based on the margin/direction
    tree_FUN <- switch(type, row = rowTree, column = colTree, stop("."))
    links_FUN <- switch(type, row = rowLinks, column = colLinks, stop("."))
    dimnames_FUN <- switch(type, row = rownames, column = colnames, stop("."))
    add_names_FUN <- switch(
        type, row = `rownames<-`, column = `colnames<-`, stop("."))
    # Check that the tree is compatible with the data, i.e., rows are matched
    # with the tree.
    links_FUN <- switch(type, row = rowLinks, column = colLinks, stop("."))
    links <- links_FUN(object)
    ind <- links[["whichTree"]] == tree_name
    if( all(!ind) ){
        stop("Tree does not have any ", type, "s to plot.", call. = FALSE)
    }
    # Get only those rows/columns that are found from the tree
    if( type == "row" ){
        object <- object[ind, ]
    } else{
        object <- object[, ind]
    }
    # Get tree and links
    tree <- tree_FUN(object, tree_name)
    links <- links_FUN(object)

    # Remove those tips that are not leaves
    tips <- sort(setdiff(tree$edge[, 2], tree$edge[, 1]))
    drop_tip <- tips[!(tips %in% unique(links$nodeNum[links$isLeaf]))]
    oldTree <- tree
    newTree <- drop.tip(oldTree, tip = drop_tip, collapse.singles = FALSE)
    # Add alias labels to tree
    track <- trackNode(oldTree)
    track <- drop.tip(track, tip = drop_tip, collapse.singles = FALSE)
    # Link tree with alias labels
    oldAlias <- links$nodeLab_alias
    newNode <- convertNode(tree = track, node = oldAlias)
    newAlias <- convertNode(tree = newTree, node = newNode)
    # Change the tree with trimmed tree and add aliases as node labels
    if( type == "row" ){
        object <- changeTree(
            x = object, rowTree = newTree, rowNodeLab = newAlias)
    } else {
        object <- changeTree(
            x = object, colTree = newTree, colNodeLab = newAlias)
    }

    # Get tree, links and row/colnames
    tree <- tree_FUN(object)
    links <- links_FUN(object)
    dimnames <- dimnames_FUN(object)
    # Get tree as table and get which node represent which row/col
    tree_data <- as_tibble(tree)
    m <- match(links$nodeNum,tree_data$node)
    node_labels <- tree_data$label[m]
    # If user wants to rename rows/cols or if some nodes cannot be found from
    # rows/cols
    if( relabel || !all(node_labels %in% dimnames) ){
        # Rename rows/cols
        new_node_labels <- getTaxonomyLabels(
            object, with_rank = TRUE, resolve_loops = TRUE)
        object <- add_names_FUN(object, new_node_labels)
    }
    # Check if there are rows/cols that are ununique. If there are, make them
    # unique.
    if( anyDuplicated(rownames(object)) ){
        warning(
            "Data includes ununique ", type, "s. Making them unique.",
            call. = FALSE)
        object <- add_names_FUN(object, make.unique(dimnames_FUN(object)))

    }
    # Rename labels of tree with row/colnames
    tree_data$label[m] <- dimnames_FUN(object)
    # Check if there are nodes that are not unique
    if( anyDuplicated(tree_data$label[-m]) ){
        warning(
            "Tree includes ununique nodes. Making them unique.", call. = FALSE)
        tree_data$label[-m] <- make.unique( tree_data$label[-m] )
    }

    # Convert tree data back to tree-format
    tree <- as.phylo(tree_data)
    # If specified, order the tree based on alphabetical order
    if(order){
        tree <- .order_tree(tree)
    }
    res <- list(object = object, tree = tree)
    return(res)
}

#' @importFrom tidytree child
.get_tree_labels_for_ordering <- function(tree_data, node){
    children <- child(tree_data, node)
    if(nrow(children) == 0L){
        return("")
    }
    labels <- children$label
    add_labels <- lapply(
        children$node,
        .get_tree_labels_for_ordering,
        tree_data = tree_data)
    unlist(
        mapply(paste,labels,add_labels,sep="__:__",SIMPLIFY = FALSE),
        use.names = FALSE)
}

#' @importFrom tidytree rootnode as_tibble
#' @importFrom ape rotateConstr
.order_tree <- function(tree){
    tree_data <- tidytree::as_tibble(tree)
    root_node <- rootnode(tree_data)
    labels <- paste0("__:__",
        .get_tree_labels_for_ordering(tree_data, root_node$node))
    tip_labels <- regmatches(labels,regexec(".*__:__(.+?)__:__$",labels))
    tip_labels <- vapply(tip_labels,"[",character(1),2L)
    o <- order(labels, decreasing = TRUE)
    contraint <- tip_labels[o]
    tree <- ape::rotateConstr(tree, rev(contraint))
    tree
}

################################################################################

.remove_taxonomic_level_from_labels <- function(labels){
    for(rank in TAXONOMY_RANKS){
        labels <- gsub(paste0(rank,":"),"",labels,ignore.case = TRUE)
    }
    labels
}

#' @importFrom tidygraph activate
#' @importFrom dplyr mutate
.add_tree_node_labels <- function(
        tree_data, show_label, remove_levels = FALSE){
    if("label" %in% colnames(tree_data)){
        tree_data <- tree_data %>% mutate(node_label = .data$label)
    }

    if(!is.logical(show_label) || length(show_label) > 1L) {
        if(is.character(show_label) && length(show_label) == nrow(tree_data)) {
            tree_data <- tree_data %>%
                mutate(node_label = show_label)
            show_label <- TRUE
        } else if(!("node_label" %in% colnames(tree_data))){
            warning("If 'show.label' is a character with length != ",
                    "number of nodes in the graph or a logical/integer ",
                    "vector, a 'label' ",
                    "column must exist in the tree data.",
                    call. = FALSE)
            show_label <- FALSE
        } else {
            if(is.numeric(show_label)){
                if(any(show_label != as.integer(show_label)) ||
                        min(show_label) < 1 ||
                        max(show_label) > nrow(tree_data)){
                    stop("If 'show.label' is numeric, values have to be whole ",
                        "numbers and must be between 1 and the number of ",
                        "nodes in the graph", call. = FALSE)
                }
                label <- rep(FALSE, nrow(tree_data))
                label[tree_data$node %in% show_label] <- TRUE
                show_label <- label
            } else if(is.character(show_label)) {
                show_label <- tree_data$node_label %in% show_label
            }
            if(is.logical(show_label) && length(show_label) != nrow(tree_data)){
                stop("If 'show.label' is logical, it must have the length as ",
                    "nodes are in the graph.", call. = FALSE)
            }
            tree_data <- tree_data %>%
                mutate(node_label = ifelse(
                    show_label, .data$node_label, NA_character_))
            show_label <- TRUE
        }
        if(all(is.na(tree_data %>% pull("node_label")))){
            show_label <- FALSE
            warning("No labels to plot.", call. = FALSE)
        }
    } else if(is.logical(show_label) && length(show_label) == 1L &&
            !show_label) {
        tree_data <- tree_data %>%
            mutate(node_label = FALSE)
    }
    if(remove_levels){
        tree_data$node_label <- .remove_taxonomic_level_from_labels(
            tree_data$node_label)
    }
    res <- list(df = tree_data, show_label = show_label)
    return(res)
}

#' @importFrom tidygraph activate
#' @importFrom dplyr mutate
.add_tree_highlights <- function(tree_data, show_highlights){
    tree_data$highlight <- FALSE

    if(!is.logical(show_highlights) || length(show_highlights) > 1L) {
        if(is.numeric(show_highlights)){
            if(any(show_highlights != as.integer(show_highlights)) ||
                    min(show_highlights) < 1 ||
                    max(show_highlights) > nrow(tree_data)){
                stop("If 'show.highlights' is numeric, values have to be ",
                    "whole numbers and must be between 1 and the number of ",
                    "nodes in the graph", call. = FALSE)
            }
            label <- rep(FALSE, nrow(tree_data))
            label[tree_data$node %in% show_highlights] <- TRUE
            show_highlights <- label
        } else if(is.character(show_highlights)) {
            show_highlights <- tree_data$label %in% show_highlights
        }
        if(is.logical(show_highlights) &&
                length(show_highlights) != nrow(tree_data)){
            stop("If 'show.highlights' is logical, it must have the length as ",
                "nodes are in the graph.", call. = FALSE)
        }
        tree_data <- tree_data %>%
            mutate(highlight = show_highlights)
        show_highlights <- TRUE
        if(!any(tree_data %>% pull("highlight"))){
            show_highlights <- FALSE
            warning("No highlights to plot.", call. = FALSE)
        }
    } else if(is.logical(show_highlights) && length(show_highlights) == 1L &&
            show_highlights){
        tree_data$highlight <- TRUE
    }
    res <- list(df = tree_data, show_highlights = show_highlights)
    return(res)
}

#' @importFrom tidygraph activate
#' @importFrom dplyr mutate
.add_tree_highlight_labels <- function(
        tree_data, show_highlight_label, remove_levels = FALSE){
    if(!any(tree_data$highlight)){
        show_highlight_label <- FALSE
        tree_data$highlight_label <- FALSE
        res <- list(df = tree_data, show_highlight_label = show_highlight_label)
        return(res)
    }

    if("label" %in% colnames(tree_data)){
        tree_data <- tree_data %>%
            mutate(highlight_label = .data$label)
    }
    if(!is.logical(show_highlight_label) || length(show_highlight_label) > 1L) {
        if(is.character(show_highlight_label) &&
                length(show_highlight_label) == nrow(tree_data)) {
            tree_data <- tree_data %>%
                mutate(highlight_label = show_highlight_label)
            show_highlight_label <- TRUE
        } else if(!("highlight_label" %in% colnames(tree_data))){
            warning("If 'show.highlight.label' is a character with length != ",
                    "number of nodes in the graph or a logical/integer ",
                    "vector, a 'label' column must exist in the tree data.",
                    call. = FALSE)
            show_highlight_label <- FALSE
        } else {
            if(is.numeric(show_highlight_label)){
                if(any(show_highlight_label !=
                            as.integer(show_highlight_label)) ||
                        min(show_highlight_label) < 1 ||
                        max(show_highlight_label) > nrow(tree_data)){
                    stop("If 'show.highlight.label' is numeric, values have ",
                        "to be whole numbers and must be between 1 and the ",
                        "number of nodes in the graph", call. = FALSE)
                }
                label <- rep(FALSE, nrow(tree_data))
                label[tree_data$node %in% show_highlight_label] <- TRUE
                show_highlight_label <- label
            } else if(is.character(show_highlight_label)) {
                show_highlight_label <-
                    tree_data$highlight_label %in% show_highlight_label
            }
            if(is.logical(show_highlight_label) &&
                    length(show_highlight_label) != nrow(tree_data)){
                stop("If 'show.highlight.label' is logical, it must have the ",
                    "length as nodes are in the graph.", call. = FALSE)
            }
            tree_data <- tree_data %>%
                mutate(highlight_label = ifelse(
                    show_highlight_label & tree_data$highlight,
                    .data$highlight_label, NA_character_))
            show_highlight_label <- TRUE
        }
        if(!any(tree_data %>% pull("highlight")) ||
                all(is.na(tree_data %>% pull("highlight_label")))){
            show_highlight_label <- FALSE
            warning("No highlights to label.", call. = FALSE)
        }
    } else if(is.logical(show_highlight_label) &&
            length(show_highlight_label) == 1L && !show_highlight_label){
        tree_data <- tree_data %>%
            mutate(highlight_label = NA_character_)
    }
    if(remove_levels){
        tree_data$highlight_label <- .remove_taxonomic_level_from_labels(
            tree_data$highlight_label)
    }
    res <- list(df = tree_data, show_highlight_label = show_highlight_label)
    return(res)
}

################################################################################

#' @importFrom tibble tibble
.get_feature_info <- function(by, se, FUN, exprs_values, var_name){
    feature_info <- try(
        FUN(se, by = by, exprs_values = exprs_values),
        silent = TRUE)
    if(is(feature_info,"try-error")){
        stop(feature_info, "for '",var_name,"'", call. = FALSE)
    }
    feature_info <- tibble(!!sym(feature_info$name) := feature_info$value)
    feature_info
}

TIP_VARIABLES <- c("tip_colour_by", "tip_shape_by", "tip_size_by")
NODE_VARIABLES <- c("node_colour_by", "node_shape_by", "node_size_by")

.get_new_var_name_value <- function(var_name_value, add){
    if(!is.null(var_name_value) && add != var_name_value){
        new_var_name_value <- paste0(
            var_name_value, ifelse(is.null(var_name_value),"", " & "), add)
    } else {
        new_var_name_value <- add
    }
    new_var_name_value
}

#' @importFrom scater retrieveFeatureInfo retrieveCellInfo
#' @importFrom dplyr bind_cols mutate relocate
#' @importFrom tibble rownames_to_column
.incorporate_tree_vis <- function(
        tree_data,
        se,
        edge_colour_by,
        edge_size_by,
        tip_colour_by,
        tip_shape_by,
        tip_size_by,
        node_colour_by,
        node_shape_by,
        node_size_by,
        colour_highlights_by,
        by_exprs_values = "counts",
        other_fields = other_fields,
        type = c("row","column")){
    type <- match.arg(type)
    type_FUN <- switch(
        type,
        row = scater::retrieveFeatureInfo,
        column = scater::retrieveCellInfo)
    variables <- c(
        edge_colour_by = edge_colour_by,
        edge_size_by = edge_size_by,
        tip_colour_by = tip_colour_by,
        tip_shape_by = tip_shape_by,
        tip_size_by = tip_size_by,
        node_colour_by = node_colour_by,
        node_shape_by = node_shape_by,
        node_size_by = node_size_by,
        colour_highlights_by = colour_highlights_by)
    edge_colour_by <- NULL
    edge_size_by <- NULL
    colour_by <- NULL
    shape_by <- NULL
    size_by <- NULL
    colour_highlights_by <- NULL
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
                colnames(tree_data) <- c(
                    DEFAULT_TREE_DATA_COLS,names(variables)[f])
                # mirror back variable name
                for(i in variables[f]){
                    var_name <- gsub("tip_|node_","",names(variables)[f][i])
                    assign(
                        var_name,
                        .get_new_var_name_value(get(var_name), variables[f][i]))
                }
                variables <- variables[!f]
            }
        }
        if(length(variables) > 0L){
            feature_info <- vector(mode = "list", length = length(variables))
            for(i in seq_along(variables)){
                # get data
                var_name <- names(variables)[i]
                feature_info[[i]] <- .get_feature_info(
                    variables[i], se = se, FUN = type_FUN,
                    exprs_values = by_exprs_values, var_name = var_name)
                # mirror back variable name, if a partial match was used
                var_name <- gsub("tip_|node_","",var_name)
                assign(var_name, .get_new_var_name_value(
                    get(var_name), colnames(feature_info[[i]])))
                # rename columns by their usage
                colnames(feature_info[[i]]) <- names(variables[i])
            }
            feature_info <- bind_cols(feature_info)
            feature_info <- feature_info %>%
                mutate(label = rownames(se)) %>%
                relocate("label")
            tree_data <- .merge_tree_vis_data(tree_data, feature_info, se)
        }
        tree_data <- .merge_tip_node_tree_data(tree_data)
    }
    if(length(other_fields) != 0L){
        for (o in other_fields) {
            other <- type_FUN(se, o, exprs_values = by_exprs_values)
            other <- other %>%
                mutate(label = rownames(se)) %>%
                relocate("label")
            tree_data <- .merge_tree_vis_data(tree_data, other, se)
        }
    }
    res <- list(
        df = tree_data,
        edge_colour_by = edge_colour_by,
        edge_size_by = edge_size_by,
        colour_by = colour_by,
        shape_by = shape_by,
        size_by = size_by,
        colour_highlights_by = colour_highlights_by)
    return(res)
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
        colour_by <- c(
            tree_data$node_colour_by[!is_leaf_o],
            tree_data$tip_colour_by[is_leaf_o])
    } else if("tip_colour_by" %in% cn) {
        colour_by <- tree_data$tip_colour_by
    } else if("node_colour_by" %in% cn) {
        colour_by <- tree_data$node_colour_by
    }
    if(all(c("tip_shape_by","node_shape_by") %in% cn)){
        shape_by <- c(
            tree_data$node_shape_by[!is_leaf_o],
            tree_data$tip_shape_by[is_leaf_o])
    } else if("tip_shape_by" %in% cn) {
        shape_by <- tree_data$tip_shape_by
    } else if("node_shape_by" %in% cn) {
        shape_by <- tree_data$node_shape_by
    }
    if(all(c("tip_size_by","node_size_by") %in% cn)){
        size_by <- c(
            tree_data$node_size_by[!is_leaf_o],
            tree_data$tip_size_by[is_leaf_o])
    } else if("tip_size_by" %in% cn) {
        size_by <- tree_data$tip_size_by
    } else if("node_size_by" %in% cn) {
        size_by <- tree_data$node_size_by
    }
    #
    tree_data <- tree_data[,cn[!grepl("tip_|node_",cn) | cn == "node_label"]]
    tree_data$colour_by <- colour_by
    tree_data$shape_by <- shape_by
    tree_data$size_by <- size_by
    # return tree_data with original ordering
    tree_data[match(bak_o,tree_data$node),]
}

.merge_tree_vis_data <- function(tree_data, feature_info, tse){
    if(anyDuplicated(tree_data$label) || anyDuplicated(feature_info$label)){
        stop("Tree is not compatible with the data.", call. = FALSE)
    }
    tree_data <- tree_data %>%
        dplyr::left_join(feature_info, by = "label")
    return(tree_data)
}

# due to a bug in ggtree/tidytree the treedata object needs to be constructed
# in a separate step
#
# also there is some data wrangling needed
#' @importFrom tidytree as.treedata
.create_treedata_for_plotting <- function(
        tree_data, tree, edge_colour_by, edge_size_by, shape_by, size_by){
    # cleanup
    if (!is.null(edge_colour_by) &&
        anyNA(tree_data$edge_colour_by) &&
        !is.numeric(tree_data$edge_colour_by)) {
        tree_data <- groupOTU(
            tree_data,
            split(tree_data$node, tree_data$edge_colour_by),
            group_name = "group")
        f_zero <- tree_data$group != 0
        f_zero <- f_zero[!is.na(f_zero)]
        tree_data$edge_colour_by[f_zero] <- as.character(
            tree_data$group[f_zero])
    }
    tree_data <- .na_replace_from_plot_data(
        tree_data, edge_size_by, shape_by, size_by)
    object <- tidytree::as.treedata(tree_data)
    # tree needs to be restored since the original leave/tip/node orientation
    # is not compatible with ladderiez = FALSE
    object@phylo <- tree
    #
    object
}

#' @importFrom ggplot2 scale_size_identity
#' @importFrom ggtree ggtree geom_tree geom_tippoint geom_nodepoint groupOTU
#'   theme_tree
.tree_plotter <- function(
        object,
        layout,
        add_legend,
        show_label,
        show_highlights,
        show_highlight_label,
        abbr_label,
        show_tips,
        show_nodes,
        edge_colour_by,
        edge_size_by,
        colour_by,
        shape_by,
        size_by,
        colour_highlights_by,
        order_tree,
        line_alpha = line.alpha,
        line.alpha = 1,
        line_width = line.width,
        line.width = NULL,
        line_width_range = line.width.range,
        line.width.range = c(0.5,3),
        point_alpha = point.alpha,
        point.alpha = 1,
        point_size = point.size,
        point.size = 2,
        point_size_range = point.size.range,
        point.size.range = c(1,4),
        label_font_size = label.font.size,
        label.font.size = 3,
        highlight_font_size = highlight.font.size,
        highlight.font.size = 3,
        branch.length = "branch.length",
        ...){
    # We capture branch.length to disable other options than plotting branches
    # as they are or with equal branch lengths. That is because user cannot
    # control the parameter as it only specifies column from the object table.
    # However, the table cannot include any other length scales than the
    # original.
    cols <- c("branch.length", "none")
    if( !(.is_a_string(branch.length) && branch.length %in% cols) ){
        stop("'branch.length' must be one of the following options: '",
            paste0(cols, collapse = "', '"), "'", call. = FALSE)
    }
    # start plotting
    plot_out <- ggtree(
        object, ladderize = !order_tree, layout = layout,
        branch.length = branch.length, ...)
    # add highlights
    plot_out <- .plot_tree_plot_highlights(
        plot_out, layout, show_highlights, show_highlight_label, abbr_label,
        colour_highlights_by, highlight_font_size = highlight_font_size)
    # add tree and adjust edges
    plot_out <- .plot_tree_edges(
        plot_out, edge_colour_by, edge_size_by, line_alpha, line_width,
        line_width_range, layout)
    # add tip and node points
    plot_out <- .plot_tree_node_points(
        plot_out, show_tips, show_nodes, colour_by, shape_by, size_by,
        point_alpha, point_size, point_size_range)
    # add tip and node labels
    plot_out <- .plot_tree_node_labels(
        plot_out, show_label, abbr_label, label_font_size)
    # add additional guides
    plot_out <- .add_extra_guide(plot_out, shape_by, size_by)
    # add abbreviation guide
    plot_out <- .add_abbr_guide(plot_out)
    # add theme
    plot_out <- .theme_plotTree(plot_out)
    # optionally hide legends
    if (!add_legend) {
        plot_out <- plot_out +
            theme(legend.position = "none")
    }
    plot_out
}

.add_label_abbreviations <- function(
        plot_out, label_col, subset = NULL){
    non_abbr_text_col <- paste0("abbr_",label_col)
    if(is.null(subset)){
        subset <- seq_len(nrow(plot_out$data))
    }
    subset <- seq_len(nrow(plot_out$data)) %in% subset
    # initialize column if not present
    if(!(non_abbr_text_col %in% colnames(plot_out$data))){
        plot_out$data[,non_abbr_text_col] <- NA_character_
    }
    #
    text <- plot_out$data[subset,label_col,drop=TRUE]
    if(length(text) > 0L){
        # save text
        bak_text <- text
        # abbreviate with unique element
        u_text <- unique(text)
        abbr <- abbreviate(
            gsub("[_]|[-][ ]","",u_text), minlength = 1, dot = TRUE)
        # reflate to original positions
        abbr <- abbr[match(text, u_text)]
        # exchange label
        plot_out$data[subset,label_col] <- abbr
        # exchange original text
        plot_out$data[subset,non_abbr_text_col] <- bak_text
    }
    plot_out
}

.get_hightlight_args <- function(nodes, colour_highlights_by){
    aes_args <- list()
    aes_args$subset <- paste0("node %in% c(",paste(nodes, collapse = ","),")")
    aes_args$extendto <- ~highlight_extendto
    if(!is.null(colour_highlights_by)) {
        aes_args$fill <- ~colour_highlights_by
    }
    new_aes <- do.call(aes_, aes_args)
    geom_args <- list(mapping = new_aes)
    geom_args$colour <- "grey20"
    if (is.null(colour_highlights_by)) {
        geom_args$fill <- "grey70"
    }
    return(list(args = geom_args))
}

.get_cladelab_args <- function(
        nodes, layout, highlight_font_size){
    aes_args <- list()
    aes_args$subset <- paste0("node %in% c(",paste(nodes, collapse = ","), ")")
    aes_args$node <- ~node
    aes_args$label <- ~highlight_label
    aes_args$offset.text <- ~highlight_offset
    new_aes <- do.call(aes_, aes_args)
    geom_args <- list(mapping = new_aes)
    if(layout %in% c("fan","circular","radial")){
        geom_args$hjust <- 0.5
        geom_args$angle <- "auto"
        geom_args$horizontal <- FALSE
    } else if(layout %in% c("inward_circular")){
        geom_args$hjust <- 0.5
        geom_args$angle <- "auto"
        geom_args$horizontal <- FALSE
    }
    geom_args$barsize <- NA
    geom_args$fontsize <- highlight_font_size
    return(list(args = geom_args))
}

#' @importFrom dplyr mutate
.calc_highlight_extendto <- function(highlight_data, layout) {
    if(layout %in% c("fan","circular","radial")){
        ans <- highlight_data %>%
            mutate(
                highlight_extendto = (max(.data$x) - .data$x) / 1.5,
                highlight_extendto = .data$highlight_extendto +
                    max(.data$x) + 0.07)
    } else if(layout %in% c("rectangular","slanted","ellipse","roundrect")){
        ans <- highlight_data %>%
            mutate(
                highlight_extendto = (max(.data$x) - .data$x) / 1.5,
                highlight_extendto = .data$highlight_extendto +
                    max(.data$x) + 0.01)
    } else if(layout %in% c("dendrogram")){
        warning("highlights with layout `dendrogram` are buggy.")
        ans <- highlight_data %>%
            mutate(
                highlight_extendto = .data$x / 1.5,
                highlight_extendto = (.data$highlight_extendto - 0.01) * -1)
    } else if(layout %in% c("inward_circular")){
        warning("highlights with layout `inward_circular` are buggy.")
        ans <- highlight_data %>%
            mutate(
                highlight_extendto = (max(.data$x) - .data$x) / 1.5,
                highlight_extendto = .data$highlight_extendto + max(.data$x) +
                    0.07,
                highlight_extendto = .data$highlight_extendto * -1)
    } else {
        ans <- highlight_data %>%
            mutate(highlight_extendto = .data$x)
    }
    ans
}

#' @importFrom dplyr mutate
.calc_highlight_label_text_offset <- function(label_data, layout){
    if(layout %in% c("fan","circular","radial")){
        ans <- label_data %>%
            mutate(highlight_offset = .data$highlight_extendto - max(.data$x) +
                    0.015 - 0.07)
    } else if(layout %in% c("rectangular","slanted","ellipse","roundrect")){
        ans <- label_data %>%
            mutate(highlight_offset = .data$highlight_extendto - max(.data$x) -
                    0.01)
    } else if(layout %in% c("dendrogram")){
        ans <- label_data %>%
            mutate(highlight_offset = .data$highlight_extendto - 0.1)
    } else if(layout %in% c("inward_circular")){
        ans <- label_data %>%
            mutate(highlight_offset = (.data$highlight_extendto *-1) -
                    max(.data$x) - 0.022)
    } else {
        ans <- label_data %>%
            mutate(highlight_offset = .data$highlight_extendto)
    }
    ans
}

#' @importFrom dplyr filter pull
#' @importFrom ggtree geom_highlight geom_cladelab
#' @importFrom ggnewscale new_scale_fill new_scale_colour
#' @importFrom tidytree rootnode
.plot_tree_plot_highlights <- function(
        plot_out, layout, show_highlights, show_highlight_label, abbr_label,
        colour_highlights_by, highlight_font_size){
    plot_out$data <- .calc_highlight_extendto(plot_out$data, layout)
    plot_out$data <- .calc_highlight_label_text_offset(plot_out$data, layout)
    if(show_highlights && nrow(plot_out$data) > 0L){
        if(layout %in% c("daylight","ape")){
            warning("highlights not supported  for layout '",layout,"'",
                    call. = FALSE)
            return(plot_out)
        }
        subset <- plot_out$data$highlight
        highlight_nodes <- plot_out$data[subset,"node",drop=TRUE]
        hl_args <- .get_hightlight_args(
            highlight_nodes, colour_highlights_by)
        plot_out <- plot_out +
            do.call(geom_highlight, hl_args$args)
        if(!is.null(colour_highlights_by)){
            plot_out <- .resolve_plot_colours(
                plot_out,
                plot_out$data[subset, "colour_highlights_by", drop=TRUE],
                colour_highlights_by, fill = TRUE, na.value = "grey70")
            plot_out <- plot_out +
                new_scale_fill() +
                new_scale_colour()
        }
        if(show_highlight_label){
            subset <- plot_out$data$highlight &
                !is.na(plot_out$data$highlight_label)
            highlight_label_nodes <- plot_out$data[subset,"node",drop=TRUE]
            if(length(highlight_label_nodes) > 0L){
                subset_abbr <- plot_out$data[,"highlight_label",drop=TRUE] %in%
                    abbr_label
                subset[!subset_abbr] <- FALSE
                plot_out <- .add_label_abbreviations(
                    plot_out, "highlight_label", which(subset))
                cl_args <- .get_cladelab_args(
                    highlight_label_nodes, layout, highlight_font_size)
                plot_out <- plot_out +
                    do.call(geom_cladelab, cl_args$args)
                ################################################################
                # fix for geom_segment getting added by geom_cladelab even
                # though barsize = NA
                plot_out$layers <- plot_out$layers[-length(plot_out$layers)]
                ################################################################
            }
        }
    }
    plot_out
}

.plot_tree_edges <- function(
        plot_out, edge_colour_by, edge_size_by, line_alpha, line_width,
        line_width_range, layout){
    # assemble arg list
    edge_out <- .get_edge_args(
        edge_colour_by, edge_size_by, alpha = line_alpha, size = line_width,
        layout = layout)
    plot_out <- plot_out +
        do.call(geom_tree, edge_out$args) +
        theme_tree()
    plot_out <- .add_extra_guide_tree(
        plot_out, edge_size_by, line_width_range)
    # adjust edge colours
    if(!is.null(edge_colour_by)){
        plot_out <- .resolve_plot_colours(
            plot_out, plot_out$data$edge_colour_by, edge_colour_by,
            na.translate = FALSE)
    }
    plot_out
}

.plot_tree_node_points <- function(
        plot_out, show_tips, show_nodes, colour_by, shape_by, size_by,
        point_alpha, point_size, point_size_range){
    point_out <- .get_point_args(
        colour_by, shape_by, size_by, alpha = point_alpha, size = point_size)
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
    if(any(c(show_tips,show_nodes)) && !is.null(size_by)){
        if(is.numeric(plot_out$data$size_by)){
            SIZEFUN <- scale_size_continuous
        } else {
            SIZEFUN <- scale_size_discrete
        }
        plot_out <- plot_out +
            SIZEFUN(range = point_size_range)
    }
    # adjust point colours
    if(!is.null(colour_by)){
        plot_out <- .resolve_plot_colours(
            plot_out, plot_out$data$colour_by, colour_by, fill = point_out$fill,
            na.translate = FALSE)
    }
    plot_out
}

#' @importFrom ggtree geom_tiplab geom_nodelab
.plot_tree_node_labels <- function(
        plot_out, show_label, abbr_label, label_font_size){
    if(show_label){
        data <- plot_out$data
        label_data <- plot_out$data %>% drop_na(.data$node_label)
        #
        f_tip <- data$node %in% label_data$node & data$isTip
        f_node <- data$node %in% label_data$node & !data$isTip
        subset <- !is.na(plot_out$data$node_label)
        subset_abbr <- plot_out$data[,"node_label",drop=TRUE] %in%
            abbr_label
        subset[!subset_abbr] <- FALSE
        plot_out <- .add_label_abbreviations(
            plot_out, "node_label",  which(subset))
        if(any(f_tip)){
            # add tip labels
            plot_out <- plot_out +
                geom_tiplab(
                    mapping = aes_string(subset = f_tip, label = "node_label"),
                    offset = 0.01, size = label_font_size)
        }
        if(any(f_node)){
            # add node labels
            plot_out <- plot_out +
                geom_nodelab(
                    mapping = aes_string(subset = f_node, label = "node_label"),
                    size = label_font_size)
        }
    }
    plot_out
}

.add_abbr_guide <- function(plot_out){
    FUN <- function(col,data){
        abbr_col <- paste0("abbr_",col)
        if(!all(c(col, abbr_col) %in% colnames(data))){
            return(NULL)
        }
        ans <- data[!is.na(data[,abbr_col,drop=TRUE]),c(col,abbr_col)]
        colnames(ans) <- c("abbr","text")
        ans
    }
    abbr <- lapply(c("node_label","highlight_label"),FUN,plot_out$data)
    abbr <- abbr[!vapply(abbr,is.null,logical(1))]
    abbr <- Reduce(rbind,abbr)
    if(!is.null(abbr) && nrow(abbr) > 0L){
        abbr <- abbr[order(abbr$text),]
        keywidth <- max(1.5,max(nchar(abbr$abbr)) * 0.2)
        guide <- guide_legend(
            title = "Abbreviations",
            keywidth = keywidth,
            keyheight = 0.75,
            label.theme = element_text(size = 8),
            override.aes = list(fill = "transparent"),
            ncol = 1)
        plot_out <- plot_out +
            scale_discrete_identity(
                aesthetics = "label",
                name = "Abbreviations:",
                breaks = abbr$abbr,
                labels = abbr$text,
                guide = guide)
    }
    plot_out
}

.theme_plotTree <- function(plot){
    plot +
        theme(
            legend.background = element_rect(fill = "transparent",colour = NA),
            legend.box.background = element_rect(
                fill = "transparent",colour = NA),
            panel.background = element_rect(fill = "transparent",colour = NA),
            plot.background = element_rect(fill = "transparent",colour = NA),
            legend.text = element_text(size = 8))
}
