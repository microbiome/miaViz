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
#' @param ... additional arguments for plotting.
#' \itemize{
#'   \item \code{layout}: layout for the plotted tree. See
#'   \code{\link[ggtree:ggtree]{ggtree}} for details.
#'
#'   \item \code{relabel.tree}: \code{Logical scalar}. Should the tip labels be
#'   relabelec using the output of
#'   \code{getTaxonomyLabels(x, with_rank = TRUE)}? (Default: \code{FALSE})
#'
#'   \item \code{order.tree}: \code{Logical scalar}. Should the tree be ordered
#'   based on alphabetic order of taxonomic levels? (Default: \code{FALSE})
#'
#'   \item \code{levels.rm}: \code{Logical scalar}. Should taxonomic level
#'   information be removed from labels? (Default: \code{FALSE})
#'
#'   \item \code{show.label}, \code{show.highlights},
#'   \code{show.highlight.label}, \code{abbr.label} \code{logical vector},
#'   \code{integer vector}. or \code{character vector}. If a \code{logical}
#'   scalar is given, should tip labels be plotted or if a
#'   logical vector is provided, which labels should be shown? If an
#'   \code{integer} or \code{character} vector is provided, it will be converted
#'   to a logical vector. The \code{integer} values must be in the range of 1
#'   and number of nodes, whereas the values of a \code{character} vector must
#'   match values of the \code{label} column in the node data. In case of a
#'   \code{character} vector only values corresponding to actual labels will be
#'   plotted and if no labels are provided no labels will be shown. (Default:
#'   \code{FALSE})
#'
#'   \item \code{add.legend}: \code{Logical scalar}. Should legends be plotted?
#'   (Default: \code{TRUE})
#'
#'   \item \code{edge.colour.by}: \code{Character scalar}. Specification of a
#'   column metadata field or a feature to colour tree edges by.
#'   (Default: \code{NULL})
#'
#'   \item \code{edge.size.by}: \code{Character scalar}. Specification of a
#'   column metadata field or a feature to size tree edges by.
#'   (Default: \code{NULL})
#'
#'   \item \code{colour.by}: \code{Character scalar}. Specification of a
#'   column metadata field or a feature to colour tree nodes by.
#'   (Default: \code{NULL})
#'
#'   \item \code{shape.by}: \code{Character scalar}. Specification of a
#'   column metadata field or a feature to shape tree nodes by.
#'   (Default: \code{NULL})
#'
#'   \item \code{size.by}: \code{Character scalar}. Specification of a
#'   column metadata field or a feature to size tree tips by.
#'   (Default: \code{NULL})
#'
#'   \item \code{show.tips}: \code{Logical scalar}. Whether to show
#'   tip points. (Default: \code{FALSE})
#'
#'   \item \code{show.nodes}: \code{Logical scalar}. Whether to show
#'   node points. (Default: \code{FALSE})
#'
#'   \item \code{colour.highlights.by}: \code{Logical scalar}. Should the
#'   highlights be colour differently? If \code{show.highlights = TRUE},
#'   \code{colour_highlights} will be set to \code{TRUE} as default.
#'   (Default: \code{FALSE})
#'
#'   \item \code{assay.type}: \code{Character scalar}. Specifies which assay to
#'   obtain expression values from, for use in point aesthetics.
#'   (Default: \code{"counts"})
#'
#'   \item \code{other.fields}: \code{Character vector}. Additional fields to
#'   include in the node information without plotting them.
#'   (Default: \code{NULL})
#' }
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
#'     altExp(GlobalPatterns,"Genus"))
#' rowData(altExp(GlobalPatterns,"Genus"))$log_mean <- log(
#'     rowData(altExp(GlobalPatterns,"Genus"))$mean)
#' rowData(altExp(GlobalPatterns,"Genus"))$detected <- rowData(
#'     altExp(GlobalPatterns,"Genus"))$detected / 100
#' top_genus <- getTop(
#'     altExp(GlobalPatterns,"Genus"),
#'     method = "mean",
#'     top = 100L,
#'     assay.type = "counts"
#' )
#' #
#' x <- altExp(GlobalPatterns,"Genus")
#' plotRowTree(
#'     x[rownames(x) %in% top_genus,],
#'     tip.colour.by = "log_mean", tip.size.by = "detected"
#' )
#'
#' # plot with tip labels
#' plotRowTree(
#'     x[rownames(x) %in% top_genus,],
#'     tip.colour.by = "log_mean",
#'     tip.size.by = "detected",
#'     show.label = TRUE
#' )
#' # plot with selected labels
#' labels <- c("Genus:Providencia", "Genus:Morganella", "0.961.60")
#' plotRowTree(
#'     x[rownames(x) %in% top_genus,],
#'     tip.colour.by = "log_mean",
#'     tip.size.by = "detected",
#'     show.label = labels,
#'     layout = "rectangular"
#' )
#'
#' # plot with labeled edges
#' plotRowTree(
#'     x[rownames(x) %in% top_genus,],
#'     edge.colour.by = "Phylum",
#'     tip.colour.by = "log_mean"
#' )
#' # if edges are sized, colours might disappear depending on plotting device
#' plotRowTree(
#'     x[rownames(x) %in% top_genus,],
#'     node.colour.by = "Phylum",
#'     edge.size.by = "detected",
#'     edge.colour.by = "log_mean"
#' )
#'
#' # aggregating data over the taxonomic levels for plotting a taxonomic tree
#' # please note that the original tree of GlobalPatterns is dropped by
#' # unsplitByRanks
#' altExps(GlobalPatterns) <- splitByRanks(GlobalPatterns)
#' top_phyla <- getTop(
#'     altExp(GlobalPatterns,"Phylum"),
#'     method = "mean",
#'     top = 10L,
#'     assay.type="counts"
#' )
#' altExps(GlobalPatterns) <- lapply(altExps(GlobalPatterns), addPerFeatureQC)
#' altExps(GlobalPatterns) <- lapply(
#'     altExps(GlobalPatterns), function(y){
#'         rowData(y)$log_mean <- log(rowData(y)$mean)
#'         rowData(y)$detected <- rowData(y)$detected / 100
#'         return(y)
#'     })
#' x <- unsplitByRanks(GlobalPatterns)
#' x <- addHierarchyTree(x)
#'
#' highlights <- c(
#'     "Phylum:Firmicutes","Phylum:Bacteroidetes",
#'     "Family:Pseudomonadaceae","Order:Bifidobacteriales")
#' plotRowTree(
#'     x[rowData(x)$Phylum %in% top_phyla,],
#'     tip.colour.by = "log_mean",
#'     node.colour.by = "log_mean",
#'     show.highlights = highlights,
#'     show.highlight.label = highlights,
#'     colour.highlights.by = "Phylum"
#' )
#'
#' # If you do not want to show internal nodes
#' plotRowTree(
#'     x[rowData(x)$Phylum %in% top_phyla,],
#'     edge.colour.by = "Phylum",
#'     edge.size.by = "detected",
#'     node.colour.by = "log_mean",
#'     show.nodes = FALSE
#' )
#'
NULL

#' @rdname plotTree
#' @export
setMethod("plotColTree", signature = c(x = "TreeSummarizedExperiment"),
    function(x, tree.name = "phylo", ...){
        p <- .plot_row_column_tree(
            x, tree.name = tree.name, type = "column", ...)
        return(p)
    }
)
#' @rdname plotTree
#' @export
setMethod("plotRowTree", signature = c(x = "TreeSummarizedExperiment"),
    function(x, tree.name = "phylo", ...){
        p <- .plot_row_column_tree(
            x, tree.name = tree.name, type = "row", ...)
        return(p)
    }
)

################################ HELP FUNCTIONS ################################

# This function is general function for creating a tree plot from TreeSE object.
# It is utilized for both row and column trees.
.plot_row_column_tree <- function(x, tree.name, type, ...){
    # We want to visualize only tips that include data. That is why we subset
    # the data.
    x <- .get_object_and_trimmed_tree(
        x, tree.name = tree.name, type = type, ...)
    # Get tree from TreeSE
    tree_FUN <- switch(type, row = rowTree, column = colTree, stop("."))
    tree <- tree_FUN(x, tree.name)
    # Get tree as a table format
    df <- .get_tree_data(tree)
    # If user wants, we add label information to the table.
    df <- .add_tree_node_labels(df, ...)
    # User can add colors to highlight specific branches. This function adds the
    # info to the table.
    df <- .add_tree_highlights(df, ...)
    # User can add labels for highlighted sectors.
    df <- .add_tree_highlight_labels(df, ...)
    # Add data from rowData, e.g., for edge colors
    args <- .incorporate_tree_vis(df = df, x = x, type = type, ...)
    # Combine point and node formatting into single column
    args <- do.call(.combine_tree_point_formatting, args)
    # Modify the argument list so that it is ready for plotting
    args[["tree"]] <- tree
    args <- do.call(.create_treedata_for_plotting, args)
    # Create a plot
    p <- do.call(.tree_plotter, args)
    return(p)
}

# This function subsets the TreeSE to
#' @importFrom ape keep.tip as.phylo drop.tip
#' @importFrom tidytree as_tibble
.get_object_and_trimmed_tree <- function(
        x,
        tree.name = "phylo",
        type = c("row", "column"),
        relabel.tree = relabel, relabel = FALSE,
        order.tree = order, order = FALSE, ...
        ){
    # Input check
    type <- match.arg(type)
    # Check that tree exists
    check_FUN <- switch(
        type, row = .check_rowTree_present, column = .check_colTree_present)
    temp <- check_FUN(tree.name, x)
    if(!.is_a_bool(relabel)){
        stop("'relabel.tree' must be either TRUE or FALSE.", call. = FALSE)
    }
    if(!.is_a_bool(order.tree)){
        stop("'order.tree' must be either TRUE or FALSE.", call. = FALSE)
    }
    #
    # Get correct functions based on the margin/direction
    tree_FUN <- switch(type, row = rowTree, column = colTree, stop("."))
    links_FUN <- switch(type, row = rowLinks, column = colLinks, stop("."))
    dimnames_FUN <- switch(type, row = rownames, column = colnames, stop("."))
    add_names_FUN <- switch(
        type, row = `rownames<-`, column = `colnames<-`, stop("."))
    # Check that the tree is compatible with the data, i.e., rows are matched
    # with the tree.
    links_FUN <- switch(type, row = rowLinks, column = colLinks, stop("."))
    links <- links_FUN(x)
    ind <- links[["whichTree"]] == tree.name
    if( all(!ind) ){
        stop("Tree does not have any ", type, "s to plot.", call. = FALSE)
    }
    # Get only those rows/columns that are found from the tree
    if( type == "row" ){
        x <- x[ind, ]
    } else{
        x <- x[, ind]
    }
    # Get tree and links
    tree <- tree_FUN(x, tree.name)
    links <- links_FUN(x)

    # Remove those tips that are not included in the data
    args <- list(x, links[["nodeLab"]], tree.name)
    names(args) <- c(
        "x",
        paste0(type, "Leaf"),
        paste0("which", .capitalize(type), "Tree"))
    x <- do.call(subsetByLeaf, args)

    # Get tree, links and row/colnames
    tree <- tree_FUN(x)
    links <- links_FUN(x)
    dimnames <- dimnames_FUN(x)
    # Get tree as table and get which node represent which row/col
    tree_data <- as_tibble(tree)
    m <- match(links$nodeNum,tree_data$node)
    node_labels <- tree_data$label[m]
    # If user wants to rename rows/cols or if some nodes cannot be found from
    # rows/cols
    if( relabel.tree || !all(node_labels %in% dimnames) ){
        # Rename rows/cols
        new_node_labels <- getTaxonomyLabels(
            x, with_rank = TRUE, resolve_loops = TRUE)
        x <- add_names_FUN(x, new_node_labels)
    }
    # Check if there are rows/cols that are ununique. If there are, make them
    # unique.
    if( anyDuplicated(rownames(x)) ){
        warning(
            "Data includes ununique ", type, "s. Making them unique.",
            call. = FALSE)
        x <- add_names_FUN(x, make.unique(dimnames_FUN(x)))

    }
    # Rename labels of tree with row/colnames
    tree_data$label[m] <- dimnames_FUN(x)
    # Check if there are nodes that are not unique
    if( anyDuplicated(tree_data$label[-m]) ){
        warning(
            "Tree includes ununique nodes. Making them unique.", call. = FALSE)
        tree_data$label[-m] <- make.unique( tree_data$label[-m] )
    }

    # Convert tree data back to tree-format
    tree <- as.phylo(tree_data)
    # If specified, order the tree based on alphabetical order
    if( order.tree ){
        tree <- .order_tree(tree)
    }

    # Add tree back to TreeSE
    args <- list(x, tree, tree.name)
    names(args) <- c(
        "x",
        ifelse(type == "row", "rowTree", "colTree"),
        ifelse(type == "row", "whichRowTree", "whichColTree")
        )
    x <- do.call(changeTree, args)
    return(x)
}

# This function sorts the tree so that tips are in alphabetical order.
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
    return(tree)
}

# This function retrieves labels from the tree for each specified node.
#' @importFrom tidytree child
.get_tree_labels_for_ordering <- function(tree_data, node){
    children <- child(tree_data, node)
    labels <- ""
    if(nrow(children) > 0L){
        labels <- children$label
        add_labels <- lapply(
            children$node,
            .get_tree_labels_for_ordering,
            tree_data = tree_data)
        labels <- unlist(
            mapply(paste,labels,add_labels,sep="__:__",SIMPLIFY = FALSE),
            use.names = FALSE)
    }
    return(labels)
}

# This function removes taxonomy ranks from labels.
.remove_taxonomic_level_from_labels <- function(labels){
    for(rank in TAXONOMY_RANKS){
        labels <- gsub(paste0(rank,":"),"",labels,ignore.case = TRUE)
    }
    return(labels)
}

# This function controls how to show nodes. User can either show all nodes or
# specify nodes to show.
#' @importFrom tidygraph activate
#' @importFrom dplyr mutate
.add_tree_node_labels <- function(
        df, show.label = show_label, show_label = FALSE,
        levels.rm = levels_rm, levels_rm = FALSE, ...){
    # To disable "no visible binding for global variable" message in cmdcheck
    node.label <- NULL
    # Input check
    if(!.is_a_bool(show.label)){
        if( (!is.logical(show.label) && !is.character(show.label) &&
                !is.numeric(show.label)) || is.null(show.label)){
            stop("'show.label' must be either TRUE or FALSE or logical, ",
                "integer or character vector. Character alues should match ",
                "the label of the tree.", call. = FALSE)
        }
    }
    if(!.is_a_bool(levels.rm)){
        stop("'level.rm' must be either TRUE or FALSE.", call. = FALSE)
    }
    #
    # Check that show.label is correct. If user specifies labels to show labels,
    # there must be labels present
    if( !(.is_a_bool(show.label) && !show.label) &&
            !"label" %in% colnames(df) &&
            !(is.character(show.label) && length(show.label) == nrow(df))) {
        warning("'show.label' is specified but no labels present in a tree.",
                call. = FALSE)
    }
    #
    # If tree includes labels, add them to new column to keep book on labels.
    # "node_label" is used in plotting and "label" is kept untouched.
    df[["node_label"]] <- NA_character_
    if( "label" %in% colnames(df)) {
        df[["node_label"]] <- df[["label"]]
    }
    # Start to remove labels based on criteria.
    # If user has specified not to show labels at all
    if( .is_a_bool(show.label) && !show.label ){
        df[["node_label"]] <- NA_character_
    }
    # If user has specified new labels for each node
    if( is.character(show.label) && length(show.label) == nrow(df) ){
        df[["node_label"]] <- node.label
    }
    # If user has specified with character vector which labels to show
    if( is.character(show.label) && length(show.label) != nrow(df) ){
        df[ !df[["node_label"]] %in% show.label, "node_label"] <- NA_character_
    }
    # If user has specified with numeric vector which labels to show
    if( is.numeric(show.label) ){
        df[ !df[["node"]] %in% show.label,  "node_label"] <- NA_character_
    }
    # If user has specified with boolean vector which labels to show
    if( is.logical(show.label) && !.is_a_bool(show.label) ){
        df[ !show.label,  "node_label"] <- NA_character_
    }
    # Check if uswer wanted to show labels but none is available
    if( all(is.na(df[["node_label"]])) &&
            .is_a_bool(show.label) && show.label ){
        warning("No labels to plot.", call. = FALSE)
    }
    # If user wants to remove taxonomy levels from labels
    if( any(!is.na(df[["node_label"]])) && levels.rm ){
        df[["node_label"]] <- .remove_taxonomic_level_from_labels(
            df[["node_label"]])
    }
    return(df)
}

# This function controls which branches are highlighted. In the plot, these
# highlights are show as a sectors behind the tree.
#' @importFrom tidygraph activate
#' @importFrom dplyr mutate
.add_tree_highlights <- function(
        df, show.highlights = show_highlights, show_highlights = FALSE, ...){
    if(!.is_a_bool(show.highlights)){
        if( (!is.logical(show.highlights) && !is.character(show.highlights) &&
                !is.numeric(show.highlights)) || is.null(show.highlights)){
            stop("'show.label' must be either TRUE or FALSE or logical, ",
                "integer or character vector. Character alues should match ",
                "the label of the tree.", call. = FALSE)
        }
    }
    df[["highlight"]] <- FALSE
    # Check show.logical is correct if it specified a boolean value
    if( is.logical(show.highlights) &&
            !(length(show.highlights) %in% c(1L, nrow(df)) ) ){
        stop("If 'show.highlights' is logical, it must specify value for each ",
            "node in a tree.", call. = FALSE)
    }
    # Numeric vector should specify integers
    if( is.numeric(show.highlights) &&
            (any(show.highlights != as.integer(show.highlights)) ||
            min(show.highlights) < 1 || max(show.highlights) > nrow(df)) ){
        stop("If 'show.highlights' is numeric, values have to be ",
            "whole numbers and must be between 1 and the number of ",
            "nodes in the tree.", call. = FALSE)
    }
    #
    # User can specify single boolean value to specify whether to highlight
    # all branches
    if( .is_a_bool(show.highlights) && show.highlights ){
        df[["highlight"]] <- TRUE
    }
    # Or the value can be a logical vector
    if( is.logical(show.highlights) && length(show.highlights) > 1L ){
        df[["highlight"]] <- show.highlights
    }
    # It can be a numeric vector specifying node numbers
    if( is.numeric(show.highlights) ){
        df[df[["node"]] %in% show.highlights, "highlight"] <- TRUE
    }
    # It can be a character vector specifying node labels
    if( is.character(show.highlights) ){
        df[df[["label"]] %in% show.highlights, "highlight"] <- TRUE
    }
    # Give warning if user wanted to highlight branches but none was found
    if( !(.is_a_bool(show.highlights) && !show.highlights) &&
            all(!df[["highlight"]]) ){
        warning("No highlights to plot.", call. = FALSE)
    }
    return(df)
}

# In addition to coloring sectors, i.e., highlighting branches, user can add
# text for these highlights
#' @importFrom tidygraph activate
#' @importFrom dplyr mutate
.add_tree_highlight_labels <- function(
        df, show.highlight.label = show_highlight_label,
        show_highlight_label = FALSE, levels.rm = levels_rm,
        levels_rm = FALSE, ...){
    # To disable "no visible binding for global variable" message in cmdcheck
    show.highlights <- NULL
    # Input check
    if(!.is_a_bool(show.highlight.label)){
        if( (!is.logical(show.highlight.label) &&
                !is.character(show.highlight.label) &&
                !is.numeric(show.highlight.label)) ||
                is.null(show.highlight.label)){
            stop("'show.highlight.label' must be either TRUE or FALSE or ",
                "logical, integer or character vector. Character alues should ",
                "match the label of the tree.", call. = FALSE)
        }
    }
    if(!.is_a_bool(levels.rm)){
        stop("'level.rm' must be either TRUE or FALSE.", call. = FALSE)
    }
    # Check show.logical is correct if it specified a boolean value
    if( is.logical(show.highlight.label) &&
            !(length(show.highlight.label) %in% c(1L, nrow(df)) ) ){
        stop("If 'show.highlight.label' is logical, it must specify value for ",
            "each node in a tree.", call. = FALSE)
    }
    # Numeric vector should specify integers
    if( is.numeric(show.highlight.label) &&
        (any(show.highlight.label != as.integer(show.highlight.label)) ||
        min(show.highlight.label) < 1 || max(show.highlight.label) > nrow(df))){
        stop("If 'show.highlight.label' is numeric, values have to be ",
            "whole numbers and must be between 1 and the number of ",
            "nodes in the tree.", call. = FALSE)
    }
    # Initialize labels with NA
    df[["highlight_label"]] <- NA_character_
    # If user has enabled highlighting
    if( any(df[["highlight"]]) && "label" %in% colnames(df) ){
        df[df[["highlight"]], "highlight_label"] <- df[
            df[["highlight"]], "label"]
        # User can specify single boolean value to specify whether to highlight
        # all branches
        if( .is_a_bool(show.highlight.label) && !show.highlights ){
            df[["highlight_label"]] <- NA_character_
        }
        # Or the value can be a logical vector
        if( is.logical(show.highlight.label) &&
                length(show.highlight.label) > 1L ){
            df[!show.highlight.label, "highlight_label"] <- NA_character_
        }
        # It can be a numeric vector specifying node numbers
        if( is.numeric(show.highlight.label) ){
            df[!df[["node"]] %in% show.highlight.label, "highlight_label"] <-
                NA_character_
        }
        # It can be a character vector specifying node labels
        if( is.character(show.highlight.label) ){
            df[!df[["highlight_label"]] %in% show.highlight.label,
                "highlight_label"] <- NA_character_
        }
    }
    # Give warning if user wanted to highlight branches but none was found
    if( !(.is_a_bool(show.highlight.label) && !show.highlight.label) &&
            all(is.na(df[["highlight_label"]])) ){
        warning("No highlights labels to plot.", call. = FALSE)
    }
    # If user wants to remove taxonomy levels from labels
    if( any(!is.na(df[["highlight_label"]])) && levels.rm ){
        df[["highlight_label"]] <- .remove_taxonomic_level_from_labels(
            df[["highlight_label"]])
    }
    return(df)
}

# This function retrieves data from rowData, colData or alternatively from
# assay. This additional data is used for example for coloring edges.
.incorporate_tree_vis <- function(
        df, x, type,
        tip.colour.by = tip.color.by, tip.color.by = tip_colour_by,
        tip_colour_by = tip_color_by, tip_color_by = NULL,
        node.colour.by = node.color.by, node.color.by = node_colour_by,
        node_colour_by = node_color_by, node_color_by = NULL,
        #
        tip.shape.by = tip_shape_by, tip_shape_by = NULL,
        node.shape.by = node_shape_by, node_shape_by = NULL,
        #
        tip.size.by = tip_size_by, tip_size_by = NULL,
        node.size.by = node_size_by, node_size_by = NULL,
        # Edge and highlights are colored separately
        edge.colour.by = edge.color.by, edge.color.by = edge_colour_by,
        edge_colour_by = edge_color_by, edge_color_by = NULL,
        edge.size.by = edge_size_by, edge_size_by = NULL,
        colour.highlights.by = color.highlights.by,
        color.highlights.by = colour_highlights_by,
        colour_highlights_by = color_highlights_by,
        color_highlights_by = NULL,
        other.fields = other_fields, other_fields = NULL,
        ...){
    # Input check
    if( !(.is_a_string(edge.colour.by) || is.null(edge.colour.by)) ){
        stop("'edge.colour.by' must be a single character value.",
            call. = FALSE)
    }
    if( !(.is_a_string(edge.size.by) || is.null(edge.size.by)) ){
        stop("'edge.size.by' must be a single character value.", call. = FALSE)
    }
    if( !(.is_a_string(colour.highlights.by) ||
            is.null(colour.highlights.by)) ){
        stop("'colour.highlights.by' must be a single character value.",
            call. = FALSE)
    }
    if( !(is.character(other.fields) || is.null(other.fields)) ){
        stop("'other.fields' must be a character value.", call. = FALSE)
    }
    # Get all the variables into single vector
    variables <- c(
        tip_colour_by = tip.colour.by,
        tip_shape_by = tip.shape.by,
        tip_size_by = tip.size.by,
        node_colour_by = node.colour.by,
        node_shape_by = node.shape.by,
        node_size_by = node.size.by,
        edge_colour_by = edge.colour.by,
        edge_size_by = edge.size.by,
        colour_highlights_by = colour.highlights.by
        )
    names(other.fields) <- other.fields
    all_variables <- c(variables, other.fields)

    # Get function for getting links
    rowlinks_FUN <- switch(
        type,
        row = rowLinks,
        column = colLinks
    )

    # Retrieve info and create a table to add to tree data
    if( !is.null(all_variables) && length(all_variables > 0L) ){
        # Get variables
        metadata_df <- lapply(all_variables, function(var){
            .retrieve_variable(x, var, type, ...)
        })
        metadata_df <- do.call(cbind.data.frame, metadata_df)
        # Combine with rowLinks so that we can add the new data to tree data
        metadata_df[["node"]] <- rowlinks_FUN(x)[["nodeNum"]]
        df <- dplyr::left_join(
            df, metadata_df, by = "node", suffix = c("", ".y"))
    }

    # Create an argument list that is passed to plotting function
    args <- c(list(df = df), variables, list(...))
    return(args)
}

# This function gets single variable as input and it tries to fetch it from
# rowData, colData or assay.
.retrieve_variable <- function(x, var, type, assay.type = "counts", ...){
    name_FUN <- switch(
        type,
        row = colnames,
        column = rownames
    )
    rowdata_FUN <- switch(
        type,
        row = rowData,
        column = colData
    )
    # Check whether the variable is available in rowData or assay
    is_rowdata <- var %in% colnames(rowdata_FUN(x))
    is_sample <- var %in% name_FUN(x)
    # Get data from rowData (or colData if colTree)
    if( is_rowdata ){
        res <- rowdata_FUN(x)[[var]]
    } else if( is_sample ){
        # Get data from assay
        .check_assay_present(assay.type, x)
        res <- assay(x, assay.type)
        if( type == "row" ){
            res <- res[, var]
        } else{
            res <- res[var, ]
        }
    } else{
        # If not found, give error
        stop("The following variable cannot be found from ",
            ifelse(type == "row", "row", "col"), "Data(x) or from ",
            ifelse(type == "row", "column", "row"),
            " names: '", var, "'", call. = FALSE)
    }
    return(res)
}

# Tree tip and node coloring, shape or size must be in single column. That is
# why we combine these columns; there should be single column for each coloring,
# shape and size.
.combine_tree_point_formatting <- function(
        df,
        tip_colour_by = NULL,
        node_colour_by = NULL,
        tip_shape_by = NULL,
        node_shape_by = NULL,
        tip_size_by = NULL,
        node_size_by = NULL,
        ...){
    # Propagate tip info to nodes
    node_var <- c("node_colour_by", "node_shape_by", "node_size_by")
    for( var in node_var ){
        if( var %in% colnames(df) &&
            anyNA(df[[var]]) && !is.numeric(df[[var]]) ){
            df <- .propagate_to_internal_nodes(df, var = var)
        }
    }
    # Check if user wanted to show nodes, but we do not have coloring,
    # shaping, or size info for them
    if( any(node_var %in% colnames(df)) ){
        df_node <- df[
            !isTip(df, df[["node"]]), colnames(df) %in% node_var, drop = FALSE]
        if( all(colSums(is.na(df_node)) == nrow(df_node)) ){
            warning("The dataset seems to include only tips and no internal ",
                    "nodes were found. That is why 'node*by' arguments are ",
                    "ignored.", call. = FALSE)
        }
    }

    # Combine node and tip formatting into single column
    df <- .combine_tip_and_node(df, "colour")
    df <- .combine_tip_and_node(df, "shape")
    df <- .combine_tip_and_node(df, "size")
    # Combine variable names into single title
    colour_by <- .combine_tip_and_node_title(tip_colour_by, node_colour_by)
    shape_by <- .combine_tip_and_node_title(tip_shape_by, node_shape_by)
    size_by <- .combine_tip_and_node_title(tip_size_by, node_size_by)

    # Based on available columns, determine whether user wanted to show nodes
    # and labels
    show_tips <- any(
        c("tip_colour_by", "tip_shape_by", "tip_size_by") %in% colnames(df))
    show_nodes <- any(
        c("node_colour_by", "node_shape_by", "node_size_by") %in% colnames(df))

    # Create a list to forward to next function
    args <- list(
        df = df,
        colour_by = colour_by,
        shape_by = shape_by,
        size_by = size_by,
        show_tips = show_tips,
        show_nodes = show_nodes
    )
    args <- c(args, list(...))
    return(args)
}

# This function combines tip and node formatting columns into one single column.
# For instance, coloring of the points are combined into single column.
.combine_tip_and_node <- function(df, var){
    tip_ind <- isTip(df, df[["node"]])
    tip_col <- paste0("tip_", var, "_by")
    node_col <- paste0("node_", var, "_by")
    final_col <- paste0(var, "_by")
    temp <- rep(NA, nrow(df))
    if( tip_col %in% colnames(df) ){
        temp[tip_ind] <- df[tip_ind, ][[tip_col]]
    }
    if( node_col %in% colnames(df) ){
        temp[!tip_ind] <- df[!tip_ind, ][[node_col]]
    }
    if( any(!is.na(temp)) ){
        df[[final_col]] <- temp
    }
    return(df)
}

# This function combines varibale names for point formatting
.combine_tip_and_node_title <- function(var1, var2){
    var <- unique(c(var1, var2))
    var <- paste0(var, collapse = " & ")
    if( var == "" ){
        var <- NULL
    }
    return(var)
}

# due to a bug in ggtree/tidytree the treedata object needs to be constructed
# in a separate step
#
# also there is some data wrangling needed
#' @importFrom tidytree as.treedata isTip
.create_treedata_for_plotting <- function(df, tree, ...){
    # We do not have info on internal nodes, info is only for rows, i.e.,
    # usually for tips. Next step propagates the info to all nodes by finding
    # the first common ancestor, and correctly grouping the edges. This
    # does not work with numeric data.
    if( !is.null(df[["edge_colour_by"]]) && anyNA(df[["edge_colour_by"]]) &&
            !is.numeric(df[["edge_colour_by"]]) ){
        df <- .propagate_to_internal_nodes(df, var = "edge_colour_by")
    }
    # Propagate also highlight color info
    if( !is.null(df[["colour_highlights_by"]]) &&
            anyNA(df[["colour_highlights_by"]]) &&
            !is.numeric(df[["colour_highlights_by"]]) ){
        df <- .propagate_to_internal_nodes(df, var = "colour_highlights_by")
    }
    # Replace NA values with default shape or size
    df <- .na_replace_from_plot_data(
        df,
        if( "edge_size_by" %in% colnames(df) ) "edge_size_by",
        if( "shape_by" %in% colnames(df) ) "shape_by",
        if( "size_by" %in% colnames(df) ) "size_by"
    )
    # Tree cannot be build with duplicated column names
    colnames(df) <- colnames(df) |> make.unique()
    # From tibble, create treedata object
    df <- as.treedata(df)
    # tree needs to be restored since the original leave/tip/node orientation
    # is not compatible with ladderize = FALSE
    df@phylo <- tree
    # Based on label and highlight availability, determine whether to show them
    res <- c(list(...), df = df)
    res[["show_label"]] <- !all(is.na(df[["node_label"]]))
    res[["show_highlights"]] <- df[["highlight"]] |> any()
    res[["show_highlight_label"]] <- any(!is.na(df[["highlight_label"]]))
    return(res)
}

# This function propagates tip information to nodes. I.e., if we have Phyla info
# it is usually only for tips as rows of TreeSE represent usually tips. This
# function propagates the information to higher level nodes.
#' @importFrom tidytree groupOTU
.propagate_to_internal_nodes <- function(df, var){
    df <- groupOTU(
        df,
        split(df[["node"]], df[[var]]),
        group_name = "group")
    f_zero <- df[["group"]] != 0
    f_zero <- f_zero[!is.na(f_zero)]
    df[[var]][f_zero] <- as.character(df$group[f_zero])
    return(df)
}

#' @importFrom ggplot2 scale_size_identity
#' @importFrom ggtree ggtree geom_tree geom_tippoint geom_nodepoint groupOTU
#'   theme_tree
.tree_plotter <- function(
        df,
        # These arguments below comes from internal function. However, it might
        # be that not all are available as user did not specify them.
        show_tips,
        show_nodes,
        show_label,
        show_highlights,
        show_highlight_label,
        edge_colour_by = NULL,
        edge_size_by = NULL,
        colour_by = NULL,
        shape_by = NULL,
        size_by = NULL,
        colour_highlights_by = NULL,
        # These parameters are for modifying visuals of the tree plot
        layout = "circular",
        add.legend = TRUE,
        abbr.label = abbr_label, abbr_label = FALSE,
        order.tree = FALSE,
        line.alpha = line_alpha,
        line_alpha = 1,
        line.width = line_width,
        line_width = NULL,
        line.width.range = line_width_range,
        line_width_range = c(0.5,3),
        point.alpha = point_alpha,
        point_alpha = 1,
        point.size = point_size,
        point_size = 2,
        point.size.range = point_size_range,
        point_size_range = c(1,4),
        label.font.size = label_font_size,
        label_font_size = 3,
        highlight.font.size = highlight_font_size,
        highlight_font_size = 3,
        branch.length = "branch.length",
        # These are just catched so that they are not fed forward
        show.label = NULL,
        show.nodes = NULL,
        relabel.tree = NULL,
        levels.rm= NULL,
        show.highlights = NULL,
        show.highlight.label = NULL,
        tip_colour_by = NULL,
        tip_shape_by = NULL,
        tip_size_by = NULL,
        node_colour_by = NULL,
        node_shape_by = NULL,
        node_size_by = NULL,
        ...){
    # Check switches
    if(!.is_a_string(layout)){
        stop("'layout' must be a single character value.", call. = FALSE)
    }
    if(!.is_a_bool(add.legend)){
        stop("'add.legend' must be either TRUE or FALSE.", call. = FALSE)
    }
    if(!.is_a_bool(abbr.label)){
        if( (!is.logical(abbr.label) && !is.character(abbr.label) &&
                !is.numeric(abbr.label)) || is.null(abbr.label)){
            stop("'abbr.label' must be either TRUE or FALSE or logical, ",
                "integer or character vector. Character alues should match ",
                "the label of the tree.", call. = FALSE)
        }
    }
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
        df, ladderize = !order.tree, layout = layout,
        branch.length = branch.length, ...)
    # add highlights
    plot_out <- .plot_tree_plot_highlights(
        plot_out, layout, show_highlights, show_highlight_label, abbr.label,
        colour_highlights_by, highlight_font_size = highlight.font.size)
    # add tree and adjust edges
    plot_out <- .plot_tree_edges(
        plot_out, edge_colour_by, edge_size_by, line.alpha, line.width,
        line.width.range, layout)
    # add tip and node points
    plot_out <- .plot_tree_node_points(
        plot_out, show_tips, show_nodes, colour_by, shape_by, size_by,
        point.alpha, point.size, point.size.range)
    # add tip and node labels
    plot_out <- .plot_tree_node_labels(
        plot_out, show_label, abbr.label, label.font.size)
    # add additional guides
    plot_out <- .add_extra_guide(plot_out, shape_by, size_by)
    # add abbreviation guide
    plot_out <- .add_abbr_guide(plot_out)
    # add theme
    plot_out <- .theme_plotTree(plot_out)
    # optionally hide legends
    if (!add.legend) {
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
        ans <- highlight_data |>
            mutate(
                highlight_extendto = (max(.data$x) - .data$x) / 1.5,
                highlight_extendto = .data$highlight_extendto +
                    max(.data$x) + 0.07)
    } else if(layout %in% c("rectangular","slanted","ellipse","roundrect")){
        ans <- highlight_data |>
            mutate(
                highlight_extendto = (max(.data$x) - .data$x) / 1.5,
                highlight_extendto = .data$highlight_extendto +
                    max(.data$x) + 0.01)
    } else if(layout %in% c("dendrogram")){
        warning("highlights with layout `dendrogram` are buggy.")
        ans <- highlight_data |>
            mutate(
                highlight_extendto = .data$x / 1.5,
                highlight_extendto = (.data$highlight_extendto - 0.01) * -1)
    } else if(layout %in% c("inward_circular")){
        warning("highlights with layout `inward_circular` are buggy.")
        ans <- highlight_data |>
            mutate(
                highlight_extendto = (max(.data$x) - .data$x) / 1.5,
                highlight_extendto = .data$highlight_extendto + max(.data$x) +
                    0.07,
                highlight_extendto = .data$highlight_extendto * -1)
    } else {
        ans <- highlight_data |>
            mutate(highlight_extendto = .data$x)
    }
    ans
}

#' @importFrom dplyr mutate
.calc_highlight_label_text_offset <- function(label_data, layout){
    if(layout %in% c("fan","circular","radial")){
        ans <- label_data |>
            mutate(highlight_offset = .data$highlight_extendto - max(.data$x) +
                    0.015 - 0.07)
    } else if(layout %in% c("rectangular","slanted","ellipse","roundrect")){
        ans <- label_data |>
            mutate(highlight_offset = .data$highlight_extendto - max(.data$x) -
                    0.01)
    } else if(layout %in% c("dendrogram")){
        ans <- label_data |>
            mutate(highlight_offset = .data$highlight_extendto - 0.1)
    } else if(layout %in% c("inward_circular")){
        ans <- label_data |>
            mutate(highlight_offset = (.data$highlight_extendto *-1) -
                    max(.data$x) - 0.022)
    } else {
        ans <- label_data |>
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
            plot_out, plot_out$data[["edge_colour_by"]], edge_colour_by,
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
        label_data <- plot_out$data |> drop_na(.data$node_label)
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
