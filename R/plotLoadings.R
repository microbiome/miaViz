#' Plot feature loadings for TreeSummarizedExperiment/SingleCellExperiment objects
#' or feature loadings numeric matrix.
#'
#' This function is used after performing a reduction method. If TSE object is
#' given it retrieves the feature loadings matrix to plot values, tree can be added to heatmap.
#' Plotting with other layouts is possible as SCE objects and numeric matrices
#' does not include a tree.
#' 
#' @param x a
#'   \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'   x.
#' 
#' @param dimred \code{Character scalar}. Name of the reduction method if there are several in reducedDim slot. 
#'  
#' @param layout \code{Character scalar}. One way to plot feature loadings of \code{c("heatmap", "barplot")}. 
#'   (Default: \code{"barplot"})
#' 
#' @param n \code{Numeric scalar}. Number of features to be plotted.
#'   (Default: \code{10})
#'   
#' @param ncomponents \code{Numeric scalar}. Number of components must be lower or equal 
#'   to the number of components choosen in the reduction method.
#'   (Default: \code{5})
#' 
#' @param tree.name \code{Character scalar}. Value specifying a rowTree from a
#'   TreeSummarizedExperiment object. 
#'   (Default: \code{"phylo"})
#'   
#' @param rank \code{Character scalar}. Value specifying a rank from taxonomyRanks
#'   (Default: \code{NULL})
#'   
#' @param add.tree \code{Logical scalar}. Whether or not to add tree to heatmap layout.
#'   (Default: \code{FALSE})
#'   
#' @param ... additional arguments for plotting. See 
#'   \code{\link{mia-plot-args}} for more details i.e. call \code{help("mia-plot-args")}     
#' 
#' @details
#' 
#' Inspired by the \code{plotASVcircular} method using phyloseq 
#' and has been converted to use TreeSummarizedExperiment/SingleCellExperiment objects.
#' TreeSummarizedExperiment/SingleCellExperiment objects are expected to have content in reducedDim slot.
#' It is impossible to add tree if only the matrix is given. Number of features must be reduced
#' before calling function or it will not be understandable. For example
#' agglomerating data by rank or prevalence (see examples).
#' 
#' @return 
#' A \code{ggplot2} object. A circular plot annotated with TreeSummarizedExperiment object.
#'
#' @name plotLoadings
#' @export
#' 
#' @author
#' Ely Seraidarian
#' Contact: \url{microbiome.github.io}
#'
#' @examples
#' 
#' # Plotting feature loadings with tree
#' library(mia)
#' library(ggtree)
#' library(scater)
#' data("GlobalPatterns", package = "mia")
#' tse <- GlobalPatterns
#' tse <- transformAssay(tse, method = "clr", pseudocount = 1)
#' tse <- agglomerateByPrevalence(tse, rank="Phylum", update.tree = TRUE)
#' tse <- runPCA(tse, ncomponents = 5, assay.type = "clr")
#' plotLoadings(tse, dimred = "PCA", layout = "heatmap", add.tree = TRUE)
#' 
#' # Plotting without tree as a barplot
#' loadings_matrix <- attr(reducedDim(tse, "PCA"), "rotation")
#' plotLoadings(loadings_matrix)
#' 
#' # Plotting more features
#' plotLoadings(loadings_matrix, n = 12)
#' 
#' # Plotting without tree as a heatmap
#' plotLoadings(loadings_matrix, layout = "heatmap")
#' 
#' # Plotting with less components
#' tse <- runPCA(tse, ncomponents = 4, assay.type = "clr")
#' loadings_matrix <- attr(reducedDim(tse, "PCA"), "rotation")
#' plotLoadings(loadings_matrix, ncomponents = 4)
#'
#' # Plotting if loadings matrix name has been changed
#' tse <- runPCA(tse, name = "myPCAmatrix", ncomponents = 5, assay.type = "clr")
#' plotLoadings(tse, dimred = "myPCAmatrix")
#' 
#' # Plotting tree with taxonomic rank classification
#' tse <- runPCA(tse, ncomponents = 5, assay.type = "clr")
#' plotLoadings(tse, dimred = "PCA", layout = "heatmap", add.tree = TRUE, rank = "Phylum")
#' 
#' # Plotting after performing LDA 
#' library(topicmodels)
#' tse <- addLDA(tse)
#' plotLoadings(tse, dimred = "LDA", ncomponents = 2)
NULL

#' @rdname plotLoadings
setGeneric("plotLoadings", signature = c("x"),
    function(x, ...) 
        standardGeneric("plotLoadings"))


#' @rdname plotLoadings
#' @export 
setMethod("plotLoadings", signature = c(x = "TreeSummarizedExperiment"),
    function(
        x, dimred, layout = "barplot", n = 10, ncomponents = 5,
        tree.name = "phylo", rank = NULL, add.tree = FALSE, ...) {
        # Check that there are reducedDim
        if( length(reducedDims(x)) == 0 ){
          stop("No reducedDims found.", call. = FALSE)
        }
        # Check dimred. It must be either string specifying the name of
        # reducedDim or an index of reducedDim.
        if( !((.is_a_string(dimred) && dimred %in% reducedDimNames(x)) ||
              .is_an_integer(dimred) && dimred > 0 &&
              dimred <= length(reducedDims(x)) ) ){
          stop("'dimred' must be a string or an integer.", call. = FALSE)
        }
        # Check add.tree
        if( !.is_a_bool(add.tree) ){
            stop("'add.tree' must be TRUE or FALSE.", call. = FALSE)
        }
        # Check that tree.name. If user wants to add tree, the tree name must
        # specify a tree
        if( !(.is_a_string(tree.name) &&
                (add.tree && tree.name %in% rowTreeNames(x))) ){
            stop(
                "'tree.name' must be a string specifying a rowTree.",
                call. = FALSE)
        }
        # Check rank
        if( !(is.null(rank) ||
                (.is_a_string(rank) && rank %in% colnames(rowData(x)))) ){
            stop(
                "'rank' must be NULL or a column from rowData(x).",
                call. = FALSE)
        }
        #
        # Get loadings matrix
        mat <- .get_loadings_matrix(x, dimred)
        if( add.tree && layout == "heatmap" ){
            # Create dataframe for tree plotting
            data_list <- .get_loadings_tree_data(mat, x, tree.name, rank)
            tree <- data_list[["tree"]]
            mat <- data_list[["loadings"]]
            # Plot tree with feature loadings
            p <- .loadings_tree_plotter(mat, tree, rank, ...)
        } else {
            # Utilize matrix method to create a plot
            p <- plotLoadings(
                loadings_matrix, layout = layout, n = n,
                ncomponents = ncomponents, ...) 
        }
    return(p)
    }
)

#' @rdname plotLoadings
#' @export 
setMethod("plotLoadings", signature = c(x = "SingleCellExperiment"),
    function(x, dimred, layout = "barplot", n = 10, ncomponents = 5, ...){
        # Check that there are reducedDim
        if( length(reducedDims(x)) == 0 ){
            stop("No reducedDims found.", call. = FALSE)
        }
        # Check dimred. It must be either string specifying the name of
        # reducedDim or an index of reducedDim.
        if( !((.is_a_string(dimred) && dimred %in% reducedDimNames(x)) ||
                .is_an_integer(dimred) && dimred > 0 &&
                dimred <= length(reducedDims(x)) ) ){
            stop("'dimred' must be a string or an integer.", call. = FALSE)
        }
        # Get loadings matrix
        mat <- .get_loadings_matrix(x, dimred)
        # Utilize matrix method to create a plot
        p <- plotLoadings(
            mat, layout = layout, n = n, ncomponents = ncomponents, ...) 
        return(p)
    }
)

#' @rdname plotLoadings
#' @export 
setMethod("plotLoadings", signature = c(x = "matrix"),
    function(x, layout = "barplot", n = 10, ncomponents = 5, ...) {
        # Input check
        .check_loadings_matrix(
            x, layout = layout, n = n, ncomponents = ncomponents, ...)
        #
        # Get data for plotting
        df <- .get_loadings_plot_data(x, layout, n, ncomponents)
        # Create a plot
        p <- .plot_loadings(df, layout = layout, ...)
        return(p)
    }
)

################################ HELP FUNCTIONS ################################

# This function fetches loadings matrix from TreeSE object. The loadings
# are fetched from attributes of reducedDim whcih means that the result must be
# first calculated with standardized method.
.get_loadings_matrix <- function(x, dimred, ...){
    # Get reducedDim
    reddim <- reducedDim(x, dimred)
    # Get loadings matrix.
    attr_names <- names(attributes(reddim))
    loading_names <- c("rotation", "loadings")
    attr_name <- attr_names[ attr_names %in% loading_names ]
    if( length(attr_name) != 1 ) {
        stop("Loadings cannot be found.", call. = FALSE)
    }
    mat <- attr(reddim, attr_name)
    # Convert to data.frame
    mat <- as.data.frame(mat)
    return(mat)
}

# This functions checks that loadings matrix is correct
.check_loadings_matrix <- function(mat, layout, n, ncomponents, ...) {
    # Check layout
    if( !(.is_a_string(layout) && layout %in% c("barplot", "heatmap")) ){
        stop("'layout' must be 'barplot' or 'heatmap',", call. = FALSE)
    }
    # Check n
    if( !(.is_an_integer(n) && n > 0 && n <= nrow(mat)) ){
        stop(
            "'n' must be a positive integer less than or equal to the total ",
            "number of features.", call. = FALSE)
    }
    # Check ncomponents
    if( !(.is_an_integer(ncomponents) && ncomponents > 0 &&
            ncomponents <= ncol(mat)) ){
        stop(
            "'ncomponents' must be a positive integer less than or equal to ",
            "the total number of components.", call. = FALSE)
    }
    return(NULL)
}

# This function manipulates the loadings data into correct format. The output
# is data.frame in long format directly usable for ggplot.
#' @importFrom dplyr %>% rownames_to_column pivot_longer
.get_loadings_plot_data <- function(df, layout, n, ncomponents) {
    # Transform into a dataframe
    df <- as.data.frame(df)
    # Keep only the number of components needed
    df <- df[ , seq_len(ncomponents), drop = FALSE]
    # If the layout is barplot, choose top features for each component
    if( layout %in% c("barplot") ){
        res <- lapply(seq_len(ncomponents), .process_component, df = df, n = n)
        # Combine to single data.frame
        res <- do.call(rbind, res)
    } else{
        # For heatmap, the whole data.frame is just converted into long format.
        components <- colnames(df)
        res <- df %>%
            rownames_to_column(var = "Feature") %>%
            pivot_longer(
                cols = components, 
                names_to = "PC", 
                values_to = "Value")
    }
    # Convert into data.frame
    res <- as.data.frame(res)
    return(res)
}

# This function subsets the data so that it selects top features that have the
# greatest loadings for single component. 
.process_component <- function(i, df, n) {
    # Get order of loadings based on absolute value
    ind <- order(-abs(df[[i]]))
    # Get top n values
    ind <- ind[seq_len(n)]
    # Get the sorted data of single PC
    df <- df[ind, i, drop = FALSE]
    # Add PC number to data.frame so that each row can be identified to belong
    # to certain PC after merging the data into single dataset
    df[["PC"]] <- colnames(df)
    # Rename so that the colnames is same for all components
    colnames(df)[[1]] <- "Value"
    # Add rownames
    df[["Feature"]] <- rownames(df)
    return(df)
}

# This functions plots a data.frame in barplot or heatmap layout.
#' @importFrom tidytext scale_y_reordered reorder_within
#' @importFrom ggplot2 geom_tile scale_fill_gradient2 geom_bar
.plot_loadings <- function(df, layout = "barplot", ...) {
    # Initialize a plot
    plot_out <- ggplot(df)
    # Either create a heatmap or barplt
    if( layout == "heatmap" ){
        plot_out <- plot_out +
            # Create a heatmap
            geom_tile(
                mapping = aes(x = PC, y = Feature, fill = Value),
                position = position_identity()
                ) +
            # Adjust color scale
            scale_fill_gradient2(
                limits = c(-1, 1),
                low = "darkslateblue", mid = "white", high = "darkred"
                )
            
    } else{
        plot_out <- plot_out +
            # Create a bar plot. Create unique facets for each PC. Each PC can
            # have unique set of features. To reorder features by each facet,
            # we use reorder_within() and scale_y_reordered().
            geom_bar(
                mapping = aes(
                    x = Value, y = reorder_within(Feature, Value, PC)),
                stat = "identity",
                width = 0.8
                ) +
            scale_y_reordered() +
            facet_wrap(~ PC, scales = "free") +
            labs(x = "Value", y = "Feature") 
        
    }
    # Adjust theme
    plot_out <- plot_out +
        theme_minimal()
    return(plot_out)
}

# This function retrieves the data for tree + heatmap plotting. The output
# is a list that includes tree and data.frame in wide format.
#' @importFrom ggtree ggtree
.get_loadings_tree_data <- function(df, x, tree.name, rank) {
    # Retrieve rowTree
    phylo <- rowTree(x, tree.name)
    # Subset data based on the tree
    ind <- rowLinks(x)[["whichTree"]] == tree.name
    if( any(!ind) ){
        warning("Data is subsetted.", call. = FALSE)
        # Subset both TreeSE and loadings dfrix
        x <- x[ind, ]
        df <- df[ind, ]
    }
    # Add feature to column
    df[["Feature"]] <- rownames(df)
    # Instead of rownames, user can also specify a column from rowData to be
    # plotted
    if( !is.null(rank) ){
        df[["Feature"]] <- rowData(x)[[rank]]
    }
    # Check that there are not too many features to plot. If there are too many
    # rank values (or rownames), it is not possible to plot.
    if( length(unique(df[["Feature"]])) > 100 ){
        stop(
            "Too many features to plot. Consider specifying 'rank'.",
            call. = FALSE)
    }
    # Add rowlinks to data
    rownames(df) <- rowLinks(x)[["nodeLab"]]
    res <- list(tree = phylo, loadings = df)
    return(res)
}

# This function is for plotting tree with heatmap. It utilizes ggtree package.
#' @importFrom ggtree ggtree gheatmap
#' @importFrom ggnewscale new_scale_fill
#' @importFrom viridis scale_fill_viridis_d
#' @importFrom ggplot2 scale_fill_gradient2
.loadings_tree_plotter <- function(
    df, tree, rank, rank.title = ifelse(!is.null(rank), rank, "Feature"),
    ...) {
    # Check rank.title
    if( !.is_a_string(rank.title) ){
        stop("'rank.title' must be a string.", call. = FALSE)
    }
    #
    # Get features
    features <- df[ , colnames(df) %in% c("Feature"), drop = FALSE]
    # Get loadings
    loadings <- df[ , !colnames(df) %in% c("Feature"), drop = FALSE]
    # Create a tree plot
    plot_out <- ggtree(tree, layout = "circular")
    # Add first inner circle (features)
    plot_out <- gheatmap(
        p = plot_out, data = features, width = 0.1, colnames_angle = 90)
    # Adjust color scale for discrete feature values
    plot_out <- plot_out +
        scale_fill_viridis_d(option = "D", name = rank.title)
    # Add outer circles (feature loadings)
    plot_out <- plot_out + new_scale_fill()
    plot_out <- gheatmap(
        p = plot_out, data = loadings, offset = 0.1, width = 0.3,
        colnames_angle = 90
        )
    # Adjust color scale in continuous scale
    plot_out <- plot_out +
        scale_fill_gradient2(
            limits = c(-1, 1),
            low = "darkslateblue", mid = "white", high = "darkred",
            name = "Value"
        )
    return(plot_out)
}
