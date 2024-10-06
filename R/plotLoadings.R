#' Plot feature loadings for TreeSummarizedExperiment 
#' objects or feature loadings numeric matrix.
#'
#' This function is used after performing a reduction method. If \code{TreeSE}
#' is given it retrieves the feature loadings matrix to plot values.
#' A tree from \code{rowTree} can be added to heatmap layout.
#' 
#' @inheritParams plotTree
#' 
#' @param dimred \code{Character scalar}.  Determines the reduced dimension to
#' plot.
#'  
#' @param layout \code{Character scalar}. Determines the layout of plot. Must be
#' either \code{"barplot"} or \code{"heatmap"}. (Default: \code{"barplot"})
#'   
#' @param ncomponents \code{Numeric scalar}. Number of components must be lower
#' or equal to the number of components chosen in the reduction method.
#' (Default: \code{5})
#' 
#' @param add.tree \code{Logical scalar}. Whether to add tree to heatmap layout.
#' (Default: \code{FALSE})
#' 
#' @param row.var \code{NULL} or \code{Character scalar}. Specifies a
#' variable from \code{rowData} to plot with tree heatmap layout.
#' (Default: \code{NULL})
#'   
#' @param ... additional parameters for plotting.
#'   \itemize{
#'   \item \code{n}: \code{Integer scalar}. Number of features to be plotted.
#'   Applicable when \code{layout="barplot"}. (Default: \code{10}))
#'   
#'   \item \code{absolute.scale}: ("barplot") \code{Logical scalar}. Specifies
#'   whether a barplot should be visualized in absoltue scale.
#'   (Default: \code{TRUE})
#' }
#' 
#' @details
#' 
#' These method visualize feature loadings of dimension reduction results.
#' Inspired by the \code{plotASVcircular} method using \code{phyloseq}.
#' \code{TreeSummarizedExperiment} object is expected to have
#' content in \code{reducedDim} slot calculated with standardized methods from
#' \code{mia} or \code{scater} package.
#' 
#' @return 
#' A \code{ggplot2} object.
#'
#' @name plotLoadings
#' @export
#'
#' @examples
#' 
#' library(mia)
#' library(scater)
#' data("GlobalPatterns", package = "mia")
#' tse <- GlobalPatterns
#' 
#' # Calculate PCA
#' tse <- agglomerateByPrevalence(tse, rank="Phylum", update.tree = TRUE)
#' tse <- transformAssay(tse, method = "clr", pseudocount = 1)
#' tse <- runPCA(tse, ncomponents = 5, assay.type = "clr")
#' 
#' #' # Plotting feature loadings with tree
#' plotLoadings(tse, dimred = "PCA", layout = "heatmap", add.tree = TRUE)
#' 
#' # Plotting matrix as a barplot
#' loadings_matrix <- attr(reducedDim(tse, "PCA"), "rotation")
#' plotLoadings(loadings_matrix)
#' 
#' # Plotting more features but less components
#' plotLoadings(tse, dimred = "PCA", ncomponents = 2, n = 12)
#' 
#' # Plotting matrix as heatmap without tree
#' plotLoadings(loadings_matrix, layout = "heatmap")
#' 
#' # Plot with less components
#' plotLoadings(tse, "PCA", layout = "heatmap", ncomponents = 2)
#' 
NULL

#' @rdname plotLoadings
setGeneric("plotLoadings", signature = c("x"),
    function(x, ...) 
        standardGeneric("plotLoadings"))


#' @rdname plotLoadings
#' @export 
setMethod("plotLoadings", signature = c(x = "TreeSummarizedExperiment"),
    function(
        x, dimred, layout = "barplot", ncomponents = 5, tree.name = "phylo",
        row.var = NULL, add.tree = FALSE, ...) {
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
        if(add.tree && !(.is_a_string(tree.name) && 
        tree.name %in% rowTreeNames(x)) ){
            stop(
                "'tree.name' must be a string specifying a rowTree.",
                call. = FALSE)
        }
        # Check row.var
        if( !(is.null(row.var) ||
                (.is_a_string(row.var) && row.var %in% colnames(rowData(x)))) ){
            stop(
                "'row.var' must be NULL or a column from rowData(x).",
                call. = FALSE)
        }
        #
        # Get loadings matrix
        mat <- .get_loadings_matrix(x, dimred)
        if( add.tree && layout == "heatmap" ){
            # Create dataframe for tree plotting
            data_list <- .get_loadings_tree_data(mat, x, tree.name, row.var)
            tree <- data_list[["tree"]]
            mat <- data_list[["loadings"]]
            # Plot tree with feature loadings
            p <- .loadings_tree_plotter(mat, tree, row.var, ...)
        } else {
            # Utilize matrix method to create a plot
            p <- plotLoadings(
                mat, layout = layout, ncomponents = ncomponents, ...) 
        }
    return(p)
    }
)

#' @rdname plotLoadings
#' @export 
setMethod("plotLoadings", signature = c(x = "SingleCellExperiment"),
    function(x, dimred, layout = "barplot", ncomponents = 5, ...){
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
            mat, layout = layout, ncomponents = ncomponents, ...) 
        return(p)
    }
)

#' @rdname plotLoadings
#' @export 
setMethod("plotLoadings", signature = c(x = "matrix"),
    function(x, layout = "barplot", ncomponents = 5, ...) {
        # Input check
        .check_loadings_matrix(
            x, layout = layout, ncomponents = ncomponents, ...)
        #
        # Get data for plotting
        df <- .get_loadings_plot_data(x, layout, ncomponents, ...)
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
    # Ensure that the type is matrix
    mat <- as.matrix(mat)
    return(mat)
}

# This functions checks that loadings matrix is correct
.check_loadings_matrix <- function(mat, layout, ncomponents, n = 10, ...) {
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
#' @importFrom tibble rownames_to_column 
#' @importFrom tidyr pivot_longer
.get_loadings_plot_data <- function(df, layout, ncomponents, n = 10, ...) {
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
    # Add column that shows the values in absolute scale, and another column
    # showing sign
    res[["Value_abs"]] <- abs(res[["Value"]])
    res[["Sign"]] <- ifelse(
        res[["Value"]] > 0, "+", ifelse(res[["Value"]] < 0, "-", ""))
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
.plot_loadings <- function(df, layout, absolute.scale = TRUE, ...) {
    #
    if( !.is_a_bool(absolute.scale) ){
        stop("'absolute.scale' must be TRUE or FALSE.", call. = FALSE)
    }
    #
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
            
    } else if( layout == "barplot" && !absolute.scale ){
        # This creates a barplot where values can be negative or positive
        # (bars can be in negative and positive side)
        plot_out <- plot_out +
            # Create a bar plot. Create unique facets for each PC. Each PC can
            # have unique set of features. To reorder features by each facet,
            # we use reorder_within() and scale_y_reordered().
            geom_bar(
                mapping = aes(
                    x = Value, y = reorder_within(Feature, Value, PC)),
                stat = "identity"
                ) +
            scale_y_reordered() +
            facet_wrap(~ PC, scales = "free") +
            labs(x = "Value", y = "Feature") 
        
    } else if( layout == "barplot" && absolute.scale ){
        # This creates a barplot where bars are in absolute scale and the sing
        # is denoted with +/-
        plot_out <- plot_out +
            # Create bars with absolute scale
            geom_bar(
                mapping = aes(
                    x = Value_abs, y = reorder_within(Feature, -Value_abs, PC)),
                stat = "identity"
            ) +
            # Add sign that tells whether the value is + or -
            geom_text(aes(
                x = max(Value_abs) + max(Value_abs)*0.1,
                y = reorder_within(Feature, Value_abs, PC),
                label = Sign,
                fontface = "bold"
                )) +
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
.get_loadings_tree_data <- function(df, x, tree.name, row.var){
    # Check that rownames of loading matrix match with rownames of TreeSE. It
    # might be that TreeSE is updated after calculating the reduced dimension
    # which is why rownames do not match.
    all_match <- all(rownames(x) %in% rownames(df)) &&
        all(rownames(df) %in% rownames(x))
    if( !all_match ){
        stop(
            "Features of loading matrix do not match with rownames(x)",
            call. = FALSE)
    }
    # Sort the loading matrix
    df <- df[match(rownames(x), rownames(df)), ]
    # Convert loadings matrix to data.frame
    df <- as.data.frame(df)
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
    if( !is.null(row.var) ){
        df[["Feature"]] <- rowData(x)[[row.var]]
    }
    # Check that there are not too many features to plot. If there are too many
    # rank values (or rownames), it is not possible to plot.
    if( length(unique(df[["Feature"]])) > 100 ){
        stop("Too many features to plot.", call. = FALSE)
    }
    # Add rowlinks to data
    rownames(df) <- rowLinks(x)[["nodeLab"]]
    res <- list(tree = phylo, loadings = df)
    return(res)
}

# This function is for plotting tree with heatmap. It utilizes ggtree package.
#' @importFrom ggtree ggtree gheatmap
#' @importFrom ggnewscale new_scale_fill
#' @importFrom ggplot2 scale_fill_gradient2 scale_fill_viridis_d
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
