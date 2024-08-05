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
    function(x,
            dimred,
            layout = "barplot",
            n = 10,
            ncomponents = 5,
            tree.name = "phylo",
            rank = NULL,
            add.tree = FALSE,
            ...) {
      
        # Making sure there is no error in parameters given by the user
        .check_parameters(x,
                        dimred = dimred,
                        layout = layout,
                        n = n,
                        ncomponents = ncomponents,
                        tree.name = tree.name,
                        rank = rank,
                        add.tree = add.tree,
                        ...)
        loading_names <- c("rotation", "loadings")
        attr_names <- names(attributes(reducedDim(x, dimred)))
        attr_name <- attr_names[ attr_names %in% loading_names ]
        if( length(attr_name) != 1 ) {
            stop("Loadings cannot be found..", call. = FALSE)
        }
        loadings_matrix <- attr(reducedDim(x, dimred), attr_name)

        
        # Checking if there are enough components in the matrix
        .check_components(loadings_matrix, ncomponents)
        
        if (add.tree && layout == "heatmap") {
            # Plot tree with feature loadings
            p <- .loadings_tree_plotter(x, loadings_matrix, ncomponents, tree.name, rank)
        } else if (layout == "heatmap") {
            # Plot features with heatmap layout
            p <- ComplexHeatmap::Heatmap(loadings_matrix, heatmap_legend_param = list(title = "Value")) 
        } else {
            # Ordering loadings and adding factor to keep the order
            df <- .get_loadings_plot_data(loadings_matrix, n, ncomponents)
            # Plot features with barplot layout
            p <- .barplot_feature_loadings(df)
        }
    return(p)
    }
)
#' @rdname plotLoadings
#' @export 
setMethod("plotLoadings", signature = c(x = "matrix"),
    function(x,
            layout = "barplot",
            n = 10,
            ncomponents = 5,
            ...) {
      
        # Making sure there is no error in parameters given by the user
        .check_parameters(x,
                        layout = layout,
                        n = n,
                        ncomponents = ncomponents,
                        ...)
                      
        # Checking if there are enough components in the matrix
        .check_components(x, ncomponents)
        
        # Plot features with heatmap layout
        if (layout == "heatmap") {
            p <- ComplexHeatmap::Heatmap(x, heatmap_legend_param = list(title = "Value"))
        } else {
            # Ordering loadings and adding factor to keep the order
            df <- .get_loadings_plot_data(x, n, ncomponents)
            # Plot features with barplot layout
            p <- .barplot_feature_loadings(df)
        }
        return(p)
    }
)

#' @rdname plotLoadings
#' @export 
setMethod("plotLoadings", signature = c(x = "SingleCellExperiment"),
    function(x,
            dimred,
            layout = "barplot",
            n = 10,
            ncomponents = 5,
            ...) {
      
        # Making sure there is no error in parameters given by the user
        .check_parameters(x,
                          dimred = dimred,
                          layout = layout,
                          n = n,
                          ncomponents = ncomponents,
                          ...)  
        
        loading_names <- c("rotation", "loadings")
        attr_names <- names(attributes(reducedDim(x, dimred)))
        attr_name <- attr_names[ attr_names %in% loading_names ]
        if( length(attr_name) != 1 ) {
            stop("Loadings cannot be found..", call. = FALSE)
        }
        loadings_matrix <- attr(reducedDim(x, dimred), attr_name)
        # Checking if there are enough components in the matrix
        .check_components(loadings_matrix, ncomponents)
        
        #Plot features with heatmap layout
        if (layout == "heatmap") {
            p <- ComplexHeatmap::Heatmap(loadings_matrix, heatmap_legend_param = list(title = "Value"))
        } else {
            # Ordering loadings and adding factor to keep the order
            df <- .get_loadings_plot_data(loadings_matrix, n, ncomponents)
            # Plot features with barplot layout
            p <- .barplot_feature_loadings(df)
        }
        return(p)
    }
)

.check_parameters <- function(x, dimred, layout, n, ncomponents, tree.name, rank, add.tree, ...) {
    
    #Check add.tree
    if (is(x, "TreeSummarizedExperiment") && !.is_a_bool(add.tree)) {
        stop("'add.tree' must be either TRUE or FALSE", call. = FALSE)
    }
    # Check tree.name
    if(is(x, "TreeSummarizedExperiment") && add.tree && !(tree.name %in% rowTreeNames(x)  && .is_a_string(tree.name))){
        stop("'tree.name' must be a single character value specifying a rowTree.", call. = FALSE)
    }
    # Checking if dimred is correct
    if( is(x, "SingleCellExperiment") && is.null(dimred) && !(dimred %in% reducedDimNames(x)  && .is_a_string(dimred))){
        stop("'dimred' must specify reducedDim.", call. = FALSE)
    }
    # Checking if layout is correct
    if ( !(layout %in% c("barplot", "heatmap") && .is_a_string(layout)) ) {
        stop("'layout' must be one of c('barplot', 'heatmap').", call. = FALSE)
    }
    # Making sure the tree is not null
    if( is(x, "TreeSummarizedExperiment") && add.tree && is.null(rowTree(x, tree.name))) {
        stop ("Tree is null.", call. = FALSE)
    }
    # Checking if n is a positive number
    if ( !(is.numeric(n) && n > 0) ) {
        stop("'n' must be a positive number.", call. = FALSE)
    }
    # Checking if ncomponents is a positive number 
    if ( !(is.numeric(ncomponents) && ncomponents > 0) ) {
        stop("'ncomponents' must be a positive number.", call. = FALSE)
    }
    # Checking if rank is correct
    if ( is(x, "TreeSummarizedExperiment") && !(is.null(rank)) &&!(any(rank %in% taxonomyRanks(x)) )) {
        stop("'rank' must be one of taxonomyRanks", call. = FALSE)
    }
    return(NULL)

}

.check_components <- function(x, ncomponents) {
    # Checking if there are enough components in the matrix
    if (ncomponents > ncol(x)) {
        stop("'ncomponents' must be lower or equal than number of components.", call. = FALSE)  
    }
    return(NULL)
}

# Function to process each component
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

#' @importFrom dplyr select
.get_loadings_plot_data <- function(df, n, ncomponents) {
  # Transform into a dataframe
  df <- as.data.frame(df)
  # Keep only the number of components needed
  df <- df[ , seq_len(ncomponents), drop = FALSE]
  # Apply the function to each component and return the list
  res <- lapply(seq_len(ncomponents), .process_component, df = df, n = n)
  # Combine to single data.frame
  res <- do.call(rbind, res)
  res <- as.data.frame(res)
  return(res)
}

#' @importFrom dplyr select
#' @importFrom ggtree ggtree gheatmap
.loadings_tree_plotter <- function(x, loadings_matrix, ncomponents, tree.name, rank) {
    # Retrieve rowTree
    phylo <- rowTree(x, tree.name)
    # Subset data based on the tree
    ind <- rowLinks(x)[["whichTree"]] == tree.name
    if( any(ind) ){
        warning("message here", call. = FALSE)
    }
    x <- x[ind, ]
    # Store plot tree
    circ <- ggtree(phylo, layout = "circular")
    df <- rowData(x)
    # Transform into a dataframe
    if (is.null(rank)) {
        df <- data.frame(Class = rownames(x))
        legend_name <- "Class"
    } else {
        df <- data.frame(Class = df[[rank]])
        legend_name <- rank
    }
    # Match labels
    rownames(df) <- rowLinks(x)$nodeLab

    # Transform into a dataframe
    df2 <- data.frame(loadings_matrix)
    
    # Match labels
    df2 <- df2[match(rownames(x), rownames(df2)), ]
    rownames(df2) <- rowLinks(x)$nodeLab
    
    
    
    # Plot tree with first inner circle (Classes)
    p <- gheatmap(
        p = circ,
        data = df,
        offset = -.1,
        width = .1,
        color = "black",
        colnames = FALSE,
        legend_title = legend_name) 
        # Plot others circles (loadings)
        for(i in 1:ncomponents){
            if(i == 1){
                p <- p +
                ggnewscale::new_scale_fill()
            }
            df3 <- select(
                df2, (all_of(i))
            )
            
            p <- gheatmap(
                p,
                df3,
                offset = i*.065,
                width = .1,
                color = "black",
                colnames = FALSE) + 
                scale_fill_gradient2(limits = c(-1,1),
                    low = "darkslateblue", mid = "white",
                    high = "darkred", name = "Value")
            p <- p +  theme(legend.key.size = unit(0.5, 'cm'),
                    plot.title = element_text(size = 18, hjust=0.5))
        }
    return(p)
}

#' @importFrom tidytext scale_y_reordered reorder_within
.barplot_feature_loadings <- function(df) {
    cnames <- unique(df$PC)
    p <- ggplot(df, aes(x = Value, y = reorder_within(Feature, Value, PC))) +
        geom_bar(stat = "identity") +
        scale_y_reordered() +
        facet_wrap(~ PC, scales = "free") +
        theme_minimal() +
        labs(x = "Value", y = "Feature") 
    return(p)
}

