#' Plot feature loadings for TreeSummarizedExperiment/SingleCellExperiment objects
#' or feature loadings numeric matrix.
#'
#' This function is used after performing a reduction method. If TSE object is
#' given it retrieves the feature loadings matrix to plot values with tree.
#' Plotting with other layouts is possible as SCE objects and numeric matrices
#' does not include a tree.
#' 
#' @param x a
#'   \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'   x.
#' 
#' @param dimred \code{Character scalar}. Dimension reduction from reducedDimNames to be plotted. 
#'   (default: \code{"PCA"})
#'  
#' @param layout One way to plot feature loadings of \code{c("heatmap", "barplot", "tree")} 
#'   (default: \code{"tree"} for TSE, \code{"heatmap"} for others)
#' 
#' @param n A numeric specifying the number of features to be plotted.
#'   (default: \code{10})
#'   
#' @param ncomponents A numeric specifying the number of components.
#'   (default: \code{5})
#' 
#' @param tree.name A single \code{character} value specifying a rowTree from a
#'   TreeSummarizedExperiment object. (By default: \code{"phylo"})
#'   
#' @param class A single \code{character} value specifying a rank from taxonomyRanks
#'   (default: \code{rownames})
#'   
#' @param ... additional arguments for plotting. See 
#'   \code{\link{mia-plot-args}} for more details i.e. call \code{help("mia-plot-args")}     
#' 
#' @details
#' 
#' Inspired by the \code{\link[diffTop:plotASVcircular]{plotASVcircular}} method using phyloseq 
#' and has been converted to use TreeSummarizedExperiment/SingleCellExperiment objects.
#' 
#' TreeSummarizedExperiment/SingleCellExperiment objects are expected to have content in reducedDim slot.
#' It is impossible to plot tree if only the matrix is given. Number of features must be reduced
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
#' plotLoadings(tse)
#' 
#' # Plotting without tree as a heatmap
#' loadings_matrix <- attr(reducedDim(tse, "PCA"), "rotation")
#' plotLoadings(loadings_matrix)
#' 
#' # Plotting without tree as a barplot
#' plotLoadings(loadings_matrix, layout = "barplot")
#' 
#' # Plotting more features
#' plotLoadings(loadings_matrix, n = 12)
#'
#' # Plotting with less components
#' tse <- runPCA(tse, ncomponents = 4, assay.type = "clr")
#' loadings_matrix <- attr(reducedDim(tse, "PCA"), "rotation")
#' plotLoadings(loadings_matrix, ncomponents = 4)
#' 
#' # Plotting if loadings matrix name has been changed
#' tse <- runPCA(tse, name = "myPCAmatrix", ncomponents = 5, assay.type = "clr")
#' plotLoadings(tse, dimred= "myPCAmatrix")
#' 
#' # Plotting tree with taxonomic rank classification
#' tse <- runPCA(tse, ncomponents = 5, assay.type = "clr")
#' plotLoadings(tse, class = "Phylum")
NULL

#' @rdname plotLoadings
setGeneric("plotLoadings", signature = c("x"),
    function(x, ...) 
        standardGeneric("plotLoadings"))


#' @rdname plotLoadings
#' @export 
setMethod("plotLoadings", signature = c(x = "TreeSummarizedExperiment"),
    function(x,
            dimred = "PCA",
            layout = "tree",
            n = 10,
            ncomponents = 5,
            tree.name = "phylo",
            class = rownames(x),
            ...) {
        
            
        # Making sure there is no error in parameters given by the user
        .check_parameters(x,
                        dimred = dimred,
                        layout = layout,
                        n = n,
                        ncomponents = ncomponents,
                        tree.name = tree.name,
                        class = class,
                        ...)
        loading_names <- c("rotation", "loadings")
        attr_names <- names(attributes(reduced_dim))
        attr_name <- attr_names[ attr_names %in% loading_names ]
        if( length(attr_name) != 1 )){
            stop("Loadings cannot be found..")
        }
        loadings_matrix <- attr(tse, attr_name)

        # Checking if there are enough components in the matrix
        .check_components(loadings_matrix, ncomponents)
        
        if (layout == "tree") {
            # Plot tree with feature loadings
            p <- .loadings_tree_plotter(x, loadings_matrix, ncomponents, tree.name, class)
        } else {
            # Ordering loadings and adding factor to keep the order
            L <- .get_loadings_plot_data(loadings_matrix, n, ncomponents)
            # Plot features with the layout selected
            p <- .plot_pca_feature_loadings(L, layout, n, ncomponents)
        }
    return(p)
    }
)
#' @rdname plotLoadings
#' @export 
setMethod("plotLoadings", signature = c(x = "matrix"),
    function(x,
            dimred = "PCA",
            layout = "heatmap",
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
                      
        # Checking if there are enough components in the matrix
        .check_components(x, ncomponents)
        # Ordering loadings and adding factor to keep the order
        df <- .get_loadings_plot_data(x, n, ncomponents)
        # Plot features with the layout selected
        p <- .plot_pca_feature_loadings(df, layout, n, ncomponents)
        
        return(p)
    }
)

#' @rdname plotLoadings
#' @export 
setMethod("plotLoadings", signature = c(x = "SingleCellExperiment"),
    function(x,
            dimred = "PCA",
            layout = "heatmap",
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
        
        loadings_matrix <- attr(reducedDim(x, dimred), "rotation")
        # Checking if there are enough components in the matrix
        .check_components(loadings_matrix, ncomponents)
            
        # Ordering loadings and adding factor to keep the order
        L <- .get_loadings_plot_data(loadings_matrix, n, ncomponents)
        # Plot features with the layout selected
        p <- .plot_pca_feature_loadings(L, layout, n, ncomponents)
        
        return(p)
    }
)

.check_parameters <- function(x, dimred, layout, n, ncomponents, tree.name, class, ...) {
    # Check tree.name
    if( is(x, "TreeSummarizedExperiment") && !(tree.name %in% rowTreeNames(x)  && .is_a_string(tree.name))){
        stop("'tree.name' must be a single character value specifying a colTree.", call. = FALSE)
    }
    # Checking if dimred is correct
    if( is(x, "TreeSummarizedExperiment") && !(dimred %in% reducedDimNames(x)  && .is_a_string(dimred))){
        stop("'dimred' must specify reducedDim.", call. = FALSE)
    }
    # Checking if layout is correct
    if ( !(layout %in% c("barplot", "tree", "heatmap") && .is_a_string(layout)) ) {
        stop("'layout' must be one of c('barplot', 'tree', 'heatmap').", call. = FALSE)
    }
    # Making sure the user doesn't try to plot the tree if he gives only the matrix
    if ( !(is(x, "TreeSummarizedExperiment")) && layout == "tree") {
        stop("TreeSummarizedExperiment object is required for the tree plotting.", call. = FALSE)
    }
    # Making sure the tree is not null
    if( layout == "tree" && !(is(x, "TreeSummarizedExperiment") && !is.null(rowTree(x, tree.name)) )) {
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
    # Checking if class is correct
    if ( is(x, "TreeSummarizedExperiment") && !(any(class %in% taxonomyRanks(x)) || identical(class, rownames(x)))) {
        stop("'class' must be one of taxonomyRanks", call. = FALSE)
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
    # Ordering loadings by absolute value
    df <- df[order(-abs(df[[i]])), ][1:n, ]
    # Ordering by actual loadings
    df <- df[order(df[[i]]), ]
    # Add factor to keep order
    df[["Feature"]] <- factor(df[["Feature"]], levels = df[["Feature"]])
    return(df)
}

#' @importFrom dplyr select
.get_loadings_plot_data <- function(x, n, ncomponents) {
  # Transform into a dataframe
  x <- as.data.frame(x)
  # Keep only the number of components needed
  names <- colnames(x)[1:ncomponents]
  df <- select(x, all_of(names))
  # Add feature labels
  x[["Feature"]] <- rownames(x)
  # Apply the function to each component and return the list
  L <- lapply(1:ncomponents, .process_component, df = x, n = n)
  return(L)
}

#' @importFrom dplyr select
#' @importFrom ggtree ggtree gheatmap
.loadings_tree_plotter <- function(x, loadings_matrix, ncomponents, tree.name, class) {
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
    if (identical(class, rownames(x))) {
        df <- data.frame(Class = class)
    } else {
        df <- data.frame(Class = df[[class]])
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
        legend_title = "Class") 
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
            p <- p + labs(title = "Tree feature loadings plot") +
              theme(legend.key.size = unit(0.5, 'cm'),
                    plot.title = element_text(size = 18, hjust=0.5))
        }
    return(p)
}

#' @importFrom dplyr arrange
.plot_pca_feature_loadings <- function(L, layout, n, ncomponents) {
    k <- seq_len(length(L))
    # Prepare the data in correct format
    df <- lapply(k, function(i){
      L[[i]][, i, drop = FALSE]
    })
    names <- unique(unlist(lapply(L, rownames)))
    cnames <- unique(unlist(lapply(L, colnames)))
    cnames <- cnames[!cnames %in% c("Feature") ]
    
    new_df <- data.frame(matrix(-1, nrow = length(names), ncol = length(cnames)), row.names = names)
    colnames(new_df) <- cnames
    for (i in seq_along(L)) {
      pc_df <- df[[i]]
      name <- colnames(pc_df)
      new_df[rownames(pc_df), name] <- pc_df[, 1]
    }
    new_df[["Feature"]] <- names
    
    
    # To long format
    df <- reshape(
      new_df, varying = colnames(new_df)[ !colnames(new_df) %in% c("Feature") ],
      v.names = "Value", 
      timevar = "PC",
      times = colnames(new_df)[ !colnames(new_df) %in% c("Feature") ],
      direction = "long")
    # Arrange dara
    df <- filter(df, Value != -1)
    if (layout == "heatmap") {
        plots <- lapply(1:ncomponents, function(i) {
            data_subset <- subset(df, PC %in% cnames[i])
            data_subset <- arrange(data_subset, Value)
            data_subset$Feature <- factor(data_subset$Feature, levels = data_subset$Feature)
            #Plot loadings
            ggplot(data_subset, aes(x = "PC", y = Feature, label = round(Value, 2)))  +
            geom_point(aes(fill = Value), size=15 - 0.4*n, shape = 22) +
            theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
                axis.text.x = element_blank()) +
            scale_fill_gradient2(limits = c(-1,1), low = "darkslateblue",
                mid = "white", high = "darkred") +
            geom_text(color="black", size=4 - 0.075 * n) +
            facet_wrap(~ PC)
      })
      p <- patchwork::wrap_plots(plots, guides = "collect")
      
    }
    
    else if (layout == "barplot") {
        plots <- lapply(1:ncomponents, function(i) {
            data_subset <- subset(df, PC %in% cnames[i])
            data_subset <- arrange(data_subset, Value)
            data_subset$Feature <- factor(data_subset$Feature, levels = data_subset$Feature)
            #Plot loadings
            ggplot(data_subset) +
                geom_bar(stat = "identity", aes(x = .data[["Value"]], y = .data[["Feature"]])) +
                xlim(-1,1) + facet_wrap(~ PC) +
                theme(axis.title.x = element_blank(), axis.title.y = element_blank())
        })
        p <- patchwork::wrap_plots(plots)
        
    }
    return(p)
}
