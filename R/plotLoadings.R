#' Plot feature loadings after performing a reducedDim 
#'
#' Inspired by the \code{\link[diffTop:plotASVcircular]{plotASVcircular}} method using phyloseq 
#' and has been converted to use TreeSummarizedExperiment objects.
#'
#' @param x a
#'   \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'   x.
#' 
#' @param dimred One dimensionality reduction method of reducedDimNames
#'   indicating which method will be used.
#'   (default: \code{dimred = "PCA"})
#'  
#' @param layout One way to plot feature loadings of \code{c("heatmap", "barplot", "screeplot", "tree")} 
#'   (default: \code{layout = "heatmap"})
#' 
#' @param n A numeric specifying the number of features to be plotted.
#'   (default: \code{n = 10})
#'   
#' @param ncomponents A numeric specifying the number of components.
#'   (default: \code{ncomponents = 5})
#' 
#' @param tree.name a single \code{character} value specifying a rowTree/colTree from
#'   \code{object}. (By default: \code{tree.name = "phylo"})
#'     
#' 
#' @details
#' It is impossible to plot tree if only the matrix is given. Number of features must be reduced
#' before calling function or it will not be understandable. 
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
#' @examples
#' 
#' 
#' # Plotting feature loadings with tree
# library(mia)
# library(ggtree)
# library(scater)
# data("GlobalPatterns", package = "mia")
# tse <- GlobalPatterns
# tse <- transformAssay(tse, method = "clr", pseudocount = 1)
# tse <- agglomerateByPrevalence(tse, rank="Phylum", update.tree = TRUE)
# tse <- runPCA(tse, ncomponents = 5, assay.type = "clr")
# plotLoadings(tse, layout = "tree")
#' 
#' # Plotting without tree as a heatmap
#' loadings_matrix <- attr(reducedDim(tse, "PCA"), "rotation")
#' plotLoadings(loadings_matrix, layout = "heatmap")
#' 
#' # Plotting without tree as a barplot
#' plotLoadings(loadings_matrix, layout = "barplot")
#' 
#' # Plotting without tree as a screeplot
#' plotLoadings(loadings_matrix, layout = "screeplot")
#' 
#' # Plotting more features
#' plotLoadings(loadings_matrix, layout = "heatmap", n = 12)
#'
#' # Plotting with less components
#' tse <- runPCA(tse, ncomponents = 4, assay.type = "clr")
#' loadings_matrix <- attr(reducedDim(tse, "PCA"), "rotation")
#' plotLoadings(loadings_matrix, layout = "heatmap", ncomponents = 4 )
#' 
#' # Plotting if loadings matrix name has been changed
#' tse <- runPCA(tse, name = "myPCAmatrix", ncomponents = 5, assay.type = "clr")
#' plotLoadings(tse, dimred= "myPCAmatrix", layout = "heatmap")
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
            layout = "heatmap",
            n = 10,
            ncomponents = 5,
            tree.name = "phylo",
            ...) {
        
            
        # Making sure there is no error in parameters given by the user
        .check_parameters(x,
                        dimred = dimred,
                        layout = layout,
                        n = n,
                        ncomponents = ncomponents,
                        tree.name = tree.name,
                        ...)  

        loadings_matrix <- attr(reducedDim(x, dimred), "rotation")
        #Checking if there are enough components in the matrix
        .check_components(loadings_matrix, ncomponents)
        # Ordering loadings and adding factor to keep the order
        L <- .get_loadings_plot_data(loadings_matrix, n, ncomponents)
        
        if (layout == "tree") {
            # Plot tree with feature loadings
            p <- .loadings_tree_plotter(x, loadings_matrix, ncomponents, tree.name)
        } else {
            # Plot features with the layout selected
            p <- .plot_pca_feature_loadings(L, layout, ncomponents)
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
            tree.name = "phylo",
            ...) {
        # Making sure there is no error in parameters given by the user
        .check_parameters(x,
                        dimred = dimred,
                        layout = layout,
                        n = n,
                        ncomponents = ncomponents,
                        tree.name = tree.name,
                        ...)
                      
        #Checking if there are enough components in the matrix
        .check_components(x, ncomponents)
        # Ordering loadings and adding factor to keep the order
        df <- .get_loadings_plot_data(x, n, ncomponents)
        # Plot features with the layout selected
        p <- .plot_pca_feature_loadings(df, layout, ncomponents)
        
        return(p)
    }
)

.check_parameters <- function(x, dimred, layout, n, ncomponents, tree.name,...) {
    # Check tree.name
    if( is(x, "TreeSummarizedExperiment") && !(tree.name %in% rowTreeNames(x)  && .is_a_string(tree.name))){
        stop("'tree.name' must be a single character value specifying a colTree.", call. = FALSE)
    }
    # Checking if dimred is correct
    if( is(x, "TreeSummarizedExperiment") && !(dimred %in% reducedDimNames(x)  && .is_a_string(dimred))){
        stop("'dimred' must specify reducedDim.", call. = FALSE)
    }
    # Checking if layout is correct
    if ( !(layout %in% c("screeplot", "barplot", "tree", "heatmap") && .is_a_string(layout)) ) {
        stop("'layout' must be one of c('screeplot', 'barplot', 'tree', 'heatmap').", call. = FALSE)
    }
    # Making sure the user doesn't try to plot the tree if he gives only the matrix
    if (is.matrix(x) && layout == "tree") {
        stop("TreeSummarizedExperiment object is required for the tree plotting.", call. = FALSE)
    }
    # Making sure the tree is not null
    if( layout == "tree" && !(is(x, "TreeSummarizedExperiment") && !is.null(rowTree(x, tree.name)) )) {
        stop ("Tree is null.", call. = FALSE)
    }
    #Checking if n is a positive number
    if ( !(is.numeric(n) && n > 0) ) {
        stop("'n' must be a positive number.", call. = FALSE)
    }
    #Checking if ncomponents is a positive number 
    if ( !(is.numeric(ncomponents) && ncomponents > 0) ) {
        stop("'ncomponents' must be a positive number.", call. = FALSE)
    }

}

.check_components <- function(x, ncomponents) {
    #Checking if there are enough components in the matrix
    if (ncomponents > length(colnames(x))) {
        stop("'ncomponents' must be lower or equal than number of components.", call. = FALSE)  
    }
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

.get_loadings_plot_data <- function(x, n, ncomponents) {
  # Transform into a dataframe
  df <- data.frame(x)
  # Keep only the number of components needed
  names <- colnames(df)[1:ncomponents]
  df <- dplyr::select(df, all_of(names))
  # Add feature labels
  df[["Feature"]] <- rownames(x)
  
  # Apply the function to each component and return the list
  L <- lapply(1:ncomponents, .process_component, df = df, n = n)
  return(L)
}


.loadings_tree_plotter <- function(x, loadings_matrix, ncomponents, tree.name) {
    # Retrieve rowTree
    phylo <- rowTree(x, tree.name)
    #Subset data based on the tree
    ind <- rowLinks(x)[["whichTree"]] == tree.name
    if( any(ind) ){
      warning("Data was subsetted")
    }
    x <- x[ind, ]
    # Store plot tree
    circ <- ggtree::ggtree(phylo, layout = "circular")
    df <- rowData(x)
    # Get distincts colors for legend
    color <- randomcoloR::distinctColorPalette(
      length(
        unique(
          df$Phylum
        )
      )
    )
    #Transform into a dataframe
    df <- data.frame(Class = df$Phylum)
    # Match labels
    rownames(df) <- phylo$tip.label
    
    #Transform into a dataframe
    df2 <- data.frame(abs(loadings_matrix))
    #Match labels
    rownames(df2) <- phylo$tip.label
    
    # Plot tree with first inner circle (Classes)
    p <- ggtree::gheatmap(
        p = circ,
        data = df,
        offset = -.1,
        width = .1,
        color = "black",
        colnames = FALSE,
        legend_title = "Class") + 
        scale_fill_manual(values = color,
            name = "Class")
        # Plot others circles (loadings)
        for(i in 1:ncomponents){
            if(i == 1){
                p <- p +
                ggnewscale::new_scale_fill()
            }
            df3 <- dplyr::select(
                df2, (all_of(i))
            )
            
            p <- ggtree::gheatmap(
                p,
                df3,
                offset = i*.065,
                width = .1,
                high = "darkslateblue",
                low = "gray98",
                color = "black",
                colnames = FALSE,
                legend_title = expression(beta[k])
            )
            p <- p + labs(title = "Tree feature loadings plot") +
              theme(legend.key.size = unit(0.5, 'cm'),
                    plot.title = element_text(size = 18, hjust=0.5))
        }
    return(p)
}

.plot_pca_feature_loadings <- function(L, layout, ncomponents) {
    if (layout == "heatmap") {
        # Transform into a dataframe and round numerics for plotting
        df <- data.frame(PC = round(L[[1]][[1]], 2), Feature = L[[1]][["Feature"]])
        # Plot first component
        p <- ggplot(df, aes(x = "PC", y = Feature, label = PC))  +
            geom_point(aes(fill = PC), size=15, shape = 22) +
            theme(axis.title.x = element_blank()) +
            scale_fill_gradient2(limits = c(-1,1), low = "darkslateblue",
                mid = "white", high = "darkred", guide = NULL) +
            geom_text(color="black", size=4) +
            labs(title="PC1")
        # Loop plotting others components
        for (i in 2:length(L)) {
            df <- data.frame(PC = round(L[[i]][[i]], 2),
                Feature = L[[i]][["Feature"]])
            # Not showing legend for each plot
            if (i != length(L)) {
                p <- p +
                    ggplot(df, aes(x = "PC", y = Feature, label = PC))  +
                    geom_point(aes(fill = PC), size=15, shape = 22) +
                    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
                    scale_fill_gradient2(limits = c(-1,1), low = "darkslateblue",
                        mid = "white", high = "darkred", guide = NULL) +
                    geom_text(color="black", size=4) +
                    labs(title=paste("PC",i,sep=""))
            # Only show legend for last one
            } else {
                p <- p +
                    ggplot(df, aes(x = "PC", y = Feature, label = PC))  +
                    geom_point(aes(fill = PC), size=15, shape = 22) +
                    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
                    scale_fill_gradient2(limits = c(-1,1), low = "darkslateblue",
                        mid = "white", high = "darkred") +
                    geom_text(color="black", size=4) +
                    labs(title=paste("PC",i,sep=""))
            }
        }
      
        return(p)
    }
    
    else if (layout == "barplot") {
        # Transform into a dataframe
        df <- data.frame(PC = L[[1]][[1]], Feature = L[[1]][["Feature"]]) 
        # Plot first component
        p <- ggplot(df, aes(x = .data[["PC"]], .data[["Feature"]])) +
            geom_bar(stat="identity") +
            theme(axis.title.x = element_blank()) +
            xlim(-1,1) +
            labs(title="PC1")
        # Loop to plot others components
        for (i in 2:length(L)) {
            df <- data.frame(PC =L[[i]][[i]],
                Feature = L[[i]][["Feature"]])
            p <- p + ggplot(df, aes(x = .data[["PC"]], .data[["Feature"]])) +
                geom_bar(stat="identity") +
                theme(axis.title.x = element_blank()) +
                xlim(-1,1) +
                labs(title=paste("PC",i,sep=""))
        }

        return(p)
    }
    
    else if (layout == "screeplot") {
        # Transform into a dataframe
        df <- data.frame(PC = L[[1]][[1]], Feature = L[[1]][["Feature"]]) 
        # Plot first component
        p <- ggplot(df, aes(x = .data[["PC"]], .data[["Feature"]])) +
            geom_bar(stat="identity") +
            theme(axis.title.x = element_blank()) +
            xlim(-1,1) +
            labs(title="PC1") +
            coord_flip()
        # Plot others components
        for (i in 2:length(L)) {
            df <- data.frame(PC =L[[i]][[i]],
                Feature = L[[i]][["Feature"]])
            p <- p + ggplot(df, aes(x = .data[["PC"]], .data[["Feature"]])) +
                geom_bar(stat="identity") +
                theme(axis.title.x = element_blank()) +
                xlim(-1,1) +
                labs(title=paste("PC",i,sep="")) +
                coord_flip()
      }
      
      return(p)
    }
}
