#' Plot feature loadings after performing a reducedDim 
#'
#' Inspired by the \code{\link[diffTop:plotASVcircular]{plotASVcircular}} method using phyloseq 
#' and has been converted to use TreeSummarizedExperiment objects.
#'
#' @param x a
#'   \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'   x.
#' 
#' @param dimred One dimensionality reduction method of \code{c("PCA", "LDA", "NMF")}
#'   indicating which method will be used.
#'   (default: \code{dimred = "PCA"})
#'  
#' @param layout One way to plot feature loadings of \code{c("heatmap", "barplot", "screeplot", "tree")} 
#'   (default: \code{layout = "heatmap"})
#' 
#' @param n A numeric specifying the number of features to be plotted 
#'   (default: \code{n = 10})
#'   
#' @details
#' # TO DO 
#' @return 
#' A \code{ggplot2} object. A circular plot annotated with TreeSummarizedExperiment object.
#'
#'
#' @examples
#' 
#' 
#' # Plotting feature loadings with tree
#' library(mia)
#' library(ggtree)
#' library(scater)
#' data("GlobalPatterns", package = "mia")
#' tse <- GlobalPatterns
#' tse <- transformAssay(tse, method = "clr", pseudocount = 1)
#' tse <- agglomerateByPrevalence(tse, rank="Phylum", update.tree = TRUE)
#' tse <- runPCA(tse, name = "PCA", ncomponents = 2, assay.type = "clr")
#' plotLoadings(tse, dimred= "PCA", layout = "tree")
#' 
#' 
#' # Plotting without tree as a heatmap
#' loadings_matrix <- attr(reducedDim(tse, "PCA"), "rotation")
#' plotLoadings(loadings_matrix, dimred= "PCA", layout = "heatmap")
#' 
#' # Plotting without tree as a barplot
#' plotLoadings(loadings_matrix, dimred= "PCA", layout = "barplot")
#' 
#' # Plotting without tree as a screeplot
#' plotLoadings(loadings_matrix, dimred= "PCA", layout = "screeplot")
#' 
#' #Plotting more features
#' plotLoadings(loadings_matrix, dimred= "PCA", layout = "heatmap", n = 12)
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
            ...) {
            
        # Making sure there is no error in parameters given by the user
        .check_parameters(x,
                        dimred = dimred,
                        layout = layout,
                        ...)   
        
        if (dimred == "PCA") {
            # Retrieving loadings matrix
            loadings_matrix <- attr(reducedDim(x, "PCA"), "rotation")
            # Ordering loadings and adding factor to keep the order
            df <- .get_loadings_plot_data(loadings_matrix, n)
        } else {
            stop("These methods are yet to be implemented", call. = FALSE)
        }
        
        if (layout == "tree") {
            # Plot tree with feature loadings
            p <- .loadings_tree_plotter(x, loadings_matrix)
        } else {
            # Plot features with the layout selected
            df1 <- data.frame(df[1])
            df2 <- data.frame(df[2])
            p <- .plot_pca_feature_loadings(df1, df2, layout)
        }
    return(p)
    }
)

setMethod("plotLoadings", signature = c(x = "matrix"),
    function(x,
            dimred = "PCA",
            layout = "heatmap",
            n = 10,
            ...) {
        # Making sure there is no error in parameters given by the user
        .check_parameters(x,
                        dimred = dimred,
                        layout = layout,
                        ...)
                      
        if (dimred == "PCA") {
            # Ordering loadings and adding factor to keep the order
            df <- .get_loadings_plot_data(x, n)
            df1 <- data.frame(df[1])
            df2 <- data.frame(df[2])
            # Plot features with the layout selected
            p <- .plot_pca_feature_loadings(df1, df2, layout)
        } else {
            stop("These methods are yet to be implemented", call. = FALSE)
        }
    return(p)
    }
)

.check_parameters <- function(x, dimred, layout, ...) {
    # Checking if dimred is correct
    if ( !(dimred %in% c("PCA", "LDA", "NMF")) ) {
        stop("'dimred' must be one of c('PCA', 'LDA', 'NMF').", call. = FALSE)
    }
    # Checking if layout is correct
    if ( !(layout %in% c("screeplot", "barplot", "tree", "heatmap")) ) {
        stop("'dimred' must be one of c('screeplot', 'barplot', 'tree', 'heatmap').", call. = FALSE)
    }
    # Making sure the user doesn't try to plot the tree if he gives only the matrix
    if (is.matrix(x) && layout == "tree") {
        stop("TreeSummarizedExperiment object is required for the tree plotting.", call. = FALSE)
    }
    # Making sure the tree is not null
    if (is(x, "TreeSummarizedExperiment")) {
        if ( is.null(rowTree(x)) && layout == "tree") {
            stop ("Tree is null.", call. = FALSE)
        }
    }
}

.get_loadings_plot_data <- function(x, n) {
    df1 <- data.frame(x)
    df2 <- data.frame(x)
    df1[["Feature"]] <- rownames(x)
    df2[["Feature"]] <- rownames(x)
    # Ordering loadings
    df1 <- df1[ order(abs(df1[["PC1"]])[1:n], decreasing = TRUE), ]
    df1 <- df1[ order(df1[["PC1"]]), ]
    df2 <- df2[ order(abs(df2[["PC2"]])[1:n], decreasing = TRUE), ]
    df2 <- df2[ order(df2[["PC2"]]), ]
    # Add factor to keep order
    df1[["Feature"]] <- factor(df1[["Feature"]], levels = df1[["Feature"]])
    df2[["Feature"]] <- factor(df2[["Feature"]], levels = df2[["Feature"]])
    list(df1, df2)
}

.loadings_tree_plotter <- function(x, loadings_matrix) {
    phylo <- rowTree(x)
    circ <- ggtree(phylo, layout = "circular")
    df <- rowData(x)
    color <- randomcoloR::distinctColorPalette(
      length(
        unique(
          df$Phylum
        )
      )
    )
    df <- data.frame(Class = df$Phylum)
    rownames(df) <- phylo$tip.label
    
    df2 <- data.frame(abs(loadings_matrix))
    rownames(df2) <- phylo$tip.label
    
    
    p <- gheatmap(
        p = circ,
        data = df,
        offset = -.1,
        width = .1,
        colnames_angle = 95,
        colnames_offset_y = .5,
        font.size = 4,
        color = "black") +
        ggplot2::scale_fill_manual(
            values = color,
            name = "Class"
        )
    
        for(i in 1:2){
            if(i == 1){
                p <- p +
                ggnewscale::new_scale_fill()
            }
            df3 <- dplyr::select(
                df2, (all_of(i))
            )
        
            p <- gheatmap(
                p,
                df3,
                offset = i*.065,
                width = .1,
                colnames_angle = 90,
                font.size = 4,
                high = "darkslateblue",
                low = "gray98",
                color = "black",
                legend_title = expression(beta[k])
            )
        }
    return(p)
}

.plot_pca_feature_loadings <- function(df1, df2, layout) {

    if (layout == "heatmap") {
        df1 <- data.frame(PC1 = round(df1$PC1,2), Feature = df1$Feature)
        df2 <- data.frame(PC2 = round(df2$PC2,2), Feature = df2$Feature)
    
        p1 <- ggplot(df1, aes(x = "PC1", y = Feature, label = PC1))  +
            geom_point(aes(fill = PC1), size=15, shape = 22) +
            theme(axis.title.x = element_blank()) +
            scale_fill_gradient2(limits = c(-1,1), low = "darkslateblue",
                mid = "white", high = "darkred", guide = NULL) +
            geom_text(color="black", size=4) +
            labs(title="PC1")
      
        p2 <- ggplot(df2, aes(x = "PC2", y = Feature, label = PC2))  +
            geom_point(aes(fill = PC2), size=15, shape = 22) + 
            theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
            scale_fill_gradient2("Loadings", limits=c(-1,1), low = "darkslateblue",
                mid = "white", high = "darkred") +
            geom_text(color="black", size=4) +
            labs(title="PC2")
      
        return(p1 + p2)
    }
    
    else if (layout == "barplot") {
    
        p1 <- ggplot(df1, aes(x = .data[["PC1"]], .data[["Feature"]])) +
            geom_bar(stat="identity") + xlim(-1,1)
        p2 <- ggplot(df2, aes(x = .data[["PC2"]], .data[["Feature"]])) +
            geom_bar(stat="identity") + xlim(-1,1)
        return(p1 + p2)
    }
    
    else if (layout == "screeplot") {
        p1 <- ggplot(df1, aes(x = .data[["PC1"]], .data[["Feature"]])) +
            geom_bar(stat="identity") + xlim(-1,1) + coord_flip()
        p2 <- ggplot(df2, aes(x = .data[["PC2"]], .data[["Feature"]])) +
            geom_bar(stat="identity") + xlim(-1,1) + coord_flip()
        return(p1 + p2)
    }
}
