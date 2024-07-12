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
#' @param n A numeric specifying the number of features to be plotted.
#'   (default: \code{n = 10})
#'   
#' @param name A single \code{character} value specifying the name of the loadings matrix if it is
#'   different from default naming. (default: \code{name = NULL})
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
#' library(mia)
#' library(ggtree)
#' library(scater)
#' data("GlobalPatterns", package = "mia")
#' tse <- GlobalPatterns
#' tse <- transformAssay(tse, method = "clr", pseudocount = 1)
#' tse <- agglomerateByPrevalence(tse, rank="Phylum", update.tree = TRUE)
#' tse <- runPCA(tse, ncomponents = 2, assay.type = "clr")
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
#' # Plotting more features
#' plotLoadings(loadings_matrix, dimred= "PCA", layout = "heatmap", n = 12)
#'
#' # Plotting with more components
#' tse <- runPCA(tse, ncomponents = 4, assay.type = "clr")
#' loadings_matrix <- attr(reducedDim(tse, "PCA"), "rotation")
#' plotLoadings(loadings_matrix, dimred= "PCA", layout = "heatmap")
#' 
#' # Plotting if loadings matrix name has been changed
#' tse <- runPCA(tse, name = "myPCAmatrix", ncomponents = 2, assay.type = "clr")
#' plotLoadings(loadings_matrix, dimred= "PCA", layout = "heatmap", name = "myPCAmatrix")
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
            name = NULL,
            ...) {
        
            
        # Making sure there is no error in parameters given by the user
        .check_parameters(x,
                        dimred = dimred,
                        layout = layout,
                        n = n,
                        ...)   
        
        if (dimred == "PCA") {
            # Retrieving loadings matrix
            if (!is.null(name)) {
                loadings_matrix <- attr(reducedDim(x, name), "rotation")
            } else {
                loadings_matrix <- attr(reducedDim(x, "PCA"), "rotation")
            }
            # Set number of components
            ncomponents <- length(colnames(loadings_matrix))
            # Ordering loadings and adding factor to keep the order
            L <- .get_loadings_plot_data(loadings_matrix, n, ncomponents)
        } else {
            stop("These methods are yet to be implemented", call. = FALSE)
        }
        
        if (layout == "tree") {
            # Plot tree with feature loadings
            p <- .loadings_tree_plotter(x, loadings_matrix, ncomponents)
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
            name = NULL,
            ...) {
        # Making sure there is no error in parameters given by the user
        .check_parameters(x,
                        dimred = dimred,
                        layout = layout,
                        n = n,
                        ...)
                      
        if (dimred == "PCA") {
            # Set number of components
            ncomponents <- length(colnames(x))
            # Ordering loadings and adding factor to keep the order
            df <- .get_loadings_plot_data(x, n, ncomponents)
            # Plot features with the layout selected
            p <- .plot_pca_feature_loadings(df, layout, ncomponents)
        } else {
            stop("These methods are yet to be implemented", call. = FALSE)
        }
    return(p)
    }
)

.check_parameters <- function(x, dimred, layout, n,...) {
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
    #Checking if n is a positive number
    if ( !(is.numeric(n) && n > 0) ) {
        stop("'n' must be a positive number.", call. = FALSE)
    }

}

.get_loadings_plot_data <- function(x, n, ncomponents) {
    L <- list()
    for (i in 1:ncomponents) {
        df <- data.frame(x)
        df[["Feature"]] <- rownames(x)
        # Ordering loadings
        df <- df[ order(abs(df[[paste("PC",as.character(i),sep="")]])[1:n], decreasing = TRUE), ]
        df <- df[ order(df[[paste("PC",as.character(i),sep="")]]), ]
        # Add factor to keep order
        df[["Feature"]] <- factor(df[["Feature"]], levels = df[["Feature"]])
        L <- append(L,df)
    }
    return(L)
}

.loadings_tree_plotter <- function(x, loadings_matrix, ncomponents) {
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
    
        for(i in 1:ncomponents){
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

.plot_pca_feature_loadings <- function(L, layout, ncomponents) {
    if (layout == "heatmap") {
        # Transform into a dataframe and round numerics for plotting
        df <- data.frame(PC = round(L[[1]], 2), Feature = L[[ncomponents + 1]])
        # Plot first component
        p <- ggplot(df, aes(x = "PC", y = Feature, label = PC))  +
            geom_point(aes(fill = PC), size=15, shape = 22) +
            theme(axis.title.x = element_blank()) +
            scale_fill_gradient2(limits = c(-1,1), low = "darkslateblue",
                mid = "white", high = "darkred", guide = NULL) +
            geom_text(color="black", size=4) +
            labs(title="PC1")
        # Loop plotting others components
        for (i in 2:(length(L)%/%(ncomponents + 1))) {
            df <- data.frame(PC = round(L[[(ncomponents + 1)*(i-1) + i]], 2),
                Feature = L[[(ncomponents + 1)*i]])
            # Not showing legend for each plot
            if (i != (length(L)%/%(ncomponents + 1))) {
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
        df <- data.frame(PC = L[[1]], Feature = L[[ncomponents + 1]]) 
        # Plot first component
        p <- ggplot(df, aes(x = .data[["PC"]], .data[["Feature"]])) +
            geom_bar(stat="identity") +
            theme(axis.title.x = element_blank()) +
            xlim(-1,1) +
            labs(title="PC1")
        # Loop to plot others components
        for (i in 2:(length(L)%/%(ncomponents + 1))) {
            df <- data.frame(PC =L[[(ncomponents + 1)*(i-1) + i]],
                Feature = L[[(ncomponents + 1)*i]])
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
        df <- data.frame(PC = L[[1]], Feature = L[[ncomponents + 1]]) 
        # Plot first component
        p <- ggplot(df, aes(x = .data[["PC"]], .data[["Feature"]])) +
            geom_bar(stat="identity") +
            theme(axis.title.x = element_blank()) +
            xlim(-1,1) +
            labs(title="PC1") +
            coord_flip()
        # Plot others components
        for (i in 2:(length(L)%/%(ncomponents + 1))) {
            df <- data.frame(PC =L[[(ncomponents + 1)*(i-1) + i]],
                Feature = L[[(ncomponents + 1)*i]])
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
