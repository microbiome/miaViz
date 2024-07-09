#' Plot tree with feature loadings after performing a reducedDim (e.g. PCA, LDA, NMF)
#'
#' Based on the \code{\link[diffTop:plotASVcircular]{plotASVcircular}} method using phyloseq 
#' and has been converted to use TreeSummarizedExperiment objects.
#'
#' @param object a
#'   \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'   object.
#' 
#' @param method One dimensionality reduction method of \code{c("PCA", "LDA", "NMF")}
#'   indicating which method will be used.
#'   (default: \code{method = "PCA"})
#'  
#' @param rank a single character defining a taxonomic rank. Must be a value of
#'   \code{taxonomyRanks()} function.
#'   
#' @param plot One way to plot feature loadings of \code{c("Heatmap", "Barplot", "Screeplot", "Tree")} 
#'   (default: \code{plot = "Heatmap"})
#'   
#' @details
#' # TO DO 
#' @return 
#' A \code{ggplot2} object. A circular plot annotated with TreeSummarizedExperiment object.
#'
#' @name plotFeatureLoadings
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
#' tse <- logNormCounts(tse)
#' tse <- agglomerateByPrevalence(tse, rank="Phylum", prevalence=0.99, update.tree = TRUE)
#' tse <- runPCA(tse, name = "PCA", ncomponents = 2)
#' plotFeatureLoadings(tse, method= "PCA", plot = "Tree")
#' 
#' 
#' # Plotting without tree as a heatmap
#' library(mia)
#' library(ggtree)
#' library(scater)
#' data("GlobalPatterns", package = "mia")
#' tse <- GlobalPatterns
#' tse <- logNormCounts(tse)
#' tse <- agglomerateByPrevalence(tse, rank="Phylum", prevalence=0.99, update.tree = TRUE)
#' tse <- runPCA(tse, name = "PCA", ncomponents = 2)
#' plotFeatureLoadings(loadings_matrix, method= "PCA", plot = "Heatmap")
#' 
#' # Plotting without tree as a barplot
#' library(mia)
#' library(ggtree)
#' library(scater)
#' data("GlobalPatterns", package = "mia")
#' tse <- GlobalPatterns
#' tse <- logNormCounts(tse)
#' tse <- agglomerateByPrevalence(tse, rank="Phylum", prevalence=0.99, update.tree = TRUE)
#' tse <- runPCA(tse, name = "PCA", ncomponents = 2)
#' plotFeatureLoadings(loadings_matrix, method= "PCA", plot = "Barplot")
#' 
#' #' # Plotting without tree as a screeplot
#' library(mia)
#' library(ggtree)
#' library(scater)
#' data("GlobalPatterns", package = "mia")
#' tse <- GlobalPatterns
#' tse <- logNormCounts(tse)
#' tse <- agglomerateByPrevalence(tse, rank="Phylum", prevalence=0.99, update.tree = TRUE)
#' tse <- runPCA(tse, name = "PCA", ncomponents = 2)
#' plotFeatureLoadings(loadings_matrix, method= "PCA", plot = "Screeplot")
NULL

#' @rdname plotFeatureLoadings
setGeneric("plotFeatureLoadings", signature = c("object"),
           function(object, ...)
             standardGeneric("plotFeatureLoadings"))


#' @rdname plotFeatureLoadings
#' @export 
setMethod("plotFeatureLoadings", signature = c(object = "TreeSummarizedExperiment"),
          function(object,
                   method = "PCA", plot = "Heatmap",
                   ...){
            
            .check_parameters(object,
                              method = method, plot = plot, ...)   
            
            if (method == "PCA") {
              loadings_matrix <- attr(reducedDim(tse, "PCA"), "rotation")
              df <- data.frame(loadings_matrix)
              df[["Feature"]] <- rownames(loadings_matrix)
              df[["Feature"]] <- factor(df[["Feature"]], levels = df[["Feature"]])
              df <- df[ order(abs(df[["PC1"]]), decreasing = TRUE)[1:10], ]
              df <- df[ order(df[["PC1"]]), ]
              df <- df[ order(abs(df[["PC2"]]), decreasing = TRUE)[1:10], ]
              df <- df[ order(df[["PC2"]]), ]
            } else {
              stop("These methods are yet to be implemented", call. = FALSE)
            }
            
            if (plot == "Tree") {
              
              phylo <- rowTree(object)
              circ <- ggtree(phylo, layout = "circular")
              df <- rowData(object)
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
                  legend_title = expression(beta[k]))
              }
              
              p
              
            } else {
              .plot_pca_feature_loadings(df, plot)
            }
          }
)

setMethod("plotFeatureLoadings", signature = c(object = "matrix"),
          function(object,
                   method = "PCA", plot = "Heatmap",
                   ...){
            
            .check_parameters(object,
                              method = method, plot = plot, ...)
                        
            if (method == "PCA") {
              df <- data.frame(object)
              df[["Feature"]] <- rownames(object)
              df[["Feature"]] <- factor(df[["Feature"]], levels = df[["Feature"]])
              df <- df[ order(abs(df[["PC1"]]), decreasing = TRUE)[1:10], ]
              df <- df[ order(df[["PC1"]]), ]
              df <- df[ order(abs(df[["PC2"]]), decreasing = TRUE)[1:10], ]
              df <- df[ order(df[["PC2"]]), ]
              .plot_pca_feature_loadings(df, plot)
            } else {
              stop("These methods are yet to be implemented", call. = FALSE)
            }
            

          }
)

.check_parameters <- function(object, method, plot, ...) {
  
  if ( !(method %in% c("PCA", "LDA", "NMF")) ) {
    stop("'method' must be one of c('PCA', 'LDA', 'NMF').", call. = FALSE)
  }
  
  if ( !(plot %in% c("Screeplot", "Barplot", "Tree", "Heatmap")) ) {
    stop("'method' must be one of c('Screeplot', 'Barplot', 'Tree', 'Heatmap').", call. = FALSE)
  }
  
  if (is.matrix(object) && plot == "Tree") {
    stop("Tree cannot be plotted for a non-TSE object.", call. = FALSE)
  }
  
  if (is(object, "TreeSummarizedExperiment")) {
    if ( is.null(rowTree(object)) && plot == "Tree") {
      stop ("Tree is null.", call. = FALSE)
    }
  }
}

.plot_pca_feature_loadings <- function(df, plot) {
  
    if (plot == "Heatmap") {
        df1 <- data.frame(PC1 = round(df$PC1,2), Feature = df$Feature)
        df2 <- data.frame(PC2 = round(df$PC2,2), Feature = df$Feature)
    
      p1 <- ggplot(df1, aes(x = "PC1", y = Feature, label = PC1))  +
        geom_point(aes(fill = PC1), size=15, shape = 22) +
        scale_fill_gradient2(limits = c(-1,1)) +
        geom_text(color="black", size=4) +
        labs(title="Feature Loadings Plot")
      
      p2 <- ggplot(df2, aes(x = "PC2", y = Feature, label = PC2))  +
        geom_point(aes(fill = PC2), size=15, shape = 22) +
        scale_fill_gradient2(limits=c(-1,1)) +
        geom_text(color="black", size=4) +
        labs(title="Feature Loadings Plot")
      
      p1 + p2
      
  }
  
  else if (plot == "Barplot") {
    
    p1 <- ggplot(df, aes(x = .data[["PC1"]], .data[["Feature"]])) +
      geom_bar(stat="identity") + xlim(-1,1)
    p2 <- ggplot(df, aes(x = .data[["PC2"]], .data[["Feature"]])) +
      geom_bar(stat="identity") + xlim(-1,1)
    p1 + p2
  }
  
  else if (plot == "Screeplot") {
    p1 <- ggplot(df, aes(x = .data[["PC1"]], .data[["Feature"]])) +
      geom_bar(stat="identity") + xlim(-1,1) + coord_flip()
    p2 <- ggplot(df, aes(x = .data[["PC2"]], .data[["Feature"]])) +
      geom_bar(stat="identity") + xlim(-1,1) + coord_flip()
    p1 + p2
    
  }
  
  
}
##############################################
#
#
# library(NMF)
# library(mia)
# library(ggtree)
# data(GlobalPatterns)
# tse <- GlobalPatterns
# tse <- transformAssay(tse, method = "relabundance")
# # Agglomerate to keep only high prevalence features
# tse <- agglomerateByPrevalence(tse, rank="Phylum", prevalence=0.99, update.tree = TRUE)
# # Perform the NMF reduction method
# x <- t(assay(tse, "counts"))
# nmf2 <- nmf(x, 2)
# H <- nmf2@fit@H
# W <- nmf2@fit@W
# # Get the feature loadings
# feature_loadings <- t(H)
# phylo <- rowTree(tse)
# circ <- ggtree(phylo, layout = "circular")
# df <- rowData(tse)
# color <- randomcoloR::distinctColorPalette(
#   length(
#     unique(
#       df$Phylum
#     )
#   )
# )
# df <- data.frame(Class = df$Phylum)
# rownames(df) <- phylo$tip.label
# 
# df2 <- data.frame(feature_loadings)
# rownames(df2) <- phylo$tip.label
# # Plot the tree and the first inner circle (Classes)
# p <- gheatmap(
#   p = circ,
#   data = df,
#   offset = -.1,
#   width = .1,
#   colnames_angle = 95,
#   colnames_offset_y = .5,
#   font.size = 4,
#   color = "black") +
#   ggplot2::scale_fill_manual(
#     values = color,
#     name = "Class"
#   )
# # Plot the feature loadings in different circles
# for(i in 1:2){
#   if(i == 1){
#     p <- p +
#       ggnewscale::new_scale_fill()
#   }
#   df3 <- dplyr::select(
#     df2, (all_of(i))
#   )
#   p <- gheatmap(
#     p,
#     df3,
#     offset = i*.065,
#     width = .1,
#     colnames_angle = 90,
#     font.size = 4,
#     high = "darkslateblue",
#     low = "gray98",
#     color = "black",
#     legend_title = expression(beta[k]))
# }
# 
# p
# 
#  
# 
############################################################################# 
# 
# 
# library(mia)
# library(topicmodels)
# library(ggtree)
# 
# data(GlobalPatterns)
# tse <- GlobalPatterns
# 
# tse <- transformAssay(tse, method = "relabundance")
# # Agglomerate to keep only high prevalence features
# tse <- agglomerateByPrevalence(tse, rank = "Class", update.tree = TRUE, prevalence = 0.99)
# # Perform LDA reduction method
# lda_model <- LDA(t(assay(tse, "counts")), k = 2)
# 
# df <- as.data.frame(t(assay(tse, "counts")))
# 
# posteriors <- topicmodels::posterior(lda_model, df)
# # Get the feature loadings
# feature_loadings <- t(as.data.frame(posteriors$terms))
# 
# 
# phylo <- rowTree(tse)
# circ <- ggtree(phylo, layout = "circular")
# df <- rowData(tse)
# color <- randomcoloR::distinctColorPalette(
#   length(
#     unique(
#       df$Phylum
#     )
#   )
# )
# df <- data.frame(Class = df$Phylum)
# rownames(df) <- phylo$tip.label
# 
# df2 <- data.frame(feature_loadings)
# rownames(df2) <- phylo$tip.label
# # Plot the tree and first inner circles (Classes)
# p <- gheatmap(
#   p = circ,
#   data = df,
#   offset = -.1,
#   width = .1,
#   colnames_angle = 95,
#   colnames_offset_y = .5,
#   font.size = 4,
#   color = "black") +
#   ggplot2::scale_fill_manual(
#     values = color,
#     name = "Class"
#   )
# # Plot the feature loadings in different circles
# for(i in 1:2){
#   if(i == 1){
#     p <- p +
#       ggnewscale::new_scale_fill()
#   }
#   df3 <- dplyr::select(
#     df2, (all_of(i))
#   )
#   p <- gheatmap(
#     p,
#     df3,
#     offset = i*.065,
#     width = .1,
#     colnames_angle = 90,
#     font.size = 4,
#     high = "darkslateblue",
#     low = "gray98",
#     color = "black",
#     legend_title = expression(beta[k]))
# }
# 
# p