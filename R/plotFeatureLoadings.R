#' Plot tree with feature loadings after performing a reducedDim (e.g. PCA, PCoA, NMDS)
#'
#' Based on the \code{\link[diffTop:plotASVcircular]{plotASVcircular}} method using phyloseq 
#' and has been converted to use TreeSummarizedExperiment objects.
#'
#' @param object a
#'   \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'   object.
#' 
#' @param method One dimensionality reduction method of \code{c("PCA", "MDS", "PCoA", "NMDS")}
#'   indicating which method will be used.
#'   (default: \code{method = "PCA"})
#'  
#' @param rank a single character defining a taxonomic rank. Must be a value of
#'   \code{taxonomyRanks()} function.
#'   
#' @param heatmap a boolean indicating whether or not to plot heatmap 
#'   (default: \code{heatmap = FALSE})
#'   
#' @details
#' # TO DO 
#' @return 
#' A \code{ggplot2} object. A circular plot annotated with TreeSummarizedExperiment object.
#'
#' @name plotFeatureLoadings
#'
#' @examples
#' library(mia)
#' library(ggtree)
#' library(scater)
#' data("GlobalPatterns", package = "mia")
#' tse <- GlobalPatterns
#' tse <- logNormCounts(tse)
#' tse <- agglomerateByPrevalence(tse, rank="Phylum", prevalence=0.99, update.tree = TRUE)
#' plotFeatureLoadings(tse, method= "PCA")
NULL

#' @rdname plotFeatureLoadings
setGeneric("plotFeatureLoadings", signature = c("object"),
           function(object, ...)
             standardGeneric("plotFeatureLoadings"))


#' @rdname plotFeatureLoadings
#' @export
setMethod("plotFeatureLoadings", signature = c(object = "TreeSummarizedExperiment"),
          function(object,
                   method = "PCA", heatmap = FALSE,
                   ...){
            
            .check_parameters(object,
                              method = method,...)   
            
            object <- .run_reduction_method(object, method = method)
            
            if (method == "PCA") {
              loadings_matrix <- attr(reducedDim(object, "PCA"), "rotation")
            } else {
              stop("These methods are yet to be implemented", call. = FALSE)
            }
            
            .plot_feature_loadings(object, loadings_matrix, heatmap)
          })

.check_parameters <- function(object, method, ...) {
  
  if ( !(method %in% c("PCA", "MDS", "NMDS", "PCoA")) ) {
    stop("'method' must be one of c('PCA', 'MDS', 'NMDS', 'PCoA').", call. = FALSE)
  }
}

.run_reduction_method <- function(object, method) {
  if (method == "PCA") {
    object <- runPCA(object, name = "PCA", ncomponents = 2)
  } else {
    stop("These methods are yet to be implemented", call. = FALSE)
  }
  object
}

.plot_feature_loadings <- function(object, loadings_matrix, heatmap) {
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
  
  if (heatmap) {
    df2 <- data.frame(loadings_matrix, phylo$tip.label)
    
    ggcorrplot(abs(loadings_matrix), 
           type = "lower", 
           method="circle", 
           title="Feature loadings",
           legend.title = " "
           ggtheme=theme_bw) + 
           scale_fill_gradient2(
             breaks=c(0, 1), 
             limit=c(0, 1))
  } else {
  
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
  }
  
  
}
# library(mia)
# library(ggtree)
# library(scater)
# data("GlobalPatterns", package = "mia")
# tse <- GlobalPatterns
# tse <- logNormCounts(tse)
# tse <- agglomerateByPrevalence(tse, rank="Phylum", prevalence=0.99, update.tree = TRUE)
# plotFeatureLoadings(tse, method= "PCA")
#
#
###########################################################
# 
# 
# library(mia)
# library(ggtree)
# library(scater)
# data("GlobalPatterns", package = "mia")
# tse <- GlobalPatterns
# tse <- logNormCounts(tse)
# # Agglomerate to keep only high prevalence features
# tse <- agglomerateByPrevalence(tse, rank="Phylum", prevalence=0.99, update.tree = TRUE)
# # Achieve PCA reduction
# tse <- runPCA(tse, name = "PCA", ncomponents = 2)
# # Get the feature loadings
# loadings_matrix <- attr(reducedDim(tse, "PCA"), "rotation")
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
# df2 <- data.frame(abs(loadings_matrix))
# rownames(df2) <- phylo$tip.label
# 
# # Plot the tree with first inner circle (Classes)
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
# 
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
###########################################################
# 
# library(mia)
# library(ggtree)
# library(scater)
# library(ggcorrplot)
# data("GlobalPatterns", package = "mia")
# tse <- GlobalPatterns
# tse <- logNormCounts(tse)
# # Agglomerate to keep only high prevalence features
# tse <- agglomerateByPrevalence(tse, rank="Phylum", prevalence=0.99, update.tree = TRUE)
# # Achieve PCA reduction
# tse <- runPCA(tse, name = "PCA", ncomponents = 2)
# # Get the feature loadings
# loadings_matrix <- attr(reducedDim(tse, "PCA"), "rotation")
# phylo <- rowTree(tse)
# circ <- ggtree(phylo, layout = "circular")
# df <- rowData(tse)
# df <- data.frame(Class = df$Phylum)
# rownames(df) <- phylo$tip.label
# df2 <- data.frame(loadings_matrix, phylo$tip.label)
# 
# # Plot feature loadings
# ggcorrplot(abs(loadings_matrix),
#            type = "lower",
#            method="circle",
#            title="Feature loadings",
#            ggtheme=theme_bw) +
#   scale_fill_gradient2(breaks=c(0, 1), limit=c(0, 1))
#
#
#
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