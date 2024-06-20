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
#'   (default: \code{method = "PCA})
#'  
#' @param rank a single character defining a taxonomic rank. Must be a value of
#'   \code{taxonomyRanks()} function.
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
#' data("GlobalPatterns", package = "mia)
#' tse <- GlobalPatterns
#' tse <- logNormCounts(tse)
#' plotFeatureLoadings(tse[rowData(tse)$Phylum %in% c("Actinobacteria", "Caldiserica", "Synergistetes", "Crenarchaeota", "Euryarchaeota"),], method= "PCA")
NULL

#' @rdname plotFeatureLoadings
setGeneric("plotFeatureLoadings", signature = c("object"),
          function(object, ...)
             standardGeneric("plotFeatureLoadings"))


#' @rdname plotTree
#' @export
setMethod("plotFeatureLoadings", signature = c(object = "TreeSummarizedExperiment"),
    function(object,
            method = "PCA", rank = "Class",
            ...){

    .check_parameters(object,
                    method = method,...)    
      
           if (rank %in% taxonomyRanks(object)) {
             object <- agglomerateByRank(object, rank = "Class", agglomerate.tree = TRUE)
           }
           object <- .run_reduction_method(object, method = method)
           
           if (method == "PCA") {
              loadings_matrix <- attr(reducedDim(object, "PCA"), "rotation")
           } else {
             stop("These methods are yet to be implemented", call. = FALSE)
           }

           .plot_feature_loadings(object,loadings_matrix, Phylum)
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

.plot_feature_loadings <- function(object,loadings_matrix, class) {
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
    
    df2 <- data.frame(loadings_matrix)
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
            df2, all_of(i)
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
  
    p
}
# 
# 
# 
# 
# 
# 
# 
# 
# library(NMF)
# library(mia)
# data(GlobalPatterns)
# tse <- GlobalPatterns
# tse <- transformAssay(tse, method = "relabundance")
# tse <- agglomerateByPrevalence(tse, rank="Genus", assay.type="relabundance",detection = 1/100, prevalence=5/100)
# tse <- tse[!(rownames(tse) %in% c("Other")),]
# x <- t(assay(tse, "counts"))
# nmf2 <- nmf(x, 2)
# H <- nmf2@fit@H
# W <- nmf2@fit@W
# feature_loadings <- t(H)
# head(feature_loadings)
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
# 
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
# 
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
# 
# 
# 
# library(mia)
# library(topicmodels)
# 
# data(GlobalPatterns)
# tse <- GlobalPatterns
# 
# tse <- transformAssay(tse, method = "relabundance")
# 
# lda_model <- LDA(t(assay(tse, "counts")), k = 2)
# 
# df <- as.data.frame(t(assay(tse, "counts")))
# 
# posteriors <- topicmodels::posterior(lda_model, df)
# 
# feature_loadings <- t(as.data.frame(posteriors$terms))
# 
# top_loadings_indices <- order(rowSums(feature_loadings), decreasing = TRUE)[1:30]
# 
# feature_loadings_reduced <- feature_loadings[top_loadings_indices, ]
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
# df2 <- data.frame(feature_loadings_reduced)
# rownames(df2) <- phylo$tip.label
# 
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
# 
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
