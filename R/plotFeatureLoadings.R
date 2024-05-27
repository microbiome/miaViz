#' Plot tree with feature loadings after performing a reducedDim (e.g. PCA, PCoA, NMDS)
#'
#' Based on the \code{\link[diffTop:plotASVcircular]{plotASVcircular}} method using phyloseq 
#' and has been converted to use TreeSummarizedExperiment objects.
#'
#' @param object a
#'   \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'   object.
#' 
#' @param nfeatures An integer to specify the number of feature loadings to plot (i.e. number of topics).
#'   (default: \code{nfeatures = 2})
#' 
#' @param taxonomyLevel 
#' # ADD ALL PARAMS
#' 
#' @details
#' # TO DO 
#' @return 
#' A \code{ggplot2} object. A circular plot annotated with TreeSummarizedExperiment object.
#'
#' @name plotFeatureLoadings
#'
#' @examples
#' # TO DO
NULL

#' @rdname plotFeatureLoadings
setGeneric("plotFeatureLoadings", signature = c("object"),
           function(object, ...)
             standardGeneric("plotFeatureLoadings"))


#' @rdname plotTree
#' @export
setMethod("plotFeatureLoadings", signature = c(object = "TreeSummarizedExperiment"),
          function(object,
                   ...){
            library(mia)
            library(ggtree)
            library(scater)
            data("GlobalPatterns", package = "mia")
            tse <- GlobalPatterns
            tse <- logNormCounts(tse)
            tse <- tse[rowData(tse)$Phylum %in% c("Actinobacteria", "Caldiserica","Synergistetes","Crenarchaeota","Euryarchaeota"), ]
            tse <- agglomerateByRank(tse, rank = "Class", agglomerate.tree = TRUE)
            tse <- runPCA(tse, name = "PCA", ncomponents = 2, ntop = 100)
            loadings_matrix <- attr(reducedDim(tse,"PCA"), "rotation")
            phylo <- rowTree(tse)
            circ <- ggtree(phylo, layout = "circular")
            df <- rowData(tse)
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
                df2, (i)
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