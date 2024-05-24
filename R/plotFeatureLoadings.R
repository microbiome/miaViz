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
            library(miaViz)
            library(scater)
            data("GlobalPatterns", package = "mia")
            tse <- GlobalPatterns
            tse <- getUniqueFeatures(tse, rank = "Phylum")
            tse <- logNormCounts(tse)
            tse <- runPCA(tse, name = "PCA", ncomponents = 2, ntop = 10)
            loadings_matrix <- attr(reducedDim(tse,"PCA"), "rotation")
            phylo <- rowTree(tse)
            circ <- ggtree(phylo, layout = "circular")
            df <- rowData(tse)
            rownames(df) <- phylo$tip.label

            color <- randomcoloR::distinctColorPalette(
              length(
                unique(
                  df$Class
                )
              )
            )
            df <- data.frame(Class = df$Class)
            p <- gheatmap(
              p = circ,
              data = df,
              offset = -.1,
              width = .1,
              colnames_angle = 95,
              colnames_offset_y = .5,
              font.size = 5) +
              ggplot2::scale_fill_manual(
                values = color,
                name = "Class"
              )
            
            for(i in 1:2){
              if(i == 1){
                p <- p +
                  ggnewscale::new_scale_fill()
              }
              df2 <- dplyr::select(
                data.frame(loadings_matrix), (i)
              )
              p <- gheatmap(
                p,
                df2,
                offset = i*.08,
                width = .1,
                colnames_angle = 90,
                colnames_offset_y = .25,
                font.size = 6,
                high = "dodgerblue",
                low = "gray98")
            }
            p <- p +
              theme_minimal(
                base_size = 3
              ) +
              theme(
                plot.title = element_text(hjust = 0.5)
              )
           
            
            
}