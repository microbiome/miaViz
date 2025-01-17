#' Visualize mediation
#'
#' \code{plotMediation()} generates a heatmap from the results of mediation
#' analysis produced with \code{mia:getMediation()} or \code{mia:addMediation()}.
#' It displays effect size and significance of the Actual Causal Mediation
#' Effect (ACME) and the Actual Direct Effect (ADE) for each mediator included
#' in the object \code{x}.
#'
#' @param x a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object or a data.frame, returned as output from
#' \code{\link[mia:getMediation]{addMediation}} or
#' \code{\link[mia:getMediation]{getMediation}}, respectively.
#'
#' @param med.name \code{Character scalar} value defining which mediation data
#' to use. (Default: \code{"mediation"})
#'
#' @param add.signif \code{Logical scalar}. Should the p-values in the plot be
#' displayed? (Default: \code{FALSE})
#'
#' @details
#' \code{plotMediation} creates a heatmap starting from the
#' output of the \code{\link[mia:getMediation]{mediation}} functions, which are
#' mia wrappers for the basic \code{\link[mediation:mediate]{mediate}} function.
#' Either a \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or a data.frame object is supported as input. When the input is a
#' SummarizedExperiment, this should contain the output of addMediation
#' in the metadata slot and the argument \code{med.name} needs to be defined.
#' When the input is a data.frame, this should be returned as output from
#' getMediation.
#'
#' @return
#' a \code{\link[ComplexHeatmap:Heatmap-class]{Heatmap}} object
#'
#' @examples
#' library(mia)
#' library(scater)
#'
#' # Load dataset
#' data(hitchip1006, package = "miaTime")
#' tse <- hitchip1006
#' 
#' # Agglomerate features by family (merely to speed up execution)
#' tse <- agglomerateByRank(tse, rank = "Phylum")
#' # Convert BMI variable to numeric
#' tse$bmi_group <- as.numeric(tse$bmi_group)
#'
#' # Apply clr transformation to counts assay
#' tse <- transformAssay(tse,
#'                       method = "clr",
#'                       pseudocount = 1)
#'
#' # Analyse mediated effect of nationality on BMI via clr-transformed features
#' # 100 permutations were done to speed up execution, but ~1000 are recommended     
#' tse <- addMediation(tse, name = "assay_mediation",
#'                     outcome = "bmi_group",
#'                     treatment = "nationality",
#'                     assay.type = "clr",
#'                     covariates = c("sex", "age"),
#'                     treat.value = "Scandinavia",
#'                     control.value = "CentralEurope",
#'                     boot = TRUE, sims = 100,
#'                     p.adj.method = "fdr")

#' plotMediation(tse, "assay_mediation")
#'
#'@name plotMediation

#' @rdname plotMediation
#' @export
setGeneric("plotMediation", function(x, ...) standardGeneric("plotMediation"))

#' @rdname plotMediation
#' @export
#' @importFrom ComplexHeatmap Heatmap
setMethod("plotMediation", signature = c(x = "data.frame"),

    function(x, add.signif = TRUE) {
      
        coef_mat <-  as.matrix(x[ , c("ACME_estimate", "ADE_estimate")])
      
        rownames(coef_mat) <- x[["Mediator"]]
        colnames(coef_mat) <- gsub("_estimate", "", colnames(coef_mat))
      
        if( add.signif ){
        
            p_mat <-  as.matrix(x[ , c("ACME_pval", "ADE_pval")])
        
            cell_fun <- function(j, i, k, y, w, h, fill) {
                if(p_mat[i, j] < 0.001) {
                    grid.text("***", k, y)
                } else if(p_mat[i, j] < 0.01) {
                    grid.text("**", k, y)
                } else if(p_mat[i, j] < 0.05) {
                    grid.text("*", k, y)
                }
            }
        
        } else {
            cell_fun <- NULL
        }
      
        p <- Heatmap(coef_mat, name = "Effect",
            cluster_rows = FALSE, cluster_columns = FALSE,
            row_names_side = "left", column_names_side = "top",
            column_names_rot = 0, rect_gp = gpar(col = "white", lwd = 2),
            cell_fun = cell_fun)
      
        return(p)
    }
)


#' @rdname plotMediation
#' @export
setMethod("plotMediation", signature = c(x = "SummarizedExperiment"),
          
    function(x, med.name = "mediation", add.signif = TRUE){
    
        med_df <- metadata(x)[[med.name]]
        p <- plotMediation(med_df, add.signif)
        
        return(p)
    }
)
