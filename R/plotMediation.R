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
#' @param name \code{Character scalar} value defining which mediation data
#' to use. (Default: \code{"mediation"})
#'
#' @param layout \code{character scalar} Determines the layout of plot. Must be
#' either "heatmap" or "forest". (Default: \code{"heatmap"})
#'
#' @param signif.threshold \code{Numeric scalar or list}. Displays significance
#'   as stars on a heatmap layout. (Default: \code{c(0.001, 0.01, 0.05)})
#'
#' @details
#' \code{plotMediation} creates a heatmap starting from the
#' output of the \code{\link[mia:getMediation]{mediation}} functions, which are
#' mia wrappers for the basic \code{\link[mediation:mediate]{mediate}} function.
#' Either a \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or a data.frame object is supported as input. When the input is a
#' SummarizedExperiment, this should contain the output of addMediation
#' in the metadata slot and the argument \code{name} needs to be defined.
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
#'
#' # Visualise results as heatmap with custom significance thresholds
#' plotMediation(tse, "assay_mediation", signif.threshold = c(0.01, 0.05))
#' 
#' # Visualise results as forest plot
#' plotMediation(tse, "assay_mediation", layout = "forest")
#'
#'@name plotMediation

#' @rdname plotMediation
#' @export
#' @importFrom tidyr pivot_longer pivot_wider
setMethod("plotMediation", signature = c(x = "data.frame"),

    function(x, layout = "heatmap", signif.threshold = c(0.001, 0.01, 0.05)) {
      
        df <- x %>%
            pivot_longer(cols = -Mediator, names_to = c("Condition", "Metric"),
                names_sep = "_") %>%
            pivot_wider(names_from = Metric, values_from = value)
      
        df <- .add_signif_threshold(df, signif.threshold)
        
        if( layout == "heatmap" ){
            p <- .plot_med_heatmap(df)
        }else if( layout == "forest" ){
            p <- .plot_med_forest(df)
        }
      
        return(p)
    }
)

#' @rdname plotMediation
#' @export
setMethod("plotMediation", signature = c(x = "SummarizedExperiment"),
          
    function(x, name = "mediation", layout = "heatmap",
        signif.threshold = c(0.001, 0.01, 0.05)){
    
        med_df <- metadata(x)[[name]]
        p <- plotMediation(med_df, layout, signif.threshold)
        
        return(p)
    }
)

.add_signif_threshold <- function(df, signif_threshold) {
  
  df[["pval_symbol"]] <- ""
  
  for( i in seq_len(length(signif_threshold)) ){
      alpha <- rev(signif_threshold)[[i]]
      df[df[["pval"]] < alpha, "pval_symbol"] <- strrep("*", i)
  }
  
  return(df)
}

#' @importFrom ggplot2 ggplot
.plot_med_forest <- function(df) {
    
    df[["Condition"]] <- factor(df[["Condition"]],
        levels = rev(unique(df[["Condition"]])))
  
    p <- ggplot(df, aes(x = estimate, y = Condition)) +
        geom_point(size = 3) +
        geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0) +
        geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
        facet_wrap(. ~ Mediator) +
        labs(x = "Estimate") +
        theme_bw() +
        theme(axis.text.y = element_text(size = 12),
            axis.title.y = element_blank())
    
    return(p)
}

#' @importFrom ggplot2 ggplot
.plot_med_heatmap <- function(df) {
  
    df[["Mediator"]] <- factor(df[["Mediator"]],
        levels = rev(unique(df[["Mediator"]])))
  
    effect_max <- max(abs(min(df[["estimate"]])), abs(max(df[["estimate"]])))
    effect_floor <- floor(10 * effect_max)

    p <- ggplot(df, aes(x = Condition, y = Mediator, fill = estimate)) +
        geom_tile(color = "white", lwd = 1.5, linetype = 1) +
        scale_fill_gradientn(colours = c("blue", "white", "red"),
            breaks = seq(-effect_floor, effect_floor) / 10,
            limits = c(-effect_max, effect_max), name = "Effect") +
        scale_x_discrete(position = "top") +
        geom_text(aes(label = pval_symbol)) +
        theme_minimal() +
        theme(axis.text.x = element_text(size = 12),
            axis.title = element_blank())
    
    return(p)
}