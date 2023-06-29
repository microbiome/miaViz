#' Plotting Dirichlet-Multinomial Mixture Model data
#'
#' To plot DMN fits generated with `mia` use \code{plotDMNFit}.
#'
#' @param x a 
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object contain the DMN data in \code{metadata}.
#'   
#' @param type the type of measure for access the goodness of fit. One of
#'   \sQuote{laplace}, \sQuote{AIC} or \sQuote{BIC}.
#'
#' @param name the name to store the result in
#'   \code{\link[SummarizedExperiment:RangedSummarizedExperiment-class]{metadata}}
#'
#' @param ... optional arguments not used.
#'
#' @return
#' \code{plotDMNFit} returns a \code{ggplot2} plot.
#'
#' @seealso
#' \code{\link[mia:calculateDMN]{calculateDMN}}
#'
#' @name plotDMN
#'
#' @examples
#' library("bluster")
#' data(dmn_se, package = "mia")
#' tse_dmm <- cluster(dmn_se, name = "DMM",
#'                    DmmParam(k = 1:4, type = "laplace"),
#'                    MARGIN = "samples", full = TRUE)
#'
#' # plot the fit
#' plotDMNFit(tse_dmm, name = "DMM", type = "laplace")
#' 
NULL

#' @rdname plotDMN
#' @export
setGeneric("plotDMNFit", signature = "x",
           function(x, name = "DMN", type = c("laplace","AIC","BIC"), ...)
               standardGeneric("plotDMNFit"))

#' @rdname plotDMN
#' @importFrom DirichletMultinomial mixture
#' @importFrom ggplot2 ggplot aes_string geom_point geom_line theme_bw labs
#' @export
setMethod("plotDMNFit", signature = c(x = "SummarizedExperiment"),
    function(x, name = "DMN", type = c("laplace","AIC","BIC")){
        #
        if (!is.null(metadata(x)[[name]]$dmm)) {
            dmn <- metadata(x)[[name]]$dmm
        } else {
            .Deprecated(old="getDMN", new="cluster", 
                    "Now runDMN and calculate DMN are deprecated. Use cluster with DMMParam parameter and full parameter set as true instead.")
            dmn <- metadata(x)$name
        }
        fit_FUN <- mia:::.get_dmn_fit_FUN(type)
        #
        k <- vapply(dmn, function(d){ncol(mixture(d))}, numeric(1))
        fit <- vapply(dmn, fit_FUN, numeric(1))
        ggplot(data.frame(k = k, fit = fit), aes_string(x = k, y = fit)) +
          geom_point() +
          geom_line() +
          theme_bw() +
          labs(x = "Number of Dirichlet Components",
               y = paste0("Model Fit (",type,")"))
    }
)
