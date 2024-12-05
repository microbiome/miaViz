#' Plotting Dirichlet-Multinomial Mixture Model data
#'
#' To plot DMN fits generated with `mia` use \code{plotDMNFit}.
#'
#' @param x a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object contain the DMN data in \code{metadata}.
#'
#' @param type \code{Character scalar}. The type of measure for access the
#' goodness of fit. One of \sQuote{laplace}, \sQuote{AIC} or \sQuote{BIC}.
#'
#' @param name \code{Character scalar}. The name to store the result in
#' \code{\link[SummarizedExperiment:RangedSummarizedExperiment-class]{metadata}}
#' (Default: \code{"DMN"})
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
#' library(mia)
#' library(bluster)
#'
#' # Get dataset
#' data("peerj13075", package = "mia")
#' tse <- peerj13075
#'
#' # Cluster the samples
#' tse <- addCluster(tse, DmmParam(k = 1:4), name = "DMM", full = TRUE)
#'
#' # Plot the fit
#' plotDMNFit(tse, name = "DMM", type = "laplace")
#'
NULL

#' @rdname plotDMN
#' @importFrom DirichletMultinomial mixture
#' @importFrom ggplot2 ggplot aes geom_point geom_line theme_bw labs
#' @export
setMethod("plotDMNFit", signature = c(x = "SummarizedExperiment"),
    function(x, name = "DMN", type = c("laplace","AIC","BIC")){
        #
        if (!is.null(metadata(x)[[name]]$dmm)) {
            dmn <- metadata(x)[[name]]$dmm
        } else {
            .Deprecated(old="getDMN", new="addCluster",
                paste0(
                    "Now runDMN and calculateDMN are deprecated. Use ",
                    "addCluster with DMMParam parameter and full ",
                    "parameter set as true instead."))
            dmn <- metadata(x)[[name]]
        }
        fit_FUN <- mia:::.get_dmn_fit_FUN(type)
        #
        k <- vapply(dmn, function(d){ncol(mixture(d))}, numeric(1))
        fit <- vapply(dmn, fit_FUN, numeric(1))
        ggplot(data.frame(k = k, fit = fit), aes(x = k, y = fit)) +
            geom_point() +
            geom_line() +
            theme_bw() +
            labs(
                x = "Number of Dirichlet Components",
                y = paste0("Model Fit (",type,")"))
    }
)
