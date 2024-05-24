#' @export
plotNMDS <- function(x, ..., ncomponents = 2){
    plotReducedDim(x, ncomponents = ncomponents, dimred = "NMDS",
                   ...)
}
