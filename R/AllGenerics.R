# All generic methods are listed here

#' @rdname getNeatOrder
setGeneric("getNeatOrder", signature = c("x"),
    function(x, centering = "mean", ...)
    standardGeneric("getNeatOrder"))

#' @rdname plotAbundance
setGeneric("plotAbundance", signature = c("x"), function(x, ...)
    standardGeneric("plotAbundance"))

#' @rdname plotAbundanceDensity
#' @export
setGeneric("plotAbundanceDensity", signature = c("x"), function(x, ...)
    standardGeneric("plotAbundanceDensity"))

#' @rdname plotCCA
#' @aliases plotRDA
#' @export
setGeneric("plotCCA", signature = c("x"), function(x, ...)
    standardGeneric("plotCCA"))

#' @rdname plotCCA
#' @aliases plotCCA
#' @export
setGeneric("plotRDA", signature = c("x"), function(x, ...)
    standardGeneric("plotRDA"))

#' @rdname plotColTile
#' @export
setGeneric("plotColTile", signature = c("object"),
    function(object, x, y, ...)
    standardGeneric("plotColTile"))

#' @rdname plotColTile
#' @export
setGeneric("plotRowTile", signature = c("object"),
    function(object, x, y, ...)
    standardGeneric("plotRowTile"))

#' @rdname plotDMN
#' @export
setGeneric("plotDMNFit", signature = "x",
    function(x, name = "DMN", type = c("laplace","AIC","BIC"), ...)
    standardGeneric("plotDMNFit"))

#' @rdname plotGraph
#' @export
setGeneric("plotColGraph", signature = c("x","y"),
    function(x, y, ...)
    standardGeneric("plotColGraph"))

#' @rdname plotGraph
#' @export
setGeneric("plotRowGraph", signature = c("x","y"),
    function(x, y, ...)
    standardGeneric("plotRowGraph"))

#' @rdname plotLoadings
setGeneric("plotLoadings", signature = c("x"),
    function(x, ...)
    standardGeneric("plotLoadings"))

#' @rdname plotMediation
#' @export
setGeneric("plotMediation", signature = c("x"),
    function(x, ...)
    standardGeneric("plotMediation"))

#' @rdname plotPrevalence
#' @export
setGeneric("plotRowPrevalence", signature = c("x"),
    function(x, ...)
    standardGeneric("plotRowPrevalence"))

#' @rdname plotPrevalence
#' @export
setGeneric("plotPrevalentAbundance", signature = c("x"),
    function(x, ...)
    standardGeneric("plotPrevalentAbundance"))

#' @rdname plotPrevalence
#' @export
setGeneric("plotPrevalence", signature = c("x"),
    function(x, ...)
    standardGeneric("plotPrevalence"))

#' @rdname plotScree
#' @export
setGeneric("plotScree", signature = c("x"),
    function(x, ...)
    standardGeneric("plotScree"))

#' @rdname plotSeries
#' @export
setGeneric("plotSeries", signature = c("object"),
    function(object, ...)
    standardGeneric("plotSeries"))

#' @rdname plotTree
setGeneric("plotRowTree", signature = c("x"),
    function(x, ...)
    standardGeneric("plotRowTree"))
#' @rdname plotTree
setGeneric("plotColTree", signature = c("x"),
    function(x, ...)
    standardGeneric("plotColTree"))

#' @rdname treeData
setGeneric("rowTreeData", signature = c("x"),
    function(x, ...)
    standardGeneric("rowTreeData"))

#' @rdname treeData
setGeneric("colTreeData", signature = c("x"),
    function(x, ...)
    standardGeneric("colTreeData"))

#' @rdname treeData
setGeneric("rowTreeData<-", signature = c("x"),
    function(x, tree.name = tree_name, tree_name = "phylo", value)
    standardGeneric("rowTreeData<-"))

#' @rdname treeData
setGeneric("colTreeData<-", signature = c("x"),
    function(x, tree.name = tree_name, tree_name = "phylo", value)
    standardGeneric("colTreeData<-"))

#' @rdname treeData
setGeneric("combineTreeData", signature = c("x"),
    function(x, other.fields = other_fields, other_fields = list())
    standardGeneric("combineTreeData"))

#' @rdname treeData
setGeneric("combineTreeData", signature = c("x"),
    function(x, other.fields = other_fields, other_fields = list())
    standardGeneric("combineTreeData"))

#' @rdname plotHistogram
#' @export
setGeneric("plotHistogram", signature = c("x"), function(x, ...)
    standardGeneric("plotHistogram"))

#' @rdname plotHistogram
#' @export
setGeneric("plotBarplot", signature = c("x"), function(x, ...)
    standardGeneric("plotBarplot"))