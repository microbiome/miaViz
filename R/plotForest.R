#' Visualize estimated results with forest plots
#' 
#' \code{plotForest()} creates a feature- or sample-wise forest plot, showing
#' estimated results from a statistical test with their confidence intervals.
#' Additionally, the plot can be enriched with the tree structure and labelled
#' with Confidence Intervals (CIs), p-values and other information from the
#' metadata.
#' 
#' @param x a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object, or a \code{data.frame} object containing statistical estimates, with
#' rownames and the variables \code{"Estimate"}, \code{"Lower"} and
#' \code{"Upper"} as described in the details.
#' 
#' @param by \code{Character scalar}. Determines whether features or samples
#' data is used for the plot. (Default: \code{"rows"})
#' 
#' @param group \code{Character scalar}. Specifies the group for plotting. Must
#' be a value from \code{names(rowData(x))} or \code{names(colData(x))}. If
#' \code{NULL}, observations are not grouped. (Default: \code{NULL})
#'
#' @param label.by \code{Character vector}. Lists the variables from
#' \code{rowData(x)} or \code{colData(x)} by which the plot should be labelled.
#' \code{"CI"} and \code{"P-Value"} are special entries with requirements
#' described in the details.
#'
#' @param show.tree \code{Logical scalar}. Should the tree structure of the data
#' be shown next to the forest plot?
#' 
#' @param ... additional parameters from \code{plotRowTree}.
#' 
#' @details
#' \code{plotForest} requires the following variables to be present in the
#' \code{rowData(x)} or \code{colData(x)}:
#' 
#' \itemize{
#'   \item \code{Estimate}: estimated results or statistics 
#'   \item \code{Lower}: lower CI boundaries
#'   \item \code{Upper}: upper CI boundaries
#'   \item \code{Pval}: associated p-values
#' }
#' 
#' When \code{label.by} includes \code{"CI"} or \code{"P-Value"}, those are
#' constructed from the four variables above.
#' 
#' @return
#' a \code{\link[ggplot2:ggplot]{ggplot}} object.
#'
#' @examples
#'
#'
#'@name plotForest
NULL

#' @rdname plotForest
#' @export
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom ggplot2 ggplot geom_text
#' @importFrom SummarizedExperiment rowData colData
setMethod("plotForest", signature = c(x = "SummarizedExperiment"),
    function(x, by = 1, group = NULL, label.by = "CI", show.tree = TRUE, ...){
    # Select data based on margin (by)
    FUN <- switch(by, rowData, colData)
    tree.FUN <- switch(by, rowTree, colTree)
    treeplot.FUN <- switch(by, plotRowTree, plotColTree)
    # Extract side information from SE
    df <- as.data.frame(FUN(x))
    tree.exists <- is(x, "TreeSummarizedExperiment") && !is.null(tree.FUN(x))
    # Check label.by
    if( !is.vector(label.by) &&
        (is.null(label.by) || is.na(label.by) || label.by == "") ){
        label.by <- c()
    }
    if( !all(label.by %in% c("CI", "P-Value") | label.by %in% names(df)) ){
        stop("All terms in 'label.by' must be in row/colData(x).",
            call. = FALSE)
    }
    # Check CI
    if( "CI" %in% label.by && !all(c("Estimate", "Lower", "Upper") %in% names(df)) ){
        stop("To show the 'CI' label, row/colData(x) must include the ",
             "variables 'Estimate', 'Lower' and 'Upper'.", call. = FALSE)
    }
    # Check P-Value
    if( "P-Value" %in% label.by && !"Pval" %in% names(df) ){
        stop("To show the 'P-Value' label, row/colData(x) must include the ",
            "variable 'Pval'.", call. = FALSE)
    }
    # Check show.tree
    if( !is.logical(show.tree) ){
        stop("'show.tree' must be TRUE or FALSE.", call. = FALSE)
    }
    if( show.tree && !tree.exists ){
        warning("'show.tree' is ignored when row/colTree(x) does not exist.",
            call. = FALSE)
    }
    # Check kwargs
    kwargs <- list(...)
    if( !show.tree && length(kwargs) != 0 ){
        warning("Arguments in '...' are ignored when 'show.tree' is FALSE.",
            call. = FALSE)
    }
    if( any(c("layout", "branch.length") %in% names(kwargs)) ){
        stop("'layout' and 'branch.length' cannot be modified for this plot.",
            call. = FALSE)
    }
    if( tree.exists ){
        # Plot tree
        tree.plot <- treeplot.FUN(
            x,
            layout = "rectangular",
            branch.length = "none",
            show.label = TRUE
        )
        # Extract tip order from tree
        tree.data <- ggplot_build(tree.plot)
        tips <- tree.data[[1]][[5]][c("y", "label")]
        tips <- tips[order(tips$y), ]
        # Order rowData based on tree tips
        tax.order <- match(tips$label, names(x))
        df <- df[tax.order, , drop = FALSE]
    }
    # Initialise plots and widths lists
    plots <- list()
    widths <- c()
    # If tree is plotted
    if( tree.exists && show.tree ){
        # Store tree plot
        plots[[1]] <- treeplot.FUN(
            x,
            layout = "rectangular",
            branch.length = "none",
            ...
        )
        # Store tree width
        widths <- c(widths, 0.5)
    }
    # Pass to PlotForest data.frame method
    plots[[length(plots) + 1]] <- plotForest(
        df,
        group = group,
        label.by = label.by
    )
    # Store widths for forest plot and label plots
    widths <- c(widths, 2 + 1/3 * length(label.by))
    # Make final plot
    p <- wrap_plots(plots) +
        plot_layout(widths = widths, guides = "collect")
    return(p)
})

setMethod("plotForest", signature = c(x = "data.frame"),
    function(x, group = NULL, label.by = "CI"){
    # Initialise plots and widths lists
    plots <- list()
    widths <- c()
    # Convert feature names to factors
    x$Feature <- factor(rownames(x), levels = rownames(x))
    rownames(x) <- NULL
    # Find x-axis limits for forest plot
    lim <- max(abs(x$Lower), x$Upper)
    # Make forest plot
    plots[[length(plots) + 1]] <- ggplot(x, aes(x = .data$Estimate, y = .data$Feature)) +
        geom_vline(xintercept = 0, linetype = "dashed", colour = "gray") +
        geom_point() +
        geom_errorbar(aes(xmin = .data$Lower, xmax = .data$Upper),
                      orientation = "y", width = 1e-2 * nrow(x)) +
        coord_cartesian(clip = "off", xlim = c(-lim, lim), ylim = c(x$Feature[1], NA)) +
        theme_bw() +
        theme(axis.title.y = element_blank(),
              panel.grid = element_blank())
    # Store forest plot width
    widths <- c(widths, 2)
    # Construct CI label
    if( "CI" %in% label.by ){
        x$CI <- paste0(round(x$Estimate, 2), " (", round(x$Lower, 2), "—", round(x$Upper, 2), ")")
    }
    # Construct P-Value label
    if( "P-Value" %in% label.by ){
        x$`P-Value` <- paste0(round(x$Pval, 3), ifelse(x$Pval < 0.05, "*", ""))
    }
    # Initialise label plot list
    label.plots <- list()
    # Iterate over label.by terms
    for( i in seq_along(label.by) ){
        
        lab <- label.by[i]
        # Set xmax for improved positioning
        if( i == 1 ){
            xmax <- 1
        }else{
            xmax <- NA
        }
        # Make plot for current label
        label.plots[[i]] <- ggplot(x, aes(x = .data$Estimate, y = .data$Feature)) +
            geom_text(x = 0, aes(label = .data[[lab]]),
                      hjust = 0, size = 150 / nrow(x)) +
            annotate("text", x = 0, y = -1.5, label = lab, hjust = 0) +
            coord_cartesian(xlim = c(0, xmax), ylim = c(x$Feature[1], NA),
                            expand = FALSE, clip = "off") +
            theme(axis.title = element_blank(),
                  axis.text = element_blank(),
                  axis.ticks = element_blank(),
                  panel.background = element_blank())
    }
    # For non-empty label plot lists
    if( length(label.plots) != 0 ){
        # Wrap and store label plots
        plots[[length(plots) + 1]] <- wrap_plots(label.plots)
        # Store label plots total width
        widths <- c(widths, 1/3 * length(label.plots))
    }
    # Make final plot
    p <- wrap_plots(plots) +
        plot_layout(widths = widths)
    return(p)
})
