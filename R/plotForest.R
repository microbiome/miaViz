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
#' @param effect.var \code{Character scalar}. Specifies the variable of x which
#' corresponds to the effects or estimated results. (Default: \code{"effect"})
#' 
#' @param id.var \code{Character scalar}. Specifies the variable of x which
#' corresponds to the effects or estimated results. (Default: \code{"rownames"})
#' 
#' @param ci.lower.var \code{Character scalar}. Specifies the variable of x
#' which corresponds to the lower CI boundaries. (Default: \code{"lower"})
#'
#' @param ci.upper.var \code{Character scalar}. Specifies the variable of x
#' which corresponds to the upper CI boundaries. (Default: \code{"upper"})
#'
#' @param pval.var \code{Character scalar}. Specifies the variable of x which
#' corrsponds to the p-values associated with \code{effect.var}.
#' (Default: \code{"pval"})
#' 
#' @param label.by \code{Character vector}. Specifies the variables of x or
#' row/colData(x) by which the plot should be labelled. \code{"CI"} and
#' \code{"P-Value"} are special entries which require either \code{effect.var},
#' \code{ci.lower.var} and \code{ci.upper.var} or \code{pval.var} to be
#' specified, respectively.
#' 
#' @param colour.by \code{Character scalar}. Specifies the group for plotting. Must
#' be a value from \code{names(rowData(x))} or \code{names(colData(x))}. If
#' \code{NULL}, observations are not grouped. (Default: \code{NULL})
#'
#' @param color.by \code{Character scalar}. Specifies the group for plotting. Must
#' be a value from \code{names(rowData(x))} or \code{names(colData(x))}. If
#' \code{NULL}, observations are not grouped. (Default: \code{NULL})
#'
#' @param show.tree \code{Logical scalar}. Should the tree structure of the data
#' be shown next to the forest plot?
#' 
#' @param ... additional parameters from \code{plotRowTree}.
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
    function(x, by = 1, effect.var = "effect", id.var = "rownames",
        ci.lower.var = "lower", ci.upper.var = "upper", pval.var = "pval",
        label.by = "CI", color.by = colour.by, colour.by = NULL,
        show.tree = TRUE, ...){
    # Check margin (by)
    by <- .check_MARGIN(by)
    # Select data based on margin (by)
    FUN <- switch(by, rowData, colData)
    tree.FUN <- switch(by, rowTree, colTree)
    treeplot.FUN <- switch(by, plotRowTree, plotColTree)
    # Extract side information from SE
    df <- as.data.frame(FUN(x))
    tree.exists <- is(x, "TreeSummarizedExperiment") && !is.null(tree.FUN(x))
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
    # Generate forest plot components
    plots <- c(plots, .plot_forest(
        x = df,
        effect.var = effect.var,
        id.var = id.var,
        ci.lower.var = ci.lower.var,
        ci.upper.var = ci.upper.var,
        pval.var = pval.var,
        label.by = label.by,
        color.by = color.by
    ))
    # Define plot widths
    widths <- c(widths, 2, 2/5 * length(label.by))
    # Make final plot
    p <- wrap_plots(plots) +
        plot_layout(widths = widths, guides = "collect")
    return(p)
})

setMethod("plotForest", signature = c(x = "data.frame"),
    function(x, effect.var = "effect", id.var = "rownames",
        ci.lower.var = "lower", ci.upper.var = "upper", pval.var = "pval",
        label.by = NULL, color.by = colour.by, colour.by = NULL){
    # Generate forest plot components
    plots <- .plot_forest(
        x = x,
        effect.var = effect.var,
        id.var = id.var,
        ci.lower.var = ci.lower.var,
        ci.upper.var = ci.upper.var,
        pval.var = pval.var,
        label.by = label.by,
        color.by = color.by
    )
    # Define plot widths
    widths <- c(2, 2/5 * length(label.by))
    # Make final plot
    p <- wrap_plots(plots) +
        plot_layout(widths = widths, guides = "collect")
    return(p)
})

.plot_forest <- function(x, effect.var, id.var, ci.lower.var,
    ci.upper.var, pval.var, label.by, color.by){
    # Check main vars
    if( !effect.var %in% names(x) ){
        stop("'effect.var' must be a variable in x.", call. = FALSE)
    }
    if( !id.var %in% c("rownames", names(x)) ){
        stop("'id.var' must be 'rownames' or a variable in x.", call. = FALSE)
    }
    if( any(!is.null(c(ci.lower.var, ci.upper.var))) &&
        any(!c(ci.lower.var, ci.upper.var) %in% names(x)) ){
        warning("'ci.lower.var' and 'ci.upper.var' are ignored when not found ",
            "in x.", call. = FALSE)
        ci.lower.var <- NULL
        ci.upper.var <- NULL
    }
    if( !is.null(pval.var) && !pval.var %in% names(x) ){
        warning("'pval.var' is ignored when not found in x.", call. = FALSE)
        pval.var <- NULL
    }
    # Check label.by
    if( !is.vector(label.by) &&
        (is.null(label.by) || is.na(label.by) || label.by == "") ){
        label.by <- c()
    }
    if( !all(label.by %in% c("CI", "P-Value") | label.by %in% names(x)) ){
        stop("All terms in 'label.by' must be variables of x.",
            call. = FALSE)
    }
    # Check CI
    ci.exists <- !is.null(ci.lower.var) && !is.null(ci.upper.var)
    if( "CI" %in% label.by && !ci.exists ){
        stop("To show the 'CI' label, x must include the variables specified ",
            "with 'effect.var', 'ci.lower.var' and 'ci.upper.var'.",
            call. = FALSE)
    }
    # Check P-Value
    if( "P-Value" %in% label.by && is.null(pval.var) ){
        stop("To show the 'P-Value' label, x must include the variable ",
            "specified with 'pval.var'.", call. = FALSE)
    }
    # Check color.by
    if( !is.null(color.by) && !color.by %in% names(x) ){
        stop("'color.by' must be a variable of x.", call. = FALSE)
    }
    # Initialise plots and widths lists
    plots <- list()
    widths <- c()
    # Convert feature names to factors
    if( id.var == "rownames" ){
        x[[id.var]] <- factor(rownames(x), levels = rownames(x))
    }else{
        x[[id.var]] <- factor(x[[id.var]], levels = unique(x[[id.var]]))
    }
    # Store first name to set y origin
    y0 <- x[[id.var]][1]
    # Remove rownames to control row ordering
    rownames(x) <- NULL
    # Find x-axis limits for forest plot
    if( ci.exists ){
        lim <- max(abs(x[[ci.lower.var]]), x[[ci.upper.var]])
    }else{
        lim <- max(abs(x[[effect.var]]))
    }
    # Account for group
    if( is.null(color.by) ){
        p <- ggplot(x, aes(x = .data[[effect.var]], y = .data[[id.var]]))
    }else{
        p <- ggplot(x, aes(x = .data[[effect.var]], y = .data[[id.var]],
            colour = .data[[color.by]]))
    }
    # Make forest plot
    p <- p + geom_vline(xintercept = 0, linetype = "dashed", colour = "gray") +
        geom_point(position = position_dodge2(width = 0.75)) +
        coord_cartesian(xlim = c(-lim, lim), ylim = c(y0, NA), clip = "off") +
        theme_bw() +
        theme(axis.title.y = element_blank())
    # If CI is defined
    if( ci.exists ){
        # Add errorbars
        p <- p + geom_errorbar(
            aes(xmin = .data[[ci.lower.var]], xmax = .data[[ci.upper.var]]),
            orientation = "y", width = 1e-2 * nrow(x),
            position = position_dodge2(width = 0.75)
        )
    }
    # Store forest plot
    plots[[length(plots) + 1]] <- p
    # Store forest plot width
    widths <- c(widths, 2)
    # Construct CI label
    if( "CI" %in% label.by ){
        x$CI <- paste0(
            round(x[[effect.var]], 2), " (",
            round(x[[ci.lower.var]], 2), "—",
            round(x[[ci.upper.var]], 2), ")"
        )
    }
    # Construct P-Value label
    if( "P-Value" %in% label.by ){
        x$`P-Value` <- paste0(
            round(x[[pval.var]], 3),
            ifelse(x[[pval.var]] < 0.05, "*", "")
        )
    }
    # Initialise label plot list
    label.plots <- list()
    label.size <- nonlinear.textsize(nrow(x))
    ann.ypos <- -nrow(x) / 30
    
    if( !is.null(color.by) && nrow(x) > length(unique(x[[id.var]])) ){
        p <- ggplot(x, aes(x = .data[[effect.var]], y = .data[[id.var]],
            colour = .data[[color.by]]))
        # Adjust annotation y-position for grouped rows
        ann.ypos <- ann.ypos / length(unique(x[[color.by]]))
    }else{
        p <- ggplot(x, aes(x = .data[[effect.var]], y = .data[[id.var]]))
    }
    # Iterate over label.by terms
    for( i in seq_along(label.by) ){
        # Retrieve current label
        lab <- label.by[i]
        # Set xmax for improved positioning
        if( i == 1 ){
            xmax <- 1
        }else{
            xmax <- NA
        }
        
        # Make plot for current label
        label.plots[[i]] <- p +
            geom_text(x = 0, aes(y = .data[[id.var]], label = .data[[lab]]),
                      hjust = 0, position = position_dodge2(width = 0.75),
                      size = label.size, show.legend = FALSE) +
            annotate("text", x = 0, y = ann.ypos, label = lab, hjust = 0) +
            scale_x_continuous(expand = expansion(mult = c(0, 0))) +
            coord_cartesian(xlim = c(0, xmax), ylim = c(y0, NA), clip = "off") +
            theme(axis.title = element_blank(),
                  axis.text = element_blank(),
                  axis.ticks = element_blank(),
                  panel.background = element_blank(),
                  panel.grid = element_blank())
    }
    # For non-empty label plot lists
    if( length(label.plots) != 0 ){
        # Wrap and store label plots
        plots[[length(plots) + 1]] <- wrap_plots(label.plots)
    }
    return(plots)
}

nonlinear.textsize <- function(n, min.size = 3, max.size = 5) {
    # Scale size inversely but bounded between min_size and max_size
    size <- max.size - (log10(n) * (max.size - min.size) / log10(100))
    # Clamp between min and max
    size <- pmax(pmin(size, max.size), min.size)
    return(size)
}

# From https://github.com/microbiome/mia/blob/devel/R/utils.R
# Check MARGIN parameters. Should be defining rows or columns.
.check_MARGIN <- function(MARGIN, name = .get_name_in_parent(MARGIN)) {
    # MARGIN must be one of the following options
    if( !(length(MARGIN) == 1L && tolower(MARGIN) %in% c(
            1, 2, "1", "2", "features", "samples", "columns", "col", "row",
            "rows", "cols")) ) {
        stop("'", name,"' must be 'rows' or 'cols'.", call. = FALSE)
    }
    # Convert MARGIN to numeric if it is not.
    MARGIN <- ifelse(tolower(MARGIN) %in% c(
        "samples", "columns", "col", 2, "cols"), 2, 1)
    return(MARGIN)
}

was_called_by <- function(caller_name) {
  calls <- sys.calls()
  # Check if caller_name appears in any of the calls
  any(sapply(calls, function(call) {
    as.character(call[[1]]) == caller_name
  }))
}
