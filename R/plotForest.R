#' Visualize estimated results with forest plots
#'
#' \code{plotForest()} creates a feature- or sample-wise forest plot, showing
#' estimated results from a statistical test with their confidence intervals.
#' Additionally, the plot can be enriched with the tree structure and labelled
#' with Confidence Intervals (CIs), p-values and other side information.
#'
#' @param x a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object, or a \code{data.frame} object containing statistical estimates.
#'
#' @param by \code{Character scalar}. Determines whether features or samples
#' data is used for the plot. (Default: \code{"rows"})
#'
#' @param effect.var \code{Character scalar}. Specifies the variable of x which
#' corresponds to the effects or estimated results. (Default: \code{"effect"})
#'
#' @param ci.lower.var \code{Character scalar}. Specifies the variable of x
#' which corresponds to the lower CI boundaries. (Default: \code{"lower"})
#'
#' @param ci.upper.var \code{Character scalar}. Specifies the variable of x
#' which corresponds to the upper CI boundaries. (Default: \code{"upper"})
#'
#' @param err.var \code{Character scalar}. Specifies the variable of x which
#' corresponds to the standard errors associated with \code{effect.var}.
#' When defined, it overwrites \code{ci.lower.var} and \code{ci.upper.var}.
#' (Default: \code{"pval"})
#'
#' @param pval.var \code{Character scalar}. Specifies the variable of x which
#' corresponds to the p-values associated with \code{effect.var}.
#' (Default: \code{"pval"})
#'
#' @param id.var \code{Character scalar}. Specifies the variable of x which
#' corresponds to the observation identifiers. When \code{"rownames"}),
#' the object rownames are used. (Default: \code{"rownames"})
#'
#' @param label.by \code{Character vector}. Specifies the variables of x or
#' row/colData(x) by which the plot should be labelled. \code{"CI"} and
#' \code{"P-Value"} are special entries which require either \code{effect.var},
#' \code{ci.lower.var} and \code{ci.upper.var} or \code{pval.var} to be
#' specified, respectively. (Default: \code{NULL})
#'
#' @param order.by \code{Character scalar}. Specifies the variable of x by which
#' observations are ordered. If \code{NULL}, the observations are ordered by
#' the tree structure if available. (Default: \code{NULL})
#'
#' @param facet.by \code{Character scalar}. Specifies the variable of x by which
#' observations are divided into horizontal facets. (Default: \code{NULL})
#'
#' @param colour.by \code{Character scalar}. Specifies the variable of x by
#' which observations are coloured. (Default: \code{NULL})
#'
#' @param color.by \code{Character scalar}. Alias to \code{colour.by}.
#'
#' @param conf.level \code{Numeric scalar}. Specifies the confidence level of
#' the interval when inferred from \code{err.var}. It is ignored when
#' \code{ci.lower.var} and \code{ci.upper.var} are defined.
#' (Default: \code{0.95})
#'
#' @param tree.name \code{Character scalar}. Specifies a row/colTree from x.
#' (Default: \code{"phylo"})
#'
#' @param show.tree \code{Logical scalar}. Should the tree structure of the data
#' be shown next to the forest plot?
#'
#' @param ... additional parameters passed to \code{\link{plotRowTree}}.
#'
#' @return
#' a \code{\link[ggplot2:ggplot]{ggplot}} object.
#'
#' @examples
#' library(mia)
#' library(maaslin3)
#'
#' # Import dataset
#' data("Tengeler2020", package = "mia")
#' tse <- Tengeler2020
#'
#' # Agglomerate by genus and subset by prevalence
#' tse <- subsetByPrevalent(tse, rank = "Genus", prevalence = 10/100)
#'
#' # Transform count assay to relative abundances
#' tse <- transformAssay(tse, assay.type = "counts", method = "relabundance")
#'
#' # Run maaslin3
#' maaslin3_out <- maaslin3(
#'     input_data = tse,
#'     output = "maaslin_results",
#'     formula = "~ patient_status",
#' )
#'
#' # Retrieve abundance results
#' maaslin3_abund <- maaslin3_out$fit_data_abundance$results
#' maaslin3_abund <- maaslin3_abund[!is.na(maaslin3_abund$coef), ]
#'
#' # Visualize abundance results
#' plotForest(
#'     maaslin3_abund,
#'     effect.var = "coef",
#'     err.var = "stderr",
#'     pval.var = "qval_joint",
#'     id.var = "feature",
#'     label.by = c("CI", "P-Value"),
#'     order.by = "coef"
#' )
#'
#' # Add abundance results to TreeSE rowData
#' rownames(maaslin3_abund) <- maaslin3_abund$feature
#' tax_order <- match(rownames(tse), rownames(maaslin3_abund))
#' rowData(tse) <- cbind(rowData(tse), maaslin3_abund[tax_order, ])
#'
#' # Visualise abundance results with tree structure
#' plotForest(
#'     tse,
#'     effect.var = "coef",
#'     err.var = "stderr",
#'     pval.var = "qval_joint",
#'     label.by = c("CI", "P-Value")
#' )
#'
#' # Retrieve prevalence results
#' maaslin3_prev <- maaslin3_out$fit_data_prevalence$results
#' maaslin3_prev <- maaslin3_prev[!is.na(maaslin3_prev$coef), ]
#'
#' # Combine abundance and prevalence results
#' maaslin3_res <- rbind(maaslin3_abund, maaslin3_prev)
#'
#' maaslin3_res$association <- c(
#'     rep("Abundance", nrow(maaslin3_abund)),
#'     rep("Prevalence", nrow(maaslin3_prev))
#' )
#'
#' # Visualize combined results
#' plotForest(
#'     maaslin3_res,
#'     effect.var = "coef",
#'     err.var = "stderr",
#'     pval.var = "qval_joint",
#'     id.var = "feature",
#'     label.by = c("CI", "P-Value"),
#'     order.by = "coef",
#'     facet.by = "association",
#'     colour.by = "association"
#' )
#'
#' @name plotForest
NULL

#' @rdname plotForest
#' @export
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom TreeSummarizedExperiment rowTree colTree rowTreeNames
#'   colTreeNames rowLinks
#' @importFrom ggplot2 ggplot_build
#' @importFrom patchwork wrap_plots plot_layout
setMethod("plotForest", signature = c(x = "TreeSummarizedExperiment"),
    function(x, by = 1L, effect.var = "effect", ci.lower.var = "lower",
        ci.upper.var = "upper", err.var = NULL, pval.var = "pval",
        id.var = "rownames", label.by = NULL, order.by = NULL, facet.by = NULL,
        color.by = colour.by, colour.by = NULL, conf.level = 0.95,
        tree.name = "phylo", show.tree = TRUE, ...){
    # Check margin (by)
    by <- .check_MARGIN(by)
    # Select data based on margin (by)
    FUN <- switch(by, rowData, colData)
    tree_FUN <- switch(by, rowTree, colTree)
    treename_FUN <- switch(by, rowTreeNames, colTreeNames)
    treeplot_FUN <- switch(by, plotRowTree, plotColTree)
    # Extract side information from SE
    df <- as.data.frame(FUN(x))
    # Check show.tree
    if( !.is_a_bool(show.tree) ){
        stop("'show.tree' must be TRUE or FALSE.", call. = FALSE)
    }
    order_tree <- is.null(order.by)
    if( !order_tree && show.tree ){
        warning("'show.tree' is ignored when 'order.by' is defined.",
            call. = FALSE)
    }
    # Check kwargs
    kwargs <- list(...)
    if( any(c("layout", "levels.rm", "branch.length") %in% names(kwargs)) ){
        stop("'layout', 'levels.rm' and 'branch.length' cannot be modified ",
            "for this plot.", call. = FALSE)
    }
    # If tree is available
    if( order_tree ){
        # Plot tree
        tree_plot <- treeplot_FUN(
            x,
            tree.name = tree.name,
            layout = "rectangular",
            branch.length = "none",
            show.label = TRUE,
            levels.rm = TRUE
        )
        # Extract tip order from tree
        tree_data <- ggplot_build(tree_plot)
        tips <- tree_data[[1]][[5]][c("y", "label")]
        tips <- tips[order(tips$y), , drop = FALSE]
        # Order rowData based on tree tips
        tax_order <- match(tips$label, rowLinks(x)$nodeLab)
        df <- df[tax_order, , drop = FALSE]
    }
    # Initialise plots and widths lists
    plots <- list()
    widths <- c()
    # If tree is plotted
    if( show.tree && order_tree ){
        # Store tree plot
        plots[[1]] <- treeplot_FUN(
            x,
            tree.name = tree.name,
            layout = "rectangular",
            branch.length = "none",
            levels.rm = TRUE,
            ...
        )
        # Store tree width
        widths <- c(widths, 0.5)
    }
    # Generate forest plot components
    plots <- c(plots, .plot_forest(
        df, effect.var, ci.lower.var, ci.upper.var, err.var, pval.var, id.var,
        label.by, order.by, facet.by, color.by, conf.level
    ))
    # Define plot widths
    widths <- c(widths, 2, 2/5 * length(label.by))
    # Make final plot
    p <- wrap_plots(plots) +
        plot_layout(widths = widths, guides = "collect")
    return(p)
})

#' @rdname plotForest
#' @export
#' @importFrom SummarizedExperiment rowData colData
setMethod("plotForest", signature = c(x = "SummarizedExperiment"),
    function(x, by = 1L, effect.var = "effect", ci.lower.var = "lower",
        ci.upper.var = "upper", err.var = NULL, pval.var = "pval",
        id.var = "rownames", label.by = NULL, order.by = NULL, facet.by = NULL,
        color.by = colour.by, colour.by = NULL, conf.level = 0.95){
    # Check margin (by)
    by <- .check_MARGIN(by)
    # Select data based on margin (by)
    FUN <- switch(by, rowData, colData)
    # Extract side information from SE
    df <- as.data.frame(FUN(x))
    # Generate forest plot components
    p <- plotForest(
        df, effect.var, ci.lower.var, ci.upper.var, err.var, pval.var, id.var,
        label.by, order.by, facet.by, color.by, conf.level
    )
    return(p)
})

#' @rdname plotForest
#' @export
#' @importFrom patchwork wrap_plots plot_layout
setMethod("plotForest", signature = c(x = "data.frame"),
    function(x, effect.var = "effect",  ci.lower.var = "lower",
        ci.upper.var = "upper", err.var = NULL, pval.var = "pval",
        id.var = "rownames", label.by = NULL, order.by = NULL, facet.by = NULL,
        color.by = colour.by, colour.by = NULL, conf.level = 0.95){
    # Check order.by
    if( !is.null(order.by) && order.by == "tree" ){
        warning("'order.by' can be set to tree only when 'x' is a ",
        "TreeSummarizedExperiment object.", call. = FALSE)
    }
    # Generate forest plot components
    plots <- .plot_forest(
        x, effect.var, ci.lower.var, ci.upper.var, err.var, pval.var, id.var,
        label.by, order.by, facet.by, color.by, conf.level
    )
    # Define plot widths
    widths <- c(2, 2/5 * length(label.by))
    # Make final plot
    p <- wrap_plots(plots) +
        plot_layout(widths = widths, guides = "collect")
    return(p)
})

#' @importFrom stats qt
.plot_forest <- function(x, effect.var, ci.lower.var, ci.upper.var, err.var,
    pval.var, id.var, label.by, order.by, facet.by, color.by, conf.level){
    # Check main vars
    if( !effect.var %in% names(x) ){
        stop("'effect.var' must be a variable in x.", call. = FALSE)
    }
    if( !is.null(err.var) && !err.var %in% names(x) ){
        warning("'err.var' is ignored when not found in x.", call. = FALSE)
        err.var <- NULL
    }
    # Derive CI from standard error
    if( !is.null(err.var) ){
        # Compute corresponding quantile
        z <- qt(1 - (1 - conf.level) / 2, df = ncol(x) - 1)
        # Compute CI boarders
        x[[ci.lower.var]] <- x[[effect.var]] - z * x[[err.var]]
        x[[ci.upper.var]] <- x[[effect.var]] + z * x[[err.var]]
    }
    if( any(!is.null(c(ci.lower.var, ci.upper.var))) &&
        any(!c(ci.lower.var, ci.upper.var) %in% names(x)) ){
        warning("'ci.lower.var' and 'ci.upper.var' are ignored when not found ",
            "in x.", call. = FALSE)
        ci.lower.var <- NULL
        ci.upper.var <- NULL
    }
    if( !is.null(err.var) &&
        (ci.lower.var != "lower" || ci.upper.var != "upper") ){
        warning("'ci.lower.var' and 'ci.upper.var' are ignored when 'err.var' ",
            "is defined.", call. = FALSE)
    }
    if( !is.null(pval.var) && !pval.var %in% names(x) ){
        warning("'pval.var' is ignored when not found in x.", call. = FALSE)
        pval.var <- NULL
    }
    if( !id.var %in% c("rownames", names(x)) ){
        stop("'id.var' must be 'rownames' or a variable in x.", call. = FALSE)
    }
    # Check label.by
    if( !is.vector(label.by) && .is_non_empty_string(label.by) ){
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
    # Check aesthetics
    aes_list <- list(
        order.by = order.by, facet.by = facet.by, color.by = color.by)
    lapply(names(aes_list), function(s.name){
        s <- aes_list[[s.name]]
        if( !is.null(s) && (length(s) != 1L || !s %in% names(x)) ){
            stop("'", s.name, "' must be a single variable of x.",
                call. = FALSE)
        }
    })
    # Check conf.level
    if( !.is_a_numeric(conf.level) || !(conf.level >= 0 & conf.level <= 1) ){
        stop("'conf.level' must be a number between 0 and 1.", call. = FALSE)
    }
    plots <- .combine_forest_components(
        x, effect.var, ci.lower.var, ci.upper.var, pval.var, id.var,
        label.by, order.by, facet.by, color.by, ci.exists
    )
    return(plots)
}

#' @importFrom patchwork wrap_plots
.combine_forest_components <- function(x, effect.var, ci.lower.var,
    ci.upper.var, pval.var, id.var, label.by, order.by, facet.by, color.by,
    ci.exists){
    # Order features by effect size
    if( !is.null(order.by) ){
        x <- x[order(x[[order.by]]), , drop = FALSE]
    }
    # Convert feature names to factors
    if( id.var == "rownames" ){
        x[[id.var]] <- factor(rownames(x), levels = rownames(x))
    }else{
        x[[id.var]] <- factor(x[[id.var]], levels = unique(x[[id.var]]))
    }
    # Remove rownames to control row ordering
    rownames(x) <- NULL
    # Initialise plots list
    plots <- list()
    # Store forest plot
    forest_plot <- .add_forest_main_plot(
        x, effect.var, ci.lower.var, ci.upper.var,
        id.var, facet.by, color.by, ci.exists
    )
    plots <- c(plots, forest_plot)
    # Construct CI label
    if( "CI" %in% label.by ){
        x$CI <- paste0(
            round(x[[effect.var]], 2), " (",
            round(x[[ci.lower.var]], 2), "\u2014",
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
    # Make label plots
    label_plots <- .add_forest_lab_plots(
        x, effect.var, id.var, label.by, color.by
    )
    # For non-empty label plot lists
    if( length(label_plots) != 0 ){
        # Wrap and add label plots
        plots <- c(plots, wrap_plots(label_plots))
    }
    return(plots)
}

#' @importFrom ggplot2 ggplot aes geom_vline geom_point geom_errorbar geom_text
#'   annotate coord_cartesian theme_bw theme element_blank scale_x_continuous
#'   expansion position_dodge position_dodge2 facet_wrap
.add_forest_main_plot <- function(x, effect.var, ci.lower.var, ci.upper.var,
    id.var, facet.by, color.by, ci.exists){
    # Find x-axis limits for forest plot
    if( ci.exists ){
        lim <- max(abs(x[[ci.lower.var]]), x[[ci.upper.var]], na.rm = TRUE)
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
        geom_point(position = position_dodge(width = 0.75)) +
        coord_cartesian(
            xlim = c(-lim, lim), ylim = c(x[[id.var]][1], NA), clip = "off") +
        theme_bw() +
        theme(axis.title.y = element_blank())
    # If CI is defined
    if( ci.exists ){
        # Add errorbars
        p <- p + geom_errorbar(
            aes(xmin = .data[[ci.lower.var]], xmax = .data[[ci.upper.var]]),
            orientation = "y", width = 1e-2 * nrow(x),
            position = position_dodge(width = 0.75)
        )
    }
    if( !is.null(facet.by )){
        p <- p + facet_wrap(as.formula(paste("~", facet.by)))
    }
    return(p)
}

#' @importFrom ggplot2 geom_text annotate expansion
.add_forest_lab_plots <- function(x, effect.var, id.var, label.by, color.by){

    # Initialise label plot list
    label_plots <- list()
    label_size <- .nonlinear_textsize(nrow(x))
    ann_ypos <- -nrow(x) / 30

    if( !is.null(color.by) && nrow(x) > length(unique(x[[id.var]])) ){
        p <- ggplot(x, aes(x = .data[[effect.var]], y = .data[[id.var]],
            colour = .data[[color.by]]))
        # Adjust annotation y-position for grouped rows
        ann_ypos <- ann_ypos / length(unique(x[[color.by]]))
    }else{
        p <- ggplot(x, aes(x = .data[[effect.var]], y = .data[[id.var]]))
    }
    # Iterate over label.by terms
    for( i in seq_along(label.by) ){
        # Retrieve current label
        lab <- label.by[i]
        # Make plot for current label
        label_plots[[i]] <- p +
            geom_text(
                x = 0, aes(y = .data[[id.var]], label = .data[[lab]]),
                hjust = 0, position = position_dodge2(width = 0.75),
                size = label_size, show.legend = FALSE) +
            annotate("text", x = 0, y = ann_ypos, label = lab, hjust = 0) +
            scale_x_continuous(expand = expansion(mult = c(0, 0))) +
            coord_cartesian(
                xlim = c(0, 1), ylim = c(x[[id.var]][1], NA), clip = "off") +
            theme(
                axis.title = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                panel.background = element_blank(),
                panel.grid = element_blank())
    }
    return(label_plots)
}

.nonlinear_textsize <- function(n, min.size = 2, max.size = 5) {
    # Scale size inversely but bounded between min_size and max_size
    size <- max.size - (log10(n) * (max.size - min.size) / log10(100))
    # Clamp between min and max
    size <- pmax(pmin(size, max.size), min.size)
    return(size)
}
