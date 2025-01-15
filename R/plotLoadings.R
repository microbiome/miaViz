#' Plot feature loadings for TreeSummarizedExperiment
#' objects or feature loadings numeric matrix.
#'
#' This function is used after performing a reduction method. If \code{TreeSE}
#' is given it retrieves the feature loadings matrix to plot values.
#' A tree from \code{rowTree} can be added to heatmap layout.
#'
#' @inheritParams plotTree
#'
#' @param dimred \code{Character scalar}.  Determines the reduced dimension to
#' plot.
#'
#' @param layout \code{Character scalar}. Determines the layout of plot. Must be
#' either \code{"barplot"}, \code{"heatmap"}, or \code{"lollipop"}.
#' (Default: \code{"barplot"})
#'
#' @param ncomponents \code{Numeric scalar}. Number of components must be lower
#' or equal to the number of components chosen in the reduction method.
#' (Default: \code{5})
#'
#' @param add.tree \code{Logical scalar}. Whether to add tree to heatmap layout.
#' (Default: \code{FALSE})
#'
#' @param row.var \code{NULL} or \code{Character scalar}. Specifies a
#' variable from \code{rowData} to plot with tree heatmap layout.
#' (Default: \code{NULL})
#'
#' @param ... additional parameters for plotting.
#' \itemize{
#'   \item \code{n}: \code{Integer scalar}. Number of features to be plotted.
#'   Applicable when \code{layout="barplot"}. (Default: \code{10}))
#'
#'   \item \code{absolute.scale}: ("barplot", "lollipop") \code{Logical scalar}.
#'   Specifies whether a barplot or a lollipop plot should be visualized in
#'   absolute scale. (Default: \code{TRUE})
#' }
#'
#' @details
#'
#' These method visualize feature loadings of dimension reduction results.
#' Inspired by the \code{plotASVcircular} method using \code{phyloseq}.
#' \code{TreeSummarizedExperiment} object is expected to have
#' content in \code{reducedDim} slot calculated with standardized methods from
#' \code{mia} or \code{scater} package.
#'
#' @return
#' A \code{ggplot2} object.
#'
#' @name plotLoadings
#' @export
#'
#' @examples
#'
#' library(mia)
#' library(scater)
#' data("GlobalPatterns", package = "mia")
#' tse <- GlobalPatterns
#'
#' # Calculate PCA
#' tse <- agglomerateByPrevalence(tse, rank="Phylum", update.tree = TRUE)
#' tse <- transformAssay(tse, method = "clr", pseudocount = 1)
#' tse <- runPCA(tse, ncomponents = 5, assay.type = "clr")
#'
#' #' # Plotting feature loadings with tree
#' plotLoadings(tse, dimred = "PCA", layout = "heatmap", add.tree = TRUE)
#'
#' # Plotting matrix as a barplot
#' loadings_matrix <- attr(reducedDim(tse, "PCA"), "rotation")
#' plotLoadings(loadings_matrix)
#'
#' # Plotting more features but less components
#' plotLoadings(tse, dimred = "PCA", ncomponents = 2, n = 12)
#'
#' # Plotting matrix as heatmap without tree
#' plotLoadings(loadings_matrix, layout = "heatmap")
#'
#' # Plot with less components
#' plotLoadings(tse, "PCA", layout = "heatmap", ncomponents = 2)
#'
NULL

#' @rdname plotLoadings
#' @export
#' @importFrom SingleCellExperiment reducedDims
setMethod("plotLoadings", signature = c(x = "TreeSummarizedExperiment"),
    function(
        x, dimred, layout = "barplot", ncomponents = 5, tree.name = "phylo",
        row.var = NULL, add.tree = FALSE, ...) {
        # Check that there are reducedDim
        if( length(reducedDims(x)) == 0 ){
            stop("No reducedDims found.", call. = FALSE)
        }
        # Check dimred. It must be either string specifying the name of
        # reducedDim or an index of reducedDim.
        if( !((.is_a_string(dimred) && dimred %in% reducedDimNames(x)) ||
            .is_an_integer(dimred) && dimred > 0 &&
            dimred <= length(reducedDims(x)) ) ){
                stop("'dimred' must be a string or an integer.", call. = FALSE)
        }
        # Check add.tree
        if( !.is_a_bool(add.tree) ){
            stop("'add.tree' must be TRUE or FALSE.", call. = FALSE)
        }
        # Check that tree.name. If user wants to add tree, the tree name must
        # specify a tree
        if(add.tree && !(.is_a_string(tree.name) &&
        tree.name %in% rowTreeNames(x)) ){
            stop(
                "'tree.name' must be a string specifying a rowTree.",
                call. = FALSE)
        }
        # Check row.var
        if( !(is.null(row.var) ||
                (.is_a_string(row.var) && row.var %in% colnames(rowData(x)))) ){
            stop(
                "'row.var' must be NULL or a column from rowData(x).",
                call. = FALSE)
        }
        #
        # Get loadings matrix
        mat <- .get_loadings_matrix(x, dimred, ...)
        if( add.tree && layout == "heatmap" ){
            # Create dataframe for tree plotting
            data_list <- .get_loadings_tree_data(mat, x, tree.name, row.var)
            tree <- data_list[["tree"]]
            mat <- data_list[["loadings"]]
            # Plot tree with feature loadings
            p <- .loadings_tree_plotter(mat, tree, row.var, ...)
        } else {
            # Utilize matrix method to create a plot
            p <- plotLoadings(
                mat, layout = layout, ncomponents = ncomponents, ...)
        }
    return(p)
    }
)

#' @rdname plotLoadings
#' @export
#' @importFrom SingleCellExperiment reducedDims
setMethod("plotLoadings", signature = c(x = "SingleCellExperiment"),
    function(x, dimred, layout = "barplot", ncomponents = 5, ...){
        # Check that there are reducedDim
        if( length(reducedDims(x)) == 0 ){
            stop("No reducedDims found.", call. = FALSE)
        }
        # Check dimred. It must be either string specifying the name of
        # reducedDim or an index of reducedDim.
        if( !((.is_a_string(dimred) && dimred %in% reducedDimNames(x)) ||
                .is_an_integer(dimred) && dimred > 0 &&
                dimred <= length(reducedDims(x)) ) ){
            stop("'dimred' must be a string or an integer.", call. = FALSE)
        }
        # Get loadings matrix
        mat <- .get_loadings_matrix(x, dimred, ...)
        # Utilize matrix method to create a plot
        p <- plotLoadings(
            mat, layout = layout, ncomponents = ncomponents, ...)
        return(p)
    }
)

#' @rdname plotLoadings
#' @export
setMethod("plotLoadings", signature = c(x = "matrix"),
    function(x, layout = "barplot", ncomponents = 5, ...) {
        # Input check
        .check_loadings_matrix(
            x, layout = layout, ncomponents = ncomponents, ...)
        #
        # Get data for plotting
        df <- .get_loadings_plot_data(x, layout, ncomponents, ...)
        # Create a plot
        p <- .plot_loadings(df, layout = layout, ...)
        return(p)
    }
)

################################ HELP FUNCTIONS ################################

# This function fetches loadings matrix from TreeSE object. The loadings
# are fetched from attributes of reducedDim whcih means that the result must be
# first calculated with standardized method.
.get_loadings_matrix <- function(
        x, dimred, loadings.name = c("rotation", "loadings", "species"), ...){
    #
    if( !is.character(loadings.name) ){
        stop("'loadings.name' must be a character value.", call. = FALSE)
    }
    #
    # Get reducedDim
    reddim <- reducedDim(x, dimred)
    # Get loadings matrix.
    attr_names <- names(attributes(reddim))
    attr_name <- attr_names[ attr_names %in% loadings.name ]
    if( length(attr_name) != 1 ) {
        stop("Loadings cannot be found.", call. = FALSE)
    }
    mat <- attr(reddim, attr_name)
    # Ensure that the type is matrix
    mat <- as.matrix(mat)
    return(mat)
}

# This functions checks that loadings matrix is correct
.check_loadings_matrix <- function(mat, layout, ncomponents, n = 10, ...) {
    # Check layout
    if( !(.is_a_string(layout) && layout %in%
            c("barplot", "heatmap", "lollipop")) ){
        stop("'layout' must be 'barplot' or 'heatmap',", call. = FALSE)
    }
    # Check n
    if( !(.is_an_integer(n) && n > 0 && n <= nrow(mat)) ){
        stop(
            "'n' must be a positive integer less than or equal to the total ",
            "number of features.", call. = FALSE)
    }
    # Check ncomponents
    if( !(.is_an_integer(ncomponents) && ncomponents > 0 &&
            ncomponents <= ncol(mat)) ){
        stop(
            "'ncomponents' must be a positive integer less than or equal to ",
            "the total number of components.", call. = FALSE)
    }
    return(NULL)
}

# This function manipulates the loadings data into correct format. The output
# is data.frame in long format directly usable for ggplot.
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
.get_loadings_plot_data <- function(df, layout, ncomponents, n = 10, ...) {
    # Transform into a dataframe
    df <- as.data.frame(df)
    # Keep only the number of components needed
    df <- df[ , seq_len(ncomponents), drop = FALSE]
    # If the layout is barplot, choose top features for each component
    if( layout %in% c("barplot", "lollipop") ){
        res <- lapply(seq_len(ncomponents), .process_component, df = df, n = n)
        # Combine to single data.frame
        res <- do.call(rbind, res)
    } else{
        # For heatmap, the whole data.frame is just converted into long format.
        components <- colnames(df)
        res <- df %>%
            rownames_to_column(var = "Feature") %>%
            pivot_longer(
                cols = components,
                names_to = "PC",
                values_to = "Value")
    }
    # Convert into data.frame
    res <- as.data.frame(res)
    # Check that values are numeric. This is the first time we test that the
    # columns were numeric. Now all the values from columns are in this column.
    if( !is.numeric(res[["Value"]]) ){
        stop("Values must be numeric.", call. = FALSE)
    }
    # Calculate max and min values along with maximum absolute value and sign
    res <- .calculate_max_and_min_for_loadings(res)
    return(res)
}

# This function subsets the data so that it selects top features that have the
# greatest loadings for single component.
.process_component <- function(i, df, n) {
    # Get order of loadings based on absolute value
    ind <- order(-abs(df[[i]]))
    # Get top n values
    ind <- ind[seq_len(n)]
    # Get the sorted data of single PC
    df <- df[ind, i, drop = FALSE]
    # Add PC number to data.frame so that each row can be identified to belong
    # to certain PC after merging the data into single dataset
    df[["PC"]] <- colnames(df)
    # Rename so that the colnames is same for all components
    colnames(df)[[1]] <- "Value"
    # Add rownames
    df[["Feature"]] <- rownames(df)
    return(df)
}

# This function calculates place for +/- sign in barplot/lollipop plot
#' @importFrom dplyr %>% group_by mutate case_when ungroup
.calculate_max_and_min_for_loadings <- function(df){
    # To disable "no visible binding for global variable" message in cmdcheck
    Value <- PC <- NULL
    # Add column that shows the values in absolute scale, and another column
    # showing sign
    df[["Value_abs"]] <- abs(df[["Value"]])
    df[["Sign"]] <- ifelse(
        df[["Value"]] > 0, "+", ifelse(df[["Value"]] < 0, "-", ""))
    # Add maximum values. This is used in scaling and placement of +/- sign
    # in barplot and lollipop plot. In absolute scale, we use the maximum
    # absolute value. In original scale, negative values gets minimum value
    # and positive values maximum. These values are for each PC.
    df <- df %>%
        group_by(PC) %>%
        mutate(
            # Calculate max of abs(Value) and add 10%
            max_scale_abs = max(abs(Value), na.rm = TRUE) +
                0.1 * max(abs(Value), na.rm = TRUE),
            # Calculate max_scale based on the sign of the Value
            max_scale = case_when(
                Value < 0 ~ min(Value, na.rm = TRUE) +
                    0.1 * min(Value, na.rm = TRUE),
                Value > 0 ~ max(Value, na.rm = TRUE) +
                    0.1 * max(Value, na.rm = TRUE),
                TRUE ~ NA_real_
            )
        ) %>%
        ungroup()
    return(df)
}

# This functions plots a data.frame in barplot or heatmap layout.
#' @importFrom ggplot2 geom_tile scale_fill_gradient2
.plot_loadings <- function(df, layout, ...) {
    # To disable "no visible binding for global variable" message in cmdcheck
    PC <- Feature <- Value <- NULL
    # Initialize a plot
    plot_out <- ggplot(df)
    # Either create a heatmap or barplot/lollipop
    if( layout == "heatmap" ){
        plot_out <- plot_out +
            # Create a heatmap
            geom_tile(
                mapping = aes(x = PC, y = Feature, fill = Value),
                position = position_identity()
                ) +
            # Adjust color scale
            scale_fill_gradient2(
                limits = c(-max(abs(df$Value)), max(abs(df$Value))),
                low = "darkblue", mid = "white", high = "darkred"
                )
    } else if( layout %in% c("barplot", "lollipop") ){
        plot_out <- .plot_bar_or_lollipop(plot_out, df, layout, ...)
    }
    # Adjust theme
    plot_out <- plot_out +
        theme_minimal()
    return(plot_out)
}

# This functions creates a barplot or lollipop plot.
#' @importFrom tidytext scale_y_reordered reorder_within
#' @importFrom ggplot2 geom_bar geom_segment geom_point geom_text
.plot_bar_or_lollipop <- function(
        plot_out, df, layout, absolute.scale = TRUE, show.color = TRUE,
        show.sign = FALSE, ...){
    # To disable "no visible binding for global variable" message in cmdcheck
    Sign <- max_scale_abs <- max_scale <- NULL
    #
    if( !.is_a_bool(absolute.scale) ){
        stop("'absolute.scale' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    if( !.is_a_bool(show.color) ){
        stop("'show.color' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    if( !.is_a_bool(show.sign) ){
        stop("'show.sign' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    # Set the variables to use for aesthetics
    value_var <- if (absolute.scale) "Value_abs" else "Value"
    # Set the y aesthetics with reorder_within, making sure 'df' is referenced
    y_aes <- reorder_within(
        df$Feature,
        # Either get values in absolute scale or not
        if(absolute.scale) -df$Value_abs else df$Value,
        df$PC
        )

    # Plot barplot or lollipop
    if (layout == "barplot") {
        # This creates a barplot
        aesthetic <- aes(
            x = !!sym(value_var),
            y = y_aes,
            # User can decide whether the bars are colored based on +/-
            fill = if(show.color) Sign else NULL
            )
        plot_out <- plot_out + geom_bar(mapping = aesthetic, stat = "identity")
    } else if (layout == "lollipop") {
        # This creates a lollipop plot
        plot_out <- plot_out +
            # Add line
            geom_segment(mapping = aes(
                x = 0, xend = !!sym(value_var),
                y = y_aes, yend = y_aes
            )) +
            # Add point at the end of the line to create "lollipop"
            geom_point(mapping = aes(
                x = !!sym(value_var),
                y = y_aes,
                # User can choose whether the point is colored based on sign
                color = if (show.color) Sign else NULL
            ))
    }

    # Add sign labels if needed
    if( show.sign ){
        plot_out <- plot_out + geom_text(aes(
            # This determines where the sign is placed, absolute scale or not
            x = if (absolute.scale) max_scale_abs else max_scale,
            y = y_aes,
            label = Sign,
            fontface = "bold"
        ))
    }

    # Customize the legend for Sign as "Effect"
    if( show.color ) {
        # Get correct function, barplot uses fill, lollipop color
        scale_FUN <- if( layout == "barplot" ) scale_fill_manual else
            scale_color_manual
        # Currently the legend has title that shows the function call and the
        # values shows + or -. Make the legend nicer.
        plot_out <- plot_out +
            scale_FUN(
                name = "Effect",
                values = c("+" = "blue", "-" = "red"),
                labels = c("+" = "positive", "-" = "negative")
            )
    }

    # Final wrangle, set facets and order the data
    plot_out <- plot_out +
        scale_y_reordered() +
        facet_wrap(~PC, scales = "free") +
        labs(x = "Value", y = "Feature")

    return(plot_out)
}

# This function retrieves the data for tree + heatmap plotting. The output
# is a list that includes tree and data.frame in wide format.
#' @importFrom ggtree ggtree
.get_loadings_tree_data <- function(df, x, tree.name, row.var){
    # Check that rownames of loading matrix match with rownames of TreeSE. It
    # might be that TreeSE is updated after calculating the reduced dimension
    # which is why rownames do not match.
    all_match <- all(rownames(x) %in% rownames(df)) &&
        all(rownames(df) %in% rownames(x))
    if( !all_match ){
        stop(
            "Features of loading matrix do not match with rownames(x)",
            call. = FALSE)
    }
    # Sort the loading matrix
    df <- df[match(rownames(x), rownames(df)), ]
    # Convert loadings matrix to data.frame
    df <- as.data.frame(df)
    # Retrieve rowTree
    phylo <- rowTree(x, tree.name)
    # Subset data based on the tree
    ind <- rowLinks(x)[["whichTree"]] == tree.name
    if( any(!ind) ){
        warning("Data is subsetted.", call. = FALSE)
        # Subset both TreeSE and loadings dfrix
        x <- x[ind, ]
        df <- df[ind, ]
    }
    # Add feature to column
    df[["Feature"]] <- rownames(df)
    # Instead of rownames, user can also specify a column from rowData to be
    # plotted
    if( !is.null(row.var) ){
        df[["Feature"]] <- rowData(x)[[row.var]]
    }
    # Check that there are not too many features to plot. If there are too many
    # rank values (or rownames), it is not possible to plot.
    if( length(unique(df[["Feature"]])) > 100 ){
        stop("Too many features to plot.", call. = FALSE)
    }
    # Add rowlinks to data
    rownames(df) <- rowLinks(x)[["nodeLab"]]
    res <- list(tree = phylo, loadings = df)
    return(res)
}

# This function is for plotting tree with heatmap. It utilizes ggtree package.
#' @importFrom ggtree ggtree gheatmap
#' @importFrom ggnewscale new_scale_fill
#' @importFrom ggplot2 scale_fill_gradient2 scale_fill_viridis_d
.loadings_tree_plotter <- function(
    df, tree, rank, rank.title = ifelse(!is.null(rank), rank, "Feature"),
    ...) {
    # Check rank.title
    if( !.is_a_string(rank.title) ){
        stop("'rank.title' must be a string.", call. = FALSE)
    }
    #
    # Get features
    features <- df[ , colnames(df) %in% c("Feature"), drop = FALSE]
    # Get loadings
    loadings <- df[ , !colnames(df) %in% c("Feature"), drop = FALSE]
    # Create a tree plot
    plot_out <- ggtree(tree, layout = "circular")
    # Add first inner circle (features)
    plot_out <- gheatmap(
        p = plot_out, data = features, width = 0.1, colnames_angle = 90)
    # Adjust color scale for discrete feature values
    plot_out <- plot_out +
        scale_fill_viridis_d(option = "D", name = rank.title)
    # Add outer circles (feature loadings)
    plot_out <- plot_out + new_scale_fill()
    plot_out <- gheatmap(
        p = plot_out, data = loadings, offset = 0.1, width = 0.3,
        colnames_angle = 90
        )
    # Adjust color scale in continuous scale
    plot_out <- plot_out +
        scale_fill_gradient2(
            limits = c(-1, 1),
            low = "darkslateblue", mid = "white", high = "darkred",
            name = "Value"
        )
    return(plot_out)
}
