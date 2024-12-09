#' Plotting abundance data
#'
#' \code{plotAbundance()} creates a barplot of feature abundances, typically
#' used to visualize the relative abundance of features at a specific taxonomy
#' rank.
#'
#' It is recommended to handle subsetting, agglomeration, and transformation
#' outside this function. However, agglomeration and relative transformation
#' can be applied using the \code{group} and \code{as.relative} parameters,
#' respectively. If one of the \code{TAXONOMY_RANKS} is selected via
#' \code{group}, \code{mia::agglomerateByRank()} is used, otherwise
#' \code{agglomerateByVariable()} is applied.
#'
#' @param x a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#'
#' @param assay.type \code{Character scalar} value defining which assay data to
#' use. (Default: \code{"relabundance"})
#'
#' @param assay_name Deprecate. Use \code{assay.type} instead.
#'
#' @param layout \code{Character scalar}. Either \dQuote{bar} or \dQuote{point}.
#'
#' @param ... additional parameters for plotting.
#' \itemize{
#'   \item \code{group}: \code{Character scalar}. Specifies the group for
#'   agglomeration. Must be a value from \code{colnames(rowData(x))}. If
#'   \code{NULL}, agglomeration is not applied. (Default: \code{NULL})
#'
#'   \item \code{as.relative}: \code{Character scalar}. Should the relative
#'   values be calculated? (Default: \code{FALSE})
#'
#'   \item \code{col.var}: \code{Character scalar}. Selects a column from
#'   \code{colData} to be plotted below the abundance plot.
#'   Continuous numeric values will be plotted as point, whereas factors and
#'   character will be plotted as colour-code bar. (Default: \code{NULL})
#'
#'   \item \code{order.row.by}: \code{Character scalar}. How to order abundance
#'   value. By name (\code{"name"}) for sorting the taxonomic labels
#'   alphabetically, by abundance (\code{"abund"}) to sort by abundance
#'   values or by a reverse order of
#'   abundance values (\code{"revabund"}). (Default: \code{"name"})
#'
#'   \item \code{row.levels}: \code{Character vector}. Specifies order of rows
#'   in a plot. Can be combined with \code{order.row.by} to control order
#'   of only certain rows. If \code{NULL}, the order follows
#'   \code{order.row.by}. (Default: \code{NULL})
#'
#'   \item \code{order.col.by}: \code{Character scalar}. from the chosen rank of
#'   abundance data or from \code{colData} to select values to order the
#'   abundance plot by. (Default: \code{NULL})
#'
#'   \item \code{col.levels}: \code{Character vector}. Specifies order of
#'   columns in a plot. Can be combined with \code{order.col.by} to control
#'   order of only certain columns. If \code{NULL}, the order follows
#'   \code{order.col.by}. (Default: \code{NULL})
#'
#'   \item \code{decreasing}: \code{Logical scalar}. If the \code{order.col.by}
#'   is defined and the values are numeric, should the values used to order in
#'   decreasing or increasing fashion? (Default: \code{FALSE})
#'
#'   \item \code{facet.rows}: \code{Logical scalar}. Should the rows in the
#'   plot be spitted into facets? (Default: \code{FALSE})
#'
#'   \item \code{facet.cols}: \code{Logical scalar}. Should the columns in the
#'   plot be spitted into facets? (Default: \code{FALSE})
#'
#'   \item \code{ncol}: \code{Numeric scalar}. if facets are applied,
#'   \code{ncol} defines many columns should be for plotting the different
#'   facets. (Default: \code{2})
#'
#'   \item \code{scales} \code{Character scalar}. Defines the behavior of the
#'   scales of each facet. The value is passed into
#'   \code{\link[ggplot2:facet_wrap]{facet_wrap}}. (Default: \code{"fixed"})
#' }
#' See \code{\link{mia-plot-args}} for more details i.e. call
#' \code{help("mia-plot-args")}
#'
#' @return
#' a \code{\link[ggplot2:ggplot]{ggplot}} object or list of two
#' \code{\link[ggplot2:ggplot]{ggplot}} objects, if `col.var` are added to
#' the plot.
#'
#' @name plotAbundance
#'
#' @examples
#' data(GlobalPatterns, package="mia")
#' tse <- GlobalPatterns
#'
#' # If rank is set to NULL (default), agglomeration is not done. However, note
#' # that there is maximum number of rows that can be plotted. That is why
#' # we take sample from the data.
#' set.seed(26348)
#' sample <- sample(rownames(tse), 20)
#' tse_sub <- tse[sample, ]
#' # Apply relative transformation
#' tse_sub <- transformAssay(tse_sub, method = "relabundance")
#' plotAbundance(tse_sub, assay.type = "relabundance")
#'
#' # Plotting counts using the first taxonomic rank as default
#' plotAbundance(
#'     tse, assay.type="counts", group = "Phylum") +
#'     labs(y="Counts")
#'
#' # Using "Phylum" as rank. Apply relative transformation to "counts" assay.
#' plotAbundance(
#'     tse, assay.type="counts", group = "Phylum", add_legend = FALSE,
#'     as.relative = TRUE)
#'
#' # Apply relative transform
#' tse <- transformAssay(tse, method = "relabundance")
#'
#' # A feature from colData or taxon from chosen rank can be used for ordering
#' # samples.
#' plotAbundance(
#'     tse, assay.type="relabundance", group = "Phylum",
#'     order.col.by = "Bacteroidetes")
#'
#' # col.var from colData can be plotted together with abundance plot.
#' # Returned object is a list that includes two plot; other visualizes
#' ## abundance other col.var.
#' plot <- plotAbundance(
#'     tse, assay.type = "relabundance", group = "Phylum",
#'     col.var = "SampleType")
#' \donttest{
#' # These two plots can be combined with wrap_plots function from patchwork
#' # package
#' library(patchwork)
#' wrap_plots(plot, ncol = 1, heights = c(0.95, 0.05))
#' }
#'
#' # Same plot as above but showing sample IDs as labels for the x axis on the
#' # top plot. Moreover, we use facets.
#' plot <- plotAbundance(
#'     tse, assay.type = "relabundance",
#'     group = "Phylum", col.var = "SampleType", add.legend = FALSE,
#'     add.x.text = TRUE, facet.cols = TRUE, scales = "free_x") +
#'     theme(axis.text.x = element_text(angle = 90))
#' plot
#'
#' # Compositional barplot with top 5 taxa and samples sorted by
#' # "Bacteroidetes"
#'
#' # Getting top taxa on a Phylum level
#' tse <- transformAssay(tse, method = "relabundance")
#' tse_phylum <- agglomerateByRank(tse, rank = "Phylum")
#' top_taxa <- getTop(tse_phylum, top = 5, assay.type = "relabundance")
#'
#' # Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
#' phylum_renamed <- lapply(rowData(tse)$Phylum, function(x){
#'     if (x %in% top_taxa) {x} else {"Other"}})
#' rowData(tse)$Phylum <- as.character(phylum_renamed)
#'
#' # Compositional barplot
#' plotAbundance(
#'     tse, assay.type="relabundance", group = "Phylum",
#'     order.row.by="abund", order.col.by = "Bacteroidetes")
NULL

#' @rdname plotAbundance
#' @importFrom ggplot2 facet_wrap
#' @export
setMethod("plotAbundance", signature = c("SummarizedExperiment"), function(
        x, assay.type = assay_name, assay_name = "counts", layout = "bar", ...){
    #
    .check_abundance_input(x, assay.type, layout, ...)
    # Get the abundance data to be plotted. Agglomerate and apply relative
    # transformation if specified.
    abund_data <- .get_abundance_data(x, assay.type, ...)
    group <- attr(abund_data, "group")
    # If the data is paired, ensure that all time points have same sample
    # set, i.e., each patient has all the time points.
    abund_data <- .add_paired_samples(abund_data, ...)
    # Order rows and columns
    abund_data <- .order_abundance_rows(abund_data, ...)
    abund_data <- .order_abundance_cols(abund_data, ...)
    # Create the main plot
    plot_out <- .abund_plotter(
        abund_data, colour_by = group, layout = layout, ...)
    # If user wants to incorporate sample information, add info as an own
    # plot or use facets
    plot_out <- .abund_plotter_incorporate_metadata(plot_out, abund_data, ...)
    return(plot_out)
    }
)

################################ HELP FUNCTIONS ################################
################################################################################
# Data handlers

.check_abundance_input <- function(
        x, assay.type, layout, order.col.by = order_sample_by,
        order_sample_by = NULL, col.var = features, features = NULL, ...){
    #
    if( nrow(x) == 0L ){
        stop("No data to plot. nrow(x) == 0L.", call. = FALSE)
    }
    .check_assay_present(assay.type, x)
    if( !(.is_a_string(layout) && layout %in% c("bar", "point")) ){
        stop("'layout' must be 'bar' or 'point'.", call. = FALSE)
    }
    if( !(is.null(order.col.by) || .is_a_string(order.col.by) ) ){
        stop("'order.col.by' must specify a single character value.",
            call. = FALSE)
    }
    if( !(is.null(col.var) || (is.character(col.var) &&
            all(col.var %in% colnames(colData(x))))) ){
        stop("'col.var' must specify a column from colData(x).",
            call. = FALSE)
    }
    # Check that all the colnames are unique in colData. The functions assume
    # that columns have unique names. Moreover, the data must not have special
    # names that we use here in these functions.
    not_allowed <- c("X", "Y", "colour_by", assay.type)
    if( any(colnames(colData(x)) %in% not_allowed) ){
        stop("colData(x) includes colnames that are not supported. Modify ",
            "the names.\nThe following names are not allowed: '",
            paste0(not_allowed, collapse = "', '"), "'", call. = FALSE)
    }
    # Leave order.col.by out of the comparison, if it specifies a feature
    order.col.by <- if( !is.null(order.col.by) && order.col.by %in%
        colnames(colData(x))) order.col.by else NULL
    all_vars <- unique(c(order.col.by, col.var))
    if( sum(colnames(colData(x)) %in% all_vars) != length(all_vars) ){
        stop("colData(x) must have unique colnames.", call. = FALSE)
    }
    return(NULL)
}

# This function wrangles the data to long format. It takes care of
# agglomeration and transformation. The outut is a single tibble with all the
# whole dataset for plotting.
.get_abundance_data <- function(
        x, assay.type, group = rank, rank = NULL,
        as.relative = use_relative, use_relative = FALSE, ...){
    # Input check
    if( !(is.null(group) || (
        .is_non_empty_string(group) && group %in% colnames(rowData(x)) )) ){
        stop("'group' must be specify a name of a column from rowData or ",
            "NULL.", call. = FALSE)
    }
    if(!.is_a_bool(as.relative)){
        stop("'as.relative' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    # Agglomerate data if user has specified
    if (!is.null(group) && group %in% taxonomyRanks(x)) {
        x <- agglomerateByRank(x, group, ...)
        # or factor that is specified by user
    } else if (!is.null(group)) {
        x <- agglomerateByVariable(x, by = "rows", f = group, ...)
    }
    # At this point, we can check how many rows there are to plot. In practice,
    # there is a limit how many rows we can plot. If there are too many, it is
    # impossible to read the plot. Moreover, the plotting takes excessive
    # amount of time. The good limit might be somewhere around 50, but it
    # might be better to have higher maximum limit so we do not limit too much.
    max_num <- 500
    if( nrow(x) > max_num ){
        stop("The data contains more than ", max_num, " rows. The abundance ",
            "plot cannot be created. Consider subsetting/agglomeration. ",
            "(Check 'group' parameter)", call. = FALSE)
    }
    # If user wants to calculate relative abundances, apply relative transform
    # and use relative assay instead of the original assay in plotting.
    if( as.relative ){
        temp_name <- "temporary_relative_abundance"
        x <- transformAssay(
            x, assay.type = assay.type, method = "relabundance",
            name = temp_name)
        assay.type <- temp_name
    }
    # Samples must have names. In theory, TreeSE can include columns without
    # names. If that is the case, add names.
    if( is.null(colnames(x)) ){
        colnames(x) <- paste0("Sample", seq_len(ncol(x)))
    }
    # Melt TreeSE
    df <- meltSE(
        x, assay.type = assay.type, row.name = "colour_by", col.name = "X",
        add.col = TRUE)
    # Add correct column name for abundance values
    colnames(df)[ colnames(df) == assay.type ] <- "Y"
    # Add group info to attributes
    attr(df, "group") <- ifelse(!is.null(group), group, "Feature")
    return(df)
}

# The paired option is for plotting the abundances so that we can compare
# time points, i.e., with facets where each row is different time point.
# However, if all the time points do not have all the samples, the patients
# are not aligned correctly (samples from certain patient are below each other).
# This function ensures that all the time points have all the patients so that
# comparison is possible.
#' @importFrom dplyr group_by summarize pull select distinct mutate
#'     row_number ungroup across
#' @importFrom tidyr complete
.add_paired_samples <- function(
        df, paired = FALSE, order.col.by = order_sample_by,
        order_sample_by = NULL, col.var = features, features = NULL, ...){
    # To disable "no visible binding for global variable" message in cmdcheck
    colour_by <- count <- X <- NULL
    #
    if(!.is_a_bool(paired)){
        stop("'paired' must be TRUE or FALSE.", call. = FALSE)
    }
    # When paired is specified, order.col.by must be a single variable name from
    # colData
    if( paired && !(.is_a_string(order.col.by) &&
            order.col.by %in% colnames(df)) ){
        stop("When 'paired=TRUE', 'order.col.by' must specify single ",
            "variable from colData(x).", call. = FALSE)
    }
    # When paired is specified, also col.data must be a single variable name
    # from colData
    if( paired && !(all(col.var %in% colnames(df))) ){
        stop("When 'paired=TRUE', 'col.var' must specify single ",
            "variable from colData(x).", call. = FALSE)
    }
    #
    # If the data is paired, and some data is missing from the repeated time
    # points, add the samples as missing.
    # Generate all combinations of sample_type and time_point
    if( paired && !is.null(order.col.by) && !is.null(col.var) ){
        # Calculate how many times each patient-time point pair is present.
        # They must be only once (or none). If they are multiple times, the
        # data is not correctly paired.
        num_pairs <- df |>
            group_by(
                across(all_of(col.var)), .data[[order.col.by]], colour_by) |>
            summarize(count = n(), .groups = "drop") |>
            pull(count)
        if (any(num_pairs > 1)) {
            stop("Data appears to contain multiple samples for some ",
                "combinations of 'col.var' and 'order.col.by'. Ensure that ",
                "each combination corresponds to a unique sample.",
                call. = FALSE)
        }
        # Get all the time point / patient combinations for each feature
        sample_pairs <- df |>
            select(all_of(col.var), all_of(order.col.by), colour_by) |>
            distinct() |>
            complete(!!!syms(col.var), .data[[order.col.by]], colour_by)
        # Join with the original data, filling missing values with NA
        df <- sample_pairs |>
            dplyr::left_join(df, by = c(order.col.by, col.var, "colour_by"))
        # Now we have a dataset that includes all patients for each timepoint.
        # Add arbitrary sample names for those samples that were added.
        df <- df |>
            group_by(colour_by) |>
            mutate(X = ifelse(is.na(as.character(X)),
                paste0("added_", row_number()), as.character(X))) |>
            ungroup() |>
            mutate(X = factor(X))
    }
    return(df)
}

# This function modifies factor of rows/features to follow the user-specified
# order.
#' @importFrom dplyr group_by summarise arrange desc distinct pull
#' @importFrom S4Vectors unfactor
.order_abundance_rows <- function(
        df, order.row.by = order_rank_by, order_rank_by = "name",
        row.levels = NULL, order.col.by = order_sample_by,
        order_sample_by = NULL, ...){
    # To disable "no visible binding for global variable" message in cmdcheck
    colour_by <- Y <- mean_abundance <- NULL
    #
    correct <- .is_a_string(order.row.by) && order.row.by %in%
        c("name", "abund", "revabund")
    if( !correct ){
        stop("'order.row.by' must be 'name', 'abund' or 'revabund'.",
            call. = FALSE)
    }
    if( !(is.null(row.levels) || is.character(row.levels)) ){
        stop("'row.levels' must include all rows.", call. = FALSE)
    }
    # The ordering factor must be found from colData or be one of the rows
    is_coldata <- .is_a_string(order.col.by) && order.col.by %in% colnames(df)
    is_feat <- .is_a_string(order.col.by) && order.col.by %in% df$colour_by
    if( !(is.null(order.col.by) || is_coldata || is_feat) ){
        stop("'order.col.by' must be a variable from colData(x) or a name ",
            "of a row.", call. = FALSE)
    }
    #
    # If user specified levels to use, we get those levels and combine them with
    # rownames so that user do not have to specify all names
    if( !is.null(row.levels) ){
        row.levels <- union(row.levels, as.character(df$colour_by))
    }
    # Order columns and rows alphabetically by default
    if( is.null(row.levels) && order.row.by == "name" ){
        row.levels <- sort(unique(unfactor(df$colour_by)))
    }
    # Get levels based on abundance
    if( is.null(row.levels) && order.row.by %in% c("abund", "revabund") ){
        row.levels <- df |>
            group_by(colour_by) |>
            summarise(mean_abundance = mean(Y, na.rm = TRUE)) |>
            # Either sort based on increasing or decreasing order
            arrange(if (order.row.by == "abund") desc(mean_abundance) else
                mean_abundance) |>
            distinct(colour_by) |>
            pull(colour_by)
    }
    # If user wants to order columns based on abundance of certain taxa, the
    # taxa will be added on top of the figure
    if( is_feat ){
        row.levels <- unique(c(order.col.by, as.character(row.levels)))
    }
    # Apply the ordering
    df$colour_by <- factor(df$colour_by, levels = row.levels)
    return(df)
}

# This function modifies factor of columns/samples to follow the user-specified
# order.
#' @importFrom dplyr filter arrange desc distinct pull
.order_abundance_cols <- function(
        df, order.col.by = order_sample_by, order_sample_by = NULL,
        col.levels = NULL, decreasing = TRUE, ...){
    # To disable "no visible binding for global variable" message in cmdcheck
    X <- Y <- colour_by <-  NULL
    # The ordering factor must be found from colData or be one of the rows
    is_coldata <- .is_a_string(order.col.by) && order.col.by %in% colnames(df)
    is_feat <- .is_a_string(order.col.by) && order.col.by %in% df$colour_by
    if( !(is.null(order.col.by) || is_coldata || is_feat) ){
        stop("'order.col.by' must be a variable from colData(x) or a name ",
            "of a row.", call. = FALSE)
    }
    if( !(is.null(col.levels) || is.character(col.levels)) ){
        stop("col.levels' must include all columns.", call. = FALSE)
    }
    if( !.is_a_bool(decreasing) ){
        stop("'decreasing' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    # If user specified levels to use, we get those levels and combine them with
    # rownames so that user do not have to specify all names
    if( !is.null(col.levels) ){
        col.levels <- union(col.levels, as.character(df$X))
    }
    # If column from colData was specified, give the order of samples based on
    # this variable
    if( is.null(col.levels) && is_coldata ){
        col.levels <- df |>
            arrange(if (decreasing) .data[[order.col.by]] else
                desc(.data[[order.col.by]]) ) |>
            distinct(X) |>
            pull(X)
    }
    # Filter for the specified feature, arrange the dataframe based on the
    # specified column and direction, and then pull unique X values
    if( is.null(col.levels) && is_feat ){
        col.levels <- df |>
            filter(colour_by == order.col.by) |>
            arrange(if (decreasing) desc(Y) else Y) |>
            distinct(X) |>
            pull(X)
    }
    # If column ordering was not specified, order alphabetically
    if( is.null(col.levels) ){
        col.levels <- sort(unique(as.character(df$X)))
    }
    # Apply the ordering
    df$X <- factor(df$X, levels = col.levels)
    return(df)
}

################################################################################
# Abundance plotters

# THis function creates the main abundance plot
#' @importFrom ggplot2 ggplot theme_classic geom_point geom_bar coord_flip
#'   scale_y_continuous
.abund_plotter <- function(
        object,
        xlab = "Samples",
        ylab = paste0(ifelse(as.relative, "Rel. ", ""),"Abundance"),
        colour_by = NULL,
        layout = "bar",
        flipped = FALSE,
        add_legend = TRUE,
        add_x_text = add.x.text,
        add.x.text = FALSE,
        add_border = add.border,
        add.border = NULL,
        bar_alpha = bar.alpha,
        bar.alpha = 0.65,
        point_alpha = point.alpha,
        point.alpha = 1,
        point_size = point.size,
        point.size = 2,
        as.relative = use_relative,
        use_relative = FALSE,
        ...
        ){
    # To disable "no visible binding for global variable" message in cmdcheck
    X <- Y <- NULL
    # Start plotting. From barplot, we exclude 0 values. As we use borders by
    # default, 0 values get also borders which looks like they have also
    # abundance.
    plot_out <- ggplot(object, aes(
        x = X,
        y = ifelse(Y == 0 & layout == "bar", NA, Y)
        )) +
        xlab(xlab) +
        ylab(ylab)
    # Either bar or point plot
    if( layout == "bar" ){
        # Get arguments for bar plot
        abund_out <- .get_bar_args(
            colour_by, alpha = bar_alpha, add_border = add_border,
            n = length(unique(object$X)))
        # Create a bar plot
        plot_out <- plot_out +
            do.call(geom_bar, c(abund_out$args, list(stat = "identity"))) +
            scale_y_continuous(expand = c(0,0))

    } else if( layout == "point" ){
        # Get arguments for point plot
        abund_out <- .get_point_args(
            colour_by, shape_by = NULL, size_by = NULL, alpha = point_alpha,
        size = point_size)
        abund_out$border <- TRUE
        # Create a point plot
        plot_out <- plot_out +
            do.call(geom_point, abund_out$args)
    }
    # Adjust point colours
    if( !is.null(colour_by) ){
        if( abund_out$border ){
            # Resolve the colour for the line colours
            plot_out <- .resolve_plot_colours(
                plot_out, object$colour_by, colour_by, fill = FALSE)
        }
        # Resolve the color for fill
        plot_out <- .resolve_plot_colours(
            plot_out, object$colour_by, colour_by, fill = TRUE)
    }
    # Adjust theme
    plot_out <- plot_out +
        theme_classic()
    # Remove legend if specified
    plot_out <- .add_legend(plot_out, add_legend)
    # Flip the plot if specified
    plot_out <- .flip_plot(plot_out, flipped, add_x_text)
    return(plot_out)
}

# This function takes care of incorporating sample metadata to plot. It is added
# either as facets or unique plot. Moreover, the function has also functinality
# to split rows to unique facets.
#' @importFrom dplyr select all_of distinct arrange select
#' @importFrom stats formula
.abund_plotter_incorporate_metadata <- function(
        plot_out, df, col.var = features, features = NULL,
        facet.cols = FALSE, facet.rows = one.facet,
        one.facet = one_facet, one_facet = FALSE, ncol = 2, scales = "fixed",
        ...){
    # To disable "no visible binding for global variable" message in cmdcheck
    X <- NULL
    #
    if( !(is.null(col.var) || (is.character(col.var) &&
            all(col.var %in% colnames(df)))) ){
        stop("'col.var' must specify columns from colData(x).", call. = FALSE)
    }
    if(!.is_a_bool(facet.cols)){
        stop("'facet.cols' must be TRUE or FALSE.", call. = FALSE)
    }
    if(!.is_a_bool(facet.rows)){
        stop("'facet.rows' must be TRUE or FALSE.", call. = FALSE)
    }
    if( sum(c(facet.rows, facet.cols)) == 2L ){
        stop("'Both 'facet.rows' and 'facet.cols' cannot be TRUE.",
            call. = FALSE)
    }
    if( !(.is_an_integer(ncol) && ncol >= 1) ){
        stop("'ncol' must be an integer value greater or equal to 1.",
            call. = FALSE)
    }
    if( !(.is_a_string(scales) && scales %in%
            c("fixed", "free", "free_x", "free_y")) ){
        stop("'scales' must be 'fixed', 'free', 'free_x' or 'free_y.",
            call. = FALSE)
    }
    #
    # facet.rows is disabled if sample metadata is plotted
    facet.rows <- if(!is.null(col.var)) FALSE else facet.rows
    # Whether to split the main plot to multiple facets by rows, i.e., each
    # taxa has unique facet.
    if( facet.rows ){
        plot_out <- plot_out +
            facet_wrap(~colour_by, ncol = ncol, scales = scales)
    }
    # Whether to facet the plot so that columns are facetted based on sample
    # metadata. This is for single metadata variable.
    if( length(col.var) == 1L && facet.cols ){
        plot_out <- plot_out + facet_wrap(
            formula(paste0("~ `", paste0(col.var, collapse = "` + `"),"`")),
            ncol = ncol,
            scales = scales)
    }
    # This is also for facetting based on sample metadata, however, this allows
    # user to facet columns based on multiple sample metadata variables, e.g.
    # time point and sample type.
    if( length(col.var) > 1L && facet.cols ){
        .require_package("ggh4x")
        plot_out <- plot_out + ggh4x::facet_nested(
            formula(paste0("~ `", paste0(col.var, collapse = "` + `"), "`")),
            scales = scales)
    }
    # If user do not want to create facets from the sample metadata, create
    # unique plots from metadata variables.
    if( !is.null(col.var) && !facet.cols ){
        # Select only sample metadata. Get it in same order that the samples
        # are. After this we have only col.var columns, and metadata includes
        # as many rows as there are samples.
        metadata <- df |>
            select(X, all_of(col.var)) |>
            distinct() |>
            arrange(X) |>
            select(-X)
        # Create a plot from metadata variables
        plot_feature_out <- .features_plotter(metadata, col.var, ...)
        # Add the metadata plots to a list with main plot
        plot_out <- c(list(abundance = plot_out), plot_feature_out)
        names(plot_out) <- c("abundance", col.var)
    }
    return(plot_out)
}

# This function takes care that each sample metadata variable is plotted as
# unique plot.
.features_plotter <- function(
        features_data,
        xlab = NULL,
        flipped = FALSE,
        add_legend = add.legend,
        add.legend = TRUE,
        add_x_text = add.x.text,
        add.x.text = FALSE,
        point_alpha = point.alpha,
        point.alpha = 1,
        point_size = point.size,
        point.size = 2,
        ...){
    # Get the name of sample metadata variables that will be plotted
    names <- colnames(features_data)
    # For each variable, create a data.frame that contains sample names,
    # variable name and values of variable
    features_data <- lapply(names, function(col){
        data.frame(
            X = factor(rownames(features_data), rownames(features_data)),
            feature_name = col,
            Y = features_data[[col]])
        })
    names(features_data) <- names
    # Loop through variables and create plot for each variable
    plots_out <- mapply(
        .feature_plotter,
        features_data,
        names(features_data),
        MoreArgs = list(
            xlab = xlab,
            flipped = flipped,
            add_legend = add_legend,
            add_x_text = add_x_text,
            point_alpha = point_alpha,
            point_size = point_size),
        SIMPLIFY = FALSE)
    names(plots_out) <-  names(features_data)
    return(plots_out)
}

# This function creates a plot from single sample metadata variable.
#' @importFrom ggplot2 ggplot aes labs geom_point geom_raster
.feature_plotter <- function(
        feature_data,
        name,
        xlab = "Samples",
        flipped,
        add_legend,
        add_x_text,
        point_alpha,
        point_size){
    # To disable "no visible binding for global variable" message in cmdcheck
    X <- Y <- NULL
    # If the values are factors, use coloring to plot them. This step is to
    # ensure that this functions works both with factors and numeric values.
    if( is.factor(feature_data$Y) ){
        feature_data$colour_by <- feature_data$Y
        feature_data$Y <- ""
        colour_by <- unique(feature_data$feature_name)
    }
    # Start plotting
    plot_out <- ggplot(
        feature_data, aes(x = X, y = Y)) +
        labs(x = xlab, y = name)
    # If there is only one value, i.e., the variable to be plotted was factor
    if( length(unique(feature_data$Y)) == 1L ){
        # Create a bar layout
        feature_out <- .get_bar_args(
            colour_by, alpha = point_alpha, add_border = FALSE)
        plot_out <- plot_out +
            do.call(geom_raster, feature_out$args) +
            scale_y_discrete(expand = c(0,0))
        # Adjust the colour scale and legend title
        plot_out <- .resolve_plot_colours(
            plot_out, feature_data$colour_by, colour_by, fill = TRUE)
        legend_pos <- "bottom"
    } else {
        # If the values are numeric, create a point layout
        feature_out <- .get_point_args(
            NULL, shape_by = NULL, size_by = NULL, alpha = point_alpha,
            size = point_size)
        plot_out <- plot_out +
            do.call(geom_point, feature_out$args)
        legend_pos <- "right"
    }
    # Adjust theme
    plot_out <- plot_out +
        theme_classic()
    # Remove legend if specified, adjust the position
    plot_out <- .add_legend(plot_out, add_legend, legend_pos)
    # Flip the plot if specified
    plot_out <- .flip_plot(plot_out, flipped, add_x_text)
    return(plot_out)
}
