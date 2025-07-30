#' @name
#' plotBoxplot
#'
#' @title
#' Create boxplot of \code{assay}, \code{rowData} or \code{colData}.
#'
#' @description
#' This methods visualizes abundances or variables from \code{rowData} or
#' \code{colData}.
#'
#' @details
#' A box plot is standard visualization technique to compare numeric values,
#' such as abundance, between categorical values, such as sample groups.
#' \code{plotBoxplot()} streamlines creation of box plots, and it offers
#' multiple options for visualization.
#'
#' @return
#' A \code{ggplot2} object.
#'
#' @param object a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#'
#' @param assay.type \code{NULL} or \code{character scalar}. Specifies the
#' abundace table to plot. (Default: \code{NULL})
#'
#' @param row.var \code{NULL} or \code{character scalar}. Specifies a variable
#' from \code{rowData(x)} to visualize. (Default: \code{NULL})
#'
#' @param col.var \code{NULL} or \code{character scalar} Specifies a variable
#' from \code{colData(x)} to visualize. (Default: \code{NULL})
#'
#' @param x \code{NULL} or \code{character vector}. Specifies a variable
#' from \code{colData(x)} or \code{rowData(x)} to visualize in x axis.
#' (Default: \code{NULL})
#'
#' @param features \code{NULL} or \code{character vector}. If \code{assay.type}
#' is specified, this specifies rows to visualize in different facets. If
#' \code{NULL}, whole data is visualized as a whole. (Default: \code{NULL})
#'
#' @param group.by \code{NULL} or \code{character vector}. Specifies a variable
#' from \code{colData(x)} or \code{rowData(x)} to group observations.
#' (Default: \code{NULL})
#'
#' @param ... Additional parameters for plotting.
#' \itemize{
#'   \item \code{point.offset}: \code{Character scalar}. Utilized method
#'   for offsetting points. The available options include:
#'   \code{"center"}, \code{"compactswarm"}, \code{"hex"}, \code{"square"},
#'   \code{"swarm"}
#'   (see \code{\link[beeswarm:beeswarm]{beeswarm::beeswarm()}} for details),
#'   \code{"frowney"}, \code{"maxout"}, \code{"minout"}, \code{"pseudorandom"},
#'   \code{"quasirandom"},  \code{"smiley"}, \code{"tukey"}, \code{"tukeyDense"}
#'   (see \code{\link[vipor:offsetSingleGroup]{vipor::offsetSingleGroup()}}
#'   for details), \code{"jitter"}, and \code{"none"},
#'   If \code{"none"}, offsetting is not applied. (Default: \code{"jitter"})
#'
#'   \item \code{colour.by}: \code{NULL} or \code{character scalar}. Specifies a
#'   variable from \code{colData(x)} or \code{rowData(x)} which is used to
#'   colour observations. (Default: \code{NULL})
#'
#'   \item \code{fill.by}: \code{NULL} or \code{character scalar}. Specifies a
#'   variable from \code{colData(x)} or \code{rowData(x)} which is used to
#'   colour observations. (Default: \code{NULL})
#'
#'   \item \code{size.by}: \code{NULL} or \code{character scalar}. Specifies a
#'   variable from \code{colData(x)} or \code{rowData(x)} which is used to scale
#'   observation points. (Default: \code{NULL})
#'
#'   \item \code{shape.by}: \code{NULL} or \code{character scalar}. Specifies a
#'   variable from \code{colData(x)} or \code{rowData(x)} which is used to shape
#'   observation points. (Default: \code{NULL})
#'
#'   \item \code{facet.by}: \code{NULL} or \code{character scalar}. Specifies a
#'   variable from \code{colData(x)} or \code{rowData(x)} which is used to facet
#'   or group observations. (Default: \code{NULL})
#'
#'   \item \code{pair.by}: \code{NULL} or \code{character scalar}. Specifies a
#'   variable from \code{colData(x)} which is used to pair observation points.
#'   (Default: \code{NULL})
#'
#'   \item \code{add.change}: \code{Logical scalar}. Whether to visualize chance
#'   of paired observations by the color of line. (Default: \code{FALSE})
#'
#'   \item \code{add.box}: \code{Logical scalar}. Whether to add a boxplot
#'   layout. (Default: \code{TRUE})
#'
#'   \item \code{add.points}: \code{Logical scalar}. Whether to add a point
#'   layout. (Default: \code{TRUE})
#'
#'   \item \code{add.proportion}: \code{Logical scalar}. Whether to add a
#'   barplot layout denoting the proportion of observations above
#'   \code{threshold}. (Default: \code{FALSE})
#'
#'   \item \code{add.threshold}: \code{Logical scalar}. Whether to add a
#'   \code{threshold} as horizontal line when \code{add.proportion = TRUE} is
#'   specified. (Default: \code{TRUE})
#'   
#'   \item \code{add.significance}: \code{Logical scalar}. Whether to
#'   automatically compute statistical significance and add p-value annotations 
#'   to the plot. Ignored if \code{p.value} is provided. (Default: \code{FALSE})
#'
#'   \item \code{mark.significance}: \code{Logical scalar}. Whether to display
#'   the computed or supplied p-values on the plot as significance asterisks. 
#'   Has no effect if neither \code{add.significance} nor #'   \code{p.value} is 
#'   used. (Default: \code{FALSE})
#'   
#'   \item \code{p.value}: \code{NULL} or \code{data.frame} with columns
#'   \code{group1}, \code{group2}, and \code{p}. If supplied, these are used
#'   as user-defined p-values for group comparisons, overriding any computed
#'   significance. (Default: \code{NULL})
#'
#'   \item \code{threshold}: \code{Numeric scalar}. Specifies threshold for the
#'   barplots. (Default: \code{0})
#'
#'   \item \code{jitter.width}: \code{Numeric scalar}. Width of jitter.
#'   (Default: \code{0.3})
#'
#'   \item \code{jitter.height}: \code{Numeric scalar}. Height of jitter.
#'   (Default: \code{0})
#'
#'   \item \code{dodge.width}: \code{Numeric scalar}. Width of dodge. How far
#'   apart the groups are plotted? (Default: \code{0})
#'
#'   \item \code{beeswarm.corral}: \code{Character scalar}. Beeswarm's "corral"
#'   method. Fed to function \code{beeswarm::beeswarm()}.
#'   (Default: \code{"none"})
#'
#'   \item \code{scales}: \code{Character scalar}. Adjust scales of facets.
#'   (Default: \code{"fixed"})
#'
#'   \item \code{box.alpha}: \code{Numeric scalar}. Transparency of the boxplot
#'   layer. (Default: \code{0.5})
#'
#'   \item \code{point.alpha}: \code{Numeric scalar}. Transparency of the point
#'   layer. (Default: \code{0.65})
#'
#'   \item \code{line.alpha}: \code{Numeric scalar}. Transparency of the line
#'   layer. (Default: \code{0.5})
#'
#'   \item \code{point.shape}: \code{Numeric scalar}. Shape of points.
#'   (Default: \code{21})
#'
#'   \item \code{point.size}: \code{Numeric scalar}. Size of points.
#'   (Default: \code{2})
#'
#'   \item \code{point.colour}: \code{Character scalar}. Colour of points.
#'   (Default: \code{"grey70"})
#'
#'   \item \code{linetype}: \code{Numeric scalar}. Type of lines.
#'   (Default: \code{1})
#'
#'   \item \code{linewidth}: \code{Numeric scalar}. Width of lines.
#'   (Default: \code{1})
#'
#'   \item \code{line.colour}: \code{Character scalar}. Colour of lines.
#'   (Default: \code{"grey70"})
#'
#'   \item \code{box.width}: \code{Numeric scalar}. Width of boxes.
#'   (Default: \code{0.75})
#'
#'   \item \code{bar.width}: \code{Numeric scalar}. Width of proportion bars.
#'   By default, it is calculated based so that the width matches with the
#'   width of boxes.
#'   
#'   \item \code{ncol}: \code{Integer scalar}. Number of columns in the 
#'   `facet_wrap` layout. (Default: \code{NULL}, which lets ggplot2 decide the 
#'   optimal layout)
#'
#'   \item \code{nrow}: \code{Integer scalar}. Number of rows in the 
#'   `facet_wrap` layout. (Default: \code{NULL}, which lets ggplot2 decide the 
#'   optimal layout)
#' }
#'
#' @examples
#' \dontrun{
#' # This example shows how to plot Boxplots from intervention studies or a 
#' # similar study design with repeated measures across multiple timepoints
#' # with two or more interventions.
#' # the Kumaraswamy2024 data has four Interventions(A-D) over three timepoints
#' # (summer-autumn-winter)
#'  
#' data(Kumaraswamy2024)
#' tse <- Kumaraswamy2024
#' 
#' # Visualize alpha diversity
#' # We already have computed 'shannon' values.
#' # Between group comparisons
#' plotBoxplot(
#'     tse, col.var = "shannon", 
#'     fill.by = "group", x = "season"
#'     )
#' # with pvalues
#' plotBoxplot(
#'     tse, col.var = "shannon", 
#'     fill.by = "group", x = "season", a
#'     dd.significance = T
#'     )
#' # Within group comparisons
#' # simple plot
#' plotBoxplot(
#'     tse, col.var = "shannon", 
#'     fill.by = "season", x = "group"
#'     )
#' # Link points across timepoints
#' plotBoxplot(
#'     tse, col.var = "shannon", 
#'     fill.by = "season", x = "group", 
#'     pair.by = "subject"
#'     )
#' # Indicate change in diversity across timepoints
#' # First filter tse to have measurements for each timepoint
#' multiple_samples_df <- as.data.frame(colData(tse)) %>% 
#'     group_by(subject) %>% 
#'     filter(n() == 3) %>% 
#'     ungroup()
#' multiple_samples <- multiple_samples_df$sample
#' # subset tse by valid samples
#' tse <- tse[, multiple_samples]
#' 
#' # Plot with change
#' plotBoxplot(
#'     tse, col.var = "shannon", 
#'     fill.by = "season", x = "group", 
#'     pair.by = "subject", add.change = TRUE
#'     )
#' # Add significance
#' plotBoxplot(
#'     tse, col.var = "shannon", 
#'     fill.by = "season", x = "group", 
#'     pair.by = "subject", add.change = TRUE, 
#'     add.significance = TRUE
#'     )
#'     
#' # Visualize relative abundance
#' # For interested features e.g :
#' interest <- c("Bifidobacterium", 
#'     "Dialister", "Veillonella", "Collinsella")
#' 
#' # Between groups
#' plotBoxplot(
#'     tse, assay.type = "relabundance", 
#'     fill.by = "group", x = "season", 
#'     features = interest, facet.by = "rownames", scales = "free"
#'     )
#' # Add significance
#' plotBoxplot(
#'     tse, assay.type = "relabundance", 
#'     fill.by = "group", x = "season", features = c("Bifidobacterium"), 
#'     facet.by = "rownames", scales = "free", 
#'     add.significance = TRUE
#'     )
#'     
#' # Within groups
#' 
#' plotBoxplot(
#'     tse, assay.type = "relabundance", 
#'     fill.by = "season", x = "group", 
#'     features = interest, facet.by = "rownames", scales = "free"
#'     )
#' # Add pvalues
#' plotBoxplot(
#'     tse, assay.type = "relabundance", 
#'     fill.by = "season", x = "group", 
#'     features = c("Bifidobacterium"), facet.by = "rownames", 
#'     add.significance = TRUE
#'     )
#' # Add change
#' plotBoxplot(
#'     tse, assay.type = "relabundance", 
#'     fill.by = "season", x = "group", 
#'     features = c("Bifidobacterium"), 
#'     facet.by = "rownames", add.significance = TRUE, 
#'     pair.by = "subject", add.change = TRUE)
#' 
#' }
#' 
#' # Other Examples
#' 
#' data("Tito2024QMP")
#' tse <- Tito2024QMP
#'
#' tse <- transformAssay(tse, method = "relabundance")
#' tse <- addAlpha(tse, index = "shannon")
#'
#' # Visualize alpha diversity
#' plotBoxplot(tse, col.var = "shannon", x = "diagnosis")
#'
#' # Visualize relative abundance of top features
#' tse <- tse[getTop(tse, 6), ]
#'
#' plotBoxplot(
#'     tse, assay.type = "relabundance",
#'     x = "diagnosis", fill.by = "diagnosis",
#'     features = rownames(tse), facet.by = "rownames"
#' )
#'
#' # Add proportion bar
#' plotBoxplot(
#'     tse, assay.type = "relabundance",
#'     x = "diagnosis", fill.by = "diagnosis",
#'     features = rownames(tse), facet.by = "rownames",
#'     add.proportion = TRUE, threshold = 0.1
#' )
#'
#' # Visualize only with beeswarm
#' plotBoxplot(
#'     tse, assay.type = "relabundance",
#'     x = "diagnosis", group.by = "diagnosis",
#'     colour.by = "colonoscopy",
#'     features = rownames(tse), facet.by = "rownames",
#'     point.offset = "swarm", add.box = FALSE
#' )
#'
#' # Do not add points
#' plotBoxplot(
#'     tse, assay.type = "relabundance",
#'     fill.by = "diagnosis",
#'     features = rownames(tse), facet.by = "rownames",
#'     add.points = FALSE
#' )
#' 
#' # With pre-calculated pvalues
#' df <- meltSE(tse, assay.type = "relabundance", add.col = TRUE)
#' 
#' # Calculate the pvalues
#' library(rstatix)
#' res <- df |> 
#'     group_by(FeatureID) |> 
#'     wilcox_test(relabundance ~ diagnosis) |> 
#'     adjust_pvalue()
#' 
#' # Plot with pvalue dataframe    
#' plotBoxplot(tse, features = rownames(tse), 
#'     x = "diagnosis", assay.type = "relabundance", 
#'     p.value = res, facet.by = "rownames")
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[scater:plotExpression]{scater::plotExpression}}
#'   \item \code{\link[scater:plotRowData]{scater::plotRowData}}
#'   \item \code{\link[scater:plotColData]{scater::plotColData}}
#' }
#'
NULL

#' @rdname plotBoxplot
#' @export
setMethod("plotBoxplot", signature = c(object = "SummarizedExperiment"),
    function(object, assay.type = NULL, row.var = NULL, col.var = NULL,
            x = NULL, features = NULL, group.by = NULL, ...){
        # Add rownames to rowData so that they are available for plotting.
        rowData(object)[["rownames"]] <- rownames(object)
        # Check input
        args <- c(list(
            tse = object, assay.type = assay.type, features = features,
            row.var = row.var, col.var = col.var, x = x, group.by = group.by),
            list(...)
        )
        temp <- do.call(.check_input_for_boxplot, args)
        # Get the data from the object
        df <- do.call(.get_data_for_boxplot, args)
        # Create a boxplot
        p <- .plot_boxplot(df, ...)
        return(p)
    }
)

################################ HELP FUNCTIONS ################################

# This function validates the input for boxplot plotter.
.check_input_for_boxplot <- function(
        tse, assay.type, features, row.var, col.var, x, group.by,
        pair.by = NULL, add.change = FALSE, colour.by = color.by,
        color.by = NULL, fill.by = NULL, size.by = NULL, shape.by = NULL,
        facet.by = NULL, add.box = TRUE, add.points = TRUE,
        add.proportion = FALSE, add.threshold = FALSE, scales = "fixed", 
        p.value = NULL, add.significance = FALSE, 
        mark.significance = FALSE, ...){
    # Either assay.type. row.var or col.var must be specified
    if( sum(c(is.null(assay.type), is.null(row.var), is.null(col.var))) != 2L ){
        stop("Please specify either 'assay.type', 'row.var', or 'col.var'.",
            call. = FALSE)
    }
    # features cannot be specified if row.var or col.var is specified
    if( is.null(assay.type) && !is.null(features) ){
        stop("'features' can only be specified when 'assay.type is ",
            "specified.", call. = FALSE)
    }
    # As features points to rownames, the TreeSE must have rownames and features
    # must match them
    if( !is.null(features) && is.null(rownames(tse)) ){
        stop("'object' must have rownames.", call. = FALSE)
    }
    if( !(is.null(features) ||
            (is.character(features) && all(features %in% rownames(tse)) )) ){
        stop("'features' must be NULL or single character value specifying ",
            " rownames.", call. = FALSE)
    }
    # If assay was specified, check that it is correct.
    if( !is.null(assay.type) ){
        .check_assay_present(assay.type, tse)
    }
    if( !.is_a_bool(add.box) ){
        stop("'add.box' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(add.points) ){
        stop("'add.points' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(add.proportion) ){
        stop("'add.proportion' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(add.threshold) ){
        stop("'add.threshold' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(add.significance) ){
        stop("'add.significance' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(mark.significance) ){
        stop("'mark.significance' must be TRUE or FALSE.", call. = FALSE)
    }
    # Check colData/rowData variables
    temp <- .check_metadata_variable(tse, row.var, row = TRUE)
    temp <- .check_metadata_variable(tse, col.var, col = TRUE)
    temp <- .check_metadata_variable(
        tse, x,
        row = length(c(row.var, assay.type))>0,
        col = length(c(col.var, assay.type))>0
    )
    temp <- .check_metadata_variable(
        tse, group.by,
        row = length(c(row.var, assay.type))>0,
        col = length(c(col.var, assay.type))>0
    )
    temp <- .check_metadata_variable(
        tse, fill.by,
        row = length(c(row.var, assay.type))>0,
        col = length(c(col.var, assay.type))>0
    )
    temp <- .check_metadata_variable(
        tse, colour.by,
        row = length(c(row.var, assay.type))>0,
        col = length(c(col.var, assay.type))>0
    )
    temp <- .check_metadata_variable(
        tse, size.by,
        row = length(c(row.var, assay.type))>0,
        col = length(c(col.var, assay.type))>0
    )
    temp <- .check_metadata_variable(
        tse, shape.by,
        row = length(c(row.var, assay.type))>0,
        col = length(c(col.var, assay.type))>0
    )
    temp <- .check_metadata_variable(
        tse, facet.by,
        row = length(c(row.var, assay.type))>0,
        col = length(c(col.var, assay.type))>0
    )
    temp <- .check_metadata_variable(tse, pair.by, col = TRUE)
    # Check that pairing variables are correct
    if( !.is_a_bool(add.change) ){
        stop("'add.change' must be TRUE or FALSE.", call. = FALSE)
    }
    if( add.change && is.null(pair.by) ){
        stop("When 'add.change' is specified, 'pair.by' must be specified.",
            call. = FALSE)
    }
    # We have to plot points in order to connect them
    if( !is.null(pair.by) && !add.points ){
        stop("When 'pair.by' is specified, 'add.points' must be enabled.",
            call. = FALSE)
    }
    
    # Check that repeated measures are balanced when 'pair.by' is used
    if (!is.null(pair.by)) {
        coldata <- as.data.frame(colData(tse))
        
        sample_counts <- coldata %>%
            group_by(.data[[pair.by]]) %>%
            summarise(n_samples = n(), .groups = "drop")
        
        if (length(unique(sample_counts$n_samples)) > 1) {
            stop(
                "When 'pair.by' is specified, all subjects must have the ",
                "same number of samples.\nRemove 'pair.by' or filter your ",
                "data to include only subjects with balanced repeats.",
                call. = FALSE
            )
        }
    }
    
    # We can either color based on variable or the difference netween paired
    # samples. There cannot be multiple coloring schemes.
    if( add.change && !is.null(colour.by) ){
        stop("Both 'add.change' and 'colour.by' cannot be specified ",
            "simultaneously.", call. = FALSE)
    }
    # x must be character or factor in box plots
    if( !is.null(x) && is.numeric(tse[[x]]) ){
        stop("'x' must specify categorical value.", call. = FALSE)
    }
    # X axes does not line correctly
    # There should be no need to facet and group based on the same variable.
    # It leads to wrong scaling in x-axis which is why we just disable it.
    # Moreover, it cannot match with x-axis value for the same reason.
    g_match <- !is.null(facet.by) && !is.null(group.by) && facet.by == group.by
    f_match <- !is.null(facet.by) && !is.null(fill.by) && facet.by == fill.by
    x_match <- !is.null(facet.by) && !is.null(x) && facet.by == x
    if( g_match || f_match || x_match ){
        stop("'facet.by' must not match with 'x', 'group.by' or 'fill.by'.",
            call. = FALSE)
    }
    # scales is passed to faceting, but as we use it also to adjust the point
    # positions, we have to test that it has correct value.
    vals <- c("fixed", "free", "free_x", "free_y")
    if( !(.is_a_string(scales) && scales %in% vals) ){
        stop("'scales' must be a single character value from the following ",
            "options: '", paste0(vals, collapse = "', '"), "'", call. = FALSE)
    }

    # p.value defines p-values to be plotted. It should be in specific format.
    # It should have p-values and groups along with faceting variable
    # if it was specified.
    if( is.null(x) && !is.null(p.value) ){
        stop("If 'x' is not specified, 'p.value' must be NULL.", call. = FALSE)
    }
    if( !(is.null(p.value) || is.data.frame(p.value)) ){
        stop("'p.value' must be NULL or data.frame.", call. = FALSE)
    }
    cols <- c("group1", "group2")
    p_cols <- c("p", "p.adj", "p.adj.signif")
    if( is.data.frame(p.value) && (
            !all(cols %in% colnames(p.value)) ||
            sum(p_cols %in% colnames(p.value)) == 0L)
        ){
        stop("'p.value' must have columns '", paste0(cols, collapse = "', '"),
            "' and one of the following columns: '",
            paste0(p_cols, collapse = "', '"), "'", call. = FALSE)
    }

    return(NULL)
}

# This function retrieves the data from TreeSE and outputs a data.frame, ready
# for the plotter function.
.get_data_for_boxplot <- function(
        tse, x = NULL, assay.type, features, row.var, col.var, group.by,
        pair.by = NULL, add.change = FALSE,
        colour.by = color.by, color.by = NULL,
        size.by = NULL, shape.by = NULL, facet.by = NULL,
        fill.by = NULL, add.proportion = FALSE, p.value = NULL,
        add.significance = FALSE,
        ...){
    # If assay.type is specified, get melted data
    all_vars <- c(x, group.by, colour.by, size.by, shape.by, facet.by, fill.by)
    if( !is.null(assay.type) ){
        # Specify whether to retrieve data from rowData or colData
        row_vars <- vapply(all_vars, function(x){
            x %in% colnames(rowData(tse))
        }, logical(1L))
        col_vars <- all_vars[ !row_vars ]
        row_vars <- all_vars[ row_vars ]
        #
        df <- meltSE(
            tse, assay.type = assay.type,
            col.var = "id",
            add.col = c(col.var, pair.by, col_vars),
            add.row = c(row.var, row_vars)
        )
    }
    # If row.var was specified, get the data from rowData
    if( !is.null(row.var) ){
        df <- rowData(tse)[, c(row.var, all_vars), drop = FALSE]
    }
    # If col.var was specified, get the data from colData
    if( !is.null(col.var) ){
        df <- colData(tse)[, c(col.var, pair.by, all_vars), drop = FALSE]
    }
    # Check that y-axis is numeric
    if( !is.numeric(df[[c(assay.type, col.var, row.var)]]) ){
        stop("Y-axis must be numeric.", call. = FALSE)
    }
    # Prevalence can be added only if values are non-negative
    is_negative <- any(!is.na(df[[c(assay.type, col.var, row.var)]]) &
        df[[c(assay.type, col.var, row.var)]]<0)
    
    # If user wants to add significance, but p-value was not defined, calculate
    # them.
    if( add.significance && is.null(p.value) ){
        p.value <- .calculate_significance(
            df, c(assay.type, col.var, row.var),
            x, facet.by, fill.by, group.by, pair.by,
            features, ...)
    }
    
    # If features were specified, subset data. The subsetting is done after
    # calculating p-values to ensure that they are calculated for whole data
    # instead of subset so that the correction is done correctly.
    if( !is.null(features) ){
        df <- df[ df[["FeatureID"]] %in% features, , drop = FALSE]
    }

    # If user specified, calculate difference
    difference <- NULL
    if( add.change ){
        df <- .calculate_paired_difference(
            df, x, c(assay.type, col.var, row.var), pair.by, group.by, fill.by,
            facet.by)
        difference <- "difference"
    }
    # If both group.by and x are specified, the groups get them x-axis
    # position based on these both variables in boxplot layer. fill.by works
    # fine without this kind of modification.
    x.box <- group.by
    if( !is.null(x) && !is.null(group.by) ){
        x.box <- "x_box"
        df[[x.box]] <- interaction(df[[x]], df[[group.by]])
    }
    # If x-axis was not specify, specify it to be 0.
    remove.x.axis <- FALSE
    if( is.null(x) ){
        x <- "x_axis"
        df[[x]] <- 0
        remove.x.axis <- TRUE
    }
    # We add jitter to points manually. The problem is that ggplot evaluates
    # jitter for each layer separately when rendering the plot. We could specify
    # seed, but the problem comes from faceting. Even though, we know the
    # points' positions before faceting, they are not the same after faceting.
    # Setting manually the jitter for points is much easier for us.
    df <- .add_fixed_jitterdodge(
        df, x, c(assay.type, col.var, row.var), group.by, fill.by, facet.by,
        ...)
    df <- df |>
        as.data.frame()

    # Add p-values placement
    if( !is.null(p.value) ){
        # Check that p.value data.frame has the correct grouping variables. If
        # group.by or fill.by are specified the comparisons are made based on
        # them. If they are not specified, the comparison are made based on 
        # x axis variable.
        .validate_pvalue(p.value, df, x, group.by, fill.by, facet.by)
        p.value <- .add_p_value_position(p.value, x, 
                                         c(assay.type, col.var, row.var), 
                                         group.by, fill.by, facet.by, df, ...)
    }

    # Add plotting options to attributes of the data.frame. Now the data.frame
    # includes all the information for plotting.
    attributes(df) <- c(
        attributes(df),
        value = c(assay.type, col.var, row.var),
        x = x,
        group.by = group.by,
        pair.by = pair.by,
        add.change = add.change,
        colour.by = colour.by,
        fill.by = fill.by,
        size.by = size.by,
        shape.by = shape.by,
        facet.by = facet.by,
        x.box = x.box,
        difference = difference,
        remove.x.axis = remove.x.axis
    )
    attributes(df)[["p.value"]] <- p.value
    return(df)
}

# This function calculates difference between paired samples.
#' @importFrom dplyr arrange across all_of group_by mutate desc
.calculate_paired_difference <- function(
        df, x, y, pair.by, group.by, fill.by, facet.by){
    # Calculate difference between paired points
    df <- df |>
        as.data.frame() |>
        arrange(across(all_of(c(pair.by, x, group.by, fill.by)))) |>
        group_by(across(all_of(c(pair.by, facet.by)))) |>
        mutate(
            difference = .data[[y]] - dplyr::lag(.data[[y]])
        ) |>
        ungroup()
    # Sort data so that last time point comes first. ggplot gets color from
    # first instance. Otherwise it would be time point 1 -> time point 2,
    # which is NA.
    df <- df |>
        arrange(desc(across(all_of(c("difference", pair.by, x)))))
    return(df)
}

# This function adjust jitter and dodging for points. Jitter means random noise
# for points' positions while dodging means that we separate groups in x-axis.
# This functions works with similar logic than position_jitterdodge. Dodge for
# boxplot is set with ggplot.
#' @importFrom dplyr ungroup
.add_fixed_jitterdodge <- function(
        df, x, y, group.by, fill.by, facet.by,
        jitter.width = 0.3, jitter.height = 0, dodge.width = 0.8,
        point.offset = "jitter", ...){
    if( !.is_a_numeric(jitter.width) ){
        stop("'jitter.width' must be numeric.", call. = FALSE)
    }
    if( !.is_a_numeric(jitter.height) ){
        stop("'jitter.height' must be numeric.", call. = FALSE)
    }
    if( !(.is_a_numeric(dodge.width) && (dodge.width >=0 && dodge.width <=1)) ){
        stop("'dodge.width' must be numeric (0,1).", call. = FALSE)
    }
    # Check that the offset method can be found from the supported methods
    beeswarm_methods <- c("swarm", "compactswarm", "center", "hex", "square")
    vipor_methods <- c(
        "quasirandom", "pseudorandom", "smiley", "maxout", "frowney", "minout",
        "tukey", "tukeyDense")
    jitter_methods <- c("jitter", "none")
    if( !(.is_a_string(point.offset) &&
          point.offset %in% c(
              beeswarm_methods, vipor_methods, jitter_methods)) ){
        stop("'point.offset' must be a single character value from the ",
            "following options: '",
            paste0(
                sort(c(beeswarm_methods, vipor_methods, jitter_methods)),
                collapse = "', '"),
            "'", call. = FALSE)
    }
    # If user do not want to offset points, we disable jitter
    if( point.offset == "none" ){
        jitter.width <- jitter.height <- 0
    }
    # Determine dodge grouping variable, if any
    dodge_var <- if (!is.null(fill.by)) fill.by else group.by
    # Convert categorical x-axis to numeric
    df <- .categorical_x_to_numeric(df, x, facet.by, ...)
    # Apply dodge
    df <- .apply_dodge(df, x, dodge_var, dodge.width)
    # Apply offset based on specified method
    if( point.offset %in% vipor_methods ){
        df <- .apply_vipor_spread(
            df, x, y, facet.by, dodge_var, dodge.width, jitter.width,
            vipor.method = point.offset, ...)
    } else if( point.offset %in% beeswarm_methods ){
        df <- .apply_beeswarm(
            df, x, y, facet.by, dodge_var, dodge.width, jitter.width,
            beeswarm.method = point.offset, ...)
    } else{
        df <- .apply_jitter(
            df, x, y, dodge_var, dodge.width, jitter.width, jitter.height)
    }
    df <- df |> ungroup()
    return(df)
}

# This function converts categorical x axis values to numeric so that we can use
# them to determine position of points
#' @importFrom dplyr group_by across all_of mutate
.categorical_x_to_numeric <- function(df, x, facet.by, scales = "fixed", ...){
    # To disable "no visible binding for global variable" message in cmdcheck
    x_point <- NULL

    # Ensure that the data is in data.frame format
    df <- df |>
        as.data.frame()
    # If there are facets, we specify jitter and dodge for each one
    # separately
    if( scales %in% c("free", "free_x") ){
        df <- df |>
            group_by(across(all_of(facet.by)))
    }
    # Add point positions
    df <- df |>
        mutate(
            # Get original x-axis position
            x_point = as.numeric(factor(.data[[x]])),
            # If x-axis was not specified, move base x-axis back to 0 as they
            # are currently starting from 1.
            x_point = if(all(.data[[x]] == 0)) x_point - 1 else x_point
        )
    return(df)
}

# If there are grouping with group.by or fill.by, we add dodge so that points
# are aligned correctly with the boxplots.
#' @importFrom dplyr mutate n_distinct
.apply_dodge <- function(df, x, dodge.var, dodge.width){
    # To disable "no visible binding for global variable" message in cmdcheck
    x_point <- NULL
    df <- df |>
        mutate(
            # Add dodge if there are groups
            x_point = if( !is.null(dodge.var) && dodge.var != x ) {
                group_index <- as.numeric(factor(.data[[dodge.var]]))
                n_groups <- n_distinct(.data[[dodge.var]])
                x_dodged <- x_point +
                    (group_index - 1 - (n_groups - 1) / 2) * dodge.width /
                    n_groups
            } else{
                x_point
            }
        )
    return(df)
}

# This function adds random jitter to points.
#' @importFrom dplyr mutate
#' @importFrom stats runif
.apply_jitter <- function(
        df, x, y, dodge.var, dodge.width, jitter.width, jitter.height){
    # To disable "no visible binding for global variable" message in cmdcheck
    x_point <-  NULL

    # Calculate spreading
    max_spread <- .get_jitter_spread(
        df, x, dodge.var, dodge.width, jitter.width)
    # Apply jitter
    df <- df |>
        mutate(
            # Add jitter for x axis
            x_point = x_point + runif(n(), -max_spread/2, max_spread/2),
            # Add jitter for y-axis
            y_point = .data[[y]] + runif(n(), -jitter.height, jitter.height)
        )
    return(df)
}

# This function calculates the jitter spread based on grouping and user-defined
# jitter width.
#' @importFrom dplyr n_distinct
.get_jitter_spread <- function(df, x, dodge.var, dodge.width, jitter.width){
    # Get the number of groups. If the coloring/grouping is the same as x axis,
    # it is not taken into account as x axis already have separate placement for
    # points.
    dodge.var <- if( !is.null(dodge.var) && dodge.var != x ) dodge.var
    n_groups <- if (is.null(dodge.var)) 1L else
        n_distinct(df[[dodge.var]])
    # We adjust jitter x axis position based on dodge and the user-
    # specified jitter width.
    max_spread <- dodge.width / n_groups
    max_spread <- jitter.width * max_spread
    return(max_spread)
}

# This function adds beeswarm to points.
#' @importFrom utils getFromNamespace
#' @importFrom dplyr group_by across all_of group_modify arrange select
#' @importFrom scales rescale
.apply_beeswarm <- function(
        df, x, y, facet.by, dodge.var, dodge.width, jitter.width,
        beeswarm.method = "swarm", beeswarm.corral = "none", ...){
    .require_package("beeswarm")
    # To suppress cmdcheck warning:
    # '::' or ':::' import not declared from: ‘beeswarm’
    point_fun <- getFromNamespace("beeswarm", "beeswarm")

    # Add row IDs to preserve original order
    df[[".row_id"]] <- seq_len(nrow(df))

    # Calculate spreading of beeswarm
    max_spread <- .get_jitter_spread(
        df, x, dodge.var, dodge.width, jitter.width)
    # We apply beeswarm for each facet, x axis variable and group
    grouping_vars <- c(facet.by, x, dodge.var) |> unique()

    # Apply beeswarm
    df <- df |>
        group_by(across(all_of(grouping_vars))) |>
        group_modify(~{
            # We calculate beeswarm with beeswarm::beeswarm()
            swarm <- point_fun(
                .x[[y]],
                method = beeswarm.method,
                corral = beeswarm.corral,
                do.plot = FALSE
            )
            # We adjust beeswarm x axis position based on dodge and the user-
            # specified jitter-width
            x_scaled <- rescale(
                swarm[["x"]], to = c(-max_spread/2, max_spread/2))
            .x[["x_point"]] <- mean(.x[["x_point"]]) + x_scaled
            .x[["y_point"]] <- swarm[["y"]]
            return(.x)
        }) |>
        arrange(.row_id) |>
        select(-.row_id)
    return(df)
}

#' @importFrom utils getFromNamespace
#' @importFrom dplyr group_by across all_of group_modify arrange select
#' @importFrom scales rescale
.apply_vipor_spread <- function(
        df, x, y, facet.by, dodge.var, dodge.width, jitter.width,
        vipor.method, ...){
    .require_package("vipor")
    # To suppress cmdcheck warning:
    # '::' or ':::' import not declared from: ‘vipor’
    point_fun <- getFromNamespace("offsetSingleGroup", "vipor")
    # vipor cannot be used with NA values as it calculates density, and density
    # cannot be calculated if there are missing values
    if( anyNA(df[[y]]) ){
        stop("Please choose another offset method. The current option, ",
            "point.offset='", vipor.method, "', cannot be used with ",
            "missing values.", call. = FALSE)
    }

    # Add row IDs to preserve original order
    df[[".row_id"]] <- seq_len(nrow(df))

    # Calculate spreading of beeswarm
    max_spread <- .get_jitter_spread(
        df, x, dodge.var, dodge.width, jitter.width)
    # We apply offset for each facet, x axis variable and group
    grouping_vars <- c(facet.by, x, dodge.var) |> unique()
    # Apply offset function from vipor package
    df <- df |>
        group_by(across(all_of(grouping_vars))) |>
        group_modify(~{
            # vipor offsets for swarm effect
            x_offsets <- point_fun(
                .x[[y]],
                method = vipor.method
            )
            # Adjust x position based on dodge + jitter width
            x_scaled <- rescale(
                x_offsets, to = c(-max_spread/2, max_spread/2))
            .x[["x_point"]] <- mean(.x[["x_point"]]) + x_scaled
            .x[["y_point"]] <- .x[[y]]
            return(.x)
        }) |>
        arrange(.row_id) |>
        select(-.row_id)
    return(df)
}

# This function calculates significance between the groups
.calculate_significance <- function(
        df, y, x, facet.by, fill.by, group.by, pair.by, features,
        paired = !is.null(pair.by),
        significance.method = "wilcox.test", p.adjust.method = "fdr",
        mark.significance = FALSE, digits = 3, ...){
    #
    if( !.is_a_bool(paired) ){
        stop("'paired' must be TRUE or FALSE.", call. = FALSE)
    }
    supported_methods <- c("wilcox.test", "wilcoxon", "t-test", "t.test")
    if( !(.is_a_string(significance.method) &&
          significance.method %in% supported_methods) ){
        stop("'significance.method' must be a single character value from the ",
            "the following options: '",
            paste0(supported_methods, collapse = "', '"), "'", call. = FALSE)
    }
    if( !.is_a_string(p.adjust.method) ){
        stop("'p.adjust.method' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(mark.significance) ){
        stop("'mark.significance' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_an_integer(digits) ){
        stop("'digits' must be a single integer value.", call. = FALSE)
    }
    .require_package("rstatix")
    #
    # Get correct test function
    FUN <- if( significance.method %in% c("wilcox.test", "wilcoxon") )
        rstatix::wilcox_test else rstatix::t_test
    # Get variables. facet.by and x will specify the grouping variables, i.e.,
    # we test the significance for these groups separately.
    grouping_vars <- c(facet.by, x)
    # fill.by. and group.by specify groups for comparison. If they are not
    # specified, we make comparison between x axis groups.
    comparison_vars <- c(fill.by, group.by) |> unique()
    if( is.null(comparison_vars) ){
        comparison_vars <- x
    }
    grouping_vars <- grouping_vars[ !grouping_vars %in% comparison_vars ]
    # Run test
    pvals <- df |>
        as.data.frame() |>
        # Create grouping to test these groups separately
        dplyr::group_split(across(all_of(grouping_vars))) |>
        purrr::map_df(function(df_group) {
            tryCatch({
                # If we calculate paired analysis, sort data so that correct
                # samples are matched with correct subjects.
                df_group <- df_group %>% arrange(across(all_of(pair.by)))
                # This runs pairwise comparisons automatically if there are more
                # than 2 groups
                res <- FUN(
                    formula = as.formula(paste0(y, " ~ ", comparison_vars)),
                    data = df_group,
                    paired = paired,
                    p.adjust.method = "none"
                )
                # Add grouping columns back - assuming grouping_vars are unique 
                # per group
                dplyr::bind_cols(
                    res,
                    df_group %>% select(all_of(grouping_vars)) %>% distinct()
                )
            }, error = function(e) NULL)
        }) |>
        rstatix::adjust_pvalue(method = p.adjust.method)
    # If user wants to mark only significance levels with asterisks
    if( mark.significance ){
        pvals <- pvals |>
            rstatix::add_significance()
    } else{
        # Otherwise we round p-values
        pvals <- pvals |>
            rstatix::p_round(digits = digits)
    }

    # If user specified features, subset the p-values
    if (!is.null(features)) {
        if (!is.null(facet.by)) {
            # Use 'rownames' for subsetting
            pvals <- pvals[pvals[["rownames"]] %in% features, , drop = FALSE]
            
            temp <- attributes(pvals)[["args"]][["data"]]
            temp <- temp[temp[["rownames"]] %in% features, , drop = FALSE]
            attributes(pvals)[["args"]][["data"]] <- temp
            
        } else if ("FeatureID" %in% colnames(pvals)) {
            # Use 'FeatureID' for subsetting
            pvals <- pvals[
                pvals$group1 %in% features & pvals$group2 %in% features,
                , drop = FALSE
            ]
            
            temp <- attributes(pvals)[["args"]][["data"]]
            temp <- temp[temp[["FeatureID"]] %in% features, , drop = FALSE]
            attributes(pvals)[["args"]][["data"]] <- temp
            
        } else {
            # Third case: explicitly do nothing
            # No filtering applied to pvals or its attributes
        }
    }
    
    return(pvals)
}

# Add positions of p-values to the data.frame
.add_p_value_position <- function(
        pvals, x, y, group.by, fill.by, facet.by, df, dodge.width = 0.8,
        step.increase = 0.12, ...
) {
    if (any(!c("y.position", "xmin", "xmax") %in% colnames(pvals))) {
        
        # Determine grouping and comparison variables
        grouping_vars <- unique(c(facet.by, x))
        comparison_vars <- unique(c(fill.by, group.by))
        if (length(comparison_vars) == 0) comparison_vars <- x
        grouping_vars <- setdiff(grouping_vars, comparison_vars)
        
        # Compute max y for each group to stagger p-value y.position
        ypos <- df |>
            group_by(across(all_of(grouping_vars))) |>
            summarise(
                y_base = max(!!rlang::sym(y), na.rm = TRUE),
                .groups = "drop"
            )
        
        # Join y_base to pvals
        # Also dplyr::left_join because left_join is deprecated in mia and 
        # breaks
        if (nrow(ypos) == 1L) {
            # Fallback: compute y_base per comparison (group1/group2)
            x_positions <- unique(c(pvals$group1, pvals$group2))
            
            ypos <- df |>
                filter(!!rlang::sym(x) %in% x_positions) |>
                group_by(across(all_of(x))) |>
                summarise(
                    y_base = max(!!rlang::sym(y), na.rm = TRUE),
                    .groups = "drop"
                ) |>
                rename(xlab = !!x)
            
            pvals <- pvals |>
                dplyr::left_join(ypos, by = c("group1" = "xlab")) |>
                rename(y1 = y_base) |>
                dplyr::left_join(ypos, by = c("group2" = "xlab")) |>
                rename(y2 = y_base) |>
                mutate(y_base = pmax(y1, y2, na.rm = TRUE)) |>
                select(-y1, -y2)
        } else {
            if ("rownames" %in% grouping_vars && "FeatureID" %in% colnames(pvals))
                pvals$rownames <- pvals$FeatureID
            pvals <- dplyr::left_join(pvals, ypos, by = grouping_vars)
        }
        
        if (x == "rownames" || !x %in% colnames(pvals) ) {
            # When x is rownames (usually unique), dodge doesn't make sense.
            # Create artificial x levels from the comparison groups
            x_levels <- unique(c(pvals$group1, pvals$group2))
            x_numeric_map <- setNames(seq_along(x_levels), x_levels)
            
            # No dodging — just place group1 and group2 at their respective 
            # x positions
            pvals$x1 <- x_numeric_map[as.character(pvals$group1)]
            pvals$x2 <- x_numeric_map[as.character(pvals$group2)]
            
            pvals <- pvals |>
                mutate(.group = group1) |>
                group_by(.group) |>
                mutate(
                    y.position = y_base + (row_number() - 1) * 
                        (step.increase * y_base)
                ) |>
                ungroup()
        } else {
            # Standard dodging logic when x is grouped
            x_levels <- sort(unique(df[[x]]))
            x_numeric_map <- setNames(seq_along(x_levels), x_levels)
            pvals$x_numeric <- x_numeric_map[as.character(pvals[[x]])]
            
            fill_levels <- if (!is.null(fill.by)) {
                sort(unique(df[[fill.by]]))
            } else {
                sort(unique(c(pvals$group1, pvals$group2)))
            }
            
            # Compute dodge offset for each group
            dodge_offsets <- seq(-dodge.width / 2, dodge.width / 2, 
                                 length.out = length(fill_levels))
            dodge_map <- setNames(dodge_offsets, fill_levels)
            
            # Apply offsets based on group1 and group2
            pvals <- pvals |>
                mutate(
                    x1 = x_numeric + dodge_map[group1],
                    x2 = x_numeric + dodge_map[group2],
                    xmin = pmin(x1, x2),
                    xmax = pmax(x1, x2),
                    .group = if (length(grouping_vars)) {
                        interaction(!!!syms(grouping_vars), drop = TRUE)
                    } else {
                        interaction(group1, drop = TRUE)
                    }
                ) |>
                group_by(.group) |>
                mutate(
                    y.position = y_base * (1 + step.increase * 
                                               (row_number() - 1))
                ) |>
                ungroup() |>
                select(-.group)
        }
        
    }
    
    return(pvals)
}

# This function is the main plotter function
.plot_boxplot <- function(
        df, add.box = TRUE, add.points = TRUE, 
        scales = "fixed", ncol = NULL, nrow = NULL,
        add.proportion = FALSE, add.threshold = FALSE, ...){
    if( !.is_a_string(scales) ){
        stop("'scales' must be a string.", call. = FALSE)
    }
    #
    # Initialize the plot
    p <- ggplot(df, aes(
        x = .data[[attributes(df)[["x"]]]],
        y = .data[[attributes(df)[["value"]]]],
        group = if(!is.null(attributes(df)[["group.by"]]))
            .data[[attributes(df)[["group.by"]]]]
    ))
    # Add boxplot
    if( add.box ){
        p <- .add_boxplot_layer(p, df, add.points, ...)
    }
    # Add optional points layer
    if( add.points ){
        p <- .add_points_layer(p, df, ...)
    }
    # Add lines connecting points
    if( !is.null(attributes(df)[["pair.by"]]) ){
        p <- .add_line_layers(p, df, ...)
    }
    # If user wants to add prevalence bar under the boxplot
    if( add.proportion ){
        p <- .add_prevalence_bar(p, df, scales, ...)
    }
    # If user wants to add horizontal line to y-axis to denote threshold
    if( add.threshold || add.proportion ){
        p <- .add_threshold_line(p, ...)
    }
    # If faceting was specified, split plot to separate panels
    if( !is.null(attributes(df)[["facet.by"]]) ){
        p <- p +
            facet_wrap(
                as.formula(paste("~", attributes(df)[["facet.by"]])),
                scales = scales, ncol = ncol, nrow = nrow
            )
    }
    # If user wants to add p values
    if( !is.null(attributes(df)[["p.value"]]) ){
        p <- .add_pvalues_to_boxplot(p, df, ...)
    }
    # Adjust themes and titles
    p <- .adjust_boxplot_theme(p, df, add.box, ...)
    return(p)
}

# This function adds boxplot layer
.add_boxplot_layer <- function(
        p, df, add.points, box.alpha = 0.5, dodge.width = 0.8, box.width = 0.75,
        point.shape = 21, ...){
    if( !.is_a_numeric(box.width) ){
        stop("'box.width' must be numeric.", call. = FALSE)
    }
    p <- p + geom_boxplot(
        mapping = aes(
            fill = if(!is.null(attributes(df)[["fill.by"]]))
                .data[[attributes(df)[["fill.by"]]]],
            group = if(!is.null(attributes(df)[["x.box"]]))
                .data[[attributes(df)[["x.box"]]]]
            ),
        # If user wants to add points, we do not add outliers as otherwise they
        # would be plotted twice.
        outlier.shape = if(add.points) NA else point.shape,
        alpha = box.alpha,
        # For boxplot, we can use ggplot's dodging functionality as these
        # positions are deterministic. For points, we set them manually with the
        # jitter.
        position = position_dodge(width = dodge.width),
        width = box.width
    )
    return(p)
}

# This function adds points to plot
.add_points_layer <- function(
        p, df, point.alpha = 0.65, point.size = 2, point.shape = 19L,
        point.colour = point.color, point.color = "grey70", ...){
    # To disable "no visible binding for global variable" message in cmdcheck
    x_point <- y_point <- NULL
    args <- list(
        mapping = aes(
            x = x_point,
            y = y_point,
            colour = if(!is.null(attributes(df)[["colour.by"]]))
                .data[[attributes(df)[["colour.by"]]]],
            shape = if(!is.null(attributes(df)[["shape.by"]]))
                .data[[attributes(df)[["shape.by"]]]],
            size = if(!is.null(attributes(df)[["size.by"]]))
                .data[[attributes(df)[["size.by"]]]],
            fill = if(!is.null(attributes(df)[["fill.by"]]))
                .data[[attributes(df)[["fill.by"]]]]
        ),
        shape = if(is.null(attributes(df)[["shape.by"]])) point.shape,
        alpha = point.alpha,
        size = if(is.null(attributes(df)[["size.by"]])) point.size,
        colour = if(is.null(attributes(df)[["colour.by"]])) point.colour
    )
    args <- args[ lengths(args) > 0 ]
    p <- p + do.call(geom_point, args)
    return(p)
}

# This function connects points with a line
.add_line_layers <- function(
        p, df, line.alpha = 0.5, linetype = 1, linewidth = 1,
        line.colour = line.color, line.color = "grey70", ...
) {
    x_point <- y_point <- NULL
    pair_by_col <- attributes(df)[["pair.by"]]
    diff_col <- attributes(df)[["difference"]]
    facet_by_col <- attributes(df)[["facet.by"]]
    
    grouping_vars <- unique(c(facet_by_col, pair_by_col))
    grouping_vars <- grouping_vars[grouping_vars %in% colnames(df)] 
    
    if (!is.null(diff_col) && diff_col %in% colnames(df)) {
        # Prepare segments for coloring by difference
        seg_df <- df %>%
            arrange(across(all_of(c(grouping_vars, "x_point")))) %>%
            group_by(across(all_of(grouping_vars))) %>%
            mutate(
                xend = dplyr::lead(x_point),
                yend = dplyr::lead(y_point),
                difference_segment = dplyr::lead(.data[[diff_col]])
            ) %>%
            filter(!is.na(xend) & !is.na(yend)) %>%
            ungroup()
        
        p <- p + geom_segment(
            data = seg_df,
            mapping = aes(
                x = x_point,
                y = y_point,
                xend = xend,
                yend = yend,
                group = .data[[pair_by_col]],
                color = .data[["difference_segment"]]
            ),
            alpha = line.alpha,
            linetype = linetype,
            linewidth = linewidth
        ) + scale_color_gradient2(
            low = "blue", mid = "white", high = "red",
            limits = c(-max(abs(seg_df[["difference_segment"]]), na.rm = TRUE),
                       max(abs(seg_df[["difference_segment"]]), na.rm = TRUE))
        )
    } else {
        # No difference: fallback to geom_path with plain color
        args <- list(
            mapping = aes(
                x = x_point,
                y = y_point,
                group = .data[[pair_by_col]],
                color = NULL
            ),
            alpha = line.alpha,
            linetype = linetype,
            linewidth = linewidth,
            colour = line.colour
        )
        args <- args[lengths(args) > 0]
        p <- p + do.call(geom_path, args)
    }
    
    return(p)
}

# This function adds bar under the boxplot to denote prevalence.
#' @importFrom dplyr group_by across all_of mutate
.add_prevalence_bar <- function(
        p, df, scales, threshold = 0, dodge.width = 0.8, ...){
    # To disable "no visible binding for global variable" message in cmdcheck
    min_val <- max_val <- x_point <- width <- y_pos <- heigth <- prevalence <-
        NULL
    if( !.is_a_numeric(threshold) ){
        stop("'threshold' must be a single numeric value.", call. = FALSE)
    }
    if( !(.is_a_numeric(dodge.width) && (dodge.width >=0 && dodge.width <=1)) ){
        stop("'dodge.width' must be numeric (0,1).", call. = FALSE)
    }
    #
    # Get bar width based on box width
    bar_width <- .get_barplot_width(df, ...)
    # Determine dodge grouping variable, if any
    dodge_var <- if (!is.null(attributes(df)[["fill.by"]]))
        attributes(df)[["fill.by"]] else attributes(df)[["group.by"]]
    # Convert categorical x-axis to numeric
    df_prev <- .categorical_x_to_numeric(
        df, attributes(df)[["x"]], attributes(df)[["facet.by"]])
    # Apply dodge
    df_prev <- .apply_dodge(
        df_prev, attributes(df)[["x"]], dodge_var, dodge.width)

    # Calculate prevalence
    grouping_var <- c(
        "x_point", attributes(df)[["facet.by"]],
        attributes(df)[["group.by"]], attributes(df)[["fill.by"]]) |> unique()
    df_prev <- df_prev |>
        group_by(across(all_of(grouping_var))) |>
        summarise(
            prevalence = mean(
                .data[[attributes(df)[["value"]]]] > threshold, na.rm = TRUE),
            min_val = min(.data[[attributes(df)[["value"]]]], na.rm = TRUE),
            max_val = max(.data[[attributes(df)[["value"]]]], na.rm = TRUE),
            .groups = "keep"
            ) |>
        ungroup()

    # Depending on the scales, bar height can also be shared by facets. If the
    # y-axis is free, we adjust the height for each facet separately.
    if( scales %in% c("free", "free_y") ){
        grouping_var <- c(attributes(df)[["facet.by"]]) |> unique()
        df_prev <- df_prev |>
            group_by(across(all_of(grouping_var)))
    }
    # Calculate bars' y-position, height and width of the bar
    df_prev <- df_prev |>
        mutate(
            y_pos = min(min_val) - (max(max_val) - min(min_val))*0.1,
            heigth = (max(max_val) - min(min_val))*0.025,
            width = bar_width
        )

    # Add prevalence bar plot
    p <- p +
        # White background bar
        geom_rect(
            data = df_prev,
            mapping = aes(
                xmin = x_point - width / 2,
                xmax = x_point + width / 2,
                ymin = y_pos - heigth / 2,
                ymax = y_pos + heigth / 2),
            fill = "white", color = "black", inherit.aes = FALSE) +
        # Filled bar
        geom_rect(
            data = df_prev,
            mapping = aes(
                xmin = x_point - width / 2,
                xmax = x_point - width / 2 + prevalence * width,
                ymin = y_pos - heigth / 2,
                ymax = y_pos + heigth / 2),
            fill = "black", inherit.aes = FALSE)
    return(p)
}

# This function calculates the width of the prevalence bar plot based on width
# of boxplot.
.get_barplot_width <- function(df, box.width = 0.75, bar.width = NULL, ...){
    if( !(is.null(bar.width) || .is_a_numeric(bar.width)) ){
        stop("'bar.width' must be numeric.", call. = FALSE)
    }
    #
    if( is.null(bar.width) ){
        # If faceting is not specified, the bar width is simply the boxplot
        # width divided by number of groups as we have to fit all boxes
        # side-by-side.
        bar.width <- (1 - (1-box.width))
        # Get number of unique groups specified by groups
        grouping_vars <- c(
            attributes(df)[["group.by"]], attributes(df)[["fill.by"]]) |>
            unique()
        grouping_vars <- grouping_vars[
            !grouping_vars %in% attributes(df)[["x"]] ]
        n_groups <- nrow(unique(df[, grouping_vars, drop = FALSE]))
        n_groups <- if( n_groups == 0L ) 1L else n_groups
        # If grouping was not specified, we have to take into account the
        # number of x axis variables.
        n_xaxis <- length(unique(df[[attributes(df)[["x"]]]]))
        n_groups <- if( n_groups == 0L ) 1/n_xaxis else n_groups
        bar.width <- bar.width / n_groups
    }
    return(bar.width)
}

# Add horizontal line to plot to denote threshold
.add_threshold_line <- function(p, threshold = 0, ...){
    if( !.is_a_numeric(threshold) ){
        stop("'threshold' must be a single numeric value.", call. = FALSE)
    }
    p <- p + geom_hline(yintercept = threshold, linetype = 2)
    return(p)
}

.validate_pvalue <- function(p.value, df, x, group.by = NULL, fill.by = NULL, 
                             facet.by = NULL) {
    if (!is.null(p.value)) {
        group_fill <- c(group.by, fill.by)
        
        # If group.by or fill.by are specified, group1/group2 must match values 
        # in those columns
        if (length(group_fill) > 0L) {
            # Get all unique values from the specified grouping columns
            group_values <- unique(unlist(lapply(group_fill, 
                                                 function(col) unique(df[[col]]))))
            if (!all(c(p.value[["group1"]], 
                       p.value[["group2"]]) %in% group_values)) {
                stop("Groups in p.value[['group1']] and p.value[['group2']] ",
                     "must match with values in 'fill.by'/'group.by'.", 
                     call. = FALSE)
            }
        }
        # If neither group.by nor fill.by specified, group1/group2 must match x
        else if (length(group_fill) == 0L &&
                 !all(c(p.value[["group1"]], 
                        p.value[["group2"]]) %in% df[[x]])) {
            stop("Groups in p.value[['group1']] and p.value[['group2']] ",
                 "must match with 'x'.", call. = FALSE)
        }
        
        # If x is present in p.value, its values must match df
        if (!is.null(x) && !is.null(p.value[[x]]) &&
            !all(p.value[[x]] %in% df[[x]])) {
            stop("Groups in p.value[['", x, "']] ",
                 "must match with 'x'.", call. = FALSE)
        }
        
        # If facet.by is specified, its values must match df
        if (!is.null(facet.by) && !is.null(p.value[[facet.by]]) &&
            !all(p.value[[facet.by]] %in% df[[facet.by]])) {
            stop("Groups in p.value[['", facet.by, "']] ",
                 "must match with 'facet.by'.", call. = FALSE)
        }
    }
}

# This function adds user-defined p-values to the plot
.add_pvalues_to_boxplot <- function(p, df, ...){
    .require_package("ggpubr")
    args <- list(...)
    # Add p-values
    if (isTRUE(args$mark.significance)) {
        p <- p + ggpubr::stat_pvalue_manual(attributes(df)[["p.value"]], 
                                            label = "p.adj.signif", ...)
    } else {
        p <- p + ggpubr::stat_pvalue_manual(attributes(df)[["p.value"]], 
                                            label = "p.adj", ...)
    }
    return(p)
}

# This function adjust the theme and titles of the plot
.adjust_boxplot_theme <- function(p, df, add.box, ...){
    p <- p + theme_classic()
    # If user did not specify x-axis, we remove all the titles and ticks from
    # x-axis.
    if( attributes(df)[["remove.x.axis"]] ){
        p <- p + theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
            )
    } else{
        # Otherwise, we add correct title from original variable name
        p <- p + labs(
            x = attributes(df)[["x"]]
        )
    }
    # If user did not want to add boxplot player, the x axis is currently
    # numeric because of point layer. Make it categorical.
    if( !add.box ){
        p <- p + scale_x_discrete(limits = levels(
            as.factor(df[[attributes(df)[["x"]]]])))
    }
    # Add correct titles for aesthetics. The titles are the original variable
    # names.
    if( !is.null(attributes(df)[["fill.by"]]) ){
        p <- p + labs(fill = attributes(df)[["fill.by"]])
    }
    if( !is.null(attributes(df)[["colour.by"]]) ){
        p <- p + labs(colour = attributes(df)[["colour.by"]])
    }
    if( attributes(df)[["add.change"]] ){
        p <- p + labs(colour = "Difference")
    }
    if( !is.null(attributes(df)[["shape.by"]]) ){
        p <- p + labs(shape = attributes(df)[["shape.by"]])
    }
    if( !is.null(attributes(df)[["size.by"]]) ){
        p <- p + labs(shape = attributes(df)[["size.by"]])
    }
    return(p)
}
