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
#'   If \code{"none"}, ofsetting is not applied. (Default: \code{"jitter"})
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
#'   \item \code{add.chance}: \code{Logical scalar}. Whether to visualize chance
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
#' }
#'
#' @examples
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
#' \dontrun{
#' library(microbiomeDataSets)
#'
#' mae <- microbiomeDataSets::peerj32()
#' tse <- getWithColData(mae, 1)
#' tse[["time_point"]] <- as.character(tse[["time"]])

#' # Create a plot showing chance between time points in abundance of
#' # Akkermansia
#' plotBoxplot(
#'     tse, x = "time_point", assay.type = "counts", fill.by = "group",
#'     features = "Akkermansia", pair.by = "subject",
#'     add.chance = TRUE, scales = "free"
#' )
#' }
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
        pair.by = NULL, add.chance = FALSE, colour.by = color.by,
        color.by = NULL, fill.by = NULL, size.by = NULL, shape.by = NULL,
        facet.by = NULL, add.box = TRUE, add.points = TRUE,
        add.proportion = FALSE, ...){
    # Either assay.type. row.var or col.var must be specified
    if( sum(c(is.null(assay.type), is.null(row.var), is.null(col.var))) != 2L ){
        stop("Please specify either 'assay.type', 'row.var', or 'col.var'.",
            call. = FALSE)
    }
    # features cannot be specified if row.var or col.var is specified
    if( is.null(assay.type) && !is.null(features) ){
        stop("'features' can be specified only when 'assay.type is ",
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
    if( !.is_a_bool(add.chance) ){
        stop("'add.chance' must be TRUE or FALSE.", call. = FALSE)
    }
    if( add.chance && is.null(pair.by) ){
        stop("When 'add.chance' is specified, 'pair.by' must be specified.",
            call. = FALSE)
    }
    # We have to plot points in order to connect them
    if( !is.null(pair.by) && !add.points ){
        stop("When 'pair.by' is specified, 'add.points' must be enabled.",
            call. = FALSE)
    }
    # We can either color based on variable or the difference netween paired
    # samples. There cannot be multiple coloring schemes.
    if( add.chance && !is.null(colour.by) ){
        stop("Both 'add.chance' and 'colour.by' cannot be specified ",
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
    return(NULL)
}

# This function retrieves the data from TreeSE and outputs a data.frame, ready
# for the plotter function.
.get_data_for_boxplot <- function(
        tse, x = NULL, assay.type, features, row.var, col.var, group.by,
        pair.by = NULL, add.chance = FALSE,
        colour.by = color.by, color.by = NULL,
        size.by = NULL, shape.by = NULL, facet.by = NULL,
        fill.by = NULL, add.proportion = FALSE,
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
    # If features were specified, subset data
    if( !is.null(features) ){
        df <- df[ df[["FeatureID"]] %in% features, , drop = FALSE]
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

    # If user specified, calculate difference
    difference <- NULL
    if( add.chance ){
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
    # seed, but the problem comes from facetting. Even though, we know the
    # points' positions before facetting, they are not the same after facetting.
    # Setting manually the jitter for points is much easier for us.
    df <- .add_fixed_jitterdodge(
        df, x, c(assay.type, col.var, row.var), group.by, fill.by, facet.by,
        ...)
    df <- df |>
        as.data.frame()

    # Add plotting options to attributes of the data.frame. Now the data.frame
    # includes all the information for plotting.
    attributes(df) <- c(
        attributes(df),
        value = c(assay.type, col.var, row.var),
        x = x,
        group.by = group.by,
        pair.by = pair.by,
        add.chance = add.chance,
        colour.by = colour.by,
        fill.by = fill.by,
        size.by = size.by,
        shape.by = shape.by,
        facet.by = facet.by,
        x.box = x.box,
        difference = difference,
        remove.x.axis = remove.x.axis
    )
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
        arrange(desc(across(all_of(c(pair.by, x)))))
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
    df <- .categorical_x_to_numeric(df, x, facet.by)
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
.categorical_x_to_numeric <- function(df, x, facet.by){
    # To disable "no visible binding for global variable" message in cmdcheck
    x_point <- NULL
    df <- df |>
        as.data.frame() |>
        # If there are facets, we specify jitter and dodge for each one
        # separately
        group_by(across(all_of(facet.by))) |>
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

# This function is the main plotter function
.plot_boxplot <- function(
        df, add.box = TRUE, add.points = TRUE, scales = "fixed",
        add.proportion = FALSE, ...){
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
    # If facetting was specified, split plot to separate panels
    if( !is.null(attributes(df)[["facet.by"]]) ){
        p <- p +
            facet_wrap(
                ~ .data[[attributes(df)[["facet.by"]]]],
                scales = scales
            )
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
        p, df, point.alpha = 0.65, point.size = 2, point.shape = 21,
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
        line.colour = line.color, line.color = "grey70", ...){
    # To disable "no visible binding for global variable" message in cmdcheck
    x_point <- y_point <- NULL
    args <- list(
        mapping = aes(
            x = x_point,
            y = y_point,
            group = .data[[attributes(df)[["pair.by"]]]],
            color = if(!is.null(attributes(df)[["difference"]]))
                .data[[attributes(df)[["difference"]]]]
        ),
        alpha = line.alpha,
        linetype = linetype,
        linewidth = linewidth,
        colour = if(is.null(attributes(df)[["difference"]])) line.colour
    )
    args <- args[ lengths(args) > 0 ]
    p <- p + do.call(geom_path, args)
    # If user wanted to also visualize difference between consecutive samples,
    # we improve the color scale to blue-white-red
    if( !is.null(attributes(df)[["difference"]]) ){
        p <- p + scale_color_gradient2(
            low = "blue", mid = "white", high = "red",
            limits = c(
                -max(abs(df[[attributes(df)[["difference"]]]])),
                max(abs(df[[attributes(df)[["difference"]]]])))
            )
    }
    return(p)
}

# This function adds bar under the boxplot to denote prevalence.
#' @importFrom dplyr group_by across all_of mutate
.add_prevalence_bar <- function(
        p, df, scales, threshold = 0, dodge.width = 0.8, add.threshold = TRUE,
        ...){
    # To disable "no visible binding for global variable" message in cmdcheck
    min_val <- max_val <- x_point <- width <- y_pos <- heigth <- prevalence <-
        NULL
    if( !.is_a_numeric(threshold) ){
        stop("'threshold' must be a single numeric value.", call. = FALSE)
    }
    if( !(.is_a_numeric(dodge.width) && (dodge.width >=0 && dodge.width <=1)) ){
        stop("'dodge.width' must be numeric (0,1).", call. = FALSE)
    }
    if( !(length(add.threshold) && is.logical(add.threshold)) ){
        stop("'add.threshold' must be TRUE or FALSE.", call. = FALSE)
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
            )
    # Calculate bars y-position. It is shared by facets.
    grouping_var <- c(attributes(df)[["facet.by"]]) |> unique()
    df_prev <- df_prev |>
        group_by(across(all_of(grouping_var))) |>
        mutate(
            y_pos = min(min_val) - (max(max_val) - min(min_val))*0.1
            )
    # Depending on the scales, bar height can also be shared by facets. If the
    # y-axis is free, we adjust the height for each facet separately.
    if( !scales %in% c("free", "free_y") ){
        df_prev <- df_prev |> ungroup()
    }
    # Calculate height and width of the bar
    df_prev <- df_prev |>
        mutate(
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

    # Finally add threshold as horizontal line
    if( add.threshold ){
        p <- p + geom_hline(yintercept = threshold, linetype = 2)
    }

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
        # If facetting is not specified, the bar width is simply the boxplot
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
    if( attributes(df)[["add.chance"]] ){
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
