#' @name
#' plotBoxplot
#'
#' @title
#' Create boxplot of \code{assay}, \code{rowData} or \code{colData}
#'
#' @description
#' This methods visualizes abundances or variables from \code{rowData} or
#' \code{colData}.
#'
#' @details
#' Add here.
#'
#' @return
#' A \code{ggplot2} object.
#'
#' @inheritParams plotAbundance
#'
#' @param assay.type \code{NULL} or \code{character scalar}. Specifies the
#' abundace table to plot. (Default: \code{NULL})
#'
#' @param features \code{NULL} or \code{character vector}. If \code{assay.type}
#' is specified, this specifies rows to visualize in different facets. If
#' \code{NULL}, whole data is visualized as a whole. (Default: \code{NULL})
#'
#' @param row.var \code{NULL} or \code{character vector}. Specifies a variable
#' from \code{rowData(x)} to visualize. (Default: \code{NULL})
#'
#' @param col.var \code{NULL} or \code{character vector} Specifies a variable
#' from \code{colData(x)} to visualize. (Default: \code{NULL})
#'
#' @param ... Additional parameters for plotting.
#' \itemize{
#'   \item \code{layout}: \code{Character scalar}. Specifies the layout of plot.
#'   Must be either \code{"histogram"} or \code{"density"}.
#'   (Default: \code{"histogram"})
#' }
#'
#' @examples
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
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
    function(object, assay.type = NULL, features = NULL, row.var = NULL,
        col.var = NULL, x = NULL, group.by = NULL, ...){
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
        tse, assay.type, features, row.var, col.var, x, group.by, pair.by = NULL,
        add.chance = FALSE, colour.by = color.by, color.by = NULL,
        fill.by = NULL, size.by = NULL, shape.by = NULL, facet.by = NULL,
        add.points = TRUE, ...){
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

# This function checks whether variable can be found from colData or rowData.
.check_metadata_variable <- function(
        tse, var, row = FALSE, col = FALSE,
        var.name = .get_name_in_parent(var)){
    # If the variable is not NULL
    if( !is.null(var) ){
        # It must be a sting and found from colData/rowData
        is_string <- .is_a_string(var)
        check_values <- c()
        check_values <- c(check_values, if(col) colnames(colData(tse)))
        check_values <- c(check_values, if(row) colnames(rowData(tse)))
        var_found <- all( var %in% check_values )
        if( !(is_string && var_found) ){
            stop("'", var.name, "' must be a single character value from the ",
                "following options: '",
                paste0(check_values, collapse = "', '"), "'", call. = FALSE)
        }
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
        fill.by = NULL,
        ...){
    # If assay.type is specified, get melted data
    all_vars <- c(x, group.by, colour.by, size.by, shape.by, facet.by, fill.by)
    if( !is.null(assay.type) ){
        # Specify whether yp retrive data from rowData or colData
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

    # subject.group <- NULL
    # if( !is.null(pair.by) && !is.null(group.by) ){
    #     subject.group <- "subject_group"
    #     df[[subject.group]] <- interaction(df[[pair.by]], df[[group.by]])
    #
    # }


    df <- df |>
        as.data.frame()
    # |>
        # arrange(across(all_of(c(pair.by, facet.by))))

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
        # subject.group = subject.group,
        difference = difference,
        remove.x.axis = remove.x.axis
    )
    return(df)
}

# This function calculates difference between paired samples.
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
.add_fixed_jitterdodge <- function(
        df, x, y, group.by, fill.by, facet.by,
        jitter.width = 0.1, jitter.height = 0, dodge.width = 0.75, ...){
    if( !.is_a_numeric(jitter.width) ){
        stop("'jitter.width' must be numeric.", call. = FALSE)
    }
    if( !.is_a_numeric(jitter.height) ){
        stop("'jitter.height' must be numeric.", call. = FALSE)
    }
    if( !.is_a_numeric(dodge.width) ){
        stop("'dodge.width' must be numeric.", call. = FALSE)
    }
    #
    # Determine dodge grouping variable, if any
    dodge.var <- if (!is.null(fill.by)) fill.by else group.by
    #
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
            x_point = if(all(.data[[x]] == 0)) x_point - 1 else x_point,
            # Add dodge if there are groups
            x_point = if (!is.null(dodge.var)) {
                group_index = as.numeric(factor(.data[[dodge.var]]))
                n_groups = n_distinct(.data[[dodge.var]])
                x_dodged = x_point +
                    (group_index - 1 - (n_groups - 1) / 2) * dodge.width /
                    n_groups
            } else{
                x_point
            },
            # Add jitter
            x_point = x_point + runif(n(), -jitter.width, jitter.width),
            # Add jitter for y-axis
            y_point = .data[[y]] + runif(n(), -jitter.height, jitter.height)
        ) |>
        ungroup()
    return(df)
}

# This function is the main plotter function
.plot_boxplot <- function(df, add.points = TRUE, scales = "fixed", ...){
    # Initialize the plot
    p <- ggplot(df, aes(
        x = .data[[attributes(df)[["x"]]]],
        y = .data[[attributes(df)[["value"]]]],
        group = if(!is.null(attributes(df)[["group.by"]])) .data[[attributes(df)[["group.by"]]]]
    ))
    # Add boxplot layer
    p <- .add_boxplot_layer(p, df, add.points, ...)
    # Add optional points layer
    if( add.points ){
        p <- .add_points_layer(p, df, ...)
    }
    # Add lines connecting points
    if( !is.null(attributes(df)[["pair.by"]]) ){
        p <- .add_line_layers(p, df, ...)
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
    p <- .adjust_boxplot_theme(p, df, ...)
    return(p)
}

# This function adds boxplot layer
.add_boxplot_layer <- function(
        p, df, add.points, box.alpha = 0.5, dodge.width = 0.8,
        point.shape = 21, ...){
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
        position = position_dodge(width = dodge.width)
    )
    return(p)
}

# This function adds points to plot
.add_points_layer <- function(
        p, df, point.alpha = 0.65, point.size = 2, point.shape = 21,
        point.colour = "grey70", ...){
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
        line.colour = "grey70",...){
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
            low="blue", mid="white", high="red",
            limits = c(
                -max(abs(df[[attributes(df)[["difference"]]]])),
                max(abs(df[[attributes(df)[["difference"]]]])))
            )
    }
    return(p)
}

# This function adjust the theme and titles of the plot
.adjust_boxplot_theme <- function(p, df, ...){
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
