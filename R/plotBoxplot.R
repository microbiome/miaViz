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

.check_input_for_boxplot <- function(
        tse, assay.type, features, row.var, col.var, x, group.by, pair.by = NULL,
        add.chance = FALSE, colour.by = color.by, color.by = NULL,
        fill.by = NULL,
        size.by = NULL, shape.by = NULL, facet.by = NULL, add.points = TRUE,
        ...){
    # Either assay.type. row.var or col.var must be specified
    if( sum(c(is.null(assay.type), is.null(row.var), is.null(col.var))) !=
        2L ){
        stop("Please specify either 'assay.type', 'row.var', or 'col.var'.",
             call. = FALSE)
    }
    # features cannot be specified if row.var or col.var is specified
    if( is.null(assay.type) && !is.null(features) ){
        stop("'features' can be specified only when 'assay.type is ",
             "specified.", call. = FALSE)
    }
    # Check that the values are correct
    if( !is.null(features) && is.null(rownames(tse)) ){
        stop("'x' must have rownames.", call. = FALSE)
    }
    if( !(is.null(assay.type) || .is_a_string(assay.type)) ){
        stop("'assay.type' must be NULL or single character value.",
             call. = FALSE)
    }
    if( !(is.null(row.var) || .is_a_string(row.var)) ){
        stop("'row.var' must be NULL or single character value.",
             call. = FALSE)
    }
    if( !(is.null(col.var) || .is_a_string(col.var)) ){
        stop("'col.var' must be NULL or single character value.",
             call. = FALSE)
    }
    if( !(is.null(features) || is.character(features)) ){
        stop("'features' must be NULL or single character value.",
             call. = FALSE)
    }
    if( !(is.null(x) || .is_a_string(x)) ){
        stop("'x' must be NULL or single character value.",
             call. = FALSE)
    }
    if( !(is.null(group.by) || .is_a_string(group.by)) ){
        stop("'group.by' must be NULL or single character value.",
             call. = FALSE)
    }
    if( !(is.null(pair.by) || .is_a_string(pair.by)) ){
        stop("'pair.by' must be NULL or single character value.",
             call. = FALSE)
    }
    if( !.is_a_bool(add.chance) ){
        stop("'add.chance' must be TRUE or FALSE.", call. = FALSE)
    }
    if( add.chance && is.null(pair.by) ){
        stop("When 'add.chance' is specified, 'pair.by' must be speicifed.", call. = FALSE)
    }
    # If parameters are specified, check that they can be found
    if( !is.null(assay.type) ){
        .check_assay_present(assay.type, tse)
    }
    if( !is.null(row.var) && !all(row.var %in% colnames(rowData(tse))) ){
        stop("'row.var' must be from the following options: '",
             paste0(colnames(rowData(tse)), collapse = "', '"), "'",
             call. = FALSE)
    }
    if( !is.null(col.var) && !all(col.var %in% colnames(colData(tse))) ){
        stop("'col.var' must be from the following options: '",
             paste0(colnames(colData(tse)), collapse = "', '"), "'",
             call. = FALSE)
    }
    if( !is.null(features) && !all(features %in% rownames(tse)) ){
        stop("'feature' must specify features from rownames(object).",
             call. = FALSE)
    }
    #


    if( !is.null(x) && !is.null(col.var) && .is_a_string(x) && !x %in% colnames(colData(tse)) ){
        stop("'x' must be from the following options: '",
             paste0(colnames(colData(tse)), collapse = "', '"), "'",
             call. = FALSE)
    }
    if( !is.null(x) && !is.null(row.var) && .is_a_string(x) && !x %in% colnames(rowData(tse)) ){
        stop("'x' must be from the following options: '",
             paste0(colnames(rowData(tse)), collapse = "', '"), "'",
             call. = FALSE)
    }
    if( !is.null(x) && !is.null(assay.type) && .is_a_string(x) && !x %in% c(colnames(colData(tse)), colnames(rowData(tse))) ){
        stop("'x' must be from the following options: '",
             paste0(c(colnames(colData(tse)), colnames(rowData(tse))), collapse = "', '"), "'",
             call. = FALSE)
    }
    if( !((!is.null(assay.type) || !is.null(col.var)) && (is.null(group.by) || all(group.by %in% colnames(colData(tse))))) ){
        stop("'group.by' must be from the following options: '",
             paste0(colnames(colData(tse)), collapse = "', '"), "'",
             call. = FALSE)
    }
    if( !is.null(row.var) && !(is.null(group.by) || all(group.by %in% colnames(rowData(tse)))) ){
        stop("'group.by' must be from the following options: '",
             paste0(colnames(rowData(tse)), collapse = "', '"), "'",
             call. = FALSE)
    }
    if( !((!is.null(assay.type) || !is.null(col.var)) && (is.null(pair.by) || all(pair.by %in% colnames(colData(tse))))) ){
        stop("'pair.by' must be from the following options: '",
             paste0(colnames(colData(tse)), collapse = "', '"), "'",
             call. = FALSE)
    }
    if( !is.null(row.var) && !(is.null(pair.by) || all(pair.by %in% colnames(rowData(tse)))) ){
        stop("'pair.by' must be from the following options: '",
             paste0(colnames(rowData(tse)), collapse = "', '"), "'",
             call. = FALSE)
    }
    if( !is.null(assay.type) ){
        cols <- c(colnames(rowData(tse)), colnames(colData(tse)))
    } else if( !is.null(col.var) ){
        cols <- colnames(colData(tse))
    } else{
        cols <- colnames(rowData(tse))
    }
    if( !is.null(colour.by) && !(.is_a_string(colour.by) && colour.by %in% cols) ){
        stop("'colour.by' must be from the following options: '",
             paste0(cols, collapse = "', '"), "'",
             call. = FALSE)
    }
    if( !is.null(fill.by) && !(.is_a_string(fill.by) && fill.by %in% cols) ){
        stop("'fill.by' must be from the following options: '",
             paste0(cols, collapse = "', '"), "'",
             call. = FALSE)
    }
    if( !is.null(size.by) && !(.is_a_string(size.by) && size.by %in% cols) ){
        stop("'size.by' must be from the following options: '",
             paste0(cols, collapse = "', '"), "'",
             call. = FALSE)
    }
    if( !is.null(shape.by) && !(.is_a_string(shape.by) && shape.by %in% cols) ){
        stop("'shape.by' must be from the following options: '",
             paste0(cols, collapse = "', '"), "'",
             call. = FALSE)
    }
    if( !is.null(facet.by) && !(.is_a_string(facet.by) && facet.by %in% cols) ){
        stop("'facet.by' must be from the following options: '",
             paste0(cols, collapse = "', '"), "'",
             call. = FALSE)
    }

    if( add.chance && !is.null(colour.by) ){
        stop("Both 'add.chance' and 'colour.by' cannot be specidied simultaneously.", call. = FALSE)
    }

    # if( !is.null(pair.by) && !(!is.null(x) && !is.null(group.by)) ){
    #     stop("'pair.by' must be specified with 'x' and 'group.by'.", call. = FALSE)
    # }

    # x must be charcter
    if( !is.null(x) && is.numeric(tse[[x]]) ){
        stop("x must be charcter.,")
    }

    if( !is.null(group.by) && !is.null(fill.by) ){
        stop("djfkfkf")
    }

    if( !is.null(pair.by) && !add.points ){
        stop("must points be there-")
    }
    return(NULL)
}

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
        row_vars <- vapply(all_vars, function(x){
            x %in% colnames(rowData(tse))
        }, logical(1L))
        col_vars <- all_vars[ !row_vars ]
        row_vars <- all_vars[ row_vars ]

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
        # colnames(df)[ colnames(df) == "FeatureID" ] <- "facet_by"
    }
    # If row.var was specified, get the data from rowData
    if( !is.null(row.var) ){
        df <- rowData(tse)[, c(row.var, all_vars), drop = FALSE]
    }
    # If col.var was specified, get the data from colData
    if( !is.null(col.var) ){
        df <- colData(tse)[, c(col.var, pair.by, all_vars), drop = FALSE]
    }

    if( !is.numeric(df[[c(assay.type, col.var, row.var)]]) ){
        stop("Y-axis must be numeric.", call. = FALSE)
    }

    if( add.chance ){
        df <- df |>
            as.data.frame() |>
            arrange(across(all_of(c(pair.by, x, group.by, fill.by)))) |>
            group_by(.data[[pair.by]]) |>
            mutate(
                difference = .data[[c(assay.type, col.var, row.var)]] - dplyr::lag(.data[[c(assay.type, col.var, row.var)]])
            ) |>
            ungroup()
        if( !is.null(x) ){
            df <- df |>
                arrange(desc(across(all_of(c(pair.by, x)))))
        }

    }

    if( !is.null(pair.by) && !is.null(group.by) ){
        df[["pair+group"]] <- interaction(df[[pair.by]], df[[group.by]])

    }

    if( !is.null(x) && !is.null(group.by) ){
        df[["x+group"]] <- interaction(df[[x]], df[[group.by]])
    }

    df <- df |> as.data.frame()

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
        facet.by = facet.by
    )

    return(df)
}

.plot_boxplot <- function(
        df, add.points = TRUE, scales = "fixed", ...){
    p <- ggplot(df, aes(
        x = if(!is.null(attributes(df)[["x"]])) .data[[attributes(df)[["x"]]]] else 0,
        y = .data[[attributes(df)[["value"]]]],
        fill = if(!is.null(attributes(df)[["fill.by"]])) .data[[attributes(df)[["fill.by"]]]],
        group = if(!is.null(attributes(df)[["group.by"]])) .data[[attributes(df)[["group.by"]]]]
    ))

    p <- .add_boxplot_layer(p, add.points, ...)

    point_position <- .get_point_position(df, add.points, ...)
    if( add.points ){
        p <- .add_points_layer(p, df, point_position, ...)
    }
    if( !is.null(attributes(df)[["pair.by"]]) ){
        p <- .add_line_layers(p, df, point_position, ...)
    }
    p <- .adjust_boxplot_theme(p, df, ...)

    if( !is.null(attributes(df)[["facet.by"]]) ){
        p <- p +
            facet_wrap(~ .data[[attributes(df)[["facet.by"]]]], scales = scales)
    }
    return(p)
}

.add_boxplot_layer <- function(p, add.points, box.alpha = 0.5,  ...){
    p <- p + geom_boxplot(
        outlier.shape = if(add.points) NA else 19,
        alpha = box.alpha
    )
    return(p)
}

.get_point_position <- function(
        df, add.points,
        add.jitter = TRUE,
        dodge.width = 0.75,
        jitter.width = if(add.jitter) 0.4 else 0,
        seed = 6785,
        ...){
    if( add.points &&  length(c(attributes(df)[["colour.by"]], attributes(df)[["group.by"]], attributes(df)[["fill.by"]], attributes(df)[["shape.by"]])) > 0L ){
        pos <- position_jitterdodge(
            dodge.width = dodge.width,
            jitter.width = jitter.width,
            seed = seed
        )
    } else if( add.points && add.jitter ){
        pos <- "jitter"
    }  else{
        pos <- "identity"
    }
    return(pos)
}

.add_points_layer <- function(p, df, point_position, ...){
    p <- p + geom_point(aes(
        colour = if(!is.null(attributes(df)[["colour.by"]])) .data[[attributes(df)[["colour.by"]]]],
        shape = if(!is.null(attributes(df)[["shape.by"]])) .data[[attributes(df)[["shape.by"]]]],
        size = if(!is.null(attributes(df)[["size.by"]])) .data[[attributes(df)[["size.by"]]]],
        fill = if(!is.null(attributes(df)[["fill.by"]])) .data[[attributes(df)[["fill.by"]]]] else "black"
    ), position = point_position)
    # if( !is.null(attributes(df)[["pair.by"]]) ){
    #     p <- p + geom_point(aes(
    #         group = if(!is.null(attributes(df)[["pair.by"]])) .data[[attributes(df)[["pair.by"]]]],
    #         colour = if(!is.null(attributes(df)[["colour.by"]])) .data[[attributes(df)[["colour.by"]]]],
    #         shape = if(!is.null(attributes(df)[["shape.by"]])) .data[[attributes(df)[["shape.by"]]]],
    #         size = if(!is.null(attributes(df)[["size.by"]])) .data[[attributes(df)[["size.by"]]]],
    #         fill = if(!is.null(attributes(df)[["fill.by"]])) .data[[attributes(df)[["fill.by"]]]] else "black"
    #     ), position = point_position)
    # } else{
    #     p <- p + geom_point(aes(
    #         colour = if(!is.null(attributes(df)[["colour.by"]])) .data[[attributes(df)[["colour.by"]]]],
    #         shape = if(!is.null(attributes(df)[["shape.by"]])) .data[[attributes(df)[["shape.by"]]]],
    #         size = if(!is.null(attributes(df)[["size.by"]])) .data[[attributes(df)[["size.by"]]]],
    #         fill = if(!is.null(attributes(df)[["fill.by"]])) .data[[attributes(df)[["fill.by"]]]] else "black"
    #     ), position = point_position)
    # }
    return(p)
}

.add_line_layers <- function(p, df, point_position, ...){

    p_data <- ggplot_build(p)

    p_data <- p_data[["data"]][[2L]]
    p_data
    p_data <- cbind(p_data, df[, colnames(df) %in% c(attributes(df)[["pair.by"]], "difference", attributes(df)[["facet.by"]]), drop = FALSE])

    p_data <- p_data[, colnames(p_data) %in% c(attributes(df)[["pair.by"]], "x", "y", "difference", attributes(df)[["facet.by"]]), drop = FALSE]
    library(dplyr)
    library(tidyr)

    p_lines <- p_data |>
        arrange(subject, x) |>
        group_by(across(all_of(c(attributes(df)[["pair.by"]], attributes(df)[["facet.by"]])))) |>
        filter(n() > 1) |>
        mutate(
            xend = lead(x),
            yend = lead(y),
            difference_lead = if ("difference" %in% names(p_data)) lead(difference) else NA_real_
        ) |>
        filter(!is.na(xend)) |>
        ungroup()

    p <- p +
        geom_segment(
            data = p_lines,
            aes(
                x = x, xend = xend, y = y, yend = yend,
                color = if(any(!is.na(.data[["difference_lead"]]))) difference_lead ),
            inherit.aes = FALSE
        )

#
#     p <- p +
#         geom_path(
#             data = df,
#             aes(
#                 group = .data[[attributes(df)[["pair.by"]]]],
#                 color = if("difference" %in% colnames(df)) difference,
#                 fill = if(!"difference" %in% colnames(df)) "black"
#             ),
#             position = point_position
#         ) |>
#         # For position_jitterdodge(), we have to specify one of the following
#         # aesthetics: "fill", "colour", "linetype", "shape", "size", or "alpha".
#         # As fill does not do nothing, we choose it. However, it gives this
#         # warning which is suppressed: "Ignoring unknown aesthetics: fill"
#         suppressWarnings()
    if( "difference" %in% colnames(df) ){
        p <- p +
            scale_color_gradient2(
                low="blue", mid="white", high="red",
                limits = c(-max(abs(df[["difference"]])), max(abs(df[["difference"]]))))
    }
    return(p)
}

.adjust_boxplot_theme <- function(p, df, ...){
    p <- p + theme_classic()
    if( is.null(attributes(df)[["x"]]) ){
        p <- p +
            theme(
                axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank()
            )
    }
    if( is.null(attributes(df)[["fill.by"]]) ){
        p <- p +
            guides(fill = "none")
    }
    if( !is.null(attributes(df)[["x"]]) ){
        p <- p + labs(
            x = attributes(df)[["x"]]
        )
    }
    if( !is.null(attributes(df)[["fill.by"]]) ){
        p <- p +
            labs(
                fill = attributes(df)[["fill.by"]]
            )
    }
    if( !is.null(attributes(df)[["colour.by"]]) ){
        p <- p +
            labs(
                colour = attributes(df)[["colour.by"]]
            )
    }
    if( attributes(df)[["add.chance"]] ){
        p <- p +
            labs(
                colour = "Difference"
            )
    }
    if( !is.null(attributes(df)[["shape.by"]]) ){
        p <- p +
            labs(
                shape = attributes(df)[["shape.by"]]
            )
    }
    if( !is.null(attributes(df)[["size.by"]]) ){
        p <- p +
            labs(
                shape = attributes(df)[["size.by"]]
            )
    }
    return(p)
}
