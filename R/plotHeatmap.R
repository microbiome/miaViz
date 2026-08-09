#' @name
#' plotHeatmap
#'
#' @title
#' Visualize assay values as a heatmap
#'
#' @description
#' Creates a heatmap from an assay stored in a
#' \code{SummarizedExperiment} or
#' \code{TreeSummarizedExperiment}. Features are shown as rows and
#' samples as columns.
#'
#' @details
#' \code{plotHeatmap} visualizes values from an assay as a heatmap. Values
#' can be optionally centered and/or scaled across samples for each feature,
#' which is useful for highlighting relative abundance patterns rather than
#' absolute abundances.
#'
#' Additional variables from \code{rowData(x)} and \code{colData(x)} can be
#' used to facet the heatmap. When multiple row or column variables are
#' provided, nested facets are created using
#' \pkg{ggh4x}.
#'
#' For \code{TreeSummarizedExperiment} objects, a row tree can be displayed
#' alongside the heatmap. When a tree is shown, only leaf nodes are plotted
#' and the heatmap rows are reordered to match the tree tip order.
#'
#' @return
#' A \code{ggplot2} object. If \code{show.tree = TRUE}, the returned object is
#' a combined \code{patchwork} object containing the row tree and the heatmap.
#'
#' @inheritParams plotAbundance
#'
#' @param row.var \code{NULL} or \code{character vector}. Variables from
#' \code{rowData(x)} used for row facetting.
#' (Default: \code{NULL})
#'
#' @param col.var \code{NULL} or \code{character vector}. Variables from
#' \code{colData(x)} used for column facetting.
#' (Default: \code{NULL})
#'
#' @param scale \code{Logical scalar}. Should assay values be scaled for each
#' feature across samples?
#' (Default: \code{FALSE})
#'
#' @param center \code{Logical scalar}. Should assay values be centered for
#' each feature across samples?
#' (Default: \code{FALSE})
#'
#' @param tree.name \code{Character scalar}. Name of the row tree to display
#' when \code{x} is a \code{TreeSummarizedExperiment}.
#' (Default: \code{"phylo"})
#'
#' @param show.tree \code{Logical scalar}. Should the row tree be displayed?
#' Only available for \code{TreeSummarizedExperiment}.
#' (Default: \code{TRUE})
#'
#' @param ... Additional parameters controlling the visualization. When
#' \code{show.tree = TRUE}, additional arguments are passed to
#' \code{\link{plotRowTree}} to control the appearance of the tree.
#' Additional parameters include, for example:
#' \itemize{
#'   \item \code{scales}: Facet scaling passed to
#'   \code{ggplot2::facet_grid()} or
#'   \code{ggh4x::facet_nested()}.
#'   (Default: \code{"free"})
#'   \item \code{tree.width}: Relative width of the tree panel when
#'   \code{show.tree = TRUE}.
#'   (Default: \code{0.2})
#' }
#'
#' @examples
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#'
#' tse <- agglomerateByPrevalence(tse, rank = "Class")
#'
#' # Plot raw counts
#' plotHeatmap(
#'     tse,
#'     assay.type = "counts"
#' )
#'
#' # Scale and center each feature
#' plotHeatmap(
#'     tse,
#'     assay.type = "counts",
#'     scale = TRUE,
#'     center = TRUE
#' )
#'
#' # Facet samples by metadata
#' plotHeatmap(
#'     tse,
#'     assay.type = "counts",
#'     col.var = "SampleType"
#' )
#'
#' # Add phylogeny
#' plotHeatmap(
#'     tse,
#'     assay.type = "counts",
#'     show.tree = TRUE
#' )
#'
#' @seealso
#' \code{\link[=plotRowTree]{plotRowTree}}
#'
NULL

#' @rdname plotHeatmap
#' @export
setMethod(
    "plotHeatmap",
    signature = c(x = "TreeSummarizedExperiment"),
    function(
        x,
        assay.type = NULL,
        row.var = NULL,
        col.var = NULL,
        scale = FALSE,
        center = FALSE,
        tree.name = "phylo",
        show.tree = FALSE,
        ...
    ) {

        # Validate tree-related arguments
        if (!.is_a_bool(show.tree)) {
            stop("'show.tree' must be TRUE or FALSE.", call. = FALSE)
        }

        # Prepare the TreeSummarizedExperiment for plotting
        if (show.tree) {

            # Row facets are incompatible with the row tree
            if (!is.null(row.var)) {
                stop(
                    "'row.var' cannot be specified when 'show.tree = TRUE'.",
                    call. = FALSE
                )
            }

            # Check that the requested row tree is available
            .check_rowTree_present(tree.name, x)

            # Create a tree plot that we will use to extract order. We have to
            # create a plot, because exact top-to-bottom ordering of tips is
            # determined when the tree is plotted. For example, branches can be
            # rotated around internal nodes without changing the tree.
            tree <- rowTree(x, tree.name)
            tree_plot <- ggtree(tree) + geom_tiplab()
            tree_data <- ggplot_build(tree_plot)
            tips <- tree_data[["data"]][[3]][, c("y", "label")]
            tips <- tips[order(tips[["y"]]), , drop = FALSE]

            tax_order <- match(
                tips[["label"]],
                rowLinks(x)[["nodeLab"]]
            )

            x <- x[tax_order, ]
            #
            # # Keep only features represented by tree leaves
            # x <- subsetByLeaf(x, whichRowTree = tree.name)
            #
            # # Reorder rows to match the tree tip order
            # links <- rowLinks(x)
            # keep <- links$whichTree == tree.name & links$isLeaf
            # links <- links[keep, "nodeNum"] |> order()
            # x <- x[links, ]
        }

        # Create the heatmap using the SummarizedExperiment method
        p <- callNextMethod(
            x = x,
            assay.type = assay.type,
            row.var = row.var,
            col.var = col.var,
            scale = scale,
            center = center,
            ...
        )

        # Add the row tree to the heatmap
        if (show.tree) {
            p <- .add_row_tree_to_heatmap(
                tse = x,
                p = p,
                tree.name = tree.name,
                ...
            )
        }

        return(p)
    }
)


#' @rdname plotHeatmap
#' @export
setMethod(
    "plotHeatmap",
    signature = c(x = "SummarizedExperiment"),
    function(
        x,
        assay.type = NULL,
        row.var = NULL,
        col.var = NULL,
        scale = FALSE,
        center = FALSE,
        ...
    ) {

        # Validate input
        .check_assay_present(assay.type, x)
        .check_metadata_variable(x, row.var, row = TRUE, multiple = TRUE)
        .check_metadata_variable(x, col.var, col = TRUE, multiple = TRUE)

        if (!.is_a_bool(scale)) {
            stop("'scale' must be TRUE or FALSE.", call. = FALSE)
        }

        if (!.is_a_bool(center)) {
            stop("'center' must be TRUE or FALSE.", call. = FALSE)
        }

        # Retrieve assay values together with the requested metadata
        df <- .get_heatmap_data(
            tse = x,
            assay.type = assay.type,
            row.var = row.var,
            col.var = col.var,
            scale = scale,
            center = center
        )

        # Update the legend title to reflect any transformation
        assay.name <- assay.type
        if (scale && center) {
            assay.type <- paste0(assay.type, " (scaled & centered)")
        } else if (scale) {
            assay.type <- paste0(assay.type, " (scaled)")
        } else if (center) {
            assay.type <- paste0(assay.type, " (centered)")
        }
        names(df)[names(df) == assay.name] <- assay.type

        # Create the heatmap
        p <- .plot_heatmap(
            df = df,
            fill = assay.type,
            scale = scale,
            ...
        )

        # Add row and/or column facets
        # No faceting requested
        if (!is.null(row.var) || !is.null(col.var)) {
            p <- .resolve_heatmap_facets(
                p = p,
                row.var = row.var,
                col.var = col.var,
                ...
            )
        }

        return(p)
    }
)

############################### HELPER FUNCTIONS ###############################

# Retrieve assay values together with the requested metadata.
#' @importFrom dplyr group_by mutate ungroup
.get_heatmap_data <- function(
        tse,
        assay.type,
        row.var = NULL,
        col.var = NULL,
        scale = FALSE,
        center = FALSE,
        ...
) {

    df <- meltSE(
        tse,
        assay.type = assay.type,
        add.row = row.var,
        add.col = col.var
    )

    # Scale and/or center each feature across samples
    if (scale || center) {
        df <- df |>
            group_by(FeatureID) |>
            mutate(
                !!assay.type := as.numeric(
                    scale(
                        .data[[assay.type]],
                        center = center,
                        scale = scale
                    )
                )
            ) |>
            ungroup()
    }

    return(df)
}


# Create the base heatmap.
#' @importFrom ggplot2 aes element_rect element_text geom_tile ggplot labs
#'     scale_fill_gradient2 theme theme_minimal
.plot_heatmap <- function(
        df,
        fill,
        scale = FALSE,
        ...
) {

    # Create the heatmap
    p <- ggplot(
        df,
        aes(
            x = SampleID,
            y = FeatureID,
            fill = .data[[fill]]
        )
    ) +
        geom_tile(
            width = 0.95,
            height = 0.95,
            colour = "white",
            linewidth = 0.2
        )

    # Apply the heatmap theme
    p <- p +
        theme_minimal() +
        theme(
            strip.background = element_rect(
                fill = "white",
                colour = "black",
                linewidth = 0.5
            ),
            strip.text = element_text(
                face = "bold",
                colour = "black"
            )
        )

    # Remove axis titles
    p <- p +
        labs(
            x = NULL,
            y = NULL
        )

    # Apply the fill scale
    p <- .resolve_heatmap_scale(
        p = p,
        values = df[[fill]],
        name = fill,
        scale = scale
    )

    return(p)
}


# Apply the fill scale.

#' @importFrom ggplot2 scale_fill_gradient2
.resolve_heatmap_scale <- function(
        p,
        values,
        name,
        scale,
        ...
) {

    if (scale) {
        p <- p +
            scale_fill_gradient2(
                low = "#2166AC",
                mid = "white",
                high = "#B2182B",
                midpoint = 0
            )
    } else {
        p <- .resolve_plot_colours(
            p,
            values,
            name,
            fill = TRUE
        )
    }

    return(p)
}


# Add row and/or column facets.
#' @importFrom ggplot2 facet_grid
.resolve_heatmap_facets <- function(
        p,
        row.var = NULL,
        col.var = NULL,
        scales = "free",
        ...
) {

    # Construct the facet formula
    rows <- if (is.null(row.var)) "." else paste(row.var, collapse = " + ")
    cols <- if (is.null(col.var)) "." else paste(col.var, collapse = " + ")

    form <- stats::as.formula(
        paste(rows, "~", cols)
    )

    # Use ggplot2 for simple facets and ggh4x for nested facets
    if (length(row.var) <= 1L && length(col.var) <= 1L) {
        p <- p +
            facet_grid(
                form,
                scales = scales
            )
    } else {
        .require_package("ggh4x")
        p <- p +
            ggh4x::facet_nested(
                form,
                scales = scales
            )
    }

    return(p)
}


# Combine the row tree and heatmap into a single plot.
#' @importFrom patchwork plot_layout
#' @importFrom ggplot2 theme
.add_row_tree_to_heatmap <- function(
        tse,
        p,
        tree.name,
        layout = "rectangular",
        tree.width = 0.2,
        ...
) {

    # Draw the row tree. Only rectangular layout is available.
    p_tree <- plotRowTree(
        tse,
        tree.name = tree.name,
        layout = "rectangular",
        ...
    )

    # Combine the tree and heatmap
    p <- (p_tree + p) +
        plot_layout(
            widths = c(tree.width, 1),
            guides = "collect"
        ) &
        theme(
            legend.position = "right"
        )

    return(p)
}
