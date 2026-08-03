#' @name
#' plotOrdination
#'
#' @title
#' Create ordination plot
#'
#' @description
#' Ordinaton plotter
#'
#' @details
#' Creates ordination plot
#'
#' @return
#' A \code{ggplot2} object.
#'
#' @param x a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#'
#' @param dimred \code{character scalar}. Name of the reduced dimension result
#' stored in \code{reducedDim(x)} to visualize.
#'
#' @param ... Additional parameters for plotting.
#' \itemize{
#'   \item \code{ncomponents}: \code{integer vector} of length 2 or
#'   \code{integer scalar}. Specifies which ordination components are plotted.
#'   If a scalar is provided, the first two components are used.
#'   (Default: \code{2L})
#'
#'   \item \code{colour.by}: \code{NULL} or \code{character scalar}. Specifies a
#'   variable from \code{colData(x)} or a feature from \code{rownames(x)} used
#'   to colour observations. Feature abundances are taken from
#'   \code{assay.type}. (Default: \code{NULL})
#'
#'   \item \code{fill.by}: \code{NULL} or \code{character scalar}. Specifies a
#'   variable from \code{colData(x)} or a feature from \code{rownames(x)} used
#'   to fill observations or ellipses. Feature abundances are taken from
#'   \code{assay.type}. Cannot be used together with
#'   \code{add.density = TRUE}. (Default: \code{NULL})
#'
#'   \item \code{shape.by}: \code{NULL} or \code{character scalar}. Specifies a
#'   categorical variable from \code{colData(x)} used for point shapes.
#'   (Default: \code{NULL})
#'
#'   \item \code{size.by}: \code{NULL} or \code{character scalar}. Specifies a
#'   variable from \code{colData(x)} used for point sizes.
#'   (Default: \code{NULL})
#'
#'   \item \code{group.by}: \code{NULL} or \code{character scalar}. Specifies a
#'   categorical variable from \code{colData(x)} used for grouping when drawing
#'   ellipses, centroids and centroid vectors. (Default: \code{NULL})
#'
#'   \item \code{linetype.by}: \code{NULL} or \code{character scalar}.
#'   Specifies a categorical variable from \code{colData(x)} used for ellipse
#'   line types. (Default: \code{NULL})
#'
#'   \item \code{pair.by}: \code{NULL} or \code{character scalar}. Specifies a
#'   variable from \code{colData(x)} identifying observations that should be
#'   connected by lines. (Default: \code{NULL})
#'
#'   \item \code{sort.by}: \code{NULL} or \code{character scalar}. Specifies a
#'   variable from \code{colData(x)} used to order observations before drawing
#'   connecting lines. (Default: \code{NULL})
#'
#'   \item \code{facet.by}: \code{NULL} or \code{character scalar}. Specifies a
#'   categorical variable from \code{colData(x)} used to split the plot into
#'   facets. (Default: \code{NULL})
#'
#'   \item \code{assay.type}: \code{character scalar}. Name of the assay used
#'   when \code{colour.by} or \code{fill.by} specifies a feature.
#'   (Default: \code{"counts"})
#'
#'   \item \code{add.points}: \code{logical scalar}. Whether to draw sample
#'   points. (Default: \code{TRUE})
#'
#'   \item \code{add.ellipse}: \code{logical scalar}. Whether to draw confidence
#'   ellipses around groups. (Default: \code{FALSE})
#'
#'   \item \code{add.density}: \code{logical scalar}. Whether to draw a
#'   two-dimensional density estimate in the background.
#'   (Default: \code{FALSE})
#'
#'   \item \code{add.centroids}: \code{logical scalar}. Whether to draw group
#'   centroids. (Default: \code{FALSE})
#'
#'   \item \code{add.centroids.lines}: \code{logical scalar}. Whether to connect
#'   observations to their group centroids. (Default: \code{FALSE})
#'
#'   \item \code{add.vectors}: \code{logical scalar}. Whether to draw vectors
#'   from the global centroid to group centroids.
#'   (Default: \code{FALSE})
#'
#'   \item \code{add.rotation}: \code{logical scalar}. Whether to draw rotation
#'   (species score) coordinates if available in the ordination result.
#'   (Default: \code{add.species})
#'
#'   \item \code{add.species}: \code{logical scalar}. Alias for
#'   \code{add.rotation}. (Default: \code{FALSE})
#'
#'   \item \code{add.expl.var}: \code{logical scalar}. Whether to append the
#'   percentage of explained variance to the axis labels when available.
#'   (Default: \code{FALSE})
#'
#'   \item \code{scales}: \code{character scalar}. Scaling used for faceted
#'   plots. Passed to \code{ggplot2::facet_wrap()}. (Default:
#'   \code{"fixed"})
#'
#'   \item \code{xlab}, \code{ylab}: \code{character scalar}. Axis labels.
#'   Defaults to the ordination component names.
#'
#'   \item \code{panel.by.eigen}: \code{logical scalar}. Whether to scale the
#'   panel aspect ratio according to the eigenvalues of the ordination when
#'   available. (Default: \code{TRUE})
#'
#'   \item \code{point.shape}: Shape used for points.
#'   (Default: \code{19})
#'
#'   \item \code{point.alpha}: \code{numeric scalar}. Transparency of points.
#'   Must be between 0 and 1. (Default: \code{0.4})
#'
#'   \item \code{ellipse.alpha}: \code{numeric scalar}. Transparency of ellipse
#'   fills. Must be between 0 and 1. (Default: \code{0.2})
#'
#'   \item \code{ellipse.linewidth}: \code{numeric scalar}. Line width of
#'   ellipse borders. (Default: \code{0.5} or \code{0} when
#'   \code{fill.by} is specified.)
#'
#'   \item \code{ellipse.linetype}: Integer specifying the ellipse line type.
#'   (Default: \code{1})
#'
#'   \item \code{confidence.level}: \code{numeric scalar}. Confidence level used
#'   for ellipse calculation. Must be between 0 and 1.
#'   (Default: \code{0.95})
#'
#'   \item \code{adjust}: \code{numeric scalar}. Multiplicative adjustment for
#'   the bandwidth used in the background density estimate.
#'   (Default: \code{1})
#' }
#'
#' @examples
#' data("Tito2024QMP")
#' tse <- Tito2024QMP
#'
#' # Compute relative abundances and an MDS ordination
#' tse <- transformAssay(tse, method = "relabundance")
#' tse <- addMDS(tse, assay.type = "relabundance", method = "bray", ncomponents = 50)
#'
#' # Basic ordination plot
#' plotOrdination(tse, "MDS")
#'
#' # Colour samples by a sample-level variable
#' plotOrdination(tse, "MDS", colour.by = "diagnosis")
#'
#' # Colour samples by the abundance of a single feature
#' plotOrdination(
#'     tse, "MDS",
#'     colour.by = rownames(tse)[1],
#'     assay.type = "relabundance")
#'
#' # Use multiple aesthetics simultaneously
#' plotOrdination(
#'     tse, "MDS",
#'     colour.by = "diagnosis",
#'     shape.by = "colonoscopy"
#' )
#'
#' # Add confidence ellipses
#' plotOrdination(
#'     tse, "MDS",
#'     colour.by = "diagnosis",
#'     group.by = "diagnosis",
#'     add.ellipse = TRUE
#' )
#'
#' # Fill ellipses instead of colouring them
#' plotOrdination(
#'     tse, "MDS",
#'     fill.by = "diagnosis",
#'     group.by = "diagnosis",
#'     add.ellipse = TRUE
#' )
#'
#' # Show explained variance in the axis labels
#' plotOrdination(
#'     tse, "MDS",
#'     colour.by = "diagnosis",
#'     add.expl.var = TRUE
#' )
#'
#' # Add group centroids
#' plotOrdination(
#'     tse, "MDS",
#'     colour.by = "diagnosis",
#'     group.by = "diagnosis",
#'     add.centroids = TRUE
#' )
#'
#' # Connect samples to their group centroids
#' plotOrdination(
#'     tse, "MDS",
#'     colour.by = "diagnosis",
#'     group.by = "diagnosis",
#'     add.centroids.lines = TRUE
#' )
#'
#' # Show vectors from the global centroid to each group centroid
#' plotOrdination(
#'     tse, "MDS",
#'     colour.by = "diagnosis",
#'     group.by = "diagnosis",
#'     add.vectors = TRUE
#' )
#'
#' # Draw a background density estimate
#' plotOrdination(
#'     tse, "MDS",
#'     add.density = TRUE
#' )
#'
#' # Split the plot into facets
#' plotOrdination(
#'     tse, "MDS",
#'     colour.by = "diagnosis",
#'     facet.by = "colonoscopy"
#' )
#'
#' # Plot different ordination components
#' plotOrdination(
#'     tse, "MDS",
#'     ncomponents = c(2, 3)
#' )
#'
#' # Customize point appearance
#' plotOrdination(
#'     tse, "MDS",
#'     colour.by = "diagnosis",
#'     point.shape = 17,
#'     point.alpha = 0.8
#' )
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[scater:plotReducedDim]{scater::plotReducedDim}}
#'   \item \code{\link[=plotCCA]{plotCCA}}
#' }
#'
NULL

#' @rdname plotOrdination
#' @export
setMethod("plotOrdination", signature = c(x = "SingleCellExperiment"),
    function(x, dimred, colour.by = color.by, color.by = NULL, ...){
        args <- .check_ordination_input(x, dimred, colour.by = colour.by, ...)
        df <- do.call(.get_ordination_data, args)
        p <- .ordination_plotter(df, ...)
        return(p)
    }
)

################################ HELP FUNCTIONS ################################

# This method checks that the input is in correct format
.check_ordination_input <- function(
        x, dimred,
        ncomponents = 2L,
        colour.by = color.by, color.by = NULL,
        fill.by = NULL,
        shape.by = NULL,
        size.by = NULL,
        group.by = NULL,
        linetype.by = NULL,
        pair.by = NULL,
        sort.by = NULL,
        facet.by = NULL,
        assay.type = "counts",
        add.points = TRUE, add.ellipse = FALSE, add.density = FALSE,
        add.centroids = FALSE, add.centroids.lines = FALSE, add.vectors = FALSE,
        add.rotation = add.species, add.species = FALSE, add.expl.var = FALSE,
        ...){
    # Check that dimred can be found
    .check_dimred_present(dimred, x)
    # Check that ncomponents is correct. We can only visualize 2 components.
    if( .is_an_integer(ncomponents) ){
        ncomponents <- seq_len(ncomponents)
    }
    if( !(.is_integer(ncomponents) && length(ncomponents) == 2L &&
            all(ncomponents > 0L & ncomponents <= ncol(reducedDim(x, dimred))))
        ){
        stop("'ncomponents' must specify columns from reducedDim(x, dimred) ",
            "with integer values.", call. = FALSE)
    }

    # Check aesthetic variables
    temp <- .check_metadata_variable(x, colour.by, FALSE, TRUE, FALSE, TRUE)
    temp <- .check_metadata_variable(x, fill.by, FALSE, TRUE, FALSE, TRUE)
    temp <- .check_metadata_variable(x, shape.by, FALSE, TRUE, FALSE, FALSE)
    temp <- .check_metadata_variable(x, size.by, FALSE, TRUE, FALSE, FALSE)
    temp <- .check_metadata_variable(x, linetype.by, FALSE, TRUE, FALSE, FALSE)
    temp <- .check_metadata_variable(x, group.by, FALSE, TRUE, FALSE, FALSE)
    temp <- .check_metadata_variable(x, pair.by, FALSE, TRUE, FALSE, FALSE)
    temp <- .check_metadata_variable(x, sort.by, FALSE, TRUE, FALSE, FALSE)
    temp <- .check_metadata_variable(x, facet.by, FALSE, TRUE, FALSE, FALSE)
    # If colour.by specifies rowname, we check assay.type as the abundance
    # values are used for coloring
    if( !is.null(colour.by) && colour.by %in% rownames(x) ){
        temp <- .check_assay_present(assay.type, x)
    }
    # Check other flags
    if( !.is_a_bool(add.points) ){
        stop("'add.points' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(add.ellipse) ){
        stop("'add.ellipse' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(add.density) ){
        stop("'add.density' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(add.centroids) ){
        stop("'add.centroids' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(add.centroids.lines) ){
        stop("'add.centroids.lines' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(add.vectors) ){
        stop("'add.vectors' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(add.rotation) ){
        stop("'add.rotation' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(add.expl.var) ){
        stop("'add.expl.var' must be TRUE or FALSE.", call. = FALSE)
    }

    # add.density cannot be specified simultaneously with fill.by as they both
    # are using fill aesthetic and we can have only on fill scale.
    if( add.density && !is.null(fill.by) ){
        stop("Both 'add.density' and 'fill.by' cannot be specified ",
            "simultaneously.", call. = FALSE)
    }
    # Put all the arguments into list that be fed to the data retrieval function
    args <- list(
        x = x, dimred = dimred,
        ncomponents = ncomponents,
        colour.by = colour.by,
        fill.by = fill.by,
        shape.by = shape.by,
        size.by = size.by,
        group.by = group.by,
        linetype.by = linetype.by,
        pair.by = pair.by,
        sort.by = sort.by,
        facet.by = facet.by,
        assay.type = assay.type,
        add.expl.var = add.expl.var
    )
    args <- c(args, list(...))
    return(args)
}

# This function retrieves the data from reducedDim
.get_ordination_data <- function(
        x, dimred, ncomponents, colour.by, fill.by, shape.by, size.by, group.by,
        linetype.by, pair.by, sort.by, facet.by, assay.type,
        add.expl.var = FALSE, ...){
    # Get data and store the original attributes that might include rotation
    # data, for instance
    df <- reducedDim(x, dimred)
    orig_attributes <- attributes(df)
    orig_attributes <- orig_attributes[
        !names(orig_attributes) %in% c("dim", "dimnames") ]
    # Subset the data to include specified columns to plot
    df <- df[, ncomponents]
    # Add colnames if they are not present
    if( is.null(colnames(df)) ){
        colnames(df) <- paste0(dimred, ncomponents)
    }
    df <- df |> as.data.frame()
    # Take names of the columns that will be plotted
    x_var <- colnames(df)[[1L]]
    y_var <- colnames(df)[[2L]]

    # Get rotation data and put it in correct format
    rotation <- NULL
    rotation_names <- c("rotation", "species")
    if( any(rotation_names %in% names(orig_attributes)) ){
        rotation_names <- rotation_names[
            rotation_names %in% names(orig_attributes)][[1L]]
        rotation <- orig_attributes[[rotation_names]] |> as.data.frame()
        rotation <- rotation[, ncomponents, drop = FALSE]
        colnames(rotation) <- colnames(df)
    }

    xlab <- x_var
    ylab <- y_var
    eigen_names <- c("eig", "percentVar")
    if( add.expl.var && any(eigen_names %in% names(orig_attributes)) ){
        eigen_names <- eigen_names[
            eigen_names %in% names(orig_attributes)][[1L]]
        eigen <- orig_attributes[[eigen_names]]
        if( eigen_names %in% c("eig") ){
            eigen <- eigen * 100
        }
        xlab <- paste0(
            xlab, " (",
            round(eigen[ncomponents][[1L]], 1),
            "%)"
        )
        ylab <- paste0(
            ylab, " (",
            round(eigen[ncomponents][[2L]], 1),
            "%)"
        )
    } else if( add.expl.var ){
        warning("No explained variance found from the data.", call. = FALSE)
    }

    # List data that is fetched from colData
    cols <- c(
        shape.by, size.by, group.by, linetype.by, pair.by, sort.by, facet.by)
    # colour.by and fill.by can also specify a feature and its abundance.
    # Either get their column name to fetch from colData or directly add them
    # to data.
    if( !is.null(colour.by) && colour.by %in% colnames(colData(x)) ){
        cols <- c(cols, colour.by)
    } else if( !is.null(colour.by) ){
        df[[colour.by]] <- assay(x, assay.type)[colour.by, ]
    }
    if( !is.null(fill.by) && fill.by %in% colnames(colData(x)) ){
        cols <- c(cols, fill.by)
    } else if( !is.null(fill.by) ){
        df[[fill.by]] <- assay(x, assay.type)[fill.by, ]
    }
    # Get colData and merge it with the data
    cols <- cols |> unique()
    cd <- colData(x)[, cols, drop = FALSE]
    df <- cbind(df, cd)

    # Check that the values are correct
    if( !is.null(group.by) && is.numeric(df[[group.by]]) ){
        stop("Values specified by 'group.by' must be categorical.",
            call. = FALSE)
    }
    if( !is.null(linetype.by) && is.numeric(df[[linetype.by]]) ){
        stop("Values specified by 'linetype.by' must be categorical.",
            call. = FALSE)
    }
    if( !is.null(facet.by) && is.numeric(df[[facet.by]]) ){
        stop("Values specified by 'facet.by' must be categorical.",
            call. = FALSE)
    }

    # Sort the data. For instance, if we want to add arrows between consecutive
    # time points, this is essential step.
    if( !is.null(sort.by) ){
        df <- df[order(df[[sort.by]]), , drop = FALSE]
    }

    # Calculate centroid data
    grouping_var <- c(
        colour.by, facet.by, group.by, fill.by) |> unique()
    df_centroids <- df |>
        group_by(across(all_of(grouping_var))) |>
        summarise(
            x_centroid = mean(.data[[x_var]], na.rm = TRUE),
            y_centroid = mean(.data[[y_var]], na.rm = TRUE),
            .groups = "drop"
        )
    # Add global mean
    df_centroids[["x_global"]] <- mean(df[[x_var]], na.rm = TRUE)
    df_centroids[["y_global"]] <- mean(df[[y_var]], na.rm = TRUE)

    # Add additional information to attributes of df.
    attributes(df) <- c(
        attributes(df),
        x = x_var,
        y = y_var,
        colour.by = colour.by,
        fill.by = fill.by,
        shape.by = shape.by,
        size.by = size.by,
        group.by = group.by,
        linetype.by = linetype.by,
        pair.by = pair.by,
        sort.by = sort.by,
        facet.by = facet.by,
        assay.type = assay.type,
        xlab = xlab,
        ylab = ylab
    )
    attr(df, "rotation") <- rotation
    attr(df, "centroids") <- df_centroids
    attr(df, "eigen") <- eigen
    return(df)
}

# This method is the main plotter function.
.ordination_plotter <- function(
        df, scales = "fixed", add.points = TRUE, add.ellipse = FALSE,
        add.density = FALSE, add.centroids = FALSE, add.centroids.lines = FALSE,
        add.vectors = FALSE, add.rotation = add.species, add.species = FALSE,
        ...){
    # Initialize the plot
    p <- ggplot(df, aes(
        x = .data[[attributes(df)[["x"]]]],
        y = .data[[attributes(df)[["y"]]]],
        colour = if( !is.null(attributes(df)[["colour.by"]]) )
            .data[[attributes(df)[["colour.by"]]]]
    ))
    grouping_var <- c(
        attributes(df)[["colour.by"]], attributes(df)[["facet.by"]],
        attributes(df)[["group.by"]], attributes(df)[["fill.by"]]) |> unique()
    # Add points, i.e., samples
    if( add.points ){
        p <- .add_ordination_points(p, df, ...)
    }
    # Add ellipses
    if( add.ellipse && length(grouping_var) > 0L ){
        p <- .add_ellipse_to_ordination(p, df, ...)
    }
    # Connect samples of subjects, groups etc
    if( !is.null(attributes(df)[["pair.by"]]) ){
        p <- .connect_ordination_points(p, df)
    }
    # Add density to background
    if( add.density ){
        # Background density
        p <- .add_point_density(p, df)
    }
    # Add group centroids
    if( add.centroids && length(grouping_var) > 0L ){
        p <- .add_centroids(p, df, ...)
    }
    # Add lines connecting points and group centroids
    if( add.centroids.lines && length(grouping_var) > 0L ){
        p <- .add_centroid_lines(p, df, grouping_var, ...)
    }
    # Add species scores
    if( add.rotation ){
        p <- .add_rotation(p, df)
    }
    # Add vectors from global centroid to group centroids
    if( add.vectors && length(grouping_var) > 0L ){
        p <- .add_centroids_vector(p, df, grouping_var, ...)
    }
    # If facetting was specified, split plot to separate panels
    if( !is.null(attributes(df)[["facet.by"]]) ){
        p <- p +
            facet_wrap(
                ~ .data[[attributes(df)[["facet.by"]]]],
                scales = scales
            )
    }
    # Adjust theme
    p <- .adjust_ordination_theme(p, df, ...)
    return(p)
}

# Add points for ordination plot
.add_ordination_points <- function(
        p, df, point.shape = 19, point.alpha = 0.4, ...){
    p <-  .add_points_layer(
        p, df,
        x = attributes(df)[["x"]],
        y = attributes(df)[["y"]],
        point.shape = point.shape,
        point.alpha = point.alpha,
        ...)
    return(p)
}

# This function adds ellipse visualization.
.add_ellipse_to_ordination <- function(
        p, df,
        ellipse.alpha = 0.2,
        ellipse.linewidth = if(is.null(attributes(df)[["fill.by"]])) 0.5 else 0,
        ellipse.linetype = 1,
        confidence.level = 0.95,
        ...){
    # To disable "no visible binding for global variable" message in cmdcheck
    color <- NULL
    #
    if( !.are_whole_numbers(ellipse.linetype) ){
        stop("'vec.linetype' must be a whole number.", call. = FALSE)
    }
    if ( !(is.numeric(ellipse.alpha) && ellipse.alpha > 0 &&
           ellipse.alpha < 1 ) ) {
        stop("'ellipse.alpha' must be a number between 0 and 1.", call. = FALSE)
    }
    if ( !(is.numeric(ellipse.linewidth) && ellipse.linewidth >= 0) ) {
        stop("'ellipse.linewidth' must be a positive number.", call. = FALSE)
    }
    if( !(is.numeric(confidence.level) && confidence.level > 0 &&
          confidence.level < 1) ) {
        stop("'confidence.level' must be a number between 0 and 1.",
             call. = FALSE)
    }
    #
    # Get all the arguments. User can fill and colour the ellipses separately.
    # However, in most of the cases that might not make sense, but it is still
    # made possible.
    args <- list(
        mapping = aes(
            group = if( !is.null(attributes(df)[["group.by"]]) )
                .data[[attributes(df)[["group.by"]]]],
            colour = if( !is.null(attributes(df)[["colour.by"]]) )
                .data[[attributes(df)[["colour.by"]]]],
            fill = if( !is.null(attributes(df)[["fill.by"]]) )
                .data[[attributes(df)[["fill.by"]]]],
            linetype = if( !is.null(attributes(df)[["linetype.by"]]) )
                .data[[attributes(df)[["linetype.by"]]]],
        ),
        geom = "polygon",
        linewidth = ellipse.linewidth,
        level = confidence.level,
        alpha = if( is.null(attributes(df)[["fill.by"]]) ) 0 else ellipse.alpha
    )
    # If user did not specify coloring, add black border to ellipses
    if( is.null(attributes(df)[["colour.by"]]) ){

    }
    if( is.null(attributes(df)[["linetype.by"]]) ){
        args[["linetype"]] = ellipse.linetype
    }
    p <- p + do.call(stat_ellipse, args)
    return(p)
}

# This method connects points with a line or arrow.
.connect_ordination_points <- function(p, df){
    # Get arguments. If user sorted the data, use directed arrow.
    args <- list(
        mapping = aes(
            x = .data[[attributes(df)[["x"]]]],
            y = .data[[attributes(df)[["y"]]]],
            group = .data[[attributes(df)[["pair.by"]]]],
        ),
        arrow = if(!is.null(attributes(df)[["sort.by"]]))
            arrow(length = unit(0.2, "cm"), type = "closed"),
        alpha = 0.4
    )
    args <- args[ lengths(args) > 0 ]
    p <- p + do.call(geom_path, args)
    return(p)
}

# This method adds group
.add_centroids <- function(p, df, ...){
    df_centroids <- attributes(df)[["centroids"]]
    # Add centroids
    p <- p +
        geom_point(
            data = df_centroids,
            aes(x = x_centroid, y = y_centroid),
            size = 5,
            shape = 4,  # cross
            stroke = 2  # thick border
        )
    return(p)
}

# This method adds group
.add_centroid_lines <- function(p, df, grouping_var, ...){
    df_centroids <- attributes(df)[["centroids"]]
    # Add centroids data to original df
    df <- df %>%
        dplyr::left_join(df_centroids, by = grouping_var)
    # Connect points with centroids
    p <- p + geom_segment(data = df, aes(
        x = x_centroid, y = y_centroid,
        xend = .data[[attributes(df)[["x"]]]],
        yend = .data[[attributes(df)[["y"]]]]
    ), alpha = 0.4)
    return(p)
}

# This methods creates vectors that start from global mean and ends to group
# centroids. This shows how the covariate correlates with the ordination.
#' @importFrom ggrepel geom_label_repel
.add_centroids_vector <- function(p, df, grouping_var, ...){
    df_centroids <- attributes(df)[["centroids"]]
    # Visualize vectors
    p <- p + geom_segment(
        data = df_centroids,
        aes(
            x = x_global, y = y_global,
            xend = x_centroid, yend = y_centroid,
        ),
        arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
        size = 1
    )
    # Add label to denote which vector belongs to which group
    p <- p + geom_label_repel(
        data = df_centroids,
        mapping = aes(
            x = x_centroid,
            y = y_centroid,
            label = .data[[grouping_var]])
    )
    return(p)
}

# This methods adds "species scores", i.e., the coordinates of species in the
# ordination.
.add_rotation <- function(p, df){
    # The points are added only if the rotation is present in the data
    if( !is.null(attributes(df)[["rotation"]]) ){
        df_species <- attributes(df)[["rotation"]]
        # Add points
        p <- p + geom_point(
            data = df_species,
            mapping = aes(
                x = .data[[attributes(df)[["x"]]]],
                y = .data[[attributes(df)[["y"]]]]
            ),
            colour = "red", size = 1)
    }
    return(p)
}

# This methods coors the backrounds based on the point density. This creates
# "landscape plots".
.add_point_density <- function(p, df, adjust = 1, ...){
    # Calculate the correct smoothing bandwidth. We could use the default
    # values, but this custom bandwidth helps us to highlight the density
    # better.
    bandwidth <- adjust * c(
        .get_bandwidth(df[[attributes(df)[["x"]]]]),
        .get_bandwidth(df[[attributes(df)[["y"]]]]))
    # Add the landscape
    p <- p +
        stat_density_2d(
            aes(fill = after_stat(density)),
            geom = "raster",
            h = bandwidth,
            contour = FALSE,
            alpha = 0.5,
        ) +
        scale_fill_gradient(name = "Density", low = "white", high = "black")
    return(p)
}

# This function calculates the smoothing bandwidth. It highlights better point-
# rich areas than the default choice.
.get_bandwidth <- function(x){
    r <- quantile(x, c(0.25, 0.75))
    4 * 1.06 * min(sd(x), (r[[2]] - r[[1]])/1.34) * length(x)^(-0.2)
}

# This function adjusts the theme of the plot
.adjust_ordination_theme <- function(
        p, df,
        xlab = attributes(df)[["xlab"]],
        ylab = attributes(df)[["ylab"]],
        coord.equal = TRUE,
        panel.by.eigen = TRUE,
        ...){
    if( !.is_a_string(xlab) ){
        stop("'xlab' must be a single character value.", call. = FALSE)
    }
    if( !.is_a_string(ylab) ){
        stop("'ylab' must be a single character value.", call. = FALSE)
    }
    if( !.is_a_bool(coord.equal) ){
        stop("'coord.equal' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(panel.by.eigen) ){
        stop("'panel.by.eigen' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    p <- p + theme_classic()
    p <- p + labs(x = xlab, y = ylab)
    if( !is.null(attributes(df)[["fill.by"]]) ){
        p <- p + labs(fill = attributes(df)[["fill.by"]])
    }
    if( !is.null(attributes(df)[["colour.by"]]) ){
        p <- p + labs(colour = attributes(df)[["colour.by"]])
    }
    if( !is.null(attributes(df)[["shape.by"]]) ){
        p <- p + labs(shape = attributes(df)[["shape.by"]])
    }
    if( !is.null(attributes(df)[["size.by"]]) ){
        p <- p + labs(shape = attributes(df)[["size.by"]])
    }
    if( !is.null(attributes(df)[["linetype.by"]]) ){
        p <- p + labs(linetype = attributes(df)[["linetype.by"]])
    }

    # Use equal scaling on both axes so that one unit on the x-axis has the same
    # physical length as one unit on the y-axis. This preserves distances and
    # angles in the ordination and avoids visual distortion.
    if( coord.equal ){
        p <- p + coord_equal()
    }
    # Make the panel dimensions proportional to the eigenvalues so that the
    # relative lengths of the axes reflect the variation explained by each
    # ordination axis.
    if( panel.by.eigen && !is.null(attributes(df)[["eigen"]]) ){
        eig <- attributes(df)[["eigen"]]
        p <- p + theme(aspect.ratio = eig[2] / eig[1])
    }

    # Adjust colors
    if( !is.null(attributes(df)[["colour.by"]]) ){
        name <- attributes(df)[["colour.by"]]
        vals <- df[[name]]
        p <- .resolve_plot_colours(
            p, vals, name,
        )
    }
    if( !is.null(attributes(df)[["fill.by"]]) ){
        name <- attributes(df)[["fill.by"]]
        vals <- df[[name]]
        p <- .resolve_plot_colours(
            p, vals, name, fill = TRUE
        )
    }

    return(p)
}
