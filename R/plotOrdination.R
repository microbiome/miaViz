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
#' @param ... Additional parameters for plotting.
#' \itemize{
#'   \item \code{colour.by}: \code{NULL} or \code{character scalar}. Specifies a
#'   variable from \code{colData(x)} or \code{rowData(x)} which is used to
#'   colour observations. (Default: \code{NULL})
#' }
#'
#' @examples
#' data("Tito2024QMP")
#' tse <- Tito2024QMP
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[scater:plotReducedDim]{scater::plotReducedDim}}
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
        add.rotation = FALSE, add.expl.var = FALSE,
        ...){
    # Check if there are any reduced dim present
    if( length(reducedDims(x)) == 0L ){
        stop("No data present in reducedDim(x).", call. = FALSE)
    }
    # Check that dimred can be found
    is_name <- .is_a_string(dimred) && dimred %in% reducedDimNames(x)
    is_index <- .is_an_integer(dimred) && dimred > 0L &&
        dimred <= length(reducedDims(x))
    if( !(is_name || is_index) ){
        stop("'dimred' must specify data from reducedDim(x). It must be one ",
            "of the following options: '",
            paste0(reducedDimNames(x), collapse = "', '"), "'", call. = FALSE)
    }
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
    temp <- .check_metadata_variable(tse, colour.by, FALSE, TRUE, FALSE, TRUE)
    temp <- .check_metadata_variable(tse, fill.by, FALSE, TRUE, FALSE, TRUE)
    temp <- .check_metadata_variable(tse, shape.by, FALSE, TRUE, FALSE, FALSE)
    temp <- .check_metadata_variable(tse, size.by, FALSE, TRUE, FALSE, FALSE)
    temp <- .check_metadata_variable(
        tse, linetype.by, FALSE, TRUE, FALSE, FALSE)
    temp <- .check_metadata_variable(tse, group.by, FALSE, TRUE, FALSE, FALSE)
    temp <- .check_metadata_variable(tse, pair.by, FALSE, TRUE, FALSE, FALSE)
    temp <- .check_metadata_variable(tse, sort.by, FALSE, TRUE, FALSE, FALSE)
    temp <- .check_metadata_variable(tse, facet.by, FALSE, TRUE, FALSE, FALSE)
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

    # Add colnames if they are not present
    if( is.null(colnames(df)) ){
        colnames(df) <- paste0(dimred, ncomponents)
    }
    df <- df |> as.data.frame()
    # Take only 2 specified columns
    df <- df[, ncomponents]
    x_var <- colnames(df)[[1L]]
    y_var <- colnames(df)[[2L]]

    # Get rotation data adnd put it in correct format
    rotation <- NULL
    rotation_names <- c("rotation")
    if( any(rotation_names %in% names(orig_attributes)) ){
        rotation_names <- rotation_names[[1L]]
        rotation <- orig_attributes[[rotation_names]] |> as.data.frame()
        rotation <- rotation[, ncomponents, drop = FALSE]
        colnames(rotation) <- colnames(df)
    }

    expl_var_name <- c("eig")
    xlab <- x_var
    ylab <- y_var
    if( add.expl.var && any(expl_var_name %in% names(orig_attributes)) ){
        eigen <- orig_attributes[expl_var_name][[1L]]
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
    return(df)
}

# This method is the main plotter function.
.ordination_plotter <- function(
        df, scales = "fixed", add.points = TRUE, add.ellipse = FALSE,
        add.density = FALSE, add.centroids = FALSE, add.centroids.lines = FALSE,
        add.vectors = FALSE, add.rotation = FALSE, ...){
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
    p <- p + ggrepel::geom_label_repel(
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
    rot_names <- c("rotation")
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
        ...){
    if( !.is_a_string(xlab) ){
        stop("'xlab' must be a single character value.", call. = FALSE)
    }
    if( !.is_a_string(ylab) ){
        stop("'ylab' must be a single character value.", call. = FALSE)
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
    return(p)
}
