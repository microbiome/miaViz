#' Plot factor data as tiles
#'
#' Relative relations of two grouping can be visualized by plotting tiles with
#' relative sizes. \code{plotColTile} and \code{plotRowTile} can be used for
#' this.
#'
#' @param object a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#'
#' @param x \code{Character scalar}. Specifies the column-level metadata field
#' to show on the x-axis.
#'   Alternatively, an \link{AsIs} vector or data.frame, see
#'   \code{?\link{retrieveFeatureInfo}} or \code{?\link{retrieveCellInfo}}. Must
#'   result in a returned \code{character} or \code{factor} vector.
#'
#' @param y \code{Character scalar}. Specifies the column-level metadata to
#' show on the y-axis.
#'   Alternatively, an \link{AsIs} vector or data.frame, see
#'   \code{?\link{retrieveFeatureInfo}} or \code{?\link{retrieveCellInfo}}. Must
#'   result in a returned \code{character} or \code{factor} vector.
#'
#' @param ... additional arguments for plotting. See
#' \code{\link{mia-plot-args}} for more details i.e. call
#' \code{help("mia-plot-args")}
#'
#' @return
#' A \code{ggplot2} object or \code{plotly} object, if more than one
#' \code{prevalences} was defined.
#'
#' @name plotColTile
#'
#' @examples
#' data(GlobalPatterns)
#' se <- GlobalPatterns
#' plotColTile(se,"SampleType","Primer")
NULL

#' @rdname plotColTile
#' @export
setMethod("plotColTile", signature = c("SummarizedExperiment"),
    function(object, x, y, ...){
        .plot_tile_data(object, type = "column", x, y, ...)
    }
)

#' @rdname plotColTile
#' @export
setMethod("plotRowTile", signature = c("SummarizedExperiment"),
    function(object, x, y, ...){
        .plot_tile_data(object, type = "row", x, y, ...)
    }
)

.get_tile_data <- function(object, type, x, y){
    retrieve_FUN <- switch(
        type,
        "row" = retrieveFeatureInfo,
        "column" = retrieveCellInfo)
    retrieve_search <- switch(
        type,
        "row" = "rowData",
        "column" = "colData")
    #
    x_by_out <- retrieve_FUN(object, x, search = retrieve_search)
    x_lab <- x_by_out$name
    y_by_out <- retrieve_FUN(object, y, search = retrieve_search)
    y_lab <- y_by_out$name
    #
    if(!is.factor(x_by_out$value) && !is.character(x_by_out$value)){
        stop("'x' must specify a factor or character vector.", call. = FALSE)
    }
    if(!is.factor(y_by_out$value) && !is.character(y_by_out$value)){
        stop("'y' must specify a factor or character vector.", call. = FALSE)
    }
    #
    list(
        data = data.frame(
            X = factor(x_by_out$value),
            Y = factor(y_by_out$value)),
        x_lab = x_lab,
        y_lab = y_lab)
}

#' @importFrom dplyr group_by mutate summarise ungroup n
.summarise_tile_data <- function(object, data, type){
    retrieve_FUN <- switch(
        type,
        "row" = retrieveFeatureInfo,
        "column" = retrieveCellInfo)
    retrieve_search <- switch(
        type,
        "row" = "rowData",
        "column" = "colData")
    x_group <- data %>%
        group_by(.data$X) %>%
        summarise(group_n = n()) %>%
        mutate(
            group_freq = .data$group_n/sum(.data$group_n),
            x = cumsum(.data$group_freq),
            xmin = c(0,.data$x[-length(.data$x)]))
    data <- data %>%
        group_by(.data$X, .data$Y) %>%
        summarise(fill_n = n(), .groups = "rowwise") %>%
        dplyr::left_join(x_group, by = "X") %>%
        mutate(fill_freq = .data$fill_n/.data$group_n) %>%
        group_by(.data$X) %>%
        mutate(
            y = cumsum(.data$fill_freq),
            ymin = c(0,.data$y[-length(.data$y)])) %>%
        ungroup()
    data
}

.plot_tile_data <- function(object, type = c("row", "column"), x, y, ...){
    type <- match.arg(type)
    tile_out <- .get_tile_data(object, type, x, y)
    tile_data <- tile_out$data
    xlab <- tile_out$x_lab
    ylab <- tile_out$y_lab
    tile_data <- .summarise_tile_data(object, tile_data, type)
    tile_data$colour_by <- tile_data$Y
    .tile_plotter(tile_data, xlab = xlab, ylab = ylab, ...)
}

.get_xcoord_mid <- function(data){
    data %>%
        ungroup() %>%
        select(.data$X, .data$x, .data$xmin) %>%
        unique() %>%
        mutate(xmid = .data$xmin + (.data$x - .data$xmin)/2 )
}

.tile_plotter <- function(
        data,
        xlab,
        ylab,
        add_legend = add.legend,
        add.legend = TRUE,
        rect_alpha = rect.alpha,
        rect.alpha = 1,
        rect_colour = rect.colour,
        rect.colour = "black",
        na.value = "grey80"){
    coord <- .get_xcoord_mid(data)
    # get plotting arguments for rect
    rect_args <- .get_rect_args(
        colour_by = ylab, alpha = rect_alpha, colour = rect_colour)
    rect_args$args$mapping$xmin <- sym("xmin")
    rect_args$args$mapping$xmax <- sym("x")
    rect_args$args$mapping$ymin <- sym("ymin")
    rect_args$args$mapping$ymax <- sym("y")
    # start plotting
    plot_out <- ggplot(data)
    plot_out <- plot_out +
        do.call(geom_rect,rect_args$args)
    # add scales
    plot_out <- plot_out +
        scale_x_continuous(
            name = paste0("Fraction (",xlab,")"),
            expand = c(0,0),
            breaks = seq(0,1,0.1),
            sec.axis = dup_axis(
                name = xlab, breaks = coord$xmid, labels = coord$X)) +
        scale_y_continuous(
            name = paste0("Fraction (",ylab,")"),
            expand = c(0,0),
            breaks = seq(0,1,0.1))
    # resolve the fill colours
    plot_out <- .resolve_plot_colours(
        plot_out, data$colour_by, ylab, fill = TRUE, na.value = na.value)
    # add legend and theme
    plot_out <- plot_out +
        theme_classic()
    if (!add_legend) {
        plot_out <- plot_out + theme(legend.position = "none")
    }
    plot_out
}
