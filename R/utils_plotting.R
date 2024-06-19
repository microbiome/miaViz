#' Additional arguments for plotting
#' 
#' To be able to fine tune plotting, several additional plotting arguments are
#' available. These are described on this page.
#' 
#' @section Tree plotting:
#' 
#' \describe{
#'   \item{\code{line_alpha}:}{Numeric scalar in \code{[0, 1]}, specifying the 
#'     transparency of the tree edges. Defaults to \code{1}.}
#'   \item{\code{line_width}:}{Numeric scalar, specifying the default width of 
#'     an edge. Defaults to NULL to use default of the \code{ggtree} package}
#'   \item{\code{line_width_range}:}{Two numeric values, the range for plotting
#'     dynamic edge widths in. Defaults to \code{c(0.5,3)}.}
#'   \item{\code{point.alpha}:}{Numeric scalar in \code{[0, 1]}, specifying the 
#'     transparency of the tips. Defaults to \code{1}.}
#'   \item{\code{point.size}:}{Numeric scalar, specifying the 
#'     default size of tips Defaults to \code{2.}.}
#'   \item{\code{point_size_range}:}{Two numeric values, the range for plotting
#'     dynamic tip sizes in. Defaults to \code{c(1,4)}.}
#'   \item{\code{label_font_size}:}{Numeric scalar, font size for the tip and 
#'     node labels. Default to \code{3}.}
#'   \item{\code{highlight_font_size}:}{Numeric scalar, font size for the 
#'     highlight labels. Default to \code{3}.}
#' }
#' 
#' @section Graph plotting:
#' 
#' \describe{
#'   \item{\code{line_alpha}:}{Numeric scalar in \code{[0, 1]}, specifying the 
#'     transparency of the tree edges. Defaults to \code{1}.}
#'   \item{\code{line_width}:}{Numeric scalar, specifying the default width of 
#'     an edge. Defaults to NULL to use default of the \code{ggraph} package}
#'   \item{\code{line_width_range}:}{Two numeric values, the range for plotting
#'     dynamic edge widths in. Defaults to \code{c(0.5,3)}.}
#'   \item{\code{point.alpha}:}{Numeric scalar in \code{[0, 1]}, specifying the 
#'     transparency of the tips. Defaults to \code{1}.}
#'   \item{\code{point.size}:}{Numeric scalar, specifying the 
#'     default size of tips Defaults to \code{2.}.}
#'   \item{\code{point_size_range}:}{Two numeric values, the range for plotting
#'     dynamic tip sizes in. Defaults to \code{c(1,4)}.}
#' }
#' 
#' @section Abundance plotting:
#' 
#' \describe{
#'   \item{\code{flipped}:}{Logical scalar. Should the plot be flipped? Defaults
#'     to \code{FALSE}.}
#'   \item{\code{add.legend}:}{Logical scalar. Should legends be plotted? 
#'     Defaults to \code{TRUE}.}
#'   \item{\code{add.x.text}:}{Logical scalar. Should x tick labels be plotted?
#'     Defaults to \code{FALSE}.}
#'   \item{\code{add_border}:}{Logical scalar. Should border of bars be plotted?
#'     Defaults to \code{FALSE}.}
#'   \item{\code{bar_alpha}:}{Numeric scalar in \code{[0, 1]}, specifying the 
#'     transparency of the bars. Defaults to \code{1}.}
#'   \item{\code{point.alpha}:}{Numeric scalar in \code{[0, 1]}, specifying the 
#'     transparency of the points. Defaults to \code{1}.}
#'   \item{\code{point.size}:}{Numeric scalar, specifying the 
#'     default size of points. Defaults to \code{2.}.}
#' }
#' 
#' @section Abundance density plotting:
#' 
#' \describe{
#'   \item{\code{add.legend}:}{Logical scalar. Should legends be plotted? 
#'     Defaults to \code{TRUE}.}
#'   \item{\code{point_shape}:}{Numeric scalar setting the shape of points. 
#'     Defaults to \code{21}.}
#'   \item{\code{point_colour}:}{Character scalar, specifying the 
#'     default colour of points. Defaults to \code{2.}.}
#'   \item{\code{point.size}:}{Numeric scalar, specifying the 
#'     default size of points. Defaults to \code{2.}.}
#'   \item{\code{point.alpha}:}{Numeric scalar in \code{[0, 1]}, specifying the 
#'     transparency of the points. Defaults to \code{1}.}
#'   \item{\code{flipped}:}{Logical scalar. Should the plot be flipped? Defaults
#'     to \code{FALSE}.}
#'   \item{\code{scales_free}:}{Logical scalar. Should \code{scales = "free"} be
#'     set for faceted plots? Defaults to \code{TRUE}.}
#'   \item{\code{angle_x_text}:}{Logical scalar. Should x tick labels be plotted?
#'     Defaults to \code{FALSE}.}
#' }
#' 
#' @section Prevalence plotting:
#' 
#' \describe{
#'   \item{\code{flipped}:}{Logical scalar, specifying whether the plot should
#'     be flipped. Defaults to \code{FALSE}.}
#'   \item{\code{add.legend}:}{Logical scalar. Should legends be plotted? 
#'     Defaults to \code{TRUE}.}
#'   \item{\code{point.alpha}:}{Numeric scalar in \code{[0, 1]}, specifying the 
#'     transparency of the tips. Defaults to \code{1}.}
#'   \item{\code{point.size}:}{Numeric scalar, specifying the 
#'     default size of tips Defaults to \code{2.}.}
#'   \item{\code{line_alpha}:}{Numeric scalar in \code{[0, 1]}, specifying the 
#'     transparency of the tree edges. Defaults to \code{1}.}
#'   \item{\code{line_type}:}{Numeric scalar, specifying the default line type.
#'     Defaults to NULL to use default of the \code{ggplot2} package}
#'   \item{\code{line_size}:}{Numeric scalar, specifying the default width of 
#'     a line. Defaults to NULL to use default of the \code{ggplot2} package}
#' }
#' 
#' @section Series plotting:
#' 
#' \describe{
#'   \item{\code{add.legend}:}{Logical scalar. Should legends be plotted? 
#'     Defaults to \code{TRUE}.}
#'   \item{\code{line_alpha}:}{Numeric scalar in \code{[0, 1]}, specifying the 
#'     transparency of the tree edges. Defaults to \code{1}.}
#'   \item{\code{line_type}:}{Numeric scalar, specifying the default line type.
#'     Defaults to NULL to use default of the \code{ggplot2} package}
#'   \item{\code{line_width}:}{Numeric scalar, specifying the default width of 
#'     a line. Defaults to NULL to use default of the \code{ggplot2} package}
#'   \item{\code{line_width_range}:}{Two numeric values, the range for plotting
#'     dynamic line widths in. Defaults to \code{c(0.5,3)}.}
#'   \item{\code{ribbon_alpha}:}{Numeric scalar in \code{[0, 1]}, specifying the 
#'     transparency of the ribbon. Defaults to \code{0.3}.}
#' }
#' 
#' @section Tile plotting:
#' 
#' \describe{
#'   \item{\code{add.legend}:}{Logical scalar. Should legends be plotted? 
#'     Defaults to \code{TRUE}.}
#'   \item{\code{rect_alpha}:}{Numeric scalar in \code{[0, 1]}, specifying the 
#'     transparency of the areas. Defaults to \code{1}.}
#'   \item{\code{rect_colour}:}{Character scalar, specifying the colour to use
#'     for colouring the borders of the areas. Defaults to \code{"black"}.}
#'   \item{\code{na.value}:}{Character scalar, specifying the colour to use
#'     for \code{NA} values. Defaults to \code{"grey80"}.}
#' }
#' 
#' @name mia-plot-args
NULL

.get_palette <- scater:::.get_palette
# Adjusted function originally developed for scater package by Aaron Lun
#' @importFrom viridis scale_fill_viridis scale_colour_viridis
#' @importFrom ggplot2 scale_fill_manual scale_colour_manual
.resolve_plot_colours <- function(plot.out = plot_out, plot_out, colour.by = colout_by, 
                                  colour_by,
                                  colour.by.name = colour_by_name,
                                  colour_by_name,
                                  fill = FALSE,
                                  type = c("normal","edges"),
                                  na.translate = TRUE,
                                  na.value = NA,
                                  ...) 
{
    if (is.null(colour.by)) {
        return(plot_out)
    }
    type <- match.arg(type)
    if(type == "normal"){
        if (fill) {
            VIRIDFUN <- scale_fill_viridis
            SCALEFUN <- scale_fill_manual
        }
        else {
            VIRIDFUN <- scale_colour_viridis
            SCALEFUN <- scale_colour_manual
        }
        option <- "D"
    } else if(type == "edges") {
        if (fill) {
            VIRIDFUN <- scale_edge_fill_viridis
            SCALEFUN <- scale_edge_fill_manual
        }
        else {
            VIRIDFUN <- scale_edge_colour_viridis
            SCALEFUN <- scale_edge_colour_manual
        }
        option <- "C"
    } else {
        stop("Unrecognized colour type")
    }
    if (is.numeric(colour.by)) {
        plot_out <- plot_out + VIRIDFUN(name = colour_by_name, option = option,
                                        na.value = na.value,
                                        ...
                                        )
    }
    else {
        nlevs_colour_by <- nlevels(as.factor(colour.by))
        if (nlevs_colour_by <= 10) {
            plot_out <- plot_out + SCALEFUN(values = .get_palette("tableau10medium"), 
                                            name = colour_by_name,
                                            na.translate = na.translate,
                                            na.value = na.value)
        }
        else {
            if (nlevs_colour_by > 10 && nlevs_colour_by <= 20) {
                plot_out <- plot_out + SCALEFUN(values = .get_palette("tableau20"), 
                                                name = colour_by_name,
                                                na.translate = na.translate,
                                                na.value = na.value)
            }
            else {
                plot_out <- plot_out + VIRIDFUN(name = colour_by_name, 
                                                discrete = TRUE,
                                                na.translate = na.translate,
                                                option = option,
                                                na.value = na.value)
            }
        }
    }
    plot_out
}

.add_extra_guide <- scater:::.add_extra_guide

.add_extra_guide_graph <- function(plot.out = plot_out, plot_out, 
    edge.width.by = edge_width_by, edge_width_by){
    guide_args <- list()
    if (!is.null(edge.width.by)) {
        guide_args$edge_width <- guide_legend(title = edge.width.by)
        plot_out <- plot_out + 
            do.call(guides, guide_args)
    }
    plot_out
}

#' @importFrom ggnewscale new_scale
.add_extra_guide_tree <- function(lot.out = plot_out, plot_out, 
    edge.size.by = edge_size_by, edge_size_by,
    line.width.range = line_width_range, line_width_range){
    if (!is.null(edge_size_by)) {
        if(is.numeric(plot_out$data$edge_size_by)){
            SIZEFUN <- scale_size_continuous
        } else {
            SIZEFUN <- scale_size_discrete
        }
        plot_out <- plot_out + 
            SIZEFUN(name = edge_size_by, range = line_width_range) +
            new_scale("size")
    }
    plot_out
}

.na_replace_from_plot_data <- function(x,
                                       edge.size.by = edge_size_by,
                                       edge_size_by = NULL,
                                       shape.by = shape_by,
                                       shape_by = NULL,
                                       size.by = size_by,
                                       size_by = NULL,
                                       default.shape = default_shape,
                                       default_shape = 21,
                                       default.size = default_size,
                                       default_size = 0,
                                       default.edge.size = default_edge_size,
                                       default_edge_size = 0){
    if(!is.null(shape.by)){
        x$shape.by[is.na(x$shape.by)] <- default_shape
    }
    if(!is.null(size.by)){
        x$size.by[is.na(x$size.by)] <- default_size
    }
    if(!is.null(edge_size_by)){
        x$edge_size_by[is.na(x$edge_size_by)] <- default_edge_size
    }
    retun(x)
}

.get_bar_args <- function (colour.by = colour_by, colour_by, alpha = 0.65, 
    add.border = add_border, add_border = NULL, n = 0) 
{
    fill_colour <- TRUE
    border <- FALSE
    aes_args <- list()

    if (!is.null(colour.by)) {
        aes_args$fill <- "colour.by"
    }
    if(!is.null(add_border) && add_border && !is.null(colour.by)){
        border <- TRUE
        aes_args$colour <- "colour.by"
    } else if(is.null(add_border) && n <= 20) {
        border <- TRUE
        aes_args$colour <- "colour.by"
    }
    aes_args <- lapply(aes_args, function(x) if (!is.null(x)) sym(x))
    new_aes <- do.call(aes, aes_args)
    geom_args <- list(mapping = new_aes, alpha = alpha)
    if (is.null(colour.by)) {
        geom_args$colour <- "grey20"
    }
    return(list(args = geom_args, fill = fill_colour, border = border))
}


# Adjusted function originally developed for scater package by Aaron Lun
.get_point_args <- function(colour.by = colour_by, colour_by, 
                            shape.by = shape_by, shape_by, 
                            size.by = size_by, size_by,
                            alpha = 0.65,
                            size = NULL,
                            shape = 21,
                            colour = "grey70") 
{
    aes_args <- list()
    fill_colour <- TRUE
    if (!is.null(shape.by)) {
        aes_args$shape <- "shape.by"
    }
    if (!is.null(colour.by)) {
        # Only shapes 21 to 25 can be filled. Filling does not work in other shapes.
        if(shape >= 21 && shape <= 25){
            aes_args$fill <- "colour.by"
        } else {
            aes_args$colour <- "colour.by"
            fill_colour <- FALSE
        }
    }
    if (!is.null(size.by)) {
        aes_args$size <- "size.by"
    }
    aes_args <- lapply(aes_args, function(x) if (!is.null(x)) sym(x))
    new_aes <- do.call(aes, aes_args)
    geom_args <- list(mapping = new_aes, alpha = alpha)
    if (is.null(colour.by)) {
        geom_args$fill <- colour
    }
    if (is.null(shape.by)) {
        geom_args$shape <- shape
    }
    if (is.null(size.by)) {
        geom_args$size <- size
    }
    return(list(args = geom_args, fill = fill_colour))
}

.get_line_args <- function(colour.by = colour_by, colour_by,
                           linetype.by = linetype_by, linetype_by,
                           size.by = size_by, size_by,
                           alpha = 0.65,
                           linetype = 1,
                           linewidth = NULL,
                           colour = "grey70") 
{
    aes_args <- list()
    if (!is.null(linetype_by)) {
        aes_args$linetype <- "linetype_by"
    }
    if (!is.null(colour.by)) {
        aes_args$colour <- "colour.by"
    }
    if (!is.null(size.by)) {
        aes_args$linewidth <- "size.by"
    }
    aes_args <- lapply(aes_args, function(x) if (!is.null(x)) sym(x))
    new_aes <- do.call(aes, aes_args)
    geom_args <- list(mapping = new_aes, alpha = alpha)
    if (is.null(colour.by)) {
        geom_args$colour <- colour
    }
    if (is.null(linetype_by)) {
        geom_args$linetype <- linetype
    }
    if (is.null(size.by)) {
        geom_args$linewidth <- linewidth
    }
    return(list(args = geom_args))
}

.get_ribbon_args <- function(colour.by = colour_by, colour_by,
                             alpha = 0.3) 
{
    aes_args <- aes(ymin = .data[["Y"]] - .data[["sd"]], ymax = .data[["Y"]] + .data[["sd"]])
    if (!is.null(colour.by)) {
        aes_args$fill <- substitute(`colour.by`)
    }
    new_aes <- do.call(aes, aes_args)
    geom_args <- list(mapping = new_aes, alpha = alpha)
    if (is.null(colour.by)) {
        geom_args$fill <- "grey70"
    }
    return(list(args = geom_args))
}

.get_edge_args <- function(edge.colour.by = edge_colour_by, edge_colour_by, 
    edge.size.by = edge_size_by, edge_size_by, alpha = 1, size = NULL,
    layout = NULL){
    aes_args <- list()
    if (!is.null(edge.colour.by)) {
        aes_args$colour <- "edge.colour.by"
    }
    if (!is.null(edge_size_by)) {
        aes_args$size <- "edge_size_by"
    }
    aes_args <- lapply(aes_args, function(x) if (!is.null(x)) sym(x))
    new_aes <- do.call(aes, aes_args)
    geom_args <- list(mapping = new_aes, alpha = alpha)
    if (is.null(edge.colour.by)) {
        geom_args$colour <- "black"
    }
    if (is.null(edge_size_by)) {
        geom_args$size <- size
    }
    # Add layout if specified
    if( !is.null(layout) ){
        geom_args$layout <- layout 
    }
    return(list(args = geom_args))
}

.get_graph_edge_args <- function(edge.colour.by = edge_colour_by, edge_colour_by,
    edge.width.by = edge_width_by, edge_width_by, alpha = 1,
    size = NULL, edge.type = edge_type, edge_type){
    edge_args <- .get_edge_args(edge.colour.by, edge.width.by, alpha, size)
    if (!is.null(edge.width.by)) {
        edge_args$args$mapping$edge_width <- sym("edge.width.by")
        edge_args$args$mapping$size <- NULL
    } else {
        edge_args$args$edge_width <- size
        edge_args$args$size <- NULL
    }
    edge_args
}

.get_rect_args <- function(colour.by = colour_by, colour_by, alpha = 1, colour = "black"){
    aes_args <- list()
    if (!is.null(colour.by)) {
        aes_args$fill <- "colour.by"
    }
    aes_args <- lapply(aes_args, function(x) if (!is.null(x)) sym(x))
    new_aes <- do.call(aes, aes_args)
    geom_args <- list(mapping = new_aes, alpha = alpha, colour = colour)
    return(list(args = geom_args))
}

.get_density_args <- function(colour.by = colour_by, colour_by, alpha = 0.65, colour = "black") {
    fill_colour <- TRUE
    aes_args <- list()
    if (!is.null(colour.by)) {
        aes_args$colour <- "colour.by"
        aes_args$fill <- "colour.by"
    }
    aes_args <- lapply(aes_args, function(x) if (!is.null(x)) sym(x))
    new_aes <- do.call(aes, aes_args)
    geom_args <- list(mapping = new_aes,
                      alpha = alpha)
    if (is.null(colour.by)) {
        geom_args$colour <- colour
        geom_args$fill <- "grey70"
    }
    return(list(args = geom_args, fill = fill_colour))
}

#' @importFrom ggplot2 coord_flip element_blank element_text
.flip_plot <- function(plot.out = plot_out, plot_out, flipped = FALSE, 
    add.x.text = add_x_text, add_x_text = FALSE,
    angle.x.text = angle_x_text, angle_x_text = TRUE){
    if (flipped) {
        plot_out <- plot_out + 
            coord_flip()
        if(!add.x.text){
            plot_out <- plot_out +
                theme(axis.text.y = element_blank(),
                      axis.ticks.y = element_blank())
        } else if(angle_x_text) {
            plot_out <- plot_out +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        }
    } else {
        if(!add.x.text){
            plot_out <- plot_out +
                theme(axis.text.x = element_blank(),
                      axis.ticks.x = element_blank())
        }
    }
    plot_out
}

#' @importFrom ggplot2 theme
.add_legend <- function(plot.out = plot_out, plot_out, 
    add.legend = add_legend, add_legend, position = c("right","bottom")){
    position <- match.arg(position)
    if(!add.legend){
        plot_out <- plot_out +
            theme(legend.position = "none")
    } else {
        if(position == "bottom"){
            plot_out <- plot_out +
                theme(legend.position = "bottom")
        }
    }
    plot_out
}
