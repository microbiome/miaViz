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
#'   \item{\code{point_alpha}:}{Numeric scalar in \code{[0, 1]}, specifying the 
#'     transparency of the tips. Defaults to \code{1}.}
#'   \item{\code{point_size}:}{Numeric scalar, specifying the 
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
#'   \item{\code{point_alpha}:}{Numeric scalar in \code{[0, 1]}, specifying the 
#'     transparency of the tips. Defaults to \code{1}.}
#'   \item{\code{point_size}:}{Numeric scalar, specifying the 
#'     default size of tips Defaults to \code{2.}.}
#'   \item{\code{point_size_range}:}{Two numeric values, the range for plotting
#'     dynamic tip sizes in. Defaults to \code{c(1,4)}.}
#' }
#' 
#' @section Abundance plotting:
#' 
#' \describe{
#'   \item{\code{flipped}:}{Logical scalar. Should the plot be flipped. Defaults
#'     to \code{FALSE}.}
#'   \item{\code{add_legend}:}{Logical scalar. Should legends be plotted? 
#'     Defaults to \code{TRUE}.}
#'   \item{\code{add_x_text}:}{Logical scalar. Should x tick labels be plotted?
#'     Defaults to \code{FALSE}.}
#'   \item{\code{add_border}:}{Logical scalar. Should border of bars be plotted?
#'     Defaults to \code{FALSE}.}
#'   \item{\code{bar_alpha}:}{Numeric scalar in \code{[0, 1]}, specifying the 
#'     transparency of the bars. Defaults to \code{1}.}
#'   \item{\code{point_alpha}:}{Numeric scalar in \code{[0, 1]}, specifying the 
#'     transparency of the tips. Defaults to \code{1}.}
#'   \item{\code{point_size}:}{Numeric scalar, specifying the 
#'     default size of tips Defaults to \code{2.}.}
#' }
#' 
#' @section Prevalence plotting:
#' 
#' \describe{
#'   \item{\code{flipped}:}{Logical scalar, specifying whether the plot should
#'     be flipped. Defaults to \code{FALSE}.}
#'   \item{\code{add_legend}:}{Logical scalar. Should legends be plotted? 
#'     Defaults to \code{TRUE}.}
#'   \item{\code{point_alpha}:}{Numeric scalar in \code{[0, 1]}, specifying the 
#'     transparency of the tips. Defaults to \code{1}.}
#'   \item{\code{point_size}:}{Numeric scalar, specifying the 
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
#'   \item{\code{add_legend}:}{Logical scalar. Should legends be plotted? 
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
#'   \item{\code{add_legend}:}{Logical scalar. Should legends be plotted? 
#'     Defaults to \code{TRUE}.}
#'   \item{\code{rect_alpha}:}{Numeric scalar in \code{[0, 1]}, specifying the 
#'     transparency of the areas. Defaults to \code{1}.}
#'   \item{\code{rect_colour}:}{Character scalar, specfiying the colour to use
#'     for colouring the borders of the areas. Defaults to \code{"black"}.}
#'   \item{\code{na.value}:}{Character scalar, specfiying the colour to use
#'     for \code{NA} values. Defaults to \code{"grey80"}.}
#' }
#' 
#' @name mia-plot-args
NULL

.get_palette <- scater:::.get_palette
# Adjusted function originally developed for scater package by Aaron Lun
#' @importFrom viridis scale_fill_viridis scale_colour_viridis
#' @importFrom ggplot2 scale_fill_manual scale_colour_manual
.resolve_plot_colours <- function(plot_out, colour_by, colour_by_name,
                                  fill = FALSE,
                                  type = c("normal","edges"),
                                  na.translate = TRUE,
                                  na.value = NA) 
{
    if (is.null(colour_by)) {
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
    if (is.numeric(colour_by)) {
        plot_out <- plot_out + VIRIDFUN(name = colour_by_name, option = option,
                                        na.value = na.value)
    }
    else {
        nlevs_colour_by <- nlevels(as.factor(colour_by))
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

.add_extra_guide_graph <- function(plot_out, edge_width_by){
    guide_args <- list()
    if (!is.null(edge_width_by)) {
        guide_args$edge_width <- guide_legend(title = edge_width_by)
        plot_out <- plot_out + 
            do.call(guides, guide_args)
    }
    plot_out
}

#' @importFrom ggnewscale new_scale
.add_extra_guide_tree <- function(plot_out, edge_size_by, line_width_range){
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

.na_replace_from_plot_data <- function(object,
                                       edge_size_by = NULL,
                                       shape_by = NULL,
                                       size_by = NULL,
                                       default_shape = 21,
                                       default_size = 0,
                                       default_edge_size = 0){
    if(!is.null(shape_by)){
        object$shape_by[is.na(object$shape_by)] <- default_shape
    }
    if(!is.null(size_by)){
        object$size_by[is.na(object$size_by)] <- default_size
    }
    if(!is.null(edge_size_by)){
        object$edge_size_by[is.na(object$edge_size_by)] <- default_edge_size
    }
    object
}

.get_bar_args <- function (colour_by, alpha = 0.65, add_border = NULL,
                           n = 0) 
{
    aes_args <- list()
    fill_colour <- TRUE
    border <- FALSE
    if (!is.null(colour_by)) {
        aes_args$fill <- "colour_by"
    }
    if(!is.null(add_border) && add_border && !is.null(colour_by)){
        border <- TRUE
        aes_args$colour <- "colour_by"
    } else if(is.null(add_border) && n <= 20) {
        border <- TRUE
        aes_args$colour <- "colour_by"
    }
    new_aes <- do.call(aes_string, aes_args)
    geom_args <- list(mapping = new_aes, alpha = alpha)
    if (is.null(colour_by)) {
        geom_args$colour <- "grey20"
    }
    return(list(args = geom_args, fill = fill_colour, border = border))
}


# Adjusted function originally developed for scater package by Aaron Lun
.get_point_args <- function(colour_by, shape_by, size_by, alpha = 0.65,
                            size = NULL, shape = 21) 
{
    aes_args <- list()
    fill_colour <- TRUE
    if (!is.null(shape_by)) {
        aes_args$shape <- "shape_by"
    }
    if (!is.null(colour_by)) {
        # Only shapes 21 to 25 can be filled. Filling does not work in other shapes.
        if(shape >= 21 && shape <= 25){
            aes_args$fill <- "colour_by"
        } else {
            aes_args$colour <- "colour_by"
            fill_colour <- FALSE
        }
    }
    if (!is.null(size_by)) {
        aes_args$size <- "size_by"
    }
    new_aes <- do.call(aes_string, aes_args)
    geom_args <- list(mapping = new_aes, alpha = alpha)
    if (is.null(colour_by)) {
        geom_args$fill <- "grey70"
    }
    if (is.null(shape_by)) {
        geom_args$shape <- shape
    }
    if (is.null(size_by)) {
        geom_args$size <- size
    }
    return(list(args = geom_args, fill = fill_colour))
}

.get_line_args <- function(colour_by, linetype_by, size_by,
                           alpha = 0.65,
                           linetype = 1,
                           size = NULL) 
{
    aes_args <- list()
    if (!is.null(linetype_by)) {
        aes_args$linetype <- "linetype_by"
    }
    if (!is.null(colour_by)) {
        aes_args$colour <- "colour_by"
    }
    if (!is.null(size_by)) {
        aes_args$size <- "size_by"
    }
    new_aes <- do.call(aes_string, aes_args)
    geom_args <- list(mapping = new_aes, alpha = alpha)
    if (is.null(colour_by)) {
        geom_args$colour <- "grey70"
    }
    if (is.null(linetype_by)) {
        geom_args$linetype <- linetype
    }
    if (is.null(size_by)) {
        geom_args$size <- size
    }
    return(list(args = geom_args))
}

.get_ribbon_args <- function(colour_by,
                             alpha = 0.3) 
{
    aes_args <- list()
    aes_args$ymin <- "Y - sd"
    aes_args$ymax <- "Y + sd"
    if (!is.null(colour_by)) {
        aes_args$fill <- "colour_by"
    }
    new_aes <- do.call(aes_string, aes_args)
    geom_args <- list(mapping = new_aes, alpha = alpha)
    if (is.null(colour_by)) {
        geom_args$fill <- "grey70"
    }
    return(list(args = geom_args))
}

.get_edge_args <- function(edge_colour_by, edge_size_by, alpha = 1, size = NULL){
    aes_args <- list()
    if (!is.null(edge_colour_by)) {
        aes_args$colour <- "edge_colour_by"
    }
    if (!is.null(edge_size_by)) {
        aes_args$size <- "edge_size_by"
    }
    new_aes <- do.call(aes_string, aes_args)
    geom_args <- list(mapping = new_aes, alpha = alpha)
    if (is.null(edge_colour_by)) {
        geom_args$colour <- "black"
    }
    if (is.null(edge_size_by)) {
        geom_args$size <- size
    }
    return(list(args = geom_args))
}

.get_graph_edge_args <- function(edge_colour_by, edge_width_by, alpha = 1,
                                 size = NULL, edge_type){
    edge_args <- .get_edge_args(edge_colour_by, edge_width_by, alpha, size)
    if (!is.null(edge_width_by)) {
        edge_args$args$mapping$edge_width <- sym("edge_width_by")
        edge_args$args$mapping$size <- NULL
    } else {
        edge_args$args$edge_width <- size
        edge_args$args$size <- NULL
    }
    edge_args
}

.get_rect_args <- function(colour_by, alpha = 1, colour = "black"){
    aes_args <- list()
    if (!is.null(colour_by)) {
        aes_args$fill <- "colour_by"
    }
    new_aes <- do.call(aes_string, aes_args)
    geom_args <- list(mapping = new_aes, alpha = alpha, colour = colour)
    return(list(args = geom_args))
}

.get_density_args <- function(colour_by, alpha = 0.65, colour = "black") {
    fill_colour <- TRUE
    aes_args <- list()
    if (!is.null(colour_by)) {
        aes_args$colour <- "colour_by"
        aes_args$fill <- "colour_by"
    }
    new_aes <- do.call(aes_string, aes_args)
    geom_args <- list(mapping = new_aes,
                      alpha = alpha)
    if (is.null(colour_by)) {
        geom_args$colour <- colour
        geom_args$fill <- "grey70"
    }
    return(list(args = geom_args, fill = fill_colour))
}

#' @importFrom ggplot2 coord_flip element_blank element_text
.flip_plot <- function(plot_out, flipped = FALSE, add_x_text = FALSE,
                       angle_x_text = TRUE){
    if (flipped) {
        plot_out <- plot_out + 
            coord_flip()
        if(!add_x_text){
            plot_out <- plot_out +
                theme(axis.text.y = element_blank(),
                      axis.ticks.y = element_blank())
        } else if(angle_x_text) {
            plot_out <- plot_out +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        }
    } else {
        if(!add_x_text){
            plot_out <- plot_out +
                theme(axis.text.x = element_blank(),
                      axis.ticks.x = element_blank())
        }
    }
    plot_out
}

#' @importFrom ggplot2 theme
.add_legend <- function(plot_out, add_legend, position = c("right","bottom")){
    position <- match.arg(position)
    if(!add_legend){
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
