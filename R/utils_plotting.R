
.resolve_plot_colours <- scater:::.resolve_plot_colours
.add_extra_guide <- scater:::.add_extra_guide

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

.get_bar_args <- function (fill_by, colour_by, alpha = 0.65, add_border = NULL,
                           n = 0) 
{
    aes_args <- list()
    fill_colour <- TRUE
    border <- FALSE
    if (!is.null(fill_by)) {
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
.get_point_args <- function (colour_by, shape_by, size_by, alpha = 0.65,
                             size = NULL) 
{
    aes_args <- list()
    fill_colour <- TRUE
    if (!is.null(shape_by)) {
        aes_args$shape <- "shape_by"
    }
    if (!is.null(colour_by)) {
        aes_args$fill <- "colour_by"
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
        geom_args$shape <- 21
    }
    if (is.null(size_by)) {
        geom_args$size <- size
    }
    return(list(args = geom_args, fill = fill_colour))
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

.get_hightlight_args <- function(highlights, colour_highlights = FALSE,
                                 extendto = 0.52){
    aes_args <- list()
    aes_args$subset <- paste0("node %in% c(",paste(highlights, collapse = ","),
                              ")")
    if (colour_highlights) {
        aes_args$fill <- ~label
    }
    new_aes <- do.call(aes_, aes_args)
    geom_args <- list(mapping = new_aes)
    geom_args$extendto = extendto
    if (!colour_highlights) {
        geom_args$fill <- "grey70"
    }
    if (colour_highlights) {
        geom_args$colour <- "grey40"
    }
    if (!colour_highlights) {
        geom_args$colour <- "grey20"
    }
    return(list(args = geom_args))
}

#' @importFrom ggplot2 coord_flip element_blank element_text
.flip_plot <- function(plot_out, flipped = FALSE, add_x_text = FALSE){
    if (flipped) {
        plot_out <- plot_out + 
            coord_flip()
        if(!add_x_text){
            plot_out <- plot_out +
                theme(axis.text.y = element_blank(),
                      axis.ticks.y = element_blank())
        }
    } else {
        if(!add_x_text){
            plot_out <- plot_out +
                theme(axis.text.x = element_blank(),
                      axis.ticks.x = element_blank())
        } else {
            plot_out <- plot_out +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        }
    }
    plot_out
}

#' @importFrom ggplot2 theme
.add_legend <- function(plot_out, add_legend){
    if(!add_legend){
        plot_out <- plot_out +
            theme(legend.position = "none")
    }
    plot_out
}

#' @importFrom cowplot plot_grid
.combine_plots <- function(plots, flipped = FALSE, ...){
    if(flipped){
        plot_out <- plot_grid(plotlist = rev(plots), nrow=1, align="h",
                              axis = "tb",
                              rel_widths = c(rep(1, length(plots) - 1L), 2))
    } else {
        plot_out <- plot_grid(plotlist = plots, ncol=1, align="v",
                              axis = "lr",
                              rel_heights = c(2, rep(1, length(plots) - 1L)))
    }
    plot_out
}
