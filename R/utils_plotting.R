
.resolve_plot_colours <- scater:::.resolve_plot_colours
.add_extra_guide <- scater:::.add_extra_guide

.na_replace_from_plot_data <- function(object,
                                       edge_size_by = NULL,
                                       shape_by = NULL,
                                       size_by = NULL,
                                       default_shape = 21,
                                       default_size = 2){
    if(!is.null(shape_by)){
        object$shape_by[is.na(object$shape_by)] <- default_shape
    }
    if(!is.null(size_by)){
        object$size_by[is.na(object$size_by)] <- default_size
    }
    object
}

# Adjusted function originally developed for scater package by Aaron Lun
.get_point_args <- function (colour_by, shape_by, size_by, alpha = 0.65, size = NULL) 
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
