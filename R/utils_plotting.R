
.coerce_to_factor <- scater:::.coerce_to_factor
.get_point_args <- scater:::.get_point_args
.resolve_plot_colours <- scater:::.resolve_plot_colours
.add_extra_guide <- scater:::.add_extra_guide

.na_replace_from_plot_data <- function(object,
                                       shape_by = NULL,
                                       size_by = NULL,
                                       default_shape = 19,
                                       default_size = 1){
    if(!is.null(shape_by)){
        object$shape_by[is.na(object$shape_by)] <- default_shape
    }
    if(!is.null(size_by)){
        object$size_by[is.na(object$size_by)] <- default_size
    }
    object
}
