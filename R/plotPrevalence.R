#' Plot prevalence information
#'
#' \code{plotPrevalence} and \code{plotRowPrevalence} visualize prevalence
#' information.
#'
#' Whereas \code{plotPrevalence} produces a line plot, \code{plotRowPrevalence}
#' returns a heatmap.
#'
#' @param x a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#'
#' @param rank,... additional arguments
#' \itemize{
#'   \item as.relative \code{Logical scalar}. Should the relative values
#'   be calculated? (Default: \code{FALSE})
#'
#'   \item ndetection \code{Integer scalar}. Determines the number of breaks
#'   calculated detection thresholds when \code{detection=NULL}. When
#'   \code{TRUE},  \code{as_relative} is then also regarded as \code{TRUE}.
#'   (Default: \code{20})
#'
#'   \item{If \code{!is.null(rank)} matching arguments are passed on to
#'     \code{\link[=agglomerate-methods]{agglomerateByRank}}. See
#'     \code{\link[=agglomerate-methods]{?agglomerateByRank}} for more details.
#'   }
#'
#'   \item{additional arguments for plotting. See
#'   \code{\link{mia-plot-args}} for more details i.e. call
#'   \code{help("mia-plot-args")}}
#' }
#'
#'
#' @param assay.type \code{Character scalar}. Defines which assay data to
#'   use. (Default: \code{"relabundance"})
#'
#' @param assay_name Deprecated. Use \code{assay.type} instead.
#'
#' @param colour.by \code{Character scalar}. Specification of a feature to
#' colour points by, see the \code{by} argument in
#' \code{\link[scater:retrieveFeatureInfo]{?retrieveFeatureInfo}} for
#' possible values. Only used with \code{layout = "point"}.
#' (Default: \code{NULL})
#'
#' @param colour_by Deprecated. Use \code{colour.by} instead.
#'
#' @param shape.by \code{Character scalar}. Specification of a feature to shape
#' points by, see the \code{by} argument in
#' \code{\link[scater:retrieveFeatureInfo]{?retrieveFeatureInfo}} for
#' possible values. Only used with \code{layout = "point"}.
#' (Default: \code{NULL})
#'
#' @param shape_by Deprecated. Use \code{shape.by} instead.
#'
#' @param size.by \code{Character scalar}. Specification of a feature to size
#' points by, see the \code{by} argument in
#' \code{\link[scater:retrieveFeatureInfo]{?retrieveFeatureInfo}} for
#' possible values. Only used with \code{layout = "point"}.
#' (Default: \code{NULL})
#'
#' @param size_by Deprecated. Use \code{size.by} instead.
#'
#' @param facet.by \code{Character scalar}. Taxonomic rank to facet the plot by.
#' Value must be of \code{taxonomyRanks(x)}
#' Argument can only be used in function plotPrevalentAbundance.
#'
#' @param facet_by Deprecated. Use \code{facet.by} instead.
#'
#' @param show.label \code{Logical scalar}, \code{character scalar} or
#' \code{integer vector} for selecting labels from the rownames of \code{x}.
#' If \code{rank} is not \code{NULL} the rownames might change.
#' (Default: \code{NULL})
#'
#' @param label Deprecated. Use \code{show.label} instead.
#'
#' @param detection \code{Numeric scalar}. Detection thresholds for
#' absence/presence. Either an absolutes value compared directly to the values
#' of \code{x} or a relative value between 0 and 1, if \code{TRUE}.
#'
#' @param detections Deprecated. Use \code{detection} instead.
#'
#' @param prevalence \code{Numeric scalar}. Prevalence thresholds (in 0 to 1).
#' The required prevalence is strictly greater by default. To include the
#' limit, set \code{include.lowest} to \code{TRUE}.
#'
#' @param prevalences Deprecated. Use \code{prevalence} instead.
#'
#' @param min.prevalence \code{Numeric scalar}. Applied as a threshold for
#'   plotting. The threshold is applied per row and column.
#'   (Default: \code{0})
#'
#' @param min_prevalence Deprecated. Use \code{min.prevalence} instead.
#'
#' @param BPPARAM A
#'   \code{\link[BiocParallel:BiocParallelParam-class]{BiocParallelParam}}
#'   object specifying whether the UniFrac calculation should be parallelized.
#'
#' @details
#' Agglomeration on different taxonomic levels is available through the
#' \code{rank} argument.
#'
#' To exclude certain taxa, preprocess \code{x} to your liking, for example
#' with subsetting via \code{getPrevalent} or
#' \code{agglomerateByPrevalence}.
#'
#' @return
#' A \code{ggplot2} object or \code{plotly} object, if more than one
#' \code{prevalence} was defined.
#'
#' @seealso
#' \code{\link[mia:getPrevalence]{getPrevalence}},
#' \code{\link[mia:getPrevalence]{agglomerateByPrevalence}},
#' \code{\link[mia:agglomerate-methods]{agglomerateByRank}}
#'
#' @name plotPrevalence
#'
#' @examples
#' data(GlobalPatterns, package = "mia")
#'
#' # Apply relative transformation
#' GlobalPatterns <- transformAssay(GlobalPatterns, method = "relabundance")
#'
#' # plotting N of prevalence exceeding taxa on the Phylum level
#' plotPrevalence(GlobalPatterns, rank = "Phylum")
#' plotPrevalence(GlobalPatterns, rank = "Phylum") + scale_x_log10()
#'
#' # plotting prevalence per taxa for different detection thresholds as heatmap
#' plotRowPrevalence(GlobalPatterns, rank = "Phylum")
#'
#' # by default a continuous scale is used for different detection levels,
#' # but this can be adjusted
#' plotRowPrevalence(
#'     GlobalPatterns, rank = "Phylum", assay.type = "relabundance",
#'     detection = c(0, 0.001, 0.01, 0.1, 0.2))
#'
#' # point layout for plotRowPrevalence can be used to visualize by additional
#' # information
#' plotPrevalentAbundance(
#'     GlobalPatterns, rank = "Family", colour.by = "Phylum") +
#'     scale_x_log10()
#'
#' # When using function plotPrevalentAbundace, it is possible to create facets
#' # with 'facet.by'.
#' plotPrevalentAbundance(
#'     GlobalPatterns, rank = "Family",
#'     colour.by = "Phylum", facet.by = "Kingdom") +
#'     scale_x_log10()
NULL

################################################################################
# plotPrevalence

#' @rdname plotPrevalence
#' @export
setMethod("plotPrevalence", signature = c(x = "SummarizedExperiment"),
    function(
            x,
            detection = detections, detections = c(0.01, 0.1, 1, 2, 5, 10, 20),
            prevalence = prevalences, prevalences = seq(0.1, 1, 0.1),
            assay.type = assay_name, assay_name = "counts",
            rank = NULL,
            BPPARAM = BiocParallel::SerialParam(),
            ...){
        # input check
        if(!all(.is_numeric_string(detection))){
            stop("'detection' must be numeric values.", call. = FALSE)
        }
        if(!all(.is_numeric_string(prevalence)) || any(prevalence < 0) ||
                any(prevalence > 1)){
            stop("'prevalence' must be numeric values between 0 and 1.",
                call. = FALSE)
        }
        .check_assay_present(assay.type, x)
        #
        # Agglomerate data if specified
        if( !is.null(rank) ){
            x <- agglomerateByRank(x, rank = rank, ...)
        }
        # Get data to plot
        plot_data <- .get_prevalence_plot_data(
            x, assay.type, detection, prevalence, BPPARAM, ...)
        plot_data$colour_by <- plot_data$colour_by * 100
        # Plot the data
        p <- .prevalence_plotter(plot_data,
                            layout = "line",
                            ylab = "N",
                            colour_by = "Prevalence [%]",
                            size_by = NULL,
                            shape_by = NULL,
                            ...)
        return(p)
    }
)

# This function returns a number which tells the number of features which exceed
# the detection and prevalence thresholds.
.get_prevalence_count <- function(d, p, mat, ...){
    length(getPrevalent(mat, detection = d, prevalence = p, ...))
}

#' @importFrom BiocParallel bpmapply bpisup bpstart bpstop SerialParam
#' @importFrom SummarizedExperiment assay
.get_prevalence_plot_data <- function(
        x, assay.type, detections, prevalences,
        BPPARAM = BiocParallel::SerialParam(), as.relative = as_relative,
        as_relative = FALSE, ...){
    # Input check
    if(!.is_a_bool(as_relative)){
        stop("'as_relative' must be TRUE or FALSE.", call. = FALSE)
    }
    if(as_relative && (any(detections < 0) || any(detections > 1))){
        stop("If 'as_relative' == TRUE, detection' must be numeric ",
            "values between 0 and 1.", call. = FALSE)
    }
    #
    # Apply relative transform if specified
    if(as_relative){
        temp_name <- "temporary_relative_abundance"
        x <- transformAssay(
            x, assay.type = assay.type, method = "relabundance",
            name = temp_name)
        assay.type <- temp_name
    }
    # Get assay
    mat <- assay(x, assay.type, withDimnames = TRUE)
    # Get all combinations of different detection/prevalence thresholds
    ans <- expand.grid(detection = detections, prevalence = prevalences)
    # Start parallelization if specified
    if (!(bpisup(BPPARAM) || is(BPPARAM, "MulticoreParam"))) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add = TRUE)
    }
    # For each detection/prevalence threshold, get how many taxa exceed the
    # limit.
    temp <- bpmapply(
        .get_prevalence_count,
        ans$detection,
        ans$prevalence,
        MoreArgs = list(mat = mat, ...),
        BPPARAM = BPPARAM,
        SIMPLIFY = FALSE)
    # Add the values to result table
    ans[["n"]] <- unlist(temp)
    colnames(ans) <- c("X","colour_by","Y")
    return(ans)
}

################################################################################
# plotPrevalentAbundance

#' @rdname plotPrevalence
#' @export
setMethod("plotPrevalentAbundance", signature = c(x = "SummarizedExperiment"),
    function(
            x,
            rank = NULL,
            assay.type = assay_name, assay_name = "counts",
            colour.by = colour_by, colour_by = NULL,
            size.by = size_by,
            size_by = NULL,
            shape.by = shape_by,
            shape_by = NULL,
            show.label = label,
            label = NULL,
            facet.by = facet_by,
            facet_by = NULL,
            ...){
        # input check
        .check_assay_present(assay.type, x)

        # Check facet.by It is FALSE by default, but user can specify it, but
        # the value must be in taxonomyRanks.
        if(!(is.null(facet.by) || facet.by %in% taxonomyRanks(x))){
            stop("'facet.by' must be in taxonomyRanks.",  call. = FALSE)
        }
        #
        # Agglomerate data if specified
        if( !is.null(rank) ){
            x <- agglomerateByRank(x, rank = rank, ...)
        }
        # Check that labels are correct (or get rownames as labels if NULL)
        label <- .norm_label(show.label, x)
        # Get prevalence data
        plot_data <- .get_prevalence_plot_point_data(
            x, assay.type, label = label, ...)
        # Get data to plot
        vis_out <- .incorporate_prevalence_vis(
            plot_data,
            se = x,
            colour_by = colour.by,
            size_by = size.by,
            shape_by = shape.by,
            label = label,
            facet_by = facet.by)
        plot_data <- vis_out$df
        colour_by <- vis_out$colour_by
        size_by <- vis_out$size_by
        shape_by <- vis_out$shape_by
        facet_by <- vis_out$facet_by
        ylab <- paste0(
            "Prevalence (", ifelse(is.null(rank), "Features", rank), ") [%]")
        # Plot the data
        plot <- .prevalence_plotter(plot_data,
                            layout = "point",
                            ylab = ylab,
                            colour_by = colour_by,
                            size_by = size_by,
                            shape_by = shape_by,
                            ...)
        # If facet.by is not NULL, user has specified it. Adds the facets to
        # the plot.
        if(!is.null(facet.by)){
            plot <- plot +
                # Create facets
                facet_wrap(vars(!!sym("facet_by")))
        }
        return(plot)
    }
)

#' @importFrom DelayedArray rowMeans
#' @importFrom SummarizedExperiment assay
.get_prevalence_plot_point_data <- function(
        x, assay.type, as_relative = TRUE, label = NULL, ...){
    # Input check
    if(!.is_a_bool(as_relative)){
        stop("'as_relative' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    # Apply relative transform if specified
    if(as_relative){
        temp_name <- "temporary_relative_abundance"
        x <- transformAssay(
            x, assay.type = assay.type, method = "relabundance",
            name = temp_name)
        assay.type <- temp_name
    }
    # Get assay
    mat <- assay(x, assay.type, withDimnames = TRUE)
    # Calculate means and prevalences
    ans <- data.frame(
        X = rowMeans(mat, na.rm = TRUE),
        Y = getPrevalence(mat, detection = 0) * 100,
        label = rownames(mat))
    # label denotes, which rownames are displayed. Remove those labels/rownames
    # that have FALSE i.e. user do not want to display them.
    if(!is.null(label)){
        ans$label[!label] <- NA
    }
    return(ans)
}

#' @importFrom scater retrieveFeatureInfo
.incorporate_prevalence_vis <- function(
        plot_data,
        se = se,
        colour_by = NULL,
        size_by = NULL,
        shape_by = NULL,
        label = NULL,
        facet_by = NULL){
    # Get user specified options for plotting. These include the information
    # on which variable is used to colour etc the plot.
    variables <- c(
        colour_by = colour_by, size_by = size_by, shape_by = shape_by,
        facet_by = facet_by)
    # If user specified some options
    if(!is.null(variables)){
        # Loop through each option
        for(i in seq_along(variables)){
            # Get data from rowData
            feature_info <- retrieveFeatureInfo(
                se, variables[i], search = "rowData")
            # Mirror back variable name, if a partial match was used
            var_name <- names(variables)[i]
            assign(
                var_name, .get_new_var_name_value(get(var_name),
                feature_info$name))
            # Add data to plot data
            plot_data[, names(feature_info$name)] <- feature_info$value
            # label is for determining which labels are shown in the final
            # plot. Remove those labels that user do not want to show if the
            # variable has classes/factors/characters that will be plotted
            # with labels. (Numbers would not be).
            if(!is.null(label) && (is.factor(feature_info$value) ||
                    is.character(feature_info$value)) ){
                plot_data[, names(feature_info$name)][!label] <- NA
            }
        }
    }
    # Create a list to be returned
    res <- list(
        df = plot_data,
        colour_by = colour_by,
        size_by = size_by,
        shape_by = shape_by,
        facet_by = facet_by)
    return(res)
}

################################################################################
# plotRowPrevalence

#' @rdname plotPrevalence
#' @export
setMethod("plotRowPrevalence", signature = c(x = "SummarizedExperiment"),
    function(
            x,
            rank = NULL,
            assay.type = assay_name, assay_name = "counts",
            detection = detections, detections = c(0.01, 0.1, 1, 2, 5, 10, 20),
            min.prevalence = min_prevalence,
            min_prevalence = 0,
            BPPARAM = BiocParallel::SerialParam(),
            ...){
        # input check
        if(!all(.is_numeric_string(detection))){
            stop("'detection' must be numeric values.", call. = FALSE)
        }
        .check_assay_present(assay.type, x)

        if(length(min.prevalence) != 1 || !.is_numeric_string(min.prevalence)){
            stop("'min.prevalence' must be single numeric values.",
                call. = FALSE)
        }
        #
        # Agglomerate data if specified
        if( !is.null(rank) ){
            x <- agglomerateByRank(x, rank = rank, ...)
        }
        # Get prevalence data
        plot_data <- .get_prevalence_plot_matrix(
            x, assay.type, detection, min.prevalence, BPPARAM, ...)
        plot_data$colour_by <- plot_data$colour_by * 100
        ylab <- ifelse(is.null(rank), "Features", rank)
        colour_by <- "Prevalence [%]"
        # Plot the data
        p <- .prevalence_plotter(plot_data,
                            layout = "heatmap",
                            ylab = ylab,
                            colour_by = colour_by,
                            size_by = NULL,
                            shape_by = NULL,
                            ...)
        return(p)
    }
)

.is_continuous <- function(i){
    i <- sort(unique(i))
    z <- round(c(i[-1L],max(i)) - i,5L)
    length(unique(z[-length(z)])) == 1L
}


#' @importFrom BiocParallel bplapply bpisup bpstart bpstop SerialParam
#' @importFrom SummarizedExperiment assay
#' @importFrom tidyr pivot_longer
#' @importFrom DelayedArray rowSums
.get_prevalence_plot_matrix <- function(
        x, assay.type, detections, min_prevalence,
        BPPARAM = BiocParallel::SerialParam(), as.relative = as_relative,
        as_relative = FALSE,
        ndetection = 20, ...){
    # To disable "no visible binding for global variable" message in cmdcheck
    ID <- NULL
    # Input check
    if(!.is_a_bool(as_relative)){
        stop("'as_relative' must be TRUE or FALSE.", call. = FALSE)
    }
    if(as_relative && (any(detections < 0) || any(detections > 1))){
        stop("If 'as_relative' == TRUE, detection' must be numeric ",
            "values between 0 and 1.", call. = FALSE)
    }
    if( !.is_an_integer(ndetection) ){
        stop("'ndetection' must be a single integer value.", call. = FALSE)
    }
    if( is.null(detections) && is.null(ndetection) ){
        stop("Either 'detection' or 'ndetection' must be specified.",
            call. = FALSE)
    }
    # If detection thesholds were not specified, calculate them based
    # on ndetection. Because the values are set to relative scale, enable
    # relative abundances.
    if( is.null(detections) ){
        detections <- seq(0, 1, length.out = ndetection + 1L)
        as_relative <- TRUE
    }
    #
    # Apply relative transform if specified
    if(as_relative){
        temp_name <- "temporary_relative_abundance"
        x <- transformAssay(
            x, assay.type = assay.type, method = "relabundance",
            name = temp_name)
        assay.type <- temp_name
    }
    #
    # Start parallelling if specified
    if (!(bpisup(BPPARAM) || is(BPPARAM, "MulticoreParam"))) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add = TRUE)
    }
    # Calculate prevalences of taxa in different detection limits and create a
    # data.frame from them
    ans <- bplapply(
        detections, function(d){
            getPrevalence(x, assay.type = assay.type, detection = d)
        }, BPPARAM = BPPARAM)
    ans <- data.frame(ans)
    colnames(ans) <- detections
    # Remove those taxa and prevalence thresholds that do not exceed the minimum
    # prevalence threshold.
    f <- ans >= min_prevalence
    ans <- ans[rowSums(f) != 0, colSums(f) != 0, drop=FALSE]
    # If there are no data to plot anymro e after subsetting, give error.
    if(any(dim(ans) == 0)){
        stop("No data left after apply threshold 'min_prevalence'.",
            call. = FALSE)
    }
    # Get the taxa order, the most abundant taxa comes first
    lvls <- rownames(ans)[order(rowSums(ans))]
    # Add rownames to table
    ans[["ID"]] <- rownames(x)[rowSums(f) != 0]
    # Convert the table to long format
    ans <- ans %>%
        pivot_longer(!ID, names_to = "detection", values_to = "prevalence")
    colnames(ans) <- c("Y","X","colour_by")
    # Round values
    ans$X <- round(as.numeric(ans$X),4) * 100
    # If the values are discrete, convert them to factors
    if( !.is_continuous(ans$X) ){
        ans$X <- factor(
            ans$X, unique(as.character(sort(as.numeric(ans$X)))))
    }
    # Convert taxa to factor. The order is based on prevalence of taxa.
    ans$Y <- factor(ans$Y,lvls)
    return(ans)
}

#' @importFrom ggplot2 ggplot geom_raster geom_point geom_line
#'   scale_fill_distiller scale_y_discrete scale_x_discrete scale_x_continuous
#'   theme_classic
.prevalence_plotter <- function(plot_data,
        layout = c("line","point","heatmap"),
        xlab = paste0(ifelse(as.relative, "Rel. ", ""),"Abundance"),
        ylab,
        colour_by,
        size_by,
        shape_by,
        flipped = FALSE,
        add_legend = add.legend,
        add.legend = TRUE,
        point_alpha = point.alpha,
        point.alpha = 1,
        point_size = point.size,
        point.size = 2,
        line_alpha = line.alpha,
        line.alpha = 1,
        line_type = line.type,
        line.type = NULL,
        line_size = line.size,
        line.size = 1,
        as.relative = as_relative,
        as_relative = FALSE,
        ...){
    # Start plotting
    plot_out <- ggplot(plot_data, aes(x = .data[["X"]], y = .data[["Y"]])) +
        labs(x = xlab, y = ylab)
    # Add diffetent layouts
    if(layout == "line"){
        # Get point and line options
        point_args <- .get_point_args(
            colour_by = colour_by, shape_by = NULL, size_by = NULL,
            alpha = point_alpha, size = point_size)
        line_args <- .get_line_args(
            colour_by = colour_by, linetype_by = NULL, size_by = NULL,
            alpha = line_alpha, linetype = line_type, linewidth = line_size)
        # Add grouping. Otherwise, the line does not follow the same value as
        # colouring.
        point_args$args$mapping$group <- sym("colour_by")
        line_args$args$mapping$group <- sym("colour_by")
        # Add point and line layouts
        plot_out <- plot_out +
            do.call(geom_point, point_args$args) +
            do.call(geom_line, line_args$args)
        # Adjust the colour scale and legend title
        plot_out <- .resolve_plot_colours(
            plot_out, plot_data$colour_by, colour_by, fill = TRUE)
        plot_out <- .resolve_plot_colours(
            plot_out, plot_data$colour_by, colour_by, fill = FALSE)
    } else if(layout == "point"){
        # Get point options
        point_args <- .get_point_args(
            colour_by = colour_by, shape_by = shape_by, size_by = size_by,
            alpha = point_alpha, size = point_size)
        # Add point layout
        plot_out <- plot_out +
            do.call(geom_point, point_args$args)
        # Adjust the colour scale and legend title
        plot_out <- .resolve_plot_colours(
            plot_out, plot_data$colour_by, colour_by, fill = TRUE,
            na.translate = FALSE)
        # Add legend for shape_by and size_by if they were used.
        plot_out <- .add_extra_guide(plot_out, shape_by, size_by)
    } else if(layout == "heatmap"){
        # Get bar options
        raster_args <- .get_bar_args(
            colour_by = colour_by, alpha = 1, add_border = FALSE)
        # Add bar layout
        plot_out <- plot_out +
            do.call(geom_raster, raster_args$args) +
            scale_fill_distiller(palette = "RdYlBu", name = colour_by) +
            scale_y_discrete(expand = c(0,0))
        # Add scale. If numeric, add continuous scaling, if discrete, add
        # discrete scaling.
        if(is.factor(plot_data$X)){
            plot_out <- plot_out +
                scale_x_discrete(expand = c(0,0))
        } else {
            plot_out <- plot_out +
                scale_x_continuous(expand = c(0,0), n.breaks = 7L)
        }
    } else {
        stop("incompatible layout.", call. = FALSE)
    }
    # Adjust theme
    plot_out <- plot_out +
        theme_classic()
    # Remove legend if specified
    plot_out <- .add_legend(plot_out, add_legend)
    # Flip plot if specified
    plot_out <- .flip_plot(
        plot_out, flipped, add_x_text = TRUE, angle_x_text = FALSE)
    return(plot_out)
}
